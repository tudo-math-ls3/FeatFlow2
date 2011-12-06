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
!# 9.) btree_insertIntoTree = btree_insertIntoTreeDble /
!#                            btree_insertIntoTreeSngl /
!#                            btree_insertIntoTreeInt
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
module binarytree
  use fsystem
  use storage
  use genoutput
  implicit none

  private
  public :: t_btree
  public :: btree_createTree
  public :: btree_releaseTree
  public :: btree_resizeTree
  public :: btree_copyToTree
  public :: btree_copyFromTreeKey
  public :: btree_copyFromTreeDble
  public :: btree_copyFromTreeSngl
  public :: btree_copyFromTreeInt
  public :: btree_insertIntoTree
  public :: btree_deleteFromTree
  public :: btree_searchInTree
  public :: btree_getItemInTree
  public :: btree_printTree
  public :: btree_getHeight
  public :: btree_infoTree
  public :: btree_duplicateTree
  public :: btree_restoreTree

!<constants>

!<constantblock description="Global flags for tree output">

  ! Tag for preorder traversal
  integer, parameter, public :: BTREE_PREORDER  = 0

  ! Tag for inorder traversal
  integer, parameter, public :: BTREE_INORDER = 1

  ! Tag for postorder traversal
  integer, parameter, public :: BTREE_POSTORDER = 2

!</constantblock>

!<constantblock description="Global flags for tree operations">

  ! Identifier for "not found in tree"
  integer, parameter, public :: BTREE_NOT_FOUND = -1

  ! Identifier for "found in tree"
  integer, parameter, public :: BTREE_FOUND = 0
  
  ! Identifier for "copy to tree"
  integer, parameter, public :: BTREE_COPY1 = 1

  ! Identifier for "copy from tree"
  integer, parameter, public :: BTREE_COPY2 = 2

!</constantblock>

!<constantblock description="Internal tags for list status">

  ! Tag for empty tree
  integer, parameter, public :: TNULL = 0

  ! Tag for root of tree
  integer, parameter, public :: TROOT = 0
  
  ! Tag for next free item
  integer, parameter :: TFREE = 0

  ! Tag for left child
  integer, parameter, public :: TLEFT = 0

  ! Tag for right child
  integer, parameter, public :: TRIGHT = 1

!</constantblock>
!</constants>

!<types>
!<typeblock>
  
  type t_btree
    ! Format-Tag:
    integer :: ctreeFormat = ST_NOHANDLE
    
    ! Current depth of the tree
    integer :: depth = 0

    ! Number of items that are currently stored in the tree
    integer :: NA = 0

    ! Total number of items that can be stored in the tree
    integer :: NNA = 0

    ! Total number of items that can initially be stored in the tree
    ! This information is needed to compute the growth of the tree
    ! after several resize operations
    integer :: NNA0 = 0

    ! Total number of resize operations
    integer :: NRESIZE = 0
    
    ! Dimension of the auxiliary Integer values to be stored
    integer :: isizeInt = 0

    ! Dimension of the auxiliary Double values to be stored
    integer :: isizeDble = 0

    ! Dimension of the auxiliary Single values to be stored
    integer :: isizeSngl = 0

    ! Factor by which the list is enlarged if new storage is allocate
    real(DP) :: dfactor = 1.5_DP

    ! Handle to the tree key data
    integer :: h_Key = ST_NOHANDLE

    ! Handle to the balance data
    integer :: h_Kbal = ST_NOHANDLE

    ! Handle to the path data
    integer :: h_Kpath = ST_NOHANDLE

    ! Handle to the children`s data
    integer :: h_Kchild = ST_NOHANDLE

    ! Handle to the tree auxiliary Integer data
    integer :: h_IData = ST_NOHANDLE

    ! Handle to the tree auxiliary Double data
    integer :: h_DData = ST_NOHANDLE

    ! Handle to the tree auxiliary Single data
    integer :: h_FData = ST_NOHANDLE

    ! Tree child structure
    ! NOTE: This array is introduced to increase performance. It
    ! should not be touched by the user. However, if the handle would
    ! be dereferenced for each operation such as search, delete,
    ! performance would be very poor.
    integer, dimension(:,:), pointer :: p_Kchild => null()

    ! Tree path structure
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:), pointer :: p_Kpath => null()

    ! Tree balance structure
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:), pointer :: p_Kbal => null()

    ! Tree key data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    real(DP), dimension(:), pointer :: p_DKey => null()

    ! Tree key data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    real(SP), dimension(:), pointer :: p_FKey => null()

    ! Tree key data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:), pointer :: p_IKey => null()

    ! Tree data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    real(DP), dimension(:,:), pointer :: p_DData => null()

    ! Tree data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    real(SP), dimension(:,:), pointer :: p_FData => null()

    ! Tree data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:,:), pointer :: p_IData => null()
  end type t_btree

!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************
      
  interface btree_copyToTree
    module procedure btree_copyToTree_handle
    module procedure btree_copyToTree_arrayDble
    module procedure btree_copyToTree_arraySngl
    module procedure btree_copyToTree_arrayInt
  end interface

  interface btree_copyFromTreeKey
    module procedure btree_copyFromTreeKey_handle
    module procedure btree_copyFromTreeKey_arrayDble
    module procedure btree_copyFromTreeKey_arraySngl
    module procedure btree_copyFromTreeKey_arrayInt
  end interface

  interface btree_copyFromTreeDble
    module procedure btree_copyFromTree_handle
    module procedure btree_copyFromTree_arrayDble
    module procedure btree_copyFromTree_arrayDble2D
  end interface

  interface btree_copyFromTreeSngl
    module procedure btree_copyFromTree_handle
    module procedure btree_copyFromTree_arraySngl
    module procedure btree_copyFromTree_arraySngl2D
  end interface

  interface btree_copyFromTreeInt
    module procedure btree_copyFromTree_handle
    module procedure btree_copyFromTree_arrayInt
    module procedure btree_copyFromTree_arrayInt2D
  end interface
  
  interface btree_insertIntoTree
    module procedure btree_insertIntoTreeDble
    module procedure btree_insertIntoTreeSngl
    module procedure btree_insertIntoTreeInt
  end interface

  interface btree_deleteFromTree
    module procedure btree_deleteFromTreeDble
    module procedure btree_deleteFromTreeSngl
    module procedure btree_deleteFromTreeInt
  end interface
  
  interface btree_searchInTree
    module procedure btree_searchInTreeDble
    module procedure btree_searchInTreeSngl
    module procedure btree_searchInTreeInt
  end interface

  interface btree_getItemInTree
    module procedure btree_getItemInTreeDble
    module procedure btree_getItemInTreeSngl
    module procedure btree_getItemInTreeInt
  end interface

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

contains

  ! ***************************************************************************

!<subroutine>

  subroutine btree_createTree(rtree, nna, ctreeFormat, isizeDble,&
                              isizeSngl, isizeInt, dfactor)

!<description>
    ! This subroutine creates a new tree
!</description>

!<input>
    ! Total number of items that can be stored in the tree
    integer, intent(in) :: nna

    ! Format-tag. Type of tree format (Double,Single,Integer)
    integer, intent(in) :: ctreeFormat

    ! Dimension of the auxiliary Double values to be stored
    integer, intent(in) :: isizeDble

    ! Dimension of the auxiliary Single values to be stored
    integer, intent(in) :: isizeSngl

    ! Dimension of the auxiliary Integer values to be stored
    integer, intent(in) :: isizeInt

    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    real(DP), intent(in), optional :: dfactor
!</input>

!<output>
    ! binary tree
    type(t_btree), intent(out) :: rtree
!</output>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize1,Isize2
    
    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1_DP) rtree%dfactor=dfactor
    end if

    ! Set tree format
    rtree%ctreeFormat=ctreeFormat

    ! Initialise tree
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
    call storage_new('btree_createTree', 'Kchild', Isize1, Isize2,&
                     ST_INT, rtree%h_Kchild, ST_NEWBLOCK_ZERO)
    call storage_getbase_int2D(rtree%h_Kchild, rtree%p_Kchild)
    
    call storage_new('btree_createTree', 'Kpath', ceiling(1.441*LOG2(nna))+1,&
                     ST_INT, rtree%h_Kpath, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rtree%h_Kpath, rtree%p_Kpath)

    call storage_new('btree_createTree', 'Kbal', nna,&
                     ST_INT, rtree%h_Kbal, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rtree%h_Kbal, rtree%p_Kbal)

    ! Allocate memory for Key
    select case(rtree%ctreeFormat)
    case (ST_DOUBLE)
      call storage_new('btree_createTree', 'Key', nna,&
                       ST_DOUBLE, rtree%h_Key, ST_NEWBLOCK_ZERO)
      call storage_getbase_double(rtree%h_Key, rtree%p_DKey)
    case (ST_SINGLE)
      call storage_new('btree_createTree', 'Key', nna,&
                       ST_SINGLE, rtree%h_Key, ST_NEWBLOCK_ZERO)
      call storage_getbase_single(rtree%h_Key, rtree%p_FKey)
    case (ST_INT)
      call storage_new('btree_createTree', 'Key', nna,&
                       ST_INT, rtree%h_Key, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(rtree%h_Key, rtree%p_IKey)
    case DEFAULT
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_createTree')
      call sys_halt()
    end select

    ! Allocate memory fo auxiliary data
    if (isizeDble > 0) then
      Isize1 = (/isizeDble, nna/)
      call storage_new('btree_createTree', 'DData', Isize1,&
                       ST_DOUBLE, rtree%h_DData, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2D(rtree%h_DData, rtree%p_DData)
    end if

    if (isizeSngl > 0) then
      Isize1 = (/isizeSngl, nna/)
      call storage_new('btree_createTree', 'FData', Isize1,&
                       ST_SINGLE, rtree%h_FData, ST_NEWBLOCK_NOINIT)
      call storage_getbase_single2D(rtree%h_FData, rtree%p_FData)
    end if

    if (isizeInt > 0) then
      Isize1 = (/isizeInt, nna/)
      call storage_new('btree_createTree', 'IData', Isize1,&
                       ST_INT, rtree%h_IData, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2D(rtree%h_IData, rtree%p_IData)
    end if

    ! Initialise tree structure
    rtree%p_Kchild(TFREE,TROOT) = 1
  end subroutine btree_createTree

  ! ***************************************************************************

!<subroutine>
  
  subroutine btree_releaseTree(rtree)

!<description>
    ! This subroutine releases an existing tree
!</description>

!<inputoutput>
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>

    ! Release memory
    if (rtree%h_Key    .ne. ST_NOHANDLE) call storage_free(rtree%h_Key)
    if (rtree%h_Kbal   .ne. ST_NOHANDLE) call storage_free(rtree%h_Kbal)
    if (rtree%h_Kchild .ne. ST_NOHANDLE) call storage_free(rtree%h_Kchild)
    if (rtree%h_Kpath  .ne. ST_NOHANDLE) call storage_free(rtree%h_Kpath)

    if (rtree%h_DData  .ne. ST_NOHANDLE) call storage_free(rtree%h_DDATA)
    if (rtree%h_FData  .ne. ST_NOHANDLE) call storage_free(rtree%h_FDATA)
    if (rtree%h_IData  .ne. ST_NOHANDLE) call storage_free(rtree%h_IDATA)

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

    nullify(rtree%p_DKey, rtree%p_FKey, rtree%p_IKey)
    nullify(rtree%p_Kbal, rtree%p_Kpath, rtree%p_Kchild)
    nullify(rtree%p_DData, rtree%p_FData, rtree%p_IData)
  end subroutine btree_releaseTree

  ! ***************************************************************************

!<subroutine>
  
  subroutine btree_resizeTree(rtree, nna)

!<description>
    ! This subroutine reallocates memory for an existing tree
!</description>

!<input>
    ! New number of total items that can be stored in the tree
    integer, intent(in) :: nna
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>

    ! Set new size
    rtree%NNA = nna
    rtree%NRESIZE = rtree%NRESIZE+1

    call storage_realloc('btree_resizeTree', nna+1,&
                         rtree%h_Kchild, ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_int2D(rtree%h_Kchild, rtree%p_Kchild)

    call storage_realloc('btree_resizeTree', nna,&
                         rtree%h_Kbal, ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_int(rtree%h_Kbal, rtree%p_Kbal)

    call storage_realloc('btree_resizeTree',&
                         ceiling(1.441*LOG2(nna))+1, rtree%h_Kpath,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_int(rtree%h_Kpath, rtree%p_Kpath)

    call storage_realloc('btree_resizeTree', nna,&
                         rtree%h_Key, ST_NEWBLOCK_ZERO, .true.)
    select case(rtree%ctreeFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(rtree%h_Key, rtree%p_DKey)
    case (ST_SINGLE)
      call storage_getbase_single(rtree%h_Key, rtree%p_FKey)
    case (ST_INT)
      call storage_getbase_int(rtree%h_Key, rtree%p_IKey)
    case DEFAULT
      call output_line('Unsupported key data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_resizeTree')
      call sys_halt()
    end select

    ! Auxiliary data
    if (rtree%isizeDble > 0) then
      call storage_realloc('btree_resizeTree', nna,&
                           rtree%h_DData, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_double2D(rtree%h_DData, rtree%p_DData)
    end if

    if (rtree%isizeSngl > 0) then
      call storage_realloc('btree_resizeTree', nna,&
                           rtree%h_FData, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_single2D(rtree%h_FData, rtree%p_FData)
    end if

    if (rtree%isizeInt > 0) then
      call storage_realloc('btree_resizeTree', nna,&
                           rtree%h_IData, ST_NEWBLOCK_NOINIT, .true.)
      call storage_getbase_int2D(rtree%h_IData, rtree%p_IData)
    end if
  end subroutine btree_resizeTree

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyToTree_handle(h_Key, rtree, h_DData, h_FData, h_IData)

!<description>
    ! This subroutine copies the content of handles to the tree.
    ! The key handle is mandatory to build up the tree data
    ! structure. The optional handles h_XData can be used to attach
    ! auxiliary data to the tree which is "carried" through all operations.
!</description>

!<input>
    ! handle to the key
    integer, intent(in) :: h_Key

    ! OPTIONAL: handle to auxiliary Double data
    integer, intent(in), optional :: h_DData

    ! OPTIONAL: handle to auxiliary Single data
    integer, intent(in), optional :: h_FData

    ! OPTIONAL: handle to auxiliary Integer data
    integer, intent(in), optional :: h_IData
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DKey
    real(SP), dimension(:), pointer :: p_FKey
    integer,  dimension(:), pointer :: p_IKey
    integer :: j

    ! Transform the content of the key handle to the tree
    select case (rtree%ctreeFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_Key, p_DKey)
      do j = 1, size(p_DKey)
        call btree_insertIntoTree(rtree, p_DKey(j))
      end do
      
    case (ST_SINGLE)
      call storage_getbase_single(h_Key, p_FKey)
      do j = 1, size(p_FKey)
        call btree_insertIntoTree(rtree, p_FKey(j))
      end do
      
    case (ST_INT)
      call storage_getbase_int(h_Key, p_IKey)
      do j = 1, size(p_IKey)
        call btree_insertIntoTree(rtree, p_IKey(j))
      end do
      
    case DEFAULT
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_handle')
      call sys_halt()
    end select
    
    ! Transform the content of the auxiliary handles to the tree.
    ! Note that the sorting key was inserted linearly. In principle,
    ! only the tree data structure needs to be generated but the
    ! actual key values are still stored in sequential order. Hence,
    ! it suffices to "copy" the handles of the auxiliary data to the
    ! handles of the tree structure
    
    if (present(h_DData) .and. rtree%isizeDble > 0) &
        call storage_copy(h_DData, rtree%h_DData)
    if (present(h_FData) .and. rtree%isizeSngl > 0) &
        call storage_copy(h_FData, rtree%h_FData)
    if (present(h_IData) .and. rtree%isizeInt  > 0) &
        call storage_copy(h_IData, rtree%h_IData)
  end subroutine btree_copyToTree_handle

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyToTree_arrayDble(p_DKey, rtree, p_DData, p_FData, p_IData)

!<description>
    ! This subroutine copies the content of arrays to the tree
    ! and makes use of a double-valued key. The optional arrays
    ! p_XData can be used to attach auxiliary data to the tree
    ! which is "carried" through all operations.
!</description>

!<input>
    ! double-valued key
    real(DP), dimension(:), intent(in) :: p_DKey

    ! OPTIONAL: auxiliary Double data
    real(DP), dimension(:,:), optional :: p_DData

    ! OPTIONAL: auxiliary Single data
    real(SP), dimension(:,:), optional :: p_FData

    ! OPTIONAL: auxiliary Integer data
    integer,  dimension(:,:), optional :: p_IData
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: j

    ! Check if tree has correct data format
    if (rtree%ctreeFormat .ne. ST_DOUBLE) then
      call output_line('Invalid data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
      call sys_halt()
    end if

    ! Build-up treee structure
    do j = 1, size(p_DKey)
      call btree_insertIntoTree(rtree, p_DKey(j))
    end do

    ! Copy auxiliary arrays to the tree.
    if (present(p_DData) .and. rtree%isizeDble > 0) then
      call storage_getsize(rtree%h_DData, Isize)
      if ((rtree%isizeDble .ne. size(p_DData,1)) .or.&
          (rtree%NA .ne. size(p_DData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arrayDble',&
                             rtree%NNA, rtree%h_DData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2D(rtree%h_DData, rtree%p_DData)
      end if
      ! Copy data
      rtree%p_DData(:,1:rtree%NA) = p_DData
    end if

    if (present(p_FData) .and. rtree%isizeSngl > 0) then
      call storage_getsize(rtree%h_FData, Isize)
      if ((rtree%isizeSngl .ne. size(p_FData,1)) .or.&
          (rtree%NA .ne. size(p_FData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arrayDble',&
                             rtree%NNA, rtree%h_FData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_single2D(rtree%h_FData, rtree%p_FData)
      end if
      ! Copy data
      rtree%p_FData(:,1:rtree%NA) = p_FData
    end if

    if (present(p_IData) .and. rtree%isizeInt > 0) then
      call storage_getsize(rtree%h_IData, Isize)
      if ((rtree%isizeInt .ne. size(p_IData,1)) .or.&
          (rtree%NA .ne. size(p_IData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arrayDble',&
                             rtree%NNA, rtree%h_IData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2D(rtree%h_IData, rtree%p_IData)
      end if
      ! Copy data
      rtree%p_IData(:,1:rtree%NA) = p_IData
    end if
    
  end subroutine btree_copyToTree_arrayDble

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyToTree_arraySngl(p_FKey, rtree, p_DData, p_FData, p_IData)

!<description>
    ! This subroutine copies the content of arrays to the tree
    ! and makes use of a single-valued key. The optional arrays
    ! p_XData can be used to attach auxiliary data to the tree
    ! which is "carried" through all operations.
!</description>

!<input>
    ! single-valued key
    real(SP), dimension(:), intent(in) :: p_FKey

    ! OPTIONAL: auxiliary Double data
    real(DP), dimension(:,:), optional :: p_DData

    ! OPTIONAL: auxiliary Single data
    real(SP), dimension(:,:), optional :: p_FData

    ! OPTIONAL: auxiliary Integer data
    integer,  dimension(:,:), optional :: p_IData
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: j

    ! Check if tree has correct data format
    if (rtree%ctreeFormat .ne. ST_SINGLE) then
      call output_line('Invalid data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
      call sys_halt()
    end if

    ! Build-up treee structure
    do j = 1, size(p_FKey)
      call btree_insertIntoTree(rtree, p_FKey(j))
    end do

    ! Copy auxiliary arrays to the tree.
    if (present(p_DData) .and. rtree%isizeDble > 0) then
      call storage_getsize(rtree%h_DData, Isize)
      if ((rtree%isizeDble .ne. size(p_DData,1)) .or.&
          (rtree%NA .ne. size(p_DData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arraySngl',&
                             rtree%NNA, rtree%h_DData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2D(rtree%h_DData, rtree%p_DData)
      end if
      ! Copy data
      rtree%p_DData(:,1:rtree%NA) = p_DData
    end if
    
    if (present(p_FData) .and. rtree%isizeSngl > 0) then
      call storage_getsize(rtree%h_FData, Isize)
      if ((rtree%isizeSngl .ne. size(p_FData,1)) .or.&
          (rtree%NA .ne. size(p_FData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arraySngl',&
                             rtree%NNA, rtree%h_FData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_single2D(rtree%h_FData, rtree%p_FData)
      end if
      ! Copy data
      rtree%p_FData(:,1:rtree%NA) = p_FData
    end if

    if (present(p_IData) .and. rtree%isizeInt > 0) then
      call storage_getsize(rtree%h_IData, Isize)
      if ((rtree%isizeInt .ne. size(p_IData,1)) .or.&
          (rtree%NA .ne. size(p_IData,2))) then
        call output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arraySngl',&
            rtree%NNA,rtree%h_IData,ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2D(rtree%h_IData,rtree%p_IData)
      end if
      ! Copy data
      rtree%p_IData(:,1:rtree%NA) = p_IData
    end if

  end subroutine btree_copyToTree_arraySngl

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyToTree_arrayInt(p_IKey, rtree, p_DData, p_FData, p_IData)

!<description>
    ! This subroutine copies the content of arrays to the tree
    ! and makes use of an integer-valued key. The optional arrays
    ! p_XData can be used to attach auxiliary data to the tree
    ! which is "carried" through all operations.
!</description>

!<input>
    ! integer-valued key
    integer, dimension(:), intent(in) :: p_IKey

    ! OPTIONAL: auxiliary Double data
    real(DP), dimension(:,:), optional :: p_DData

    ! OPTIONAL: auxiliary Single data
    real(SP), dimension(:,:), optional :: p_FData

    ! OPTIONAL: auxiliary Integer data
    integer,  dimension(:,:), optional :: p_IData
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: j

    ! Check if tree has correct data format
    if (rtree%ctreeFormat .ne. ST_INT) then
      call output_line('Invalid data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
      call sys_halt()
    end if

    ! Build-up treee structure
    do j = 1, size(p_IKey)
      call btree_insertIntoTree(rtree, p_IKey(j))
    end do

    ! Copy auxiliary arrays to the tree.
    if (present(p_DData) .and. rtree%isizeDble > 0) then
      call storage_getsize(rtree%h_DData, Isize)
      if ((rtree%isizeDble .ne. size(p_DData,1)) .or.&
          (rtree%NA .ne. size(p_DData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arrayInt',&
                             rtree%NNA, rtree%h_DData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2D(rtree%h_DData, rtree%p_DData)
      end if
      ! Copy data
      rtree%p_DData(:,1:rtree%NA) = p_DData
    end if
    
    if (present(p_FData) .and. rtree%isizeSngl > 0) then
      call storage_getsize(rtree%h_FData, Isize)
      if ((rtree%isizeSngl .ne. size(p_FData,1)) .or.&
          (rtree%NA .ne. size(p_FData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arrayInt',&
                             rtree%NNA, rtree%h_FData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_single2D(rtree%h_FData, rtree%p_FData)
      end if
      ! Copy data
      rtree%p_FData(:,1:rtree%NA) = p_FData
    end if

    if (present(p_IData) .and. rtree%isizeInt > 0) then
      call storage_getsize(rtree%h_IData, Isize)
      if ((rtree%isizeInt .ne. size(p_IData,1)) .or.&
          (rtree%NA .ne. size(p_IData,2))) then
        call output_line('Invalid auxiliary data!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
        call sys_halt()
      end if
      ! Do we have to reallocate memory?
      if (Isize(2) .ne. rtree%NNA) then
        call storage_realloc('btree_copyToTree_arrayInt',&
                             rtree%NNA, rtree%h_IData, ST_NEWBLOCK_NOINIT)
        call storage_getbase_int2D(rtree%h_IData, rtree%p_IData)
      end if
      ! Copy data
      rtree%p_IData(:,1:rtree%NA) = p_IData
    end if

  end subroutine btree_copyToTree_arrayInt

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTreeKey_handle(rtree, h_Key, corder)

!<description>
    ! This subroutine copies the key of the tree to the handle h_Key
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! handle to the key
    integer, intent(inout) :: h_Key
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DKey
    real(SP), dimension(:), pointer :: p_FKey
    integer,  dimension(:), pointer :: p_IKey
    integer :: isize
    integer :: j,iorder

    ! Check if handle needs to be (re-)allocated
    if (h_Key .eq. ST_NOHANDLE) then
      call storage_new('btree_copyFromTreeKey_handle', 'Key', rtree%NA,&
                       rtree%ctreeFormat, h_Key, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_Key,isize)
      if (isize < rtree%NA) then
        call storage_realloc('btree_copyFromTreeKey_handle', rtree%NA,&
                             h_Key, ST_NEWBLOCK_ZERO, .false.)
      end if
    end if

    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    ! What kind of key is used?
    select case(rtree%ctreeFormat)
    case (ST_DOUBLE)
      call storage_getbase_double(h_Key, p_DKey); j=0
      select case(iorder)
      case (BTREE_PREORDER)
        call preorderKeyDble(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        call inorderKeyDble(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        call postorderKeyDble(rtree%p_Kchild(TRIGHT, TROOT))
      end select
      
    case (ST_SINGLE)
      call storage_getbase_single(h_Key, p_FKey); j=0
      select case(iorder)
      case (BTREE_PREORDER)
        call preorderKeySngl(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        call inorderKeySngl(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        call postorderKeySngl(rtree%p_Kchild(TRIGHT, TROOT))
    end select

    case (ST_INT)
      call storage_getbase_int(h_Key, p_IKey); j=0
      select case(iorder)
      case (BTREE_PREORDER)
        call preorderKeyInt(rtree%p_Kchild(TRIGHT,TROOT))
      case (BTREE_INORDER)
        call inorderKeyInt(rtree%p_Kchild(TRIGHT,TROOT))
      case (BTREE_POSTORDER)
        call postorderKeyInt(rtree%p_Kchild(TRIGHT,TROOT))
      end select
      
    case DEFAULT
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Copy the content of the Double key to an array
    
    recursive subroutine preorderKeyDble(i)
      integer, intent(in) :: i
      
            j=j+1; p_DKey(j)=rtree%p_DKey(i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderKeyDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderKeyDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderKeyDble

    !**************************************************************
    ! Copy the content of the Double key to an array
    
    recursive subroutine inorderKeyDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderKeyDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DKey(j)=rtree%p_DKey(i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderKeyDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderKeyDble

    !**************************************************************
    ! Copy the content of the Double key to an array
    
    recursive subroutine postorderKeyDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderKeyDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderKeyDble(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_DKey(j)=rtree%p_DKey(i)
    end subroutine postorderKeyDble
       
    !**************************************************************
    ! Copy the content of the Single key to an array

    recursive subroutine preorderKeySngl(i)
      integer, intent(in) :: i
      
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderKeySngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderKeySngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderKeySngl

    !**************************************************************
    ! Copy the content of the Single key to an array

    recursive subroutine inorderKeySngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderKeySngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderKeySngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderKeySngl

    !**************************************************************
    ! Copy the content of the Single key to an array

    recursive subroutine postorderKeySngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderKeySngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderKeySngl(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
    end subroutine postorderKeySngl

    !**************************************************************
    ! Copy the content of the Integer key to an array
    
    recursive subroutine preorderKeyInt(i)
      integer, intent(in) :: i
      
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderKeyInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderKeyInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderKeyInt
    
    !**************************************************************
    ! Copy the content of the Integer key to an array

    recursive subroutine inorderKeyInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderKeyInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderKeyInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderKeyInt
    
    !**************************************************************
    ! Copy the content of the Integer key to an array
    
    recursive subroutine postorderKeyInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderKeyInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderKeyInt(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
    end subroutine postorderKeyInt
  end subroutine btree_copyFromTreeKey_handle

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTreeKey_arrayDble(rtree, p_DKey, corder)

!<description>
    ! This subroutine copies the key of the tree to the double-valued array p_DKey
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! double-valued array
    real(DP), dimension(:), intent(inout) :: p_DKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: j,iorder

    ! Check data format
    if (rtree%ctreeFormat .ne. ST_DOUBLE) then
      call output_line('Invalid data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayDble')
      call sys_halt()
    end if

    ! Check size of array
    if (size(p_DKey) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayDble')
      call sys_halt()
    end if
    
    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder
    
    select case(iorder)
    case (BTREE_PREORDER)
      j=0; call preorderKeyDble(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_INORDER)
      j=0; call inorderKeyDble(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_POSTORDER)
      j=0; call postorderKeyDble(rtree%p_Kchild(TRIGHT, TROOT))
    end select

  contains
    
    ! Here, the real working routine follows
    
    !**************************************************************
    ! Copy the content of the Double key to an array
    
    recursive subroutine preorderKeyDble(i)
      integer, intent(in) :: i
      
      j=j+1; p_DKey(j)=rtree%p_DKey(i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderKeyDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderKeyDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderKeyDble

    !**************************************************************
    ! Copy the content of the Double key to an array
    
    recursive subroutine inorderKeyDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderKeyDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DKey(j)=rtree%p_DKey(i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderKeyDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderKeyDble

    !**************************************************************
    ! Copy the content of the Double key to an array
    
    recursive subroutine postorderKeyDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderKeyDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderKeyDble(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_DKey(j)=rtree%p_DKey(i)
    end subroutine postorderKeyDble
  end subroutine btree_copyFromTreeKey_arrayDble

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTreeKey_arraySngl(rtree, p_FKey, corder)

!<description>
    ! This subroutine copies the key of the tree to the single-valued array p_FKey
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! single-valued array
    real(SP), dimension(:), intent(inout) :: p_FKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: j,iorder

    ! Check data format
    if (rtree%ctreeFormat .ne. ST_SINGLE) then
      call output_line('Invalid data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arraySngl')
      call sys_halt()
    end if

    ! Check size of array
    if (size(p_FKey) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arraySngl')
      call sys_halt()
    end if
    
    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    select case(iorder)
    case(BTREE_PREORDER)
      j=0; call preorderKeySngl(rtree%p_Kchild(TRIGHT, TROOT))
    case(BTREE_INORDER)
      j=0; call inorderKeySngl(rtree%p_Kchild(TRIGHT, TROOT))
    case(BTREE_POSTORDER)
      j=0; call postorderKeySngl(rtree%p_Kchild(TRIGHT, TROOT))
    end select

  contains

    ! Here, the real working routine follows

    !**************************************************************
    ! Copy the content of the Single key to an array

    recursive subroutine preorderKeySngl(i)
      integer, intent(in) :: i
      
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderKeySngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderKeySngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderKeySngl

    !**************************************************************
    ! Copy the content of the Single key to an array

    recursive subroutine inorderKeySngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderKeySngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderKeySngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderKeySngl

    !**************************************************************
    ! Copy the content of the Single key to an array

    recursive subroutine postorderKeySngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderKeySngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderKeySngl(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
    end subroutine postorderKeySngl
  end subroutine btree_copyFromTreeKey_arraySngl

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTreeKey_arrayInt(rtree, p_IKey, corder)

!<description>
    ! This subroutine copies the key of the tree to the integer-valued array p_IKey
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! integer-valued array
    integer, dimension(:), intent(inout) :: p_IKey
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: j,iorder

    ! Check data format
    if (rtree%ctreeFormat .ne. ST_INT) then
      call output_line('Invalid data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayInt')
      call sys_halt()
    end if

    ! Check size of array
    if (size(p_IKey) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayInt')
      call sys_halt()
    end if
    
    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    select case(iorder)
    case (BTREE_PREORDER)
      j=0; call preorderKeyInt(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_INORDER)
      j=0; call inorderKeyInt(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_POSTORDER)
      j=0; call postorderKeyInt(rtree%p_Kchild(TRIGHT, TROOT))
    end select

  contains
    
    ! Here, the real working routine follows

    !**************************************************************
    ! Copy the content of the Integer key to an array
    
    recursive subroutine preorderKeyInt(i)
      integer, intent(in) :: i
      
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderKeyInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderKeyInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderKeyInt

    !**************************************************************
    ! Copy the content of the Integer key to an array
    
    recursive subroutine inorderKeyInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderKeyInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderKeyInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderKeyInt

    !**************************************************************
    ! Copy the content of the Integer key to an array
    
    recursive subroutine postorderKeyInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderKeyInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderKeyInt(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
    end subroutine postorderKeyInt
  end subroutine btree_copyFromTreeKey_arrayInt

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTree_handle(rtree, ctype, h_Data, mask, corder)

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
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! type of data
    integer, intent(in) :: ctype

    ! OPTIONAL: mask of components to be copied
    integer, dimension(:), intent(in), optional :: mask

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! handle to the data
    integer, intent(inout) :: h_Data
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DData
    real(SP), dimension(:), pointer :: p_FData
    integer,  dimension(:), pointer :: p_IData
    real(DP), dimension(:,:), pointer :: p_DData2D
    real(SP), dimension(:,:), pointer :: p_FData2D
    integer,  dimension(:,:), pointer :: p_IData2D
    integer, dimension(2) :: Isize
    integer :: iisize
    integer :: idimension,idatatype
    
    ! Check if handle is already associated
    if (h_Data .eq. ST_NOHANDLE) then
      
      ! Get second working dimension
      Isize(2)=rtree%NA

      ! What kind of data should be copied?
      select case(ctype)
      case (ST_DOUBLE)

        ! Get first working dimension
        Isize(1)=rtree%isizeDble
        if (present(mask)) Isize(1) = size(mask)
        
        if (Isize(1) > 1) then
          call storage_new('btree_copyFromTree_handle', 'DData', Isize,&
                           ST_DOUBLE,h_Data,ST_NEWBLOCK_NOINIT)
          call storage_getbase_double2D(h_Data, p_DData2D)
          call btree_copyFromTree_arrayDble2D(rtree, p_DData2D, mask, corder)
        elseif (present(mask)) then
          call storage_new('btree_copyFromTree_handle', 'DData', Isize(2),&
                           ST_DOUBLE, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_double(h_Data, p_DData)
          call btree_copyFromTree_arrayDble(rtree, p_DData, mask(1), corder)
        else
          call storage_new('btree_copyFromTree_handle', 'DData', Isize(2),&
                           ST_DOUBLE, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_double(h_Data,p_DData)
          call btree_copyFromTree_arrayDble(rtree, p_DData, 1, corder)
        end if

      case (ST_SINGLE)
        
        ! Get first working dimension
        Isize(1)=rtree%isizeSngl
        if (present(mask)) Isize(1) = size(mask)

        if (Isize(1) > 1) then
          call storage_new('btree_copyFromTree_handle', 'FData', Isize,&
                           ST_SINGLE, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_single2D(h_Data, p_FData2D)
          call btree_copyFromTree_arraySngl2D(rtree, p_FData2D, mask, corder)
        elseif (present(mask)) then
          call storage_new('btree_copyFromTree_handle', 'FData', Isize(2),&
                           ST_SINGLE, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_single(h_Data, p_FData)
          call btree_copyFromTree_arraySngl(rtree, p_FData, mask(1), corder)
        else
          call storage_new('btree_copyFromTree_handle', 'FData', Isize(2),&
                           ST_SINGLE, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_single(h_Data, p_FData)
          call btree_copyFromTree_arraySngl(rtree, p_FData, 1, corder)
        end if

      case (ST_INT)
        
        ! Get first working dimension
        Isize(1)=rtree%isizeInt
        if (present(mask)) Isize(1) = size(mask)

        if (Isize(1) > 1) then
          call storage_new('btree_copyFromTree_handle', 'IData', Isize,&
                           ST_INT, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int2D(h_Data, p_IData2D)
          call btree_copyFromTree_arrayInt2D(rtree, p_IData2D, mask, corder)
        elseif (present(mask)) then
          call storage_new('btree_copyFromTree_handle', 'IData', Isize(2),&
                           ST_INT, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int(h_Data, p_IData)
          call btree_copyFromTree_arrayInt(rtree, p_IData, mask(1), corder)
        else
          call storage_new('btree_copyFromTree_handle', 'IData', Isize(2),&
                           ST_INT, h_Data, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int(h_Data, p_IData)
          call btree_copyFromTree_arrayInt(rtree, p_IData, 1, corder)
        end if
        
      case DEFAULT
        call output_line('Unsupported data format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
        call sys_halt()
      end select

    else   ! The handle is already associated

      ! Check if data type is valid
      call storage_getdatatype(h_Data, idatatype)
      if (idatatype .ne. ctype) then
        call output_line('Data type mismatch!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
        call sys_halt()
      end if

      ! Are we 1D- or 2D-array?
      call storage_getdimension(h_Data, idimension)
      select case(idimension)

      case (1)   !!! 1D-ARRAY  !!!

        ! Do we have to reallocate the array?
        call storage_getsize(h_Data, iisize)
        if (iisize < rtree%NA)&
            call storage_realloc('btree_copyFromTree_handle', rtree%NA,&
                                 h_Data, ST_NEWBLOCK_NOINIT, .false.)

        ! What kind of data type are we?
        select case(ctype)
        case (ST_DOUBLE)
          call storage_getbase_double(h_Data, p_DData)

          ! Which component should be copied
          if (present(mask)) then
            if (size(mask) > 1) then
              call output_line('For 1D-array mask can only have one entry!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
              call sys_halt()
            end if
            call btree_copyFromTree_arrayDble(rtree, p_DData, mask(1), corder)
          elseif (rtree%isizeDble .eq. 1) then
            call btree_copyFromTree_arrayDble(rtree, p_DData, 1, corder)
          else
            call output_line('A 1D-array was given but there are more than ' // &
                             'one components in the tree!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
            call sys_halt()
          end if
                    
        case (ST_SINGLE)
          call storage_getbase_single(h_Data, p_FData)

          ! Which component should be copied
          if (present(mask)) then
            if (size(mask) > 1) then
              call output_line('For 1D-array mask can only have one entry!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
              call sys_halt()
            end if
            call btree_copyFromTree_arraySngl(rtree, p_FData, mask(1), corder)
          elseif (rtree%isizeSngl .eq. 1) then
            call btree_copyFromTree_arraySngl(rtree, p_FData, 1, corder)
          else
            call output_line('A 1D-array was given but there are more than ' // &
                             'one components in the tree!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
            call sys_halt()
          end if
          
        case (ST_INT)
          call storage_getbase_int(h_Data, p_IData)

          ! Which component should be copied
          if (present(mask)) then
            if (size(mask) > 1) then
              call output_line('For 1D-array mask can only have one entry!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
              call sys_halt()
            end if
            call btree_copyFromTree_arrayInt(rtree, p_IData, mask(1), corder)
          elseif (rtree%isizeInt .eq. 1) then
            call btree_copyFromTree_arrayInt(rtree, p_IData, 1, corder)
          else
            call output_line('A 1D-array was given but there are more than ' // &
                             'one components in the tree!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
            call sys_halt()
          end if

        case DEFAULT
          call output_line('Unsupported data type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
          call sys_halt()
        end select

      case (2)   !!! 2D-ARRAY !!!

        ! Do we have to reallocate the array?
        call storage_getsize(h_Data, Isize)
        if (Isize(2) < rtree%NA)&
            call storage_realloc('btree_copyFromTree_handle', rtree%NA,&
                                 h_Data, ST_NEWBLOCK_NOINIT, .false.)

        ! What kind of data type are we
        select case(ctype)
        case (ST_DOUBLE)
          call storage_getbase_double2D(h_Data, p_DData2D)
          call btree_copyFromTree_arrayDble2D(rtree, p_DData2D, mask, corder)

        case (ST_SINGLE)
          call storage_getbase_single2D(h_Data, p_FData2D)
          call btree_copyFromTree_arraySngl2D(rtree, p_FData2D, mask, corder)

        case (ST_INT)
          call storage_getbase_int2D(h_Data, p_IData2D)
          call btree_copyFromTree_arrayInt2D(rtree, p_IData2D, mask, corder)

        case DEFAULT
          call output_line('Unsupported data type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
          call sys_halt()
        end select

      case DEFAULT
        call output_line('Unsupported data dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
        call sys_halt()
      end select
    end if
  end subroutine btree_copyFromTree_handle

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTree_arrayDble(rtree, p_Ddata, mask, corder)

!<description>
    ! This subroutine copies one component of the auxiliary double data
    ! attached to the tree to the given array. The component must be specified
    ! by the mask, e.g. mask=4 copies the fourth component of the data.
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! mask of component to be copied
    integer, intent(in) :: mask
    
    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! double data array
    real(DP), dimension(:),intent(inout) :: p_DData
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: j,iorder

    ! Check if auxiliary double data is available
    if (rtree%isizeDble .eq. 0) then
      call output_line('No double data available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble')
      call sys_halt()
    end if

    ! Check if given array is large enough in its second dimension
    if (size(p_DData) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble')
      call sys_halt()
    end if

    ! Check if mask is valid
    if (mask < 1 .or. mask > rtree%isizeDble) then
      call output_line('Invalid mask!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble')
      call sys_halt()
    end if
    
    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder
    
    select case(iorder)
    case (BTREE_PREORDER)
      j=0; call preorderDataDble(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_INORDER)
      j=0; call inorderDataDble(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_POSTORDER)
      j=0; call postorderDataDble(rtree%p_Kchild(TRIGHT, TROOT))
    end select

  contains

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Double data to an array
    
    recursive subroutine preorderDataDble(i)
      integer, intent(in) :: i

      j=j+1; p_DData(j)=rtree%p_DData(mask,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataDble

    !**************************************************************
    ! Copy the content of the Double data to an array
    
    recursive subroutine inorderDataDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DData(j)=rtree%p_DData(mask,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataDble

    !**************************************************************
    ! Copy the content of the Double data to an array
    
    recursive subroutine postorderDataDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataDble(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_DData(j)=rtree%p_DData(mask,i)
    end subroutine postorderDataDble

  end subroutine btree_copyFromTree_arrayDble

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTree_arrayDble2D(rtree, p_Ddata, mask, corder)

!<description>
    ! This subroutine copies some (part of) the auxiliary double data
    ! attached to the tree to the given array. The part of the data to
    ! be copied can be specified by the optional mask, e.g.
    ! mask=(/1,2,7/) only copies the first, second and seventh
    ! components of the auxiliary data.
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! OPTIONAL: mask of components to be copied
    integer, dimension(:), intent(in), optional :: mask

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! double data array
    real(DP), dimension(:,:),intent(inout) :: p_DData
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: j,iorder

    ! Check if auxiliary double data is available
    if (rtree%isizeDble .eq. 0) then
      call output_line('No double data available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
      call sys_halt()
    end if

    ! Check if given array is large enough in its second dimension
    Isize = shape(p_Ddata)
    if (Isize(2) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
      call sys_halt()
    end if

    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    ! Copy the content to the array
    if (present(mask)) then
      ! Check if given array has correct size in its first dimension
      if (Isize(1) .ne. size(mask)) then
        call output_line('Array dimensions do not match mask!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
        call sys_halt()
      end if

      ! Check if mask is valid
      if (any(mask < 1) .or. any(mask > rtree%isizeDble)) then
        call output_line('Invalid mask!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
        call sys_halt()
      end if
      select case(iorder)
      case (BTREE_PREORDER)
        j=0; call preorderDataDbleMask(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        j=0; call inorderDataDbleMask(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        j=0; call postorderDataDbleMask(rtree%p_Kchild(TRIGHT, TROOT))
      end select
    else
      ! Check if given array has correct size in its first dimension
      if (Isize(1) .ne. rtree%isizeDble) then
        call output_line('Array too small!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
        call sys_halt()
      end if
      select case(iorder)
      case (BTREE_PREORDER)
        j=0; call preorderDataDble(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        j=0; call inorderDataDble(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        j=0; call postorderDataDble(rtree%p_Kchild(TRIGHT, TROOT))
      end select
    end if

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Double data to an array

    recursive subroutine preorderDataDble(i)
      integer, intent(in) :: i
      
      j=j+1; p_DData(:,j)=rtree%p_DData(:,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataDble

    !**************************************************************
    ! Copy the content of the Double data to an array

    recursive subroutine inorderDataDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DData(:,j)=rtree%p_DData(:,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataDble

    !**************************************************************
    ! Copy the content of the Double data to an array

    recursive subroutine postorderDataDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataDble(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_DData(:,j)=rtree%p_DData(:,i)
    end subroutine postorderDataDble

    !**************************************************************
    ! Copy the content of the Double data to an array (masked)

    recursive subroutine preorderDataDbleMask(i)
      integer, intent(in) :: i
      
      j=j+1; p_DData(:,j)=rtree%p_DData(mask,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataDbleMask(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataDbleMask(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataDbleMask

    !**************************************************************
    ! Copy the content of the Double data to an array (masked)

    recursive subroutine inorderDataDbleMask(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataDbleMask(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DData(:,j)=rtree%p_DData(mask,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataDbleMask(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataDbleMask

    !**************************************************************
    ! Copy the content of the Double data to an array (masked)

    recursive subroutine postorderDataDbleMask(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataDbleMask(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataDbleMask(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_DData(:,j)=rtree%p_DData(mask,i)
    end subroutine postorderDataDbleMask

  end subroutine btree_copyFromTree_arrayDble2D

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTree_arraySngl(rtree, p_FData, mask, corder)

!<description>
    ! This subroutine copies one component of the auxiliary single data
    ! attached to the tree to the given array. The component must be specified
    ! by the mask, e.g. mask=4 copies the fourth component of the data.
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! mask of component to be copied
    integer, intent(in) :: mask

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! single data array
    real(SP), dimension(:),intent(inout) :: p_FData
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: j,iorder

    ! Check if auxiliary single data is available
    if (rtree%isizeSngl .eq. 0) then
      call output_line('No single data available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl')
      call sys_halt()
    end if

    ! Check if given array is large enough in its second dimension
    if (size(p_FData) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl')
      call sys_halt()
    end if

    ! Check if mask is valid
    if (mask < 1 .or. mask > rtree%isizeSngl) then
      call output_line('Invalid mask!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl')
      call sys_halt()
    end if
    
     ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    select case(iorder)
    case (BTREE_PREORDER)
      j=0; call preorderDataSngl(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_INORDER)
      j=0; call inorderDataSngl(rtree%p_Kchild(TRIGHT, TROOT))
    case (BTREE_POSTORDER)
      j=0; call postorderDataSngl(rtree%p_Kchild(TRIGHT, TROOT))
    end select

  contains

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Single data to an array
    
    recursive subroutine preorderDataSngl(i)
      integer, intent(in) :: i
      
      j=j+1; p_FData(j)=rtree%p_FData(mask,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataSngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataSngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataSngl

    !**************************************************************
    ! Copy the content of the Single data to an array
    
    recursive subroutine inorderDataSngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataSngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FData(j)=rtree%p_FData(mask,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataSngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataSngl

    !**************************************************************
    ! Copy the content of the Single data to an array
    
    recursive subroutine postorderDataSngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataSngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataSngl(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_FData(j)=rtree%p_FData(mask,i)
    end subroutine postorderDataSngl

  end subroutine btree_copyFromTree_arraySngl

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTree_arraySngl2D(rtree, p_Fdata, mask, corder)

!<description>
    ! This subroutine copies some (part of) the auxiliary single data
    ! attached to the tree to the given array. The part of the data to
    ! be copied can be specified by the optional mask, e.g.
    ! mask=(/1,2,7/) only copies the first, second and seventh
    ! components of the auxiliary data.
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! OPTIONAL: mask of components to be copied
    integer, dimension(:), intent(in), optional :: mask

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! single data array
    real(SP), dimension(:,:),intent(inout) :: p_FData
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: j,iorder

    ! Check if auxiliary single data is available
    if (rtree%isizeSngl .eq. 0) then
      call output_line('No single data available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
      call sys_halt()
    end if
    
    ! Check if given array is large enough in its second dimension
    Isize = shape(p_Fdata)
    if (Isize(2) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
      call sys_halt()
    end if

    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    ! Copy the content to the array
    if (present(mask)) then
      ! Check if given array has correct size in its first dimension
      if (Isize(1) .ne. size(mask)) then
        call output_line('Array dimensions do not match mask!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
        call sys_halt()
      end if

      ! Check if mask is valid
      if (any(mask < 1) .or. any(mask > rtree%isizeSngl)) then
        call output_line('Invalid mask!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
        call sys_halt()
      end if
      select case(iorder)
      case (BTREE_PREORDER)
        j=0; call preorderDataSnglMask(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        j=0; call inorderDataSnglMask(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        j=0; call postorderDataSnglMask(rtree%p_Kchild(TRIGHT, TROOT))
      end select
    else
      ! Check if given array has correct size in its first dimension
      if (Isize(1) .ne. rtree%isizeSngl) then
        call output_line('Array too small!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
        call sys_halt()
      end if
      select case(iorder)
      case (BTREE_PREORDER)
        j=0; call preorderDataSngl(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        j=0; call inorderDataSngl(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        j=0; call postorderDataSngl(rtree%p_Kchild(TRIGHT, TROOT))
      end select
    end if

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Single data to an array

    recursive subroutine preorderDataSngl(i)
      integer, intent(in) :: i
      
      j=j+1; p_FData(:,j)=rtree%p_FData(:,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataSngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataSngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataSngl

    !**************************************************************
    ! Copy the content of the Single data to an array

    recursive subroutine inorderDataSngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataSngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FData(:,j)=rtree%p_FData(:,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataSngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataSngl

    !**************************************************************
    ! Copy the content of the Single data to an array

    recursive subroutine postorderDataSngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataSngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataSngl(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_FData(:,j)=rtree%p_FData(:,i)
    end subroutine postorderDataSngl

    !**************************************************************
    ! Copy the content of the Single data to an array (masked)

    recursive subroutine preorderDataSnglMask(i)
      integer, intent(in) :: i
      
      j=j+1; p_FData(:,j)=rtree%p_FData(mask,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataSnglMask(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataSnglMask(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataSnglMask

    !**************************************************************
    ! Copy the content of the Single data to an array (masked)

    recursive subroutine inorderDataSnglMask(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataSnglMask(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FData(:,j)=rtree%p_FData(mask,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataSnglMask(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataSnglMask

    !**************************************************************
    ! Copy the content of the Single data to an array (masked)

    recursive subroutine postorderDataSnglMask(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataSnglMask(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataSnglMask(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_FData(:,j)=rtree%p_FData(mask,i)
    end subroutine postorderDataSnglMask

  end subroutine btree_copyFromTree_arraySngl2D

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTree_arrayInt(rtree, p_Idata, mask, corder)

!<description>
    ! This subroutine copies one component of the auxiliary integer data
    ! attached to the tree to the given array. The component must be specified
    ! by the mask, e.g. mask=4 copies the fourth component of the data.
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! mask of component to be copied
    integer, intent(in) :: mask

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! integer data array
    integer, dimension(:),intent(inout) :: p_IData
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: j,iorder

    ! Check if auxiliary integer data is available
    if (rtree%isizeInt .eq. 0) then
      call output_line('No integer data available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt')
      call sys_halt()
    end if

    ! Check if given array is large enough in its second dimension
    if (size(p_IData) < rtree%NA) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt')
      call sys_halt()
    end if

    ! Check if mask is valid
    if (mask < 1 .or. mask > rtree%isizeInt) then
      call output_line('Invalid mask!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt')
      call sys_halt()
    end if
    
    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    select case(iorder)
    case(BTREE_PREORDER)
      j=0; call preorderDataInt(rtree%p_Kchild(TRIGHT, TROOT))
    case(BTREE_INORDER)
      j=0; call inorderDataInt(rtree%p_Kchild(TRIGHT, TROOT))
    case(BTREE_POSTORDER)
      j=0; call postorderDataInt(rtree%p_Kchild(TRIGHT, TROOT))
    end select

  contains

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Integer data to an array
    
    recursive subroutine preorderDataInt(i)
      integer, intent(in) :: i
      
      j=j+1; p_IData(j)=rtree%p_IData(mask,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataInt

    !**************************************************************
    ! Copy the content of the Integer data to an array
    
    recursive subroutine inorderDataInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IData(j)=rtree%p_IData(mask,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataInt

    !**************************************************************
    ! Copy the content of the Integer data to an array
    
    recursive subroutine postorderDataInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataInt(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_IData(j)=rtree%p_IData(mask,i)
    end subroutine postorderDataInt

  end subroutine btree_copyFromTree_arrayInt

  ! ***************************************************************************

!<subroutine>

  subroutine btree_copyFromTree_arrayInt2D(rtree, p_Idata, mask, corder)

!<description>
    ! This subroutine copies some (part of) the auxiliary integer data
    ! attached to the tree to the given array. The part of the data to
    ! be copied can be specified by the optional mask, e.g.
    ! mask=(/1,2,7/) only copies the first, second and seventh
    ! components of the auxiliary data.
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
    
    ! OPTIONAL: mask of components to be copied
    integer, dimension(:), intent(in), optional :: mask

    ! OPTIONAL: ordering strategy: BTREE_xxxORDER
    integer, intent(in), optional :: corder
!</input>

!<inputoutput>
    ! integer data array
    integer, dimension(:,:),intent(inout) :: p_IData
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: j,iorder

    ! Check if auxiliary integer data is available
    if (rtree%isizeInt .eq. 0) then
      call output_line('No integer data available!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
      call sys_halt()
    end if

    ! Check if given array is large enough in its second dimension
    Isize = shape(p_IData)
    if (Isize(2) < rtree%NA) then
      call output_line('Invalid mask!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
      call sys_halt()
    end if

    ! Get ordering strategy
    iorder = BTREE_INORDER
    if (present(corder)) iorder = corder

    ! Copy the content to the array
    if (present(mask)) then
      ! Check if given array has correct size in its first dimension
      if (Isize(1) .ne. size(mask)) then
        call output_line('Array dimensions do not match mask!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
        call sys_halt()
      end if

      ! Check if mask is valid
      if (any(mask < 1) .or. any(mask > rtree%isizeInt)) then
        call output_line('Invalid mask!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
        call sys_halt()
      end if
      select case(iorder)
      case (BTREE_PREORDER)
        j=0; call preorderDataIntMask(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        j=0; call inorderDataIntMask(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        j=0; call postorderDataIntMask(rtree%p_Kchild(TRIGHT, TROOT))
      end select
    else
      ! Check if given array has correct size in its first dimension
      if (Isize(1) .ne. rtree%isizeInt) then
        call output_line('Array too small!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
        call sys_halt()
      end if
      select case(iorder)
      case (BTREE_PREORDER)
        j=0; call preorderDataInt(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_INORDER)
        j=0; call inorderDataInt(rtree%p_Kchild(TRIGHT, TROOT))
      case (BTREE_POSTORDER)
        j=0; call postorderDataInt(rtree%p_Kchild(TRIGHT, TROOT))
      end select
    end if

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Integer data to an array

    recursive subroutine preorderDataInt(i)
      integer, intent(in) :: i
      
      j=j+1; p_IData(:,j)=rtree%p_IData(:,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataInt

    !**************************************************************
    ! Copy the content of the Integer data to an array

    recursive subroutine inorderDataInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IData(:,j)=rtree%p_IData(:,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataInt

    !**************************************************************
    ! Copy the content of the Integer data to an array

    recursive subroutine postorderDataInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataInt(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_IData(:,j)=rtree%p_IData(:,i)
    end subroutine postorderDataInt

    !**************************************************************
    ! Copy the content of the Integer data to an array (masked)

    recursive subroutine preorderDataIntMask(i)
      integer, intent(in) :: i
      
      j=j+1; p_IData(:,j)=rtree%p_IData(mask,i)
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDataIntMask(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDataIntMask(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDataIntMask

    !**************************************************************
    ! Copy the content of the Integer data to an array (masked)

    recursive subroutine inorderDataIntMask(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDataIntMask(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IData(:,j)=rtree%p_IData(mask,i)
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDataIntMask(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDataIntMask

    !**************************************************************
    ! Copy the content of the Integer data to an array (masked)

    recursive subroutine postorderDataIntMask(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDataIntMask(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDataIntMask(rtree%p_Kchild(TRIGHT,i))
      j=j+1; p_IData(:,j)=rtree%p_IData(mask,i)
    end subroutine postorderDataIntMask
    
  end subroutine btree_copyFromTree_arrayInt2D

  ! ***************************************************************************

!<subroutine>

  subroutine btree_insertIntoTreeDble(rtree, dkey, DData, FData, IData, iposOpt)

!<description>
    ! This subroutine inserts a new Double key into the tree and
    ! (possibly) some auxiliary data
!</description>

!<input>
    ! key
    real(DP), intent(in) :: dkey
    
    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: FData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<output>
    ! OPTIONAL: Position of the new entry. If ipos < 0 then
    ! the entry already exists and ABS(ipos) is its position
    integer, intent(out), optional :: iposOpt
!</output>
!</subroutine>
    
    ! local variables
    integer :: ipos,jpos
    
    ! Check if tree format is ok
    if (rtree%ctreeFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_insertIntoTreeDble')
      call sys_halt()
    end if

    ! Check if key is already stored in tree
    if (btree_searchInTree(rtree,dkey,jpos) .eq. BTREE_NOT_FOUND) then
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%p_Kchild(TFREE,TROOT)

      if (present(iposOpt)) iposOpt=abs(ipos)

      ! Check if memory needs to be expanded
      if (abs(ipos) > rtree%nna) &
          call btree_resizeTree(rtree,ceiling(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      if (ipos > 0) then
        rtree%p_Kchild(TFREE,TROOT) = ipos+1
      else
        ipos                      = abs(ipos)
        rtree%p_Kchild(TFREE,TROOT) = rtree%p_Kchild(TFREE,ipos)
      end if
      
      ! Store given data in memory
      rtree%p_DKey(ipos)     = dkey
      rtree%p_Kbal(ipos)     = 0
      rtree%p_Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      if ((rtree%isizeInt > 0) .and. &
          present(IData)) rtree%p_IData(:,ipos) = IData

      if ((rtree%isizeDble > 0) .and. &
          present(DData)) rtree%p_DData(:,ipos) = DData

      if ((rtree%isizeSngl > 0) .and. &
          present(FData)) rtree%p_FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%p_Kchild(merge(TLEFT,TRIGHT,jpos < 0),abs(jpos)) = ipos
      rtree%p_Kpath(rtree%depth+1) = merge(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      if (rtree%depth > 0) call rebalanceAfterInsert(rtree,rtree%depth)

    elseif(present(iposOpt)) then
      iposOpt=-rtree%p_Kchild(merge(TLEFT,TRIGHT,jpos < 0),abs(jpos))
    end if
  end subroutine btree_insertIntoTreeDble
  
  ! ***************************************************************************

!<subroutine>

  subroutine btree_insertIntoTreeSngl(rtree, skey, DData, FData, IData, iposOpt)

!<description>
    ! This subroutine inserts a new Single key into the tree and
    ! (possibly) some auxiliary data
!</description>

!<input>
    ! key
    real(SP), intent(in) :: skey
    
    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: FData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<output>
    ! OPTIONAL: Position of the new entry. If ipos < 0 then
    ! the entry already exists and ABS(ipos) is its position
    integer, intent(out), optional :: iposOpt
!</output>
!</subroutine>
    
    ! local variables
    integer :: ipos,jpos
    
    ! Check if tree format is ok
    if (rtree%ctreeFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_insertIntoTreeSngl')
      call sys_halt()
    end if

    ! Check if key is already stored in tree
    if (btree_searchInTree(rtree,skey,jpos) .eq. BTREE_NOT_FOUND) then
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%p_Kchild(TFREE,TROOT)

      if (present(iposOpt)) iposOpt=abs(ipos)

      ! Check if memory needs to be expanded
      if (abs(ipos) > rtree%nna) &
          call btree_resizeTree(rtree,ceiling(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      if (ipos > 0) then
        rtree%p_Kchild(TFREE,TROOT) = ipos+1
      else
        ipos                      = abs(ipos)
        rtree%p_Kchild(TFREE,TROOT) = rtree%p_Kchild(TFREE,ipos)
      end if
      
      ! Store given data in memory
      rtree%p_FKey(ipos)     = skey
      rtree%p_Kbal(ipos)     = 0
      rtree%p_Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      if ((rtree%isizeInt > 0) .and. &
          present(IData)) rtree%p_IData(:,ipos) = IData

      if ((rtree%isizeDble > 0) .and. &
          present(DData)) rtree%p_DData(:,ipos) = DData

      if ((rtree%isizeSngl > 0) .and. &
          present(FData)) rtree%p_FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%p_Kchild(merge(TLEFT,TRIGHT,jpos < 0),abs(jpos)) = ipos
      rtree%p_Kpath(rtree%depth+1) = merge(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      if (rtree%depth > 0) call rebalanceAfterInsert(rtree,rtree%depth)

    elseif(present(iposOpt)) then
      iposOpt=-rtree%p_Kchild(merge(TLEFT,TRIGHT,jpos < 0),abs(jpos))
    end if
  end subroutine btree_insertIntoTreeSngl

  ! ***************************************************************************

!<subroutine>

  subroutine btree_insertIntoTreeInt(rtree, ikey, DData, FData, IData, iposOpt)

!<description>
    ! This subroutine inserts a new Integer key into the tree and
    ! (possibly) some auxiliary data
!</description>

!<input>
    ! key
    integer, intent(in) :: ikey
    
    ! OPTIONAL: Double data
    real(DP), dimension(:), intent(in), optional :: DData

    ! OPTIONAL: Single data
    real(SP), dimension(:), intent(in), optional :: FData

    ! OPTIONAL: Integer data
    integer, dimension(:), intent(in), optional :: IData
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<output>
    ! OPTIONAL: Position of the new entry. If ipos < 0 then
    ! the entry already exists and ABS(ipos) is its position
    integer, intent(out), optional :: iposOpt
!</output>
!</subroutine>
    
    ! local variables
    integer :: ipos,jpos
    
    ! Check if tree format is ok
    if (rtree%ctreeFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_insertIntoTreeInt')
      call sys_halt()
    end if

    ! Check if key is already stored in tree
    if (btree_searchInTree(rtree,ikey,jpos) .eq. BTREE_NOT_FOUND) then
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%p_Kchild(TFREE,TROOT)

      if (present(iposOpt)) iposOpt=abs(ipos)
      
      ! Check if memory needs to be expanded
      if (abs(ipos) > rtree%nna) &
          call btree_resizeTree(rtree,ceiling(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      if (ipos > 0) then
        rtree%p_Kchild(TFREE,TROOT) = ipos+1
      else
        ipos                      = abs(ipos)
        rtree%p_Kchild(TFREE,TROOT) = rtree%p_Kchild(TFREE,ipos)
      end if
      
      ! Store given data in memory
      rtree%p_IKey(ipos)     = ikey
      rtree%p_Kbal(ipos)     = 0
      rtree%p_Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      if ((rtree%isizeInt > 0) .and. &
          present(IData)) rtree%p_IData(:,ipos) = IData

      if ((rtree%isizeDble > 0) .and. &
          present(DData)) rtree%p_DData(:,ipos) = DData

      if ((rtree%isizeSngl > 0) .and. &
          present(FData)) rtree%p_FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%p_Kchild(merge(TLEFT,TRIGHT,jpos < 0),abs(jpos)) = ipos
      rtree%p_Kpath(rtree%depth+1) = merge(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      if (rtree%depth > 0) call rebalanceAfterInsert(rtree,rtree%depth)
      
    elseif(present(iposOpt)) then
      iposOpt=-rtree%p_Kchild(merge(TLEFT,TRIGHT,jpos < 0),abs(jpos))
    end if
  end subroutine btree_insertIntoTreeInt

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine rebalanceAfterInsert(rtree, i)

!<description>
    ! This subroutine rebalances the AVL tree after insertion
!</description>

!<input>
    ! starting node
    integer, intent(in) :: i
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: v,x,w
    integer :: dir
    
    v   = rtree%p_Kpath(i)
    dir = merge(TLEFT,TRIGHT,v < 0)
    v   = abs(v)
    
    ! Which rotation needs to be performed?
    select case (dir)
      
    case (TLEFT)
      ! Node V has been reached from the left
      select case(rtree%p_Kbal(v))
        
      case (1) ! bal(v)=1
        rtree%p_Kbal(v) = 0
        
      case (0) ! bal(v)=0
        rtree%p_Kbal(v) = -1
        if (v .ne. rtree%p_Kchild(TRIGHT,TROOT)) &
            call rebalanceAfterInsert(rtree,i-1)
        
      case (-1) ! bal(v)=-1
        x = rtree%p_Kchild(TLEFT,v)
        w = rtree%p_Kchild(TRIGHT,x)
        
        select case(rtree%p_Kbal(x))
          
        case (-1,0) ! bal(x)=-1 or bal(x)=0
          ! Single-Right-rotation
          rtree%p_Kchild(TLEFT,v)  = rtree%p_Kchild(TRIGHT,x)
          rtree%p_Kchild(TRIGHT,x) = v
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)+1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = x
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = x
          end if
          
        case (1) ! bal(x)=1
          ! Double-Left-Right-rotation
          rtree%p_Kchild(TLEFT,v)  = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TRIGHT,x) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,w)  = x
          rtree%p_Kchild(TRIGHT,w) = v
          rtree%p_Kbal(v)          = -min(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -max(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = w
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = w
          end if
          
        end select
        
      end select
      
    case (TRIGHT)
      ! Node V has been reached from the right
      select case(rtree%p_Kbal(v))
        
      case (-1) ! bal(v)=-1
        rtree%p_Kbal(v) = 0
        
      case (0) ! bal(v)=0
        rtree%p_Kbal(v) = 1
        if (v .ne. rtree%p_Kchild(TRIGHT,TROOT))&
            call rebalanceAfterInsert(rtree,i-1)
        
      case (1) ! bal(v)=1
        x = rtree%p_Kchild(TRIGHT,v)
        w = rtree%p_Kchild(TLEFT,x)
        
        select case(rtree%p_Kbal(x))
          
        case (0,1) ! bal(x)=0 or bal(x)=1
          ! Single-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,x)
          rtree%p_Kchild(TLEFT,x)  = v
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)-1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = x
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = x
          end if
          
        case (-1) ! bal(x)=-1
          ! Double-Right-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,x)  = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TLEFT,w)  = v
          rtree%p_Kchild(TRIGHT,w) = x
          rtree%p_Kbal(v)          = -max(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -min(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = w
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = w
          end if
          
        end select
        
      end select
      
    end select
  end subroutine rebalanceAfterInsert

  ! ***************************************************************************

!<function>

  function btree_deleteFromTreeDble(rtree, dkey) result(f)

!<description>
    ! This functions deletes a Double key from the tree
!</description>

!<input>
    ! Key
    real(DP), intent(in) :: dkey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<result>
    ! Result of the deletion BTREE_NOT_FOUND / BTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos,jpos

    ! Check if tree format is ok
    if (rtree%ctreeFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_deleteFromTreeDble')
      call sys_halt()
    end if

    ! Search for key
    f=btree_searchInTree(rtree,dkey,ipred)
    if (f .eq. BTREE_NOT_FOUND) return

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))

    ! Check if node to be deleted has two children.
    if ((rtree%p_Kchild(TLEFT,ipos) .ne. TNULL) .and. &
        (rtree%p_Kchild(TRIGHT,ipos) .ne. TNULL)) then
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%p_Kchild(TLEFT,ipos)
      ipred = -ipos
      do
        rtree%depth              = rtree%depth+1
        rtree%p_Kpath(rtree%depth) = ipred
        
        if (rtree%p_Kchild(TRIGHT,jpos) .eq. TNULL) then
          ! Change key
          rtree%p_DKey(ipos) = rtree%p_DKey(jpos)

          ! Change auxiliary data
          if (rtree%isizeDble > 0) rtree%p_DData(:,ipos) = rtree%p_DData(:,jpos)
          if (rtree%isizeSngl > 0) rtree%p_FData(:,ipos) = rtree%p_FData(:,jpos)
          if (rtree%isizeInt  > 0) rtree%p_IData(:,ipos) = rtree%p_IData(:,jpos)
          exit
        end if
        ipred = jpos
        jpos = rtree%p_Kchild(TRIGHT,jpos)
      end do
    end if

    ! Node to be deleted has less than two childen
    ipos = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
    
    if (rtree%p_Kchild(TLEFT,ipos) .eq. TNULL) then
      if (rtree%p_Kchild(TRIGHT,ipos) .eq. TNULL) then
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) = TNULL
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      else
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) =&
            rtree%p_Kchild(TRIGHT,ipos)
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      end if
    else
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) =&
          rtree%p_Kchild(TLEFT,ipos)
      rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
      rtree%p_Kchild(TFREE,TROOT) = -ipos
    end if
    
    ! Invoke rebalance procedure
    if (rtree%depth > 0) call rebalanceAfterDeletion(rtree, rtree%depth)
  end function btree_deleteFromTreeDble

  ! ***************************************************************************

!<function>

  function btree_deleteFromTreeSngl(rtree, skey) result(f)

!<description>
    ! This functions deletes a Single key from the tree
!</description>

!<input>
    ! Key
    real(SP), intent(in) :: skey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<result>
    ! Result of the deletion BTREE_NOT_FOUND / BTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos,jpos

    ! Check if tree format is ok
    if (rtree%ctreeFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_deleteFromTreeSngl')
      call sys_halt()
    end if

    ! Search for key
    f=btree_searchInTree(rtree,skey,ipred)
    if (f .eq. BTREE_NOT_FOUND) return

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))

    ! Check if node to be deleted has two children.
    if ((rtree%p_Kchild(TLEFT,ipos) .ne. TNULL) .and. &
        (rtree%p_Kchild(TRIGHT,ipos) .ne. TNULL)) then
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%p_Kchild(TLEFT,ipos)
      ipred = -ipos
      do
        rtree%depth              = rtree%depth+1
        rtree%p_Kpath(rtree%depth) = ipred
        
        if (rtree%p_Kchild(TRIGHT,jpos) .eq. TNULL) then
          ! Change key
          rtree%p_FKey(ipos) = rtree%p_FKey(jpos)

          ! Change auxiliary data
          if (rtree%isizeDble > 0) rtree%p_DData(:,ipos) = rtree%p_DData(:,jpos)
          if (rtree%isizeSngl > 0) rtree%p_FData(:,ipos) = rtree%p_FData(:,jpos)
          if (rtree%isizeInt  > 0) rtree%p_IData(:,ipos) = rtree%p_IData(:,jpos)
          exit
        end if
        ipred = jpos
        jpos = rtree%p_Kchild(TRIGHT,jpos)
      end do
    end if

    ! Node to be deleted has less than two childen
    ipos = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
    
    if (rtree%p_Kchild(TLEFT,ipos) .eq. TNULL) then
      if (rtree%p_Kchild(TRIGHT,ipos) .eq. TNULL) then
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) = TNULL
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      else
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) =&
            rtree%p_Kchild(TRIGHT,ipos)
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      end if
    else
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) =&
          rtree%p_Kchild(TLEFT,ipos)
      rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
      rtree%p_Kchild(TFREE,TROOT) = -ipos
    end if
    
    ! Invoke rebalance procedure
    if (rtree%depth > 0) call rebalanceAfterDeletion(rtree, rtree%depth)
  end function btree_deleteFromTreeSngl

  ! ***************************************************************************

!<function>

  function btree_deleteFromTreeInt(rtree, ikey) result(f)

!<description>
    ! This functions deletes a Integer key from the tree
!</description>

!<input>
    ! Key
    integer, intent(in) :: ikey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<result>
    ! Result of the deletion BTREE_NOT_FOUND / BTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: ipred,ipos,jpos

    ! Check if tree format is ok
    if (rtree%ctreeFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_deleteFromTreeInt')
      call sys_halt()
    end if

    ! Search for key
    f=btree_searchInTree(rtree,ikey,ipred)
    if (f .eq. BTREE_NOT_FOUND) return

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))

    ! Check if node to be deleted has two children.
    if ((rtree%p_Kchild(TLEFT,ipos) .ne. TNULL) .and. &
        (rtree%p_Kchild(TRIGHT,ipos) .ne. TNULL)) then
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%p_Kchild(TLEFT,ipos)
      ipred = -ipos
      do
        rtree%depth              = rtree%depth+1
        rtree%p_Kpath(rtree%depth) = ipred
        
        if (rtree%p_Kchild(TRIGHT,jpos) .eq. TNULL) then
          ! Change key
          rtree%p_IKey(ipos) = rtree%p_IKey(jpos)

          ! Change auxiliary data
          if (rtree%isizeDble > 0) rtree%p_DData(:,ipos) = rtree%p_DData(:,jpos)
          if (rtree%isizeSngl > 0) rtree%p_FData(:,ipos) = rtree%p_FData(:,jpos)
          if (rtree%isizeInt  > 0) rtree%p_IData(:,ipos) = rtree%p_IData(:,jpos)
          exit
        end if
        ipred = jpos
        jpos = rtree%p_Kchild(TRIGHT,jpos)
      end do
    end if

    ! Node to be deleted has less than two childen
    ipos = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
    
    if (rtree%p_Kchild(TLEFT,ipos) .eq. TNULL) then
      if (rtree%p_Kchild(TRIGHT,ipos) .eq. TNULL) then
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) = TNULL
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      else
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) =&
            rtree%p_Kchild(TRIGHT,ipos)
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      end if
    else
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred)) =&
          rtree%p_Kchild(TLEFT,ipos)
      rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
      rtree%p_Kchild(TFREE,TROOT) = -ipos
    end if
    
    ! Invoke rebalance procedure
    if (rtree%depth > 0) call rebalanceAfterDeletion(rtree, rtree%depth)
  end function btree_deleteFromTreeInt

  ! ***************************************************************************

!<subroutine>
  
  recursive subroutine rebalanceAfterDeletion(rtree, i)

!<description>
    ! This subroutine rebalances the AVL tree after deletion
!</description>

!<input>
    ! starting node
    integer, intent(in) :: i
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: v,x,w
    integer :: dir,xbal
    
    v   = rtree%p_Kpath(i)
    dir = merge(TLEFT,TRIGHT,v < 0)
    v   = abs(v)
    
      ! Which rotation needs to be performed?
    select case (dir)
      
    case (TLEFT)
      ! Node V has been reached from the left
      select case (rtree%p_Kbal(v))
        
      case (0) ! bal(v)=0
        rtree%p_Kbal(v) = 1
        
      case (-1) ! bal(v)=-1
        rtree%p_Kbal(v) = 0
        if (v .ne. rtree%p_Kchild(TRIGHT,TROOT))&
            call rebalanceAfterDeletion(rtree,i-1)
        
      case (1) ! bal(v)=1
        x = rtree%p_Kchild(TRIGHT,v)
        
        select case(rtree%p_Kbal(x))
          
        case (-1) ! bal(x)=-1
          w = rtree%p_Kchild(TLEFT,x)
          ! Double-Right-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,x)  = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TLEFT,w)  = v
          rtree%p_Kchild(TRIGHT,w) = x
          rtree%p_Kbal(v)          = -max(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -min(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = w
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = w
            call rebalanceAfterDeletion(rtree,i-1)
          end if
          
        case DEFAULT ! bal(x)=0 or bal(x)=1
          ! Single-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,x)
          rtree%p_Kchild(TLEFT,x)  = v
          xbal                   = rtree%p_Kbal(x)
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)-1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = x
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = x
            if (xbal .eq. 1) call rebalanceAfterDeletion(rtree,i-1)
          end if
          
        end select
        
      end select
      
    case (TRIGHT)
      ! Node V has been reached from the right
      select case (rtree%p_Kbal(v))
        
      case (0) ! bal(v)=0
        rtree%p_Kbal(v) = -1
        
      case (1) ! bal(v)=1
        rtree%p_Kbal(v) = 0
        if (v .ne. rtree%p_Kchild(TRIGHT,TROOT))&
            call rebalanceAfterDeletion(rtree,i-1)
        
      case (-1) ! bal(v)=-1
        x=rtree%p_Kchild(TLEFT,v)
        
        select case(rtree%p_Kbal(x))
          
        case (1) ! bal(x)=1
          w = rtree%p_Kchild(TRIGHT,x)
          ! Double-Left-Right-rotation
          rtree%p_Kchild(TLEFT,v ) = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TRIGHT,x) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,w)  = x
          rtree%p_Kchild(TRIGHT,w) = v
          rtree%p_Kbal(v)          = -min(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -max(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = w
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = w
            call rebalanceAfterDeletion(rtree,i-1)
          end if
          
        case DEFAULT ! bal(x)=0 or bal(x)=-1
          ! Single-Right-rotation
          rtree%p_Kchild(TLEFT,v)  = rtree%p_Kchild(TRIGHT,x)
          rtree%p_Kchild(TRIGHT,x) = v
          xbal                   = rtree%p_Kbal(x)
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)+1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          if (v .eq. rtree%p_Kchild(TRIGHT,TROOT)) then
            rtree%p_Kchild(TRIGHT,TROOT) = x
          else
            rtree%p_Kchild(merge(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                abs(rtree%p_Kpath(i-1))) = x
            if (xbal .eq. -1) call rebalanceAfterDeletion(rtree,i-1)
          end if
          
        end select
        
      end select
      
    end select
  end subroutine rebalanceAfterDeletion

  ! ***************************************************************************
  
!<function>

  function btree_searchInTreeDble(rtree, dkey, ipos) result(f)
!<description>
    ! This subroutine searches for a given Double key in the tree and
    ! returns the position of its predecessor.
!</description>

!<input>
    ! Key
    real(DP), intent(in) :: dkey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<output>
    ! Position of the predecessor
    integer, intent(out) :: ipos
!</output>

!<result>
    ! Result of the search BTREE_NOT_FOUND / BTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: jpos
    integer :: dir

    ! Check if list format is ok
    if (rtree%ctreeFormat .ne. ST_DOUBLE) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_searchInTreeDble')
      call sys_halt()
    end if
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%p_Kchild(dir,ipos)
    rtree%depth = 0
    
    search: do
      if (jpos .eq. TNULL) then
        ipos = merge(-ipos,ipos,dir .eq. TLEFT)
        exit search
      end if
      
      if (rtree%p_DKey(jpos) .eq. dkey) then
        f    = BTREE_FOUND
        ipos = merge(-ipos,ipos,dir .eq. TLEFT)
        exit search
      end if
      
      ipos                     = jpos
      dir                      = merge(TLEFT,TRIGHT,rtree%p_DKey(ipos) > dkey)
      jpos                     = rtree%p_Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%p_Kpath(rtree%depth) = merge(-ipos,ipos,dir .eq. TLEFT)
    end do search
  end function btree_searchInTreeDble

  ! ***************************************************************************
  
!<function>

  function btree_searchInTreeSngl(rtree, skey, ipos) result(f)
!<description>
    ! This subroutine searches for a given Single key in the tree and
    ! returns the position of its predecessor.
!</description>

!<input>
    ! Key
    real(SP), intent(in) :: skey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<output>
    ! Position of the predecessor
    integer, intent(out) :: ipos
!</output>

!<result>
    ! Result of the search BTREE_NOT_FOUND / BTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: jpos
    integer :: dir

    ! Check if list format is ok
    if (rtree%ctreeFormat .ne. ST_SINGLE) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_searchInTreeSngl')
      call sys_halt()
    end if
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%p_Kchild(dir,ipos)
    rtree%depth = 0
    
    search: do
      if (jpos .eq. TNULL) then
        ipos = merge(-ipos,ipos,dir .eq. TLEFT)
        exit search
      end if
      
      if (rtree%p_FKey(jpos) .eq. skey) then
        f    = BTREE_FOUND
        ipos = merge(-ipos,ipos,dir .eq. TLEFT)
        exit search
      end if
      
      ipos                     = jpos
      dir                      = merge(TLEFT,TRIGHT,rtree%p_FKey(ipos) > skey)
      jpos                     = rtree%p_Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%p_Kpath(rtree%depth) = merge(-ipos,ipos,dir .eq. TLEFT)
    end do search
  end function btree_searchInTreeSngl

  ! ***************************************************************************
  
!<function>

  function btree_searchInTreeInt(rtree, ikey, ipos) result(f)
!<description>
    ! This subroutine searches for a given Integer key in the tree and
    ! returns the position of its predecessor.
!</description>

!<input>
    ! Key
    integer, intent(in) :: ikey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<output>
    ! Position of the predecessor
    integer, intent(out) :: ipos
!</output>

!<result>
    ! Result of the search BTREE_NOT_FOUND / BTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    integer :: jpos
    integer :: dir

    ! Check if list format is ok
    if (rtree%ctreeFormat .ne. ST_INT) then
      call output_line('Unsupported data format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_searchInTreeInt')
      call sys_halt()
    end if
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%p_Kchild(dir,ipos)
    rtree%depth = 0
    
    search: do
      if (jpos .eq. TNULL) then
        ipos = merge(-ipos,ipos,dir .eq. TLEFT)
        exit search
      end if
      
      if (rtree%p_IKey(jpos) .eq. ikey) then
        f    = BTREE_FOUND
        ipos = merge(-ipos,ipos,dir .eq. TLEFT)
        exit search
      end if
      
      ipos                     = jpos
      dir                      = merge(TLEFT,TRIGHT,rtree%p_IKey(ipos) > ikey)
      jpos                     = rtree%p_Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%p_Kpath(rtree%depth) = merge(-ipos,ipos,dir .eq. TLEFT)
    end do search
  end function btree_searchInTreeInt

  ! ***************************************************************************

!<function>

  function btree_getItemInTreeDble(rtree, dkey) result(ipos)

!<description>
    ! This subroutine searches for a given Double key in the tree and
    ! returns its position. If the item cannot not be found than
    ! program execution is terminated.
!</description>

!<input>
    ! Key
    real(DP), intent(in) :: dkey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<result>
    ! Position of the item
    integer :: ipos
!</result>
!</function>
  
    ! local variables
    integer :: ipred
    
    if (btree_searchInTree(rtree,dkey,ipred) .ne. BTREE_FOUND) then
      call output_line('Unable to find item in tree!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_getItemInTreeDble')
      call sys_halt()
    end if
    ipos = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
  end function btree_getItemInTreeDble

  ! ***************************************************************************

!<function>

  function btree_getItemInTreeSngl(rtree, skey) result(ipos)

!<description>
    ! This subroutine searches for a given Single key in the tree and
    ! returns its position. If the item cannot not be found than
    ! program execution is terminated.
!</description>

!<input>
    ! Key
    real(SP), intent(in) :: skey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<result>
    ! Position of the item
    integer :: ipos
!</result>
!</function>
  
    ! local variables
    integer :: ipred
    
    if (btree_searchInTree(rtree,skey,ipred) .ne. BTREE_FOUND) then
      call output_line('Unable to find item in tree!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_getItemInTreeSngl')
      call sys_halt()
    end if
    ipos = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
  end function btree_getItemInTreeSngl
  
  ! ***************************************************************************

!<function>

  function btree_getItemInTreeInt(rtree, ikey) result(ipos)

!<description>
    ! This subroutine searches for a given Integer key in the tree and
    ! returns its position. If the item cannot not be found than
    ! program execution is terminated.
!</description>

!<input>
    ! Key
    integer, intent(in) :: ikey
!</input>

!<inputoutput>
    ! binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>

!<result>
    ! Position of the item
    integer :: ipos
!</result>
!</function>
  
    ! local variables
    integer :: ipred
    
    if (btree_searchInTree(rtree,ikey,ipred) .ne. BTREE_FOUND) then
      call output_line('Unable to find item in tree!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_getItemInTreeInt')
      call sys_halt()
    end if
    ipos = rtree%p_Kchild(merge(TLEFT,TRIGHT,ipred < 0),abs(ipred))
  end function btree_getItemInTreeInt

  ! ***************************************************************************

!<subroutine>
  
  subroutine btree_printTree(rtree, cordering)

!<description>
    ! This subroutine prints the content of the tree
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree

    ! type of traversal: BTREE_xxORDER
    integer, intent(in) :: cordering
!</input>
!</subroutine>

    ! Which kind of traversal should be applied
    select case (cordering)
      
    case (BTREE_PREORDER)
      if (rtree%p_Kchild(TRIGHT,TROOT) .ne. TNULL) then
        select case(rtree%ctreeFormat)
        case (ST_DOUBLE)
          call preorderDble(rtree%p_Kchild(TRIGHT,TROOT))
        case (ST_SINGLE)
          call preorderSngl(rtree%p_Kchild(TRIGHT,TROOT))
        case (ST_INT)
          call preorderInt(rtree%p_Kchild(TRIGHT,TROOT))
        case DEFAULT
          call output_line('Unsupported data format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'btree_printTree')
          call sys_halt()
        end select
      end if
      
    case (BTREE_INORDER)
      if (rtree%p_Kchild(TRIGHT,TROOT) .ne. TNULL) then
        select case(rtree%ctreeFormat)
        case (ST_DOUBLE)
          call inorderDble(rtree%p_Kchild(TRIGHT,TROOT))
        case (ST_SINGLE)
          call inorderSngl(rtree%p_Kchild(TRIGHT,TROOT))
        case (ST_INT)
          call inorderInt(rtree%p_Kchild(TRIGHT,TROOT))
        case DEFAULT
          call output_line('Unsupported data format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'btree_printTree')
          call sys_halt()
        end select
      end if
      
    case (BTREE_POSTORDER)
      if (rtree%p_Kchild(TRIGHT,TROOT) .ne. TNULL) then
        select case(rtree%ctreeFormat)
        case (ST_DOUBLE)
          call postorderDble(rtree%p_Kchild(TRIGHT,TROOT))
        case (ST_SINGLE)
          call postorderSngl(rtree%p_Kchild(TRIGHT,TROOT))
        case (ST_INT)
          call postorderInt(rtree%p_Kchild(TRIGHT,TROOT))
        case DEFAULT
          call output_line('Unsupported data format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'btree_printTree')
          call sys_halt()
        end select
      end if
    end select
  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! Print the content of the tree in preorder for Double key
    
    recursive subroutine preorderDble(i)
      integer, intent(in) :: i
      
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_sdEL(rtree%p_DKey(i),16))//','
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderDble

    !**************************************************************
    ! Print the content of the tree in preorder for Single key
    
    recursive subroutine preorderSngl(i)
      integer, intent(in) :: i
      
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_sdEL(real(rtree%p_FKey(i),DP),16))//','
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderSngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderSngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderSngl

    !**************************************************************
    ! Print the content of the tree in preorder for Integer key
    
    recursive subroutine preorderInt(i)
      integer, intent(in) :: i
      
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_siL(rtree%p_IKey(i),16))//','
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call preorderInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call preorderInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine preorderInt

    !**************************************************************
    ! Print the content of the tree in postorder for Double key
    
    recursive subroutine postorderDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderDble(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderDble(rtree%p_Kchild(TRIGHT,i))
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_sdEL(rtree%p_DKey(i),16))//','
    end subroutine postorderDble

    !**************************************************************
    ! Print the content of the tree in postorder for Single key
    
    recursive subroutine postorderSngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderSngl(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderSngl(rtree%p_Kchild(TRIGHT,i))
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_sdEL(real(rtree%p_FKey(i),DP),16))//','
    end subroutine postorderSngl

    !**************************************************************
    ! Print the content of the tree in postorder for Integer key
    
    recursive subroutine postorderInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call postorderInt(rtree%p_Kchild(TLEFT,i))
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call postorderInt(rtree%p_Kchild(TRIGHT,i))
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_siL(rtree%p_IKey(i),16))//','
    end subroutine postorderInt

    !**************************************************************
    ! Print the content of the tree in inorder for Double key
    
    recursive subroutine inorderDble(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderDble(rtree%p_Kchild(TLEFT,i))
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_sdEL(rtree%p_DKey(i),16))//','
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderDble(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderDble

    !**************************************************************
    ! Print the content of the tree in inorder for Single key
    
    recursive subroutine inorderSngl(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderSngl(rtree%p_Kchild(TLEFT,i))
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_sdEL(real(rtree%p_FKey(i),DP),16))//','
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderSngl(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderSngl

    !**************************************************************
    ! Print the content of the tree in inorder for Integer key
    
    recursive subroutine inorderInt(i)
      integer, intent(in) :: i
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL)&
          call inorderInt(rtree%p_Kchild(TLEFT,i))
      write(*,FMT='(A)',ADVANCE='NO')&
          trim(sys_siL(rtree%p_IKey(i),16))//','
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL)&
          call inorderInt(rtree%p_Kchild(TRIGHT,i))
    end subroutine inorderInt
  end subroutine btree_printTree

  ! ***************************************************************************
  
!<function>
  
  pure function btree_getHeight(rtree) result(h)

!<description>
    ! This function computes the height of a given tree
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
!</input>

!<result>
    ! height of the tree
    integer :: h
!</result>
!</function>
    
    h = height(rtree%p_Kchild(TRIGHT,TROOT))
    
  contains
    
    pure recursive function height(i) result(h)
      integer, intent(in) :: i
      integer :: h,hl,hr
      
      if (rtree%p_Kchild(TLEFT,i) .ne. TNULL) then
        hl = height(rtree%p_Kchild(TLEFT,i))
      else
        hl = 0
      end if
      
      if (rtree%p_Kchild(TRIGHT,i) .ne. TNULL) then
        hr = height(rtree%p_Kchild(TRIGHT,i))
      else
        hr = 0
      end if
      
      h = max(hl,hr)+1
    end function height
  end function btree_getHeight

  ! ***************************************************************************

!<subroutine>
  
  subroutine btree_infoTree(rtree)

!<description>
    ! This subroutine outputs statistical info about the tree
!</description>

!<input>
    ! binary tree
    type(t_btree), intent(in) :: rtree
!</input>
!</subroutine>

    call output_line('Tree statistics:')
    call output_line('----------------')
    call output_line('ctreeFormat: '//trim(sys_siL(rtree%ctreeFormat,5)))
    call output_line('NA:          '//trim(sys_siL(rtree%NA,15)))
    call output_line('NNA:         '//trim(sys_siL(rtree%NNA,15)))
    call output_line('NNA0:        '//trim(sys_siL(rtree%NNA0,15)))
    call output_line('NRESIZE:     '//trim(sys_siL(rtree%NRESIZE,15)))
    call output_line('h_Key:       '//trim(sys_siL(rtree%h_Key,15)))
    call output_line('h_Kbal:      '//trim(sys_siL(rtree%h_Kbal,15)))
    call output_line('h_Kpath:     '//trim(sys_siL(rtree%h_Kpath,15)))
    call output_line('h_Kchild:    '//trim(sys_siL(rtree%h_Kchild,15)))
    call output_line('h_DData:     '//trim(sys_siL(rtree%h_DData,15)))
    call output_line('h_FData:     '//trim(sys_siL(rtree%h_FData,15)))
    call output_line('h_IData:     '//trim(sys_siL(rtree%h_IData,15)))
    call output_lbrk()
    call output_line('Current data  memory usage: '//&
        trim(sys_sdL(100*rtree%NA/real(rtree%NNA,DP),2))//'%')

  end subroutine btree_infoTree

  ! ***************************************************************************

!<subroutine>

  subroutine btree_duplicateTree(rtree, rtreeBackup)

!<description>
    ! This subroutine makes a copy of a binary tree in memory.
    ! It does not make sense to share some information between binary trees,
    ! so each vectors is physically copied from the source tree to the destination tree.
!</description>

!<input>
    ! Source tree
    type(t_btree), intent(in) :: rtree
!</input>

!<inputoutput>
    ! Destination tree
    type(t_btree), intent(inout) :: rtreeBackup
!</inputoutput>
!</subroutine>

    ! Release backup tree
    call btree_releaseTree(rtreeBackup)

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
    if (rtree%h_Key .ne. ST_NOHANDLE) then
      call storage_copy(rtree%h_Key,rtreeBackup%h_Key)
      select case(rtreeBackup%ctreeFormat)
      case(ST_DOUBLE)
        call storage_getbase_double(rtreeBackup%h_Key,rtreeBackup%p_DKey)
      case(ST_SINGLE)
        call storage_getbase_single(rtreeBackup%h_Key,rtreeBackup%p_FKey)
      case(ST_INT)
        call storage_getbase_int(rtreeBackup%h_Key,rtreeBackup%p_IKey)
      case DEFAULT
        call output_line('Unsupported data format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'btree_duplicateTree')
        call sys_halt()
      end select
    end if
          
    if (rtree%h_Kbal .ne. ST_NOHANDLE) then
      call storage_copy(rtree%h_Kbal,rtreeBackup%h_Kbal)
      call storage_getbase_int(rtreeBackup%h_Kbal,rtreeBackup%p_Kbal)
    end if

    if (rtree%h_Kpath .ne. ST_NOHANDLE) then
      call storage_copy(rtree%h_Kpath,rtreeBackup%h_Kpath)
      call storage_getbase_int(rtreeBackup%h_Kpath,rtreeBackup%p_Kpath)
    end if

    if (rtree%h_Kchild .ne. ST_NOHANDLE) then
      call storage_copy(rtree%h_Kchild,rtreeBackup%h_Kchild)
      call storage_getbase_int2D(rtreeBackup%h_Kchild,rtreeBackup%p_Kchild)
    end if

    if (rtree%h_DData .ne. ST_NOHANDLE) then
      call storage_copy(rtree%h_DData,rtreeBackup%h_DData)
      call storage_getbase_double2D(rtreeBackup%h_DData,rtreeBackup%p_DData)
    end if

    if (rtree%h_FData .ne. ST_NOHANDLE) then
      call storage_copy(rtree%h_FData,rtreeBackup%h_Fdata)
      call storage_getbase_single2D(rtreeBackup%h_FData,rtreeBackup%p_FData)
    end if

    if (rtree%h_IData .ne. ST_NOHANDLE) then
      call storage_copy(rtree%h_IData,rtreeBackup%h_IData)
      call storage_getbase_int2D(rtreeBackup%h_IData,rtreeBackup%p_IData)
    end if
  end subroutine btree_duplicateTree

  ! ***************************************************************************

!<subroutine>

  subroutine btree_restoreTree(rtreeBackup, rtree)

!<description>
    ! This subroutine restores a binary tree from a previous backup.
    ! The format of both trees must be the same.
!</description>

!<input>
    ! Backup of a binary tree
    type(t_btree), intent(in) :: rtreeBackup
!</input>

!<inputoutput>
    ! Destination binary tree
    type(t_btree), intent(inout) :: rtree
!</inputoutput>
!</subroutine>

    ! Check that both trees are compatible
    if (rtree%ctreeFormat .ne. rtreeBackup%ctreeFormat .or.&
        rtree%isizeDble   .ne. rtreeBackup%isizeDble   .or.&
        rtree%isizeSngl   .ne. rtreeBackup%isizeSngl   .or.&
        rtree%isizeInt    .ne. rtreeBackup%isizeInt) then
      call output_line('Incompatible binary trees!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'btree_restoreTree')
      call sys_halt()
    end if

    ! Release binary tree
    call btree_releaseTree(rtree)

    ! Duplicate the backup
    call btree_duplicateTree(rtreeBackup,rtree)
  end subroutine btree_restoreTree

  ! ***************************************************************************
  
!<function>
  
  elemental function LOG2(i)

!<description>
    ! Compute log_2(i)
!</description>

!<input>
    ! integer value
    integer, intent(in) :: i
!</input>

!<result>
    ! log_2(i)
    real(DP) :: log2
!</result>
!</function>
        
    LOG2=log(real(i,DP))/log(2._DP)
  end function LOG2
end module binarytree

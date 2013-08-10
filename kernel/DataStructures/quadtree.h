#ifndef _QUADTREE_H_
#define _QUADTREE_H_

!##############################################################################
!# ****************************************************************************
!# <name> FEAT2_PP_TEMPLATE_T(quadtree,T) </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module header file implements a (linear) quadtree. The
!# implementation is based on the description of quadtrees by
!#
!# R. Lohner, Applied CFD Techniques. An Introduction based on
!#            Finite Element Methods, Wiley, 2008
!#
!# The following routines are available:
!#
!# 1.) qtree_create
!#     -> Create a new quadtree structure
!#
!# 2.) qtree_release
!#     -> Release an existing quadtree
!#
!# 3.) qtree_copy
!#     -> Copy data to/from the quadtree
!#
!# 4.) qtree_insert
!#     -> Insert data into quadtree
!#
!# 5.) qtree_delete
!#     -> Delete data from quadtree
!#
!# 6.) qtree_find
!#     -> Search data in quadtree
!#
!# 7.) qtree_print
!#     -> Write quadtree to file
!#
!# 8.) qtree_info
!#     -> Output info about quadtree
!#
!# 9.) qtree_getDirection
!#     -> Get direction for the next node
!#
!# 10.) qtree_getSize
!#      -> Return number of vertices in quadtree
!#
!# 11.) qtree_getBoundingBox
!#      -> Return the outer bounding box
!#
!# 12.) qtree_getX
!#      -> Return the X-value at a given position
!#
!# 13.) qtree_getY
!#      -> Return the Y-value at a given position
!#
!# 14.) qtree_duplicate
!#      -> Create a duplicate / backup of a quadtree
!#
!# 15.) qtree_restore
!#      -> Restore a quadtree from a previous backup
!#
!# 16.) qtree_rebuild
!#      -> Rebuilds the structure of a quadtree
!#
!# 17.) qtree_reposition
!#      -> Reposition item in quadtree
!#
!# 18.) qtree_cast
!#      -> Casts a quadtree to a generic object
!#
!# 19.) qtree_uncast
!#      -> Casts a generic object to a quadtree
!#
!# For the internal use the following routines are available:
!#
!# 1.) resizeNVT
!#     -> Resize the number of stored values
!#
!# 2.) resizeNNODE
!#     -> Resize the number of stored nodes
!#
!# </purpose>
!##############################################################################

#include "kernel/template.h"

  implicit none

  private
  public :: FEAT2_PP_TEMPLATE_T(t_quadtree,T)
  public :: qtree_create
  public :: qtree_release
  public :: qtree_copy
  public :: qtree_insert
  public :: qtree_delete
  public :: qtree_find
  public :: qtree_print
  public :: qtree_info
  public :: qtree_getDirection
  public :: qtree_getSize
  public :: qtree_getBoundingBox
  public :: qtree_getX
  public :: qtree_getY
  public :: qtree_duplicate
  public :: qtree_restore
  public :: qtree_rebuild
  public :: qtree_reposition
  public :: qtree_cast
  public :: qtree_uncast

  interface qtree_create
    module procedure FEAT2_PP_TEMPLATE_T(qtree_create,T)
  end interface

  interface qtree_release
    module procedure FEAT2_PP_TEMPLATE_T(qtree_release,T)
  end interface

  interface qtree_copy
    module procedure FEAT2_PP_TEMPLATE_T(qtree_cpy1,T)
    module procedure FEAT2_PP_TEMPLATE_T(qtree_cpy2,T)
    module procedure FEAT2_PP_TEMPLATE_T(qtree_cpy3,T)
    module procedure FEAT2_PP_TEMPLATE_T(qtree_cpy4,T)
  end interface

  interface qtree_insert
    module procedure FEAT2_PP_TEMPLATE_T(qtree_insert,T)
  end interface

  interface qtree_delete
    module procedure FEAT2_PP_TEMPLATE_T(qtree_delete1,T)
    module procedure FEAT2_PP_TEMPLATE_T(qtree_delete2,T)
  end interface

  interface qtree_find
    module procedure FEAT2_PP_TEMPLATE_T(qtree_find,T)
  end interface

  interface qtree_print
    module procedure FEAT2_PP_TEMPLATE_T(qtree_print,T)
  end interface

  interface qtree_info
    module procedure FEAT2_PP_TEMPLATE_T(qtree_info,T)
  end interface

  interface qtree_getDirection
    module procedure FEAT2_PP_TEMPLATE_T(qtree_getDirection,T)
  end interface

  interface qtree_getSize
    module procedure FEAT2_PP_TEMPLATE_T(qtree_getSize,T)
  end interface

  interface qtree_getBoundingBox
    module procedure FEAT2_PP_TEMPLATE_T(qtree_getBoundingBox,T)
  end interface

  interface qtree_getX
    module procedure FEAT2_PP_TEMPLATE_T(qtree_getX,T)
  end interface

  interface qtree_getY
    module procedure FEAT2_PP_TEMPLATE_T(qtree_getY,T)
  end interface

  interface qtree_duplicate
    module procedure FEAT2_PP_TEMPLATE_T(qtree_duplicate,T)
  end interface

  interface qtree_restore
    module procedure FEAT2_PP_TEMPLATE_T(qtree_restore,T)
  end interface

  interface qtree_rebuild
    module procedure FEAT2_PP_TEMPLATE_T(qtree_rebuild,T)
  end interface

  interface qtree_reposition
    module procedure FEAT2_PP_TEMPLATE_T(qtree_reposition,T)
  end interface

  interface qtree_cast
    module procedure FEAT2_PP_TEMPLATE_T(qtree_cast,T)
  end interface

  interface qtree_uncast
    module procedure FEAT2_PP_TEMPLATE_T(qtree_uncast,T)
  end interface

!<types>

!<typeblock>

  ! A linear quadtree implemented as array
  type FEAT2_PP_TEMPLATE_T(t_quadtree,T)
    private

    ! Number of next free node
    integer :: NFREE = 0

    ! Number of data items stored per node
    ! Default value: QTREE_MAX
    integer :: NDATA = QTREE_MAX

    ! Number of vertices currently stored in the quadtree
    integer :: NVT = 0

    ! Total number of vertices that can be stored  in the quadtree
    integer :: NNVT = 0

    ! Number of nodes currently store in the quadtree
    integer :: NNODE = 0

    ! Total number of nodes that can be stored in the quadtree
    integer :: NNNODE = 0

    ! Total number of resize operations
    integer :: NRESIZE = 0

    ! Factor by which the quadtree is enlarged if new storage is allocate
    real(DP) :: dfactor = 1.5_DP

#ifdef T_STORAGE
    ! Handle to data vector
    integer :: h_Data = ST_NOHANDLE

    ! Handle to bounding box
    integer :: h_BdBox = ST_NOHANDLE
#endif

    ! Pointer to data vector
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_Data => null()

    ! Pointer to bounding box
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_BdBox => null()

    ! Handle to quadtree structure
    !   KNODE(QTREE_STATUS,INODE) : < 0, the node has been subdivided
    !                               = 0, the ndoe is empty
    !                               > 0, the number of points stored in the node
    !   KNODE(QTREE_PARENT,INODE) : > 0, the node the present node came from
    !                               < 0, position of the next free node which has been deleted
    !   KNODE(QTREE_PARPOS,INODE) : > 0, the position in the node the present node came from
    !   KNODE(1:NDATA,INODE)      : for KNODE(QTREE_STATUS,INODE) > 0 : the points stored in the node
    !                               for KNODE(QTREE_STATUS,INODE) < 0 : the node into which the present node was subdivided
    integer :: h_Knode = ST_NOHANDLE

    ! Pointer to quadtree structure
    integer, dimension(:,:), pointer :: p_Knode => null()

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_create,T)(rquadtree, nnvt, nnnode,&
                                        xmin, ymin, xmax, ymax, dfactor, ndata)

!<description>
    ! This subroutine creates a new quadtree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the quadtree
    integer, intent(in) :: nnvt

    ! Total number of nodes that should be stored in the quadtree
    integer, intent(in) :: nnnode

    ! Dimensions of the initial bounding box
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: xmin,ymin,xmax,ymax

    ! OPTIONAL: Factor by which the quadtree should be enlarged if
    ! new storage has to be allocated
    real(DP), intent(in), optional :: dfactor

    ! OPTIONAL: Number of data items stored per node
    integer, optional :: ndata
!</input>

!<output>
    ! Quadtree structure
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(out) :: rquadtree
!</output>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize,Ilbound,Iubound

    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1.0_DP) rquadtree%dfactor = dfactor
    end if

    ! Set number of data items
    if (present(ndata)) then
      if (ndata .gt. QTREE_MAX) rquadtree%NDATA = ndata
    end if

    ! Set values
    rquadtree%NNNODE  = nnnode
    rquadtree%NNVT    = nnvt
    rquadtree%NNODE   = 0
    rquadtree%NVT     = 0
    rquadtree%NRESIZE = 0
    rquadtree%NFREE   = 0

    ! Allocate memory and associate pointers
    Ilbound = (/QTREE_PARPOS,1/); Iubound = (/rquadtree%NDATA,nnnode/)
    call storage_new('qtree_create', 'p_Knode', Ilbound,&
                     Iubound, ST_INT, rquadtree%h_Knode, ST_NEWBLOCK_ZERO)
    call storage_getbase(rquadtree%h_Knode, rquadtree%p_Knode)

#ifdef T_STORAGE
    Isize = (/2, nnvt/)
    call storage_new('qtree_create', 'p_Data',&
                     Isize, T_STORAGE, rquadtree%h_Data, ST_NEWBLOCK_ZERO)
    call storage_getbase(rquadtree%h_Data, rquadtree%p_Data)

    Isize = (/4, nnnode/)
    call storage_new('qtree_create', 'p_BdBox',&
                     Isize, T_STORAGE, rquadtree%h_BdBox, ST_NEWBLOCK_ZERO)
    call storage_getbase(rquadtree%h_BdBox, rquadtree%p_BdBox)
#else
    allocate(rquadtree%p_Data(2,nnvt))
    allocate(rquadtree%p_BdBox(4,nnnode))
#endif

    ! Initialise first node
    rquadtree%nnode = 1
    rquadtree%p_Knode(QTREE_STATUS,1) = QTREE_EMPTY
    rquadtree%p_Knode(QTREE_PARENT,1) = QTREE_EMPTY
    rquadtree%p_Knode(QTREE_PARPOS,1) = QTREE_EMPTY
    rquadtree%p_Knode(1:rquadtree%ndata, 1) = 0

    ! Set co-ordinates of the bounding box
    rquadtree%p_BdBox(QTREE_XMIN,1) = xmin
    rquadtree%p_BdBox(QTREE_YMIN,1) = ymin
    rquadtree%p_BdBox(QTREE_XMAX,1) = xmax
    rquadtree%p_BdBox(QTREE_YMAX,1) = ymax

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_release,T)(rquadtree)

!<description>
    ! This subroutine releases an existing quadtree.
!</description>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    ! Free memory
    if (rquadtree%h_Knode .ne. ST_NOHANDLE) call storage_free(rquadtree%h_Knode)
#ifdef T_STORAGE
    if (rquadtree%h_Data .ne. ST_NOHANDLE) call storage_free(rquadtree%h_Data)
    if (rquadtree%h_BdBox .ne. ST_NOHANDLE) call storage_free(rquadtree%h_BdBox)
#else
    deallocate(rquadtree%p_BdBox, rquadtree%p_Data)
#endif
    nullify(rquadtree%p_Knode, rquadtree%p_BdBox, rquadtree%p_Data)

    ! Reset values
    rquadtree%NNNODE  = 0
    rquadtree%NNVT    = 0
    rquadtree%NNODE   = 0
    rquadtree%NVT     = 0
    rquadtree%NRESIZE = 0
    rquadtree%NFREE   = 0

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_cpy3,T)(rquadtree, h_DataDest)

!<description>
    ! This subroutine copies the content of the quadtree to a handle.
    ! If the handle is not associated, then a new handle with correct
    ! size is allocated. Otherwise, the handle is reallocated if
    ! it does not provide enough memory.
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree
!</input>

!<inputoutput>
    ! Handle to the coordinate vector
    integer, intent(inout) :: h_DataDest
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variables
    T_TYPE, dimension(:,:), pointer :: p_DataDest
    integer, dimension(2) :: Isize

    ! Check if handle is associated
    if (h_DataDest .eq. ST_NOHANDLE) then
      Isize = (/2, rquadtree%NVT/)
      call storage_new('qtree_cpy3', 'p_DataDest',&
          Isize, T_STORAGE, h_DataDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_DataDest, Isize)
      if (Isize(2) < rquadtree%NVT) then
        call storage_realloc('qtree_cpy3',&
            rquadtree%NVT, h_DataDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! Set pointers
    call storage_getbase(h_DataDest, p_DataDest)

    ! Call copy routine
    call qtree_copy(rquadtree, p_DataDest)
#else
    call output_line('Quadtree does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'qtree_cpy3')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_cpy4,T)(rquadtree, DataDest)

!<description>
    ! This subroutine copies the content of the quadtree to an array.
    ! If the array is too small, then this subroutines terminates.
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree
!</input>

!<inputoutput>
    ! Coordinate vector
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), intent(inout) :: DataDest
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize
    integer :: i

    ! Check size of array
    Isize = shape(DataDest)
    if (Isize(1) .ne. 2 .or. Isize(2) < rquadtree%NVT) then
      call output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_cpy4')
      call sys_halt()
    end if

    ! Copy data
    do i = 1, rquadtree%NVT
      DataDest(:,i) = rquadtree%p_Data(:,i)
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_cpy1,T)(h_DataSrc, rquadtree)

!<description>
    ! This subroutine copies the content of a handle to the quadtree.
!</description>

!<input>
    ! Handle to the coordinate vector
    integer, intent(in) :: h_DataSrc
!</input>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_DataSrc

    ! Set pointer
    call storage_getbase(h_DataSrc, p_DataSrc)

    ! Call copy routine
    call qtree_copy(p_DataSrc, rquadtree)
#else
    call output_line('Quadtree does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'qtree_cpy1')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_cpy2,T)(DataSrc, rquadtree)

!<description>
    ! This subroutine copies the content of an array to the quadtree.
!</description>

!<input>
    ! Coordinate vector
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), intent(in) :: DataSrc
!</input>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: i

    if (size(DataSrc,1) .ne. 2) then
      call output_line('First dimension of array must be 2!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_cpy2')
      call sys_halt()
    end if

    ! Adjust dimension of data array
    rquadtree%NVT = size(DataSrc,2)
    if (rquadtree%NNVT .lt. rquadtree%NVT) then
      call resizeNVT(rquadtree, rquadtree%NVT)
    end if

    ! Copy data array
    do i = 1, rquadtree%NVT
      rquadtree%p_Data(:,i) = DataSrc(:,i)
    end do

    ! Rebuild structure
    call qtree_rebuild(rquadtree)

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(qtree_insert,T)(rquadtree, data, ivt, inode,&
      fcb_isEqual) result(iresult)

!<description>
    ! This function inserts a new coordinate item to the quadtree. The
    ! new position IVT is returned. The optional value INODE serves
    ! as starting node in the quadtree.  If there is no space left in
    ! this node, then it is subdivided into four leaves and the
    ! insertion procedure continues recursively.  If this optional
    ! value is not present, then the starting node in the quadtree is
    ! searched for internally.
!</description>

!<input>
    ! Coordinates of the new vertex
    FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data

    ! OPTIONAL: Number of the node to which vertex should be inserted.
    ! If there is no space left, then the next free position will be used
    integer, intent(in), optional :: inode

    ! OPTIONAL: callback function to overwrite the default isEqual function
    interface
      pure logical function fcb_isEqual(data1, data2)
        use fsystem
        FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
      end function fcb_isEqual
    end interface
    optional :: fcb_isEqual
!</input>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>

!<output>
    ! OPTIONAL: Number of the inserted vertex
    integer, intent(out), optional :: ivt
!</output>

!<result>
    ! Result of the insertion:
    !   QTREE_FAILED
    !   QTREE_FOUND
    !   QTREE_INSERTED
    integer :: iresult
!</result>
!</function>

    ! local variables
    integer :: isize,jnode,jpos,jvt

    if (present(inode)) then
      jnode = inode
    else
      ! Search potential candidate for insertion
      if (present(fcb_isEqual)) then
        iresult = qtree_find(rquadtree, data, jnode, jpos, jvt, fcb_isEqual)
      else
        iresult = qtree_find(rquadtree, data, jnode, jpos, jvt, isEqual)
      end if
      if (iresult .eq. QTREE_FOUND) then
        if (present(ivt)) ivt = jvt
        return
      end if
    end if

    ! Check if there is enough space left in the nodal component of the quadtree
    if (rquadtree%NVT .eq. rquadtree%NNVT) then
      isize = max(rquadtree%NNVT+1, ceiling(rquadtree%dfactor*rquadtree%NNVT))
      call resizeNVT(rquadtree, isize)
    end if

    ! Update values
    rquadtree%NVT = rquadtree%NVT+1
    jvt = rquadtree%NVT
    rquadtree%p_Data(:,jvt) = data

    ! Insert entry recursively
    if (present(fcb_isEqual)) then
      iresult = insert(jvt, jnode, fcb_isEqual)
    else
      iresult = insert(jvt, jnode, isEqual)
    end if

    ! Check success
    if (iresult .eq. QTREE_FAILED) then
      rquadtree%NVT = rquadtree%NVT-1
      rquadtree%p_Data(:,jvt) = 0
      jvt = 0
    end if

    if (present(ivt)) ivt = jvt

  contains

    !**************************************************************
    ! Here, the recursive insertion routine follows

    recursive function insert(ivt, inode, fcb_isEqual) result(iresult)

      integer, intent(in) :: ivt,inode
      integer :: iresult

      interface
        pure logical function fcb_isEqual(data1, data2)
          use fsystem
          FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
        end function fcb_isEqual
      end interface

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE) :: xmin,ymin,xmax,ymax,xmid,ymid
      integer :: i,jvt,jnode,nnode,isize,jresult


      ! Check status of current node
      if (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. rquadtree%ndata) then

        ! Node is full

        ! Check if there are deleted nodes
        if (rquadtree%NFREE .ne. QTREE_EMPTY) then

          ! Reuse memory of the four nodes which have been deleted lately
          nnode = rquadtree%NFREE

          ! Update pointer to next free nodes
          rquadtree%NFREE = -rquadtree%p_Knode(QTREE_PARENT, nnode)

          ! Decrease starting position of first node by one
          nnode = nnode-1

        else

          ! Otherwise, create new memory if required
          if (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode) then
            isize = max(rquadtree%nnode+QTREE_MAX,&
                        ceiling(rquadtree%dfactor*rquadtree%NNNODE))
            call resizeNNODE(rquadtree, isize)
          end if

          ! New nodes are stored after all existing nodes
          nnode = rquadtree%NNODE

        end if

        ! Increase number of nodes
        rquadtree%NNODE = rquadtree%NNODE+QTREE_MAX

        ! Compute spatial coordinates of bounding boxes
        xmin = rquadtree%p_BdBox(QTREE_XMIN,inode)
        ymin = rquadtree%p_BdBox(QTREE_YMIN,inode)
        xmax = rquadtree%p_BdBox(QTREE_XMAX,inode)
        ymax = rquadtree%p_BdBox(QTREE_YMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)

        ! NW-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NW) = QTREE_NW
        rquadtree%p_BdBox(:,nnode+QTREE_NW)            = (/xmin,ymid,xmid,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NW) = 0

        ! SW-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_BdBox(:,nnode+QTREE_SW)            = (/xmin,ymin,xmid,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SW) = 0

        ! SE-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_BdBox(:,nnode+QTREE_SE)            = (/xmid,ymin,xmax,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SE) = 0

        ! NE-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_BdBox(:,nnode+QTREE_NE)            = (/xmid,ymid,xmax,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NE) = 0

        ! Add the data values from INODE to the four new nodes
        do i = 1, rquadtree%ndata
          jvt = rquadtree%p_Knode(i, inode)
          jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Data(:,jvt), inode)
          jresult = insert(jvt, jnode, fcb_isEqual)

          if (jresult .eq. QTREE_FAILED) then
            call output_line('Internal error in insertion!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'qtree_insert')
            call sys_halt()
          end if
        end do

        ! Mark the current node as subdivided and set pointers to its four children
        rquadtree%p_Knode(QTREE_STATUS, inode) =     QTREE_SUBDIV
        rquadtree%p_Knode(1:QTREE_MAX, inode)  = - (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                    nnode+QTREE_SE, nnode+QTREE_NE/)

        ! Add the new entry to the next position recursively
        jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Data(:,ivt), inode)
        iresult = insert(ivt, jnode, fcb_isEqual)

      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .ge. QTREE_EMPTY) then

        ! There is still some space in the node
        rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)+1
        rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = ivt

        ! Set success flag
        iresult = QTREE_INSERTED

      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. QTREE_SUBDIV) then

        ! Proceed to correcponding sub-tree
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, data, inode), inode)
        iresult = insert(ivt, jnode, fcb_isEqual)

      else

        ! Set failure flag
        iresult = QTREE_FAILED

      end if

    end function insert

  end function

  ! ***************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(qtree_delete1,T)(rquadtree, data, ivt,&
      fcb_isEqual) result(iresult)

!<description>
    ! This function deletes an item from the quadtree.
    ! The value IVT returns the number of the item which is
    ! moved to the position of the deleted item. If the deleted
    ! item was the last one, then IVT=NVT is returned.
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data

    ! OPTIONAL: callback function to overwrite the default isEqual function
    interface
      pure logical function fcb_isEqual(data1, data2)
        use fsystem
        FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
      end function fcb_isEqual
    end interface
    optional :: fcb_isEqual
!</input>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>

!<output>
    ! OPTIONAL: Number of the vertex that is deleted
    integer, intent(out), optional :: ivt
!</output>

!<result>
    ! Result of the deletion:
    !   QTREE_FAILED
    !   QTREE_DELETED
    integer :: iresult
!</result>
!</function>

    ! local variable
    integer :: jvt

    ! Delete item starting at root
    if (present(fcb_isEqual)) then
      iresult = delete(1, jvt, fcb_isEqual)
    else
      iresult = delete(1, jvt, isEqual)
    end if

    if (present(ivt)) ivt=jvt

  contains

    !**************************************************************
    ! Here, the recursive deletion routine follows

    recursive function delete(inode, ivt, fcb_isEqual) result(iresult)

      integer, intent(in) :: inode
      integer, intent(inout) :: ivt
      integer :: iresult

      interface
        pure logical function fcb_isEqual(data1, data2)
          use fsystem
          FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
        end function fcb_isEqual
      end interface

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE), dimension(2) :: DataTmp
      integer, dimension(QTREE_MAX) :: Knode
      integer :: i,jvt,jnode,ipos,jpos,nemptyChildren


      ! Check status of current node
      select case(rquadtree%p_Knode(QTREE_STATUS, inode))

      case (QTREE_SUBDIV)   ! Node is subdivided

        ! Compute child INODE which to look recursively.
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, data, inode), inode)
        iresult = delete(jnode, ivt, fcb_isEqual)

        ! Save values from current node
        Knode = rquadtree%p_Knode(1:QTREE_MAX, inode)

        ! Check if three or more children are empty
        nemptyChildren = 0

        do i = 1, QTREE_MAX
          jnode = -Knode(i)
          if (rquadtree%p_Knode(QTREE_STATUS, jnode) .eq. QTREE_EMPTY) then
            nemptyChildren = nemptyChildren+1
          elseif (rquadtree%p_Knode(QTREE_STATUS, jnode) .eq. QTREE_SUBDIV) then
            ! If the child is not a leaf, then do not compress this sub-tree
            return
          end if
        end do

        if (nemptyChildren .ge. QTREE_MAX-1) then

          ! Mark node as empty
          rquadtree%p_Knode(QTREE_STATUS, inode) = QTREE_EMPTY

          ! Copy data from non-empty child (if any) and mark nodes as deleted
          do i = 1, QTREE_MAX
            jnode = -Knode(i)
            if (rquadtree%p_Knode(QTREE_STATUS, jnode) .gt. QTREE_EMPTY) then
              ! Copy status of node
              rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, jnode)

              ! Copy data
              rquadtree%p_Knode(1:rquadtree%NDATA, inode) = rquadtree%p_Knode(1:rquadtree%NDATA, jnode)
            end if

            ! Mark node as deleted
            rquadtree%p_Knode(QTREE_STATUS, jnode) = QTREE_DEL

            ! Set pointer to next free position
            rquadtree%p_Knode(QTREE_PARENT, jnode) = -rquadtree%NFREE

          end do

          ! Update pointer to next free position
          rquadtree%NFREE = -Knode(1)

          ! Reduce number of nodes
          rquadtree%NNODE = rquadtree%NNODE-QTREE_MAX

        end if


      case (QTREE_EMPTY)   ! Node is empty so it cannot contain the item

        iresult = QTREE_FAILED
        ivt = 0


      case (QTREE_DEL)   ! Node is deleted -> serious error

        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_delete1')
        call sys_halt()


      case default   ! Node is a non-empty leaf

        ! Search for (x,y) in current node
        do ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)

          ! Get vertex number
          ivt = rquadtree%p_Knode(ipos, inode)

          if (fcb_isEqual(rquadtree%p_Data(:, ivt), data)) then

            ! Physically remove the item IVT from node INODE
            jpos = rquadtree%p_Knode(QTREE_STATUS, inode)

            rquadtree%p_Knode(ipos, inode) = rquadtree%p_Knode(jpos, inode)
            rquadtree%p_Knode(jpos, inode) = 0
            rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)-1

            ! If IVT is not last item then find item with largest number JVT
            if (ivt .ne. rquadtree%NVT) then
              DataTmp(:) = rquadtree%p_Data(:,rquadtree%NVT)
              if (qtree_find(rquadtree, DataTmp(:),&
                             jnode, jpos, jvt, fcb_isEqual) .eq. QTREE_FOUND) then

                ! Move last item JVT to position IVT
                rquadtree%p_Data(:, ivt) = rquadtree%p_Data(:, jvt)
                rquadtree%p_Knode(jpos, jnode) = ivt
              else
                call output_line('Internal error in deletion!',&
                                 OU_CLASS_ERROR,OU_MODE_STD,'qtree_delete')
                call sys_halt()
              end if

              ! Set number of removed vertex
              ivt = rquadtree%NVT
            end if

            ! Decrease number of vertices
            rquadtree%NVT = rquadtree%NVT-1

            ! We have found the item IVT in node INODE
            iresult = QTREE_DELETED

            ! That is it
            return
          end if
        end do

        ! We have not found the item
        iresult = QTREE_FAILED
        ivt = 0

      end select

    end function delete

  end function

  ! ***************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(qtree_delete2,T)(rquadtree, ivt, ivtReplace,&
      fcb_isEqual) result(iresult)

!<description>
    ! This function deletes vertex with number IVT from the quadtree.
!</description>

!<input>
    ! Number of the vertex to be deleted
    integer, intent(in) :: ivt

    ! OPTIONAL: callback function to overwrite the default isEqual function
    interface
      pure logical function fcb_isEqual(data1, data2)
        use fsystem
        FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
      end function fcb_isEqual
    end interface
    optional :: fcb_isEqual
!</input>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>

!<output>
    ! OPTIONAL: Number of the vertex that replaces the deleted vertex
    integer, intent(out), optional :: ivtReplace
!</output>

!<result>
    ! Result of the deletion:
    !   QTREE_FAILED
    !   QTREE_DELETED
    integer :: iresult
!</result>
!</function>

    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(2) :: data

    if (ivt .le. rquadtree%NVT) then
      ! Get coordinates and invoke deletion routine
      data    = rquadtree%p_Data(:,ivt)
      iresult = qtree_delete(rquadtree, data, ivtReplace, fcb_isEqual)
    else
      iresult = QTREE_FAILED
    end if

  end function

  ! ***************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(qtree_find,T)(rquadtree, data, inode, ipos, ivt,&
      fcb_isEqual) result(iresult)

!<description>
    ! This subroutine searches for given coordinates in the quadtree.
    ! The result of the search operation is returned by the value IRESULT.
    ! If the item was found than INODE is the number of the node, IPOS
    ! is the position of the item in the node and IVT is number of the
    ! item in the data array. Otherwise, INODE is the number of the leaf,
    ! where the item would be placed in case of insertion.
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree

    ! Coordinates that should be searched
    FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data

    ! OPTIONAL: callback function to overwrite the default isEqual function
    interface
      pure logical function fcb_isEqual(data1, data2)
        use fsystem
        FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
      end function fcb_isEqual
    end interface
    optional :: fcb_isEqual
!</input>

!<output>
    ! OPTIONAL: Number of the node in which the given coordinates are
    integer, intent(out), optional  :: inode

    ! OPTIONAL: Position of the coordinates in the node
    integer, intent(out), optional :: ipos

    ! OPTIONAL: Number of the vertex the coordinates correspond to
    integer, intent(out), optional :: ivt
!</output>

!<result>
    ! Result of the searching:
    !   QTREE_NOT_FOUND
    !   QTREE_FOUND
    integer :: iresult
!</result>
!</function>

    ! local variables
    integer :: jnode,jpos,jvt

    ! Initialise
    jnode = 1; jpos = 1; jvt = 1

    ! Search for item
    if (present(fcb_isEqual)) then
      iresult = search(jnode, jpos, jvt, fcb_isEqual)
    else
      iresult = search(jnode, jpos, jvt, isEqual)
    end if

    if (present(inode)) inode=jnode
    if (present(ipos)) ipos=jpos
    if (present(ivt)) ivt=jvt

  contains

    !**************************************************************
    ! Here, the recursive searching routine follows

    recursive function search(inode, ipos, ivt, fcb_isEqual) result(iresult)

      integer, intent(inout) :: inode,ipos,ivt
      integer :: iresult

      interface
        pure logical function fcb_isEqual(data1, data2)
          use fsystem
          FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
        end function fcb_isEqual
      end interface


      ! Check status of current node
      select case(rquadtree%p_Knode(QTREE_STATUS, inode))

      case (QTREE_SUBDIV)   ! Node is subdivided

        ! Compute child INODE which to look recursively.
        inode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, data, inode), inode)
        iresult = search(inode, ipos, ivt, fcb_isEqual)


      case (QTREE_EMPTY)   ! Node is empty so it cannot contain the item

        iresult = QTREE_NOT_FOUND


      case (QTREE_DEL)   ! Node is deleted -> serious error

        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_find')
        call sys_halt()


      case default   ! Node is a non-empty leaf

        ! Search for (x,y) in current node
        do ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)

          ! Get vertex number
          ivt = rquadtree%p_Knode(ipos, inode)

          if (fcb_isEqual(rquadtree%p_Data(:, ivt), data)) then

            ! We have found the item IVT in node INODE
            iresult = QTREE_FOUND

            ! That is it
            return
          end if
        end do

        ! We have not found the item
        iresult = QTREE_NOT_FOUND

      end select

    end function search

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(qtree_getDirection,T)(rquadtree, data, inode) result(idirection)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Data.
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree

    ! Coordinates
    FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data

    ! Number of node
    integer, intent(in) :: inode
!</input>

!<result>
    ! Further search direction
    integer :: idirection
!</result>
!</function>

    ! local variables
    FEAT2_PP_TTYPE(T_TYPE) :: xmid,ymid

    ! Compute midpoint of current node
    xmid = 0.5*(rquadtree%p_BdBox(QTREE_XMIN, inode)+&
                rquadtree%p_BdBox(QTREE_XMAX, inode))
    ymid = 0.5*(rquadtree%p_BdBox(QTREE_YMIN, inode)+&
                rquadtree%p_BdBox(QTREE_YMAX, inode))

    if (data(1) > xmid) then
      if (data(2) > ymid) then
        idirection = QTREE_NE
      else
        idirection = QTREE_SE
      end if
    else
      if (data(2) > ymid) then
        idirection = QTREE_NW
      else
        idirection = QTREE_SW
      end if
    end if

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_print,T)(rquadtree, cfilename)

!<description>
    ! This subroutine writes the content of the quadtree to a file
    ! which can be visualised by means of Matlab
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree

    ! filename of the output file
    character(LEN=*), intent(in) :: cfilename
!</input>
!</subroutine>

    ! local variables
    FEAT2_PP_TTYPE(T_TYPE) :: xmin,xmax,ymin,ymax
    integer :: iunit

    iunit = sys_getFreeUnit()
    open(UNIT=iunit, FILE=trim(adjustl(cfilename)))
    xmin = rquadtree%p_BdBox(QTREE_XMIN, 1)
    ymin = rquadtree%p_BdBox(QTREE_YMIN, 1)
    xmax = rquadtree%p_BdBox(QTREE_XMAX, 1)
    ymax = rquadtree%p_BdBox(QTREE_YMAX, 1)
    call print(xmin, ymin, xmax, ymax, 1)
    close(UNIT=iunit)

  contains

    !**************************************************************
    ! Here, the recursive print routine follows

    recursive subroutine print(xmin, ymin, xmax, ymax, inode)
      FEAT2_PP_TTYPE(T_TYPE), intent(in) :: xmin,ymin,xmax,ymax
      integer, intent(in) :: inode
      FEAT2_PP_TTYPE(T_TYPE) :: xmid,ymid
      integer :: i,j

      write(UNIT=iunit,fmt=*) 'rect'
      write(UNIT=iunit,fmt=*) xmin, ymin, xmax, ymax, inode

      if (rquadtree%p_Knode(QTREE_STATUS,inode) .eq. QTREE_SUBDIV) then

        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)

        call print(xmin, ymid, xmid, ymax, -rquadtree%p_Knode(QTREE_NW, inode))
        call print(xmin, ymin, xmid, ymid, -rquadtree%p_Knode(QTREE_SW, inode))
        call print(xmid, ymin, xmax, ymid, -rquadtree%p_Knode(QTREE_SE, inode))
        call print(xmid, ymid, xmax, ymax, -rquadtree%p_Knode(QTREE_NE, inode))

      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) > QTREE_EMPTY) then

        do i = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          write(UNIT=iunit,fmt=*) 'node'
          j = rquadtree%p_Knode(i, inode)
          write(UNIT=iunit,fmt=*) rquadtree%p_Data(1,j), rquadtree%p_Data(2,j), j
        end do

      end if

    end subroutine print

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_info,T)(rquadtree)

!<description>
    ! This subroutine outputs statistical info about the quadtree
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree
!</input>
!</subroutine>

    call output_line('Quadtree:')
    call output_line('---------')
    call output_line('NVT:     '//trim(sys_siL(rquadtree%NVT,15)))
    call output_line('NNVT:    '//trim(sys_siL(rquadtree%NNVT,15)))
    call output_line('NNODE:   '//trim(sys_siL(rquadtree%NNODE,15)))
    call output_line('NNNODE:  '//trim(sys_siL(rquadtree%NNNODE,15)))
    call output_line('NDATA:   '//trim(sys_siL(rquadtree%NDATA,15)))
    call output_line('NRESIZE: '//trim(sys_siL(rquadtree%NRESIZE,5)))
    call output_line('dfactor: '//trim(sys_sdL(rquadtree%dfactor,2)))
#ifdef T_STORAGE
    call output_line('h_Data: '//trim(sys_siL(rquadtree%h_Data,15)))
    call output_line('h_BdBox: '//trim(sys_siL(rquadtree%h_BdBox,15)))
#endif
    call output_line('h_Knode: '//trim(sys_siL(rquadtree%h_Knode,15)))
    call output_lbrk()
    call output_line('Current data memory usage: '//&
        trim(sys_sdL(100*rquadtree%NVT/real(rquadtree%NNVT,DP),2))//'%')
    call output_line('Current node memory usage: '//&
        trim(sys_sdL(100*rquadtree%NNODE/real(rquadtree%NNNODE,DP),2))//'%')

  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(qtree_getSize,T)(rquadtree) result(nvt)

!<description>
    ! This function returns the number of vertices stored in the quadtree
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree
!</input>

!<result>
    ! number of vertices in quadtree
    integer :: nvt
!</result>
!</function>

    nvt = rquadtree%NVT

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(qtree_getBoundingBox,T)(rquadtree, inode) result(Bdbox)

!<description>
    ! This function returns the bounding box of the specified node.
    ! If no node number is given, then the outer bounding box is returned.
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree

    ! OPTIONAL: number of node for which bounding box should be returned
    integer, intent(in), optional :: inode
!</input>

!<result>
    ! bounding box
    FEAT2_PP_TTYPE(T_TYPE), dimension(4) :: BdBox
!</result>
!</function>

    if (present(inode)) then
      if (inode > rquadtree%NVT) then
        call output_line('Node number exceeds quadtree dimension',&
            OU_CLASS_ERROR,OU_MODE_STD,'qtree_getBoundingBox')
        call sys_halt()
      end if
      BdBox = rquadtree%p_BdBox(QTREE_XMIN:QTREE_YMAX, inode)
    else
      BdBox = rquadtree%p_BdBox(QTREE_XMIN:QTREE_YMAX, 1)
    end if

  end function

  !************************************************************************

!<function>

  elemental function FEAT2_PP_TEMPLATE_T(qtree_getX,T)(rquadtree, ivt) result(x)

!<description>
    ! This function returns the X-value at the given position.
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree

    ! position in the quadtree
    integer, intent(in) :: ivt
!</input>

!<result>
    FEAT2_PP_TTYPE(T_TYPE) :: x
!</result>
!</function>

    x = rquadtree%p_Data(1, ivt)

  end function

  !************************************************************************

!<function>

  elemental function FEAT2_PP_TEMPLATE_T(qtree_getY,T)(rquadtree, ivt) result(y)

!<description>
    ! This function returns the Y-value at the given position.
!</description>

!<input>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree

    ! position in the quadtree
    integer, intent(in) :: ivt
!</input>

!<result>
    FEAT2_PP_TTYPE(T_TYPE) :: y
!</result>
!</function>

    y = rquadtree%p_Data(2, ivt)

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_duplicate,T)(rquadtree, rquadtreeBackup)

!<description>
    ! This subroutine makes a copy of a quadtree in memory.
    ! It does not make sense to share some information between quadtrees,
    ! so each vectors is physically copied from the source quadtree
    ! to the destination quadtree.
!</description>

!<input>
    ! Source quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtree
!</input>

!<inputoutput>
    ! Destination quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtreeBackup
!</inputoutput>
!</subroutine>

    ! Release backup quadtree
    call qtree_release(rquadtreeBackup)

    ! Copy all data
    rquadtreeBackup = rquadtree

    ! Reset handles
    rquadtreeBackup%h_Data  = ST_NOHANDLE
    rquadtreeBackup%h_BdBox = ST_NOHANDLE
    rquadtreeBackup%h_Knode = ST_NOHANDLE

    ! Copy storage blocks
    if (rquadtree%h_Knode .ne. ST_NOHANDLE) then
      call storage_copy(rquadtree%h_Knode, rquadtreeBackup%h_Knode)
      call storage_getbase(rquadtreeBackup%h_Knode,&
                           rquadtreeBackup%p_Knode)
    end if

#ifdef T_STORAGE
    if (rquadtree%h_Data .ne. ST_NOHANDLE) then
      call storage_copy(rquadtree%h_Data, rquadtreeBackup%h_Data)
      call storage_getbase(rquadtreeBackup%h_Data,&
                           rquadtreeBackup%p_Data)
    end if
#else
    if (associated(rquadtree%p_Data)) then
      allocate(rquadtreeBackup%h_Data(size(rquadtree%p_Data,1),&
                                      size(rquadtree%p_Data,2))
      rquadtreeBackup%p_Data = rquadtree%p_Data
    end if
#endif

#ifdef T_STORAGE
    if (rquadtree%h_BdBox .ne. ST_NOHANDLE) then
      call storage_copy(rquadtree%h_BdBox, rquadtreeBackup%h_BdBox)
      call storage_getbase(rquadtreeBackup%h_BdBox,&
                           rquadtreeBackup%p_BdBox)
    end if
#else
    if (associated(rquadtree%p_BdBox)) then
      allocate(rquadtreeBackup%h_Data(size(rquadtree%p_BdBox,1),&
                                      size(rquadtree%p_BdBox,2))
      rquadtreeBackup%p_BdBox = rquadtree%p_BdBox
    end if
#endif 

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_restore,T)(rquadtreeBackup, rquadtree)

!<description>
    ! This subroutine restores a quadtree from a previous backup.
!</description>

!<input>
    ! Backup of an quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in) :: rquadtreeBackup
!</input>

!<inputoutput>
    ! Destination quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    ! Release quadtree
    call qtree_release(rquadtree)

    ! Duplicate the backup
    call qtree_duplicate(rquadtreeBackup, rquadtree)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_rebuild,T)(rquadtree)

!<description>
    ! This subroutine rebuilds the structure of a quadtree
!</description>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Inodes
    integer :: h_Inodes

    ! Initialization
    rquadtree%NNODE = 1
    rquadtree%NFREE = 0

    rquadtree%NNNODE = size(rquadtree%p_Knode,2)
    rquadtree%NNVT   = size(rquadtree%p_Data,2)

    ! Check if quadtree contains data
    if (rquadtree%NVT .eq. 0) then

      ! Initialise first node
      rquadtree%nnode = 1
      rquadtree%p_Knode(QTREE_STATUS, 1) = QTREE_EMPTY
      rquadtree%p_Knode(QTREE_PARENT, 1) = QTREE_EMPTY
      rquadtree%p_Knode(QTREE_PARPOS, 1) = QTREE_EMPTY
      rquadtree%p_Knode(1:rquadtree%NDATA, 1) = 0

    else

      ! Create temporary storage
      h_Inodes = ST_NOHANDLE
      call storage_new('qtree_rebuild', 'Inodes', rquadtree%NVT,&
                       ST_INT, h_Inodes, ST_NEWBLOCK_ORDERED)
      call storage_getbase_int(h_Inodes, p_Inodes)

      ! Create quadtree top-down
      call rebuild(1, 1, rquadtree%NVT)

      ! Release temporary storage
      call storage_free(h_Inodes)

    end if

  contains

    !**************************************************************
    ! Here, the recursive top-down rebuild routine follows

    recursive subroutine rebuild(inode, istart, iend)

      integer, intent(in) :: inode, istart, iend

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE) :: xmin,ymin,xmax,ymax,xmid,ymid
      integer :: i,isize,jnode,nnode,imid1,imid2,imid3


      ! Check if istart > iend then the node is empty
      if (istart .gt. iend) return

      ! Check if current partition can be stored in a leave
      if (iend-istart+1 .le. rquadtree%NDATA) then

        ! Store vertex numbers
        do i = istart, iend
          rquadtree%p_Knode(i-istart+1,inode) = p_Inodes(i)
        end do

        ! Store number of vertices
        rquadtree%p_Knode(QTREE_STATUS, inode) = iend-istart+1

      else

        ! Otherwise, the current partition is subdivided into four subpartitions

        if (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode) then
          isize = max(rquadtree%nnode+QTREE_MAX,&
                      ceiling(rquadtree%dfactor*rquadtree%NNNODE))
          call resizeNNODE(rquadtree, isize)
        end if

        ! Store the number of node
        nnode = rquadtree%NNODE
        rquadtree%NNODE = nnode+QTREE_MAX

        ! Mark the current node as subdivided and set pointers to its four children
        rquadtree%p_Knode(QTREE_STATUS,inode) =     QTREE_SUBDIV
        rquadtree%p_Knode(1:QTREE_MAX,inode)  = - (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                    nnode+QTREE_SE, nnode+QTREE_NE/)

        ! Determine coordinates for bounding boxes
        xmin = rquadtree%p_BdBox(QTREE_XMIN,inode)
        ymin = rquadtree%p_BdBox(QTREE_YMIN,inode)
        xmax = rquadtree%p_BdBox(QTREE_XMAX,inode)
        ymax = rquadtree%p_BdBox(QTREE_YMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)

        ! NW-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NW) = QTREE_NW
        rquadtree%p_BdBox(:,nnode+QTREE_NW)            = (/xmin,ymid,xmid,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NW) = 0

        ! SW-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_BdBox(:,nnode+QTREE_SW)            = (/xmin,ymin,xmid,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SW) = 0

        ! SE-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_BdBox(:,nnode+QTREE_SE)            = (/xmid,ymin,xmax,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SE) = 0

        ! NE-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_BdBox(:,nnode+QTREE_NE)            = (/xmid,ymid,xmax,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NE) = 0


        ! Swap nodes with respect to x-axis
        imid1 = partition(istart, iend, 1, xmid)

        ! Swap nodes with respect to y-axis
        imid2 = partition(istart, imid1-1, 2, ymid)
        imid3 = partition(imid1, iend, 2, ymid)

        if (istart .lt. imid2) then
          jnode = -rquadtree%p_Knode(QTREE_SW, inode)
          call rebuild(jnode, istart, imid2-1)
        end if

        if (imid2 .lt. imid1) then
          jnode = -rquadtree%p_Knode(QTREE_NW, inode)
          call rebuild(jnode, imid2, imid1-1)
        end if

        if (imid1 .lt. imid3) then
          jnode = -rquadtree%p_Knode(QTREE_SE, inode)
          call rebuild(jnode, imid1, imid3-1)
        end if

        if (imid3 .le. iend) then
          jnode = -rquadtree%p_Knode(QTREE_NE, inode)
          call rebuild(jnode, imid3, iend)
        end if

      end if

    end subroutine rebuild

    !**************************************************************
    ! Here, the partitioning routine follows

    function partition(istart, iend, idim, dmid) result(imid)

      FEAT2_PP_TTYPE(T_TYPE), intent(in) :: dmid
      integer, intent(in) :: istart, iend, idim
      integer :: imid

      ! local variables
      integer :: i


      if (istart .gt. iend) then

        ! If istart > iend then there are no items in this partition
        imid = iend+1

      else if (istart .eq. iend) then

        ! If istart = iend then there is one item in this partition
        imid = merge(iend, iend+1, rquadtree%p_Data(idim, p_Inodes(iend)) .gt. dmid)

      else

        ! Otherwise, sort items in partition
        call quicksort(istart, iend, idim)

        ! Find location of first item belonging to "IS GREATER" partition
        do i = istart, iend
          if (rquadtree%p_Data(idim, p_Inodes(i)) .gt. dmid) then
            imid = i
            return
          end if
        end do

        ! If we end up here, then all items belong to "NOT GREATER" partition
        imid = iend+1

      end if

    end function partition

    !**************************************************************
    ! Here, the quicksort routine follows

    recursive subroutine quicksort(istart, iend, idim)

      integer, intent(in) :: istart, iend, idim

      ! local variables
      integer :: isplit

      if (istart .lt. iend) then
        isplit = split(istart, iend, idim)
        call quicksort(istart, isplit-1, idim)
        call quicksort(isplit+1, iend, idim)
      end if

    end subroutine quicksort

    !**************************************************************
    ! Here, the splitting routine follows

    function split(istart, iend, idim) result(isplit)

      integer, intent(in) :: istart, iend, idim
      integer :: isplit

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE) :: dpivot
      integer :: i, j, iaux

      ! Initialization
      i = istart
      j = iend-1
      dpivot = rquadtree%p_Data(idim, p_Inodes(iend))

      do
        do while((i .lt. iend) .and. (rquadtree%p_Data(idim, p_Inodes(i)) .le. dpivot))
          i = i+1
        end do

        do while((j .gt. istart) .and. (rquadtree%p_Data(idim, p_Inodes(j)) .ge. dpivot))
          j = j-1
        end do

        ! Swap entries if needed
        if (i .lt. j) then
          iaux        = p_Inodes(i)
          p_Inodes(i) = p_Inodes(j)
          p_Inodes(j) = iaux
        end if
        if (i .ge. j) exit
      end do

      ! Swap entry i with last entry
      iaux           = p_Inodes(i)
      p_Inodes(i)    = p_Inodes(iend)
      p_Inodes(iend) = iaux

      isplit = i

    end function split

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(qtree_reposition,T)(rquadtree, DataOld, DataNew,&
      fcb_isEqual) result(iresult)

!<description>
    ! This function modifies the coordinates of an item in the quadtree.
    ! First, the item with coordinates DDATA is searched. If the new
    ! coordinates belong to the same node (e.g. they are comprised in the
    ! Nodes bounding box), then the coordinates are simply updated.
    ! Otherwise, the vertex is deleted from the quadtree and re-inserted
    ! with the new coordinates, whereby the vertex number does not change.
!</description>

!<input>
    ! Old coordinates of the vertex
    FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: DataOld

    ! New coordinates of the vertex
    FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: DataNew

    ! OPTIONAL: callback function to overwrite the default isEqual function
    interface
      pure logical function fcb_isEqual(data1, data2)
        use fsystem
        FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
      end function fcb_isEqual
    end interface
    optional :: fcb_isEqual
!</input>

!<inputoutput>
    ! Quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>

!<result>
    ! Result of the insertion:
    !   QTREE_FAILED
    !   QTREE_REPOSITIONED
    integer :: iresult
!</result>
!</function>

    ! local variables
    integer :: ivt,jvt,jnode,jpos

    ! Move item starting at root
    if (present(fcb_isEqual)) then
      iresult = move(1, ivt, fcb_isEqual)
    else
      iresult = move(1, ivt, isEqual)
    end if

    if (iresult .eq. QTREE_DELETED) then
      ! Search potential candidate for insertion
      if (present(fcb_isEqual)) then
        iresult = qtree_find(rquadtree, DataNew, jnode, jpos, jvt, fcb_isEqual)
      else
        iresult = qtree_find(rquadtree, DataNew, jnode, jpos, jvt, isEqual)
      end if
      if (iresult .eq. QTREE_FOUND) then
        call output_line('Duplicate entry in quadtree!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_reposition')
        call sys_halt()
      end if

      ! Update values
      rquadtree%p_Data(:,ivt) = DataNew

      ! Insert entry recursively
      iresult = insert(ivt, jnode)
    end if

  contains

    !**************************************************************
    ! Here, the recursive move routine follows

    recursive function move(inode, ivt, fcb_isEqual) result(iresult)

      integer, intent(in) :: inode
      integer, intent(inout) :: ivt
      integer :: iresult

      interface
        pure logical function fcb_isEqual(data1, data2)
          use fsystem
          FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
        end function fcb_isEqual
      end interface

      ! local variables
      integer, dimension(QTREE_MAX) :: Knode
      integer :: i,jnode,ipos,jpos,nemptyChildren


      ! Check status of current node
      select case(rquadtree%p_Knode(QTREE_STATUS, inode))

      case (QTREE_SUBDIV)   ! Node is subdivided

        ! Compute child INODE which to look recursively.
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, DataOld, inode), inode)
        iresult = move(jnode, ivt, fcb_isEqual)

        ! Save values from current node
        Knode = rquadtree%p_Knode(1:QTREE_MAX, inode)

        ! Check if three or more children are empty
        nemptyChildren = 0

        do i = 1, QTREE_MAX
          jnode = -Knode(i)
          if (rquadtree%p_Knode(QTREE_STATUS, jnode) .eq. QTREE_EMPTY) then
            nemptyChildren = nemptyChildren+1
          elseif (rquadtree%p_Knode(QTREE_STATUS, jnode) .eq. QTREE_SUBDIV) then
            ! If the child is not a leaf, then do not compress this sub-tree
            return
          end if
        end do

        if (nemptyChildren .ge. QTREE_MAX-1) then

          ! Mark node as empty
          rquadtree%p_Knode(QTREE_STATUS, inode) = QTREE_EMPTY

          ! Copy data from non-empty child (if any) and mark nodes as deleted
          do i = 1, QTREE_MAX
            jnode = -Knode(i)
            if (rquadtree%p_Knode(QTREE_STATUS, jnode) .gt. QTREE_EMPTY) then
              ! Copy status of node
              rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, jnode)

              ! Copy data
              rquadtree%p_Knode(1:rquadtree%NDATA, inode) = rquadtree%p_Knode(1:rquadtree%NDATA, jnode)
            end if

            ! Mark node as deleted
            rquadtree%p_Knode(QTREE_STATUS, jnode) = QTREE_DEL

            ! Set pointer to next free position
            rquadtree%p_Knode(QTREE_PARENT, jnode) = -rquadtree%NFREE

          end do

          ! Update pointer to next free position
          rquadtree%NFREE = -Knode(1)

          ! Reduce number of nodes
          rquadtree%NNODE = rquadtree%NNODE-QTREE_MAX

        end if


      case (QTREE_EMPTY)   ! Node is empty so it cannot contain the item

        iresult = QTREE_FAILED
        ivt = 0


      case (QTREE_DEL)   ! Node is deleted -> serious error

        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_reposition')
        call sys_halt()


      case default   ! Node is a non-empty leaf

        ! Search for (x,y) in current node
        do ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)

          ! Get vertex number
          ivt = rquadtree%p_Knode(ipos, inode)

          if (fcb_isEqual(rquadtree%p_Data(:, ivt), DataOld)) then

            ! Check if the new coordinates can be stored in the same node
            if ((rquadtree%p_BdBox(QTREE_XMIN,inode) .le. DataNew(1)) .and.&
                (rquadtree%p_BdBox(QTREE_XMAX,inode) .ge. DataNew(1)) .and.&
                (rquadtree%p_BdBox(QTREE_YMIN,inode) .le. DataNew(2)) .and.&
                (rquadtree%p_BdBox(QTREE_YMAX,inode) .ge. DataNew(2))) then

              ! Just update coordinates
              rquadtree%p_Data(:,ivt) = DataNew

              ! We have updated the item IVT in node INODE
              iresult = QTREE_REPOSITIONED

              ! That is it
              return
            end if

            ! Physically remove the item IVT from node INODE
            jpos = rquadtree%p_Knode(QTREE_STATUS, inode)

            rquadtree%p_Knode(ipos, inode) = rquadtree%p_Knode(jpos, inode)
            rquadtree%p_Knode(jpos, inode) = 0
            rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)-1

            ! We have found the item IVT in node INODE
            iresult = QTREE_DELETED

            ! That is it
            return
          end if
        end do

        ! We have not found the item
        iresult = QTREE_FAILED
        ivt = 0

      end select

    end function move

    !**************************************************************
    ! Here, the recursive insertion routine follows

    recursive function insert(ivt, inode) result(iresult)

      integer, intent(in) :: ivt,inode
      integer :: iresult

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE) :: xmin,ymin,xmax,ymax,xmid,ymid
      integer :: i,jvt,jnode,nnode,isize,jresult


      ! Check status of current node
      if (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. rquadtree%ndata) then

        ! Node is full

        ! Check if there are deleted nodes
        if (rquadtree%NFREE .ne. QTREE_EMPTY) then

          ! Reuse memory of the four nodes which have been deleted lately
          nnode = rquadtree%NFREE

          ! Update pointer to next free nodes
          rquadtree%NFREE = -rquadtree%p_Knode(QTREE_PARENT, nnode)

          ! Decrease starting position of first node by one
          nnode = nnode-1

        else

          ! Otherwise, create new memory if required
          if (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode) then
            isize = max(rquadtree%nnode+QTREE_MAX,&
                        ceiling(rquadtree%dfactor*rquadtree%NNNODE))
            call resizeNNODE(rquadtree, isize)
          end if

          ! New nodes are stored after all existing nodes
          nnode = rquadtree%NNODE

        end if

        ! Increase number of nodes
        rquadtree%NNODE = rquadtree%NNODE+QTREE_MAX

        ! Compute spatial coordinates of bounding boxes
        xmin = rquadtree%p_BdBox(QTREE_XMIN,inode)
        ymin = rquadtree%p_BdBox(QTREE_YMIN,inode)
        xmax = rquadtree%p_BdBox(QTREE_XMAX,inode)
        ymax = rquadtree%p_BdBox(QTREE_YMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)

        ! NW-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NW) = QTREE_NW
        rquadtree%p_BdBox(:,nnode+QTREE_NW)            = (/xmin,ymid,xmid,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NW) = 0

        ! SW-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_BdBox(:,nnode+QTREE_SW)            = (/xmin,ymin,xmid,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SW) = 0

        ! SE-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_BdBox(:,nnode+QTREE_SE)            = (/xmid,ymin,xmax,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SE) = 0

        ! NE-node
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_BdBox(:,nnode+QTREE_NE)            = (/xmid,ymid,xmax,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NE) = 0

        ! Add the data values from INODE to the four new nodes
        do i = 1, rquadtree%ndata
          jvt = rquadtree%p_Knode(i, inode)
          jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Data(:,jvt), inode)
          jresult = insert(jvt, jnode)

          if (jresult .eq. QTREE_FAILED) then
            call output_line('Internal error in insertion!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'qtree_reposition')
            call sys_halt()
          end if
        end do

        ! Mark the current node as subdivided and set pointers to its four children
        rquadtree%p_Knode(QTREE_STATUS, inode) =     QTREE_SUBDIV
        rquadtree%p_Knode(1:QTREE_MAX, inode)  = - (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                    nnode+QTREE_SE, nnode+QTREE_NE/)

        ! Add the new entry to the next position recursively
        jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Data(:,ivt), inode)
        iresult = insert(ivt, inode)

      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .ge. QTREE_EMPTY) then

        ! There is still some space in the node
        rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)+1
        rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = ivt

        ! Set success flag
        iresult = QTREE_INSERTED

      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. QTREE_SUBDIV) then

        ! Proceed to correcponding sub-tree
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, DataNew, inode), inode)
        iresult = insert(ivt, jnode)

      else

        ! Set failure flag
        iresult = QTREE_FAILED

      end if

    end function insert

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_cast,T)(rquadtree, rgenericObject)

!<description>
    ! This subroutine casts the given quadtree to a generic object.
!</description>

!<input>
    ! The quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(in), target :: rquadtree
!</input>

!<output>
    ! The generic object
    type(t_genericObject), intent(out) :: rgenericObject
!</output>
!</subroutine>

    ! Internal data structure
    type t_void_ptr
      type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), pointer :: p_robj => null()
    end type t_void_ptr
    
    ! Internal variables
    type(t_void_ptr) :: rptr

    ! Wrap quadtree by void pointer structure
    rptr%p_robj => rquadtree
    
    ! Transfer the void pointer structure to the generic object
    rgenericObject = transfer(rptr, rgenericObject)
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(qtree_uncast,T)(rgenericObject, p_rquadtree)

!<description>
    ! This subroutine casts the given generic object into a quadtree.
!</description>

!<input>
    ! The generic object
    type(t_genericObject), intent(in) :: rgenericObject
!</input>

!<output>
    ! The quadtree
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), pointer :: p_rquadtree
!</output>
!</function>

    ! Internal data structure
    type t_void_ptr
      type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), pointer :: p_robj => null()
    end type

    ! Internal variables
    type(t_void_ptr) :: rptr

    ! Transfer the generic object to the void pointer structure
    rptr = transfer(rgenericObject, rptr)
    
    ! Unwrap quadtree from void pointer structure
    p_rquadtree => rptr%p_robj

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine resizeNVT(rquadtree, nnvt)

!<description>
    ! This subroutine reallocates memory for an existing quadtree
!</description>

!<input>
    ! New number of vertices that should be stored in the quadtree
    integer, intent(in) :: nnvt
!</input>

!<inputoutput>
    ! Quadtree that should be resized
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

#ifndef T_STORAGE
    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_Data
#endif

#ifdef T_STORAGE
    call storage_realloc('resizeNVT', nnvt, rquadtree%h_Data,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase(rquadtree%h_Data, rquadtree%p_Data)
#else
    allocate(p_Ddata(size(rquadtree%p_Data,1),size(rquadtree%p_Data,2))
    p_Ddata = rquadtree%p_Data
    deallocate(rquadtree%p_Data)
    allocate(rquadtree%p_Data(size(p_Ddata,1),nnvt))
    rquadtree%p_Data(:,1:rquadtree%NNVT) = p_Data
    deallocate(p_Data)
#endif

    rquadtree%NNVT    = nnvt
    rquadtree%NRESIZE = rquadtree%NRESIZE+1

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine resizeNNODE(rquadtree, nnnode)

!<description>
    ! This subroutine reallocates memory for an existing quadtree
!</description>

!<input>
    ! New number of nodes that should be stored in the quadtree
    integer, intent(in) :: nnnode
!</input>

!<inputoutput>
    ! Quadtree that should be resized
    type(FEAT2_PP_TEMPLATE_T(t_quadtree,T)), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

#ifndef T_STORAGE
    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_BdBox
#endif

    call storage_realloc('resizeNNODE', nnnode, rquadtree%h_BdBox,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase(rquadtree%h_BdBox, rquadtree%p_BdBox)

#ifdef T_STORAGE
    call storage_realloc('resizeNNODE', nnnode, rquadtree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase(rquadtree%h_Knode, rquadtree%p_Knode)
#else
    allocate(p_Ddata(size(rquadtree%p_BdBox,1),size(rquadtree%p_BdBox,2))
    p_BdBox = rquadtree%p_BdBox
    deallocate(rquadtree%p_BdBox)
    allocate(rquadtree%p_BdBox(size(p_BdBox,1),nnnode))
    rquadtree%p_BdBox(:,1:rquadtree%NNNODE) = p_BdBox
    deallocate(p_BdBox)
#endif

    rquadtree%NNNODE  = nnnode
    rquadtree%NRESIZE = rquadtree%NRESIZE+1

  end subroutine

  !************************************************************************

!<function>

  pure function isEqual(data1,data2) result(bisEqual)

!<description>
    ! This function checks if both data items are equal. This auxiliary
    ! function will be used in all comparisons throughout this module
    ! unless a user-defined isEqual function is provided.
!</description>

!<input>
    ! Coordinates of two vertices to be compared for equality
    FEAT2_PP_TTYPE(T_TYPE), dimension(2), intent(in) :: data1,data2
!</input>

!<result>
    ! Logical switch indicating if both vertices are equal
    logical :: bisEqual
!</result>

#ifdef T_STORAGE
    bisEqual = (maxval(abs(data1-data2)) .lt. FEAT2_PP_TEMPLATE_T(SYS_EPSREAL_,T))
#else
    bisEqual = .false.
#endif

  end function

#endif

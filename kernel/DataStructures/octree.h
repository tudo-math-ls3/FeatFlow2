#ifndef _OCTREE_H_
#define _OCTREE_H_

!##############################################################################
!# ****************************************************************************
!# <name> FEAT2_PP_TEMPLATE_T(octree,T) </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module header file implements a (linear) octree. The
!# implementation is based on the description of octrees by
!#
!# R. Lohner, Applied CFD Techniques. An Introduction based on
!#            Finite Element Methods, Wiley, 2008
!#
!# The following routines are available:
!#
!# 1.) otree_create
!#     -> Create a new octree structure
!#
!# 2.) otree_release
!#     -> Release an existing octree
!#
!# 3.) otree_copy
!#     -> Copy data to/from the octree
!#
!# 4.) otree_insert
!#     -> Insert data into octree
!#
!# 5.) otree_delete
!#     -> Delete data from octree
!#
!# 6.) otree_find
!#     -> Search data in octree
!#
!# 7.) otree_print
!#      -> Write octree to file
!#
!# 8.) otree_info
!#     -> Output info about octree
!#
!# 9.) otree_getDirection
!#     -> Get direction for the next node
!#
!# 10.) otree_getSize
!#      -> Return number of vertices in octree
!#
!# 11.) otree_getBoundingBox
!#      -> Return the outer bounding box
!#
!# 12.) otree_getX
!#      -> Return the X-value at a given position
!#
!# 13.) otree_getY
!#      -> Return the Y-value at a given position
!#
!# 14.) otree_getZ
!#      -> Return the Z-value at a given position
!#
!# 15.) otree_duplicate
!#      -> Create a duplicate / backup of an octree
!#
!# 16.) otree_restore
!#      -> Restore an octree from a previous backup
!#
!# 17.) otree_rebuild
!#      -> Rebuilds the structure of an octree
!#
!# 18.) otree_reposition
!#      -> Reposition item in octree
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

#include "../template.h"

  implicit none

  private
  public :: FEAT2_PP_TEMPLATE_T(t_octree,T)
  public :: otree_create
  public :: otree_release
  public :: otree_copy
  public :: otree_insert
  public :: otree_delete
  public :: otree_find
  public :: otree_print
  public :: otree_info
  public :: otree_getDirection
  public :: otree_getSize
  public :: otree_getBoundingBox
  public :: otree_getX
  public :: otree_getY
  public :: otree_getZ
  public :: otree_duplicate
  public :: otree_restore
  public :: otree_rebuild

  interface otree_create
    module procedure FEAT2_PP_TEMPLATE_T(otree_create,T)
  end interface

  interface otree_release
    module procedure FEAT2_PP_TEMPLATE_T(otree_release,T)
  end interface

  interface otree_copy
    module procedure FEAT2_PP_TEMPLATE_T(otree_cpy1,T)
    module procedure FEAT2_PP_TEMPLATE_T(otree_cpy2,T)
    module procedure FEAT2_PP_TEMPLATE_T(otree_cpy3,T)
    module procedure FEAT2_PP_TEMPLATE_T(otree_cpy4,T)
  end interface

  interface otree_insert
    module procedure FEAT2_PP_TEMPLATE_T(otree_insert,T)
  end interface

  interface otree_delete
    module procedure FEAT2_PP_TEMPLATE_T(otree_delete1,T)
    module procedure FEAT2_PP_TEMPLATE_T(otree_delete2,T)
  end interface

  interface otree_find
    module procedure FEAT2_PP_TEMPLATE_T(otree_find,T)
  end interface

  interface otree_print
    module procedure FEAT2_PP_TEMPLATE_T(otree_print,T)
  end interface

  interface otree_info
    module procedure FEAT2_PP_TEMPLATE_T(otree_info,T)
  end interface

  interface otree_getDirection
    module procedure FEAT2_PP_TEMPLATE_T(otree_getDirection,T)
  end interface

  interface otree_getSize
    module procedure FEAT2_PP_TEMPLATE_T(otree_getSize,T)
  end interface

  interface otree_getBoundingBox
    module procedure FEAT2_PP_TEMPLATE_T(otree_getBoundingBox,T)
  end interface

  interface otree_getX
    module procedure FEAT2_PP_TEMPLATE_T(otree_getX,T)
  end interface

  interface otree_getY
    module procedure FEAT2_PP_TEMPLATE_T(otree_getY,T)
  end interface

  interface otree_getZ
    module procedure FEAT2_PP_TEMPLATE_T(otree_getZ,T)
  end interface

  interface otree_duplicate
    module procedure FEAT2_PP_TEMPLATE_T(otree_duplicate,T)
  end interface

  interface otree_restore
    module procedure FEAT2_PP_TEMPLATE_T(otree_restore,T)
  end interface

  interface otree_rebuild
    module procedure FEAT2_PP_TEMPLATE_T(otree_rebuild,T)
  end interface

!<types>

!<typeblock>

  ! A linear octree implemented as array
  type FEAT2_PP_TEMPLATE_T(t_octree,T)
    private

    ! Number of next free node
    integer :: NFREE = 0

    ! Number of data items stored per node
    ! Default value: OTREE_MAX
    integer :: NDATA = OTREE_MAX

    ! Number of vertices currently stored in the octree
    integer :: NVT = 0

    ! Total number of vertices that can be stored  in the octree
    integer :: NNVT = 0

    ! Number of nodes currently store in the octree
    integer :: NNODE = 0

    ! Total number of nodes that can be stored in the octree
    integer :: NNNODE = 0

    ! Total number of resize operations
    integer :: NRESIZE = 0

    ! Factor by which the octree is enlarged if new storage is allocate
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

    ! Handle to octree structure
    !   KNODE(OTREE_STATUS,INODE) : < 0, the node has been subdivided
    !                               = 0, the node is empty
    !                               > 0, the number of points stored in the node
    !   KNODE(OTREE_PARENT,INODE) : > 0, the node the present node came from
    !                               < 0, position of the next free node which has been deleted
    !   KNODE(OTREE_PARPOS,INODE) : > 0, the position in the node the present node came from
    !   KNODE(1:NDATA,INODE)      : for KNODE(OTREE_STATUS,INODE) > 0 : the points stored in the node
    !                               for KNODE(OTREE_STATUS,INODE) < 0 : the nodes into which the present node was subdivided
    integer :: h_Knode = ST_NOHANDLE

    ! Pointer to octree structure
    integer, dimension(:,:), pointer :: p_Knode => null()

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_create,T)(roctree, nnvt, nnnode,&
                                        xmin, ymin, zmin, xmax, ymax, zmax,&
                                        dfactor, ndata)

!<description>
    ! This subroutine creates a new octree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the octree
    integer, intent(in) :: nnvt

    ! Total number of nodes that should be stored in the octree
    integer, intent(in) :: nnnode

    ! Dimensions of the initial bounding box
    FEAT2_PP_TTYPE(T_TYPE), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax

    ! OPTIONAL: Factor by which the octree should be enlarged if
    ! new storage has to be allocated
    real(DP), intent(in), optional :: dfactor

    ! OPTIONAL: Number of data items stored per node
    integer, optional :: ndata
!</input>

!<output>
    ! Octree structure
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(out) :: roctree
!</output>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize, Ilbound, Iubound

    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1.0_DP) roctree%dfactor=dfactor
    end if

    ! Set number of data items
    if (present(ndata)) then
      if (ndata .gt. OTREE_MAX) roctree%NDATA = ndata
    end if

    ! Set values
    roctree%NNNODE  = nnnode
    roctree%NNVT    = nnvt
    roctree%NNODE   = 0
    roctree%NVT     = 0
    roctree%NRESIZE = 0
    roctree%NFREE   = 0

    ! Allocate memory and associate pointers
    Ilbound = (/OTREE_PARPOS,1/); Iubound = (/roctree%NDATA,nnnode/)
    call storage_new('otree_create', 'p_Knode', Ilbound,&
                     Iubound, ST_INT, roctree%h_Knode, ST_NEWBLOCK_ZERO)
    call storage_getbase(roctree%h_Knode, roctree%p_Knode)

#ifdef T_STORAGE
    Isize = (/3, nnvt/)
    call storage_new('otree_create', 'p_Data',&
                     Isize, ST_DOUBLE, roctree%h_Data, ST_NEWBLOCK_ZERO)
    call storage_getbase(roctree%h_Data, roctree%p_Data)

    Isize = (/6, nnnode/)
    call storage_new('otree_create', 'p_BdBox',&
                     Isize, ST_DOUBLE, roctree%h_BdBox, ST_NEWBLOCK_ZERO)
    call storage_getbase(roctree%h_BdBox, roctree%p_BdBox)
#else
    allocate(roctree%p_Data(2,nnvt))
    allocate(roctree%p_BdBox(4,nnnode))
#endif

    ! Initialise first node
    roctree%nnode = 1
    roctree%p_Knode(OTREE_STATUS,1) = OTREE_EMPTY
    roctree%p_Knode(OTREE_PARENT,1) = OTREE_EMPTY
    roctree%p_Knode(OTREE_PARPOS,1) = OTREE_EMPTY
    roctree%p_Knode(1:roctree%ndata, 1) = 0

    ! Set co-ordinates of the bounding box
    roctree%p_BdBox(OTREE_XMIN,1) = xmin
    roctree%p_BdBox(OTREE_YMIN,1) = ymin
    roctree%p_BdBox(OTREE_ZMIN,1) = zmin
    roctree%p_BdBox(OTREE_XMAX,1) = xmax
    roctree%p_BdBox(OTREE_YMAX,1) = ymax
    roctree%p_BdBox(OTREE_ZMAX,1) = zmax

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_release,T)(roctree)

!<description>
    ! This subroutine releases an existing octree
!</description>

!<inputoutput>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    ! Free memory
    if (roctree%h_Knode .ne. ST_NOHANDLE) call storage_free(roctree%h_Knode)
#ifdef T_STORAGE
    if (roctree%h_Data .ne. ST_NOHANDLE) call storage_free(roctree%h_Data)
    if (roctree%h_BdBox .ne. ST_NOHANDLE) call storage_free(roctree%h_BdBox)
#else
    deallocate(roctree%p_BdBox, roctree%p_Data)
#endif
    nullify(roctree%p_Knode, roctree%p_BdBox, roctree%p_Data)

    ! Reset values
    roctree%NNNODE  = 0
    roctree%NNVT    = 0
    roctree%NNODE   = 0
    roctree%NVT     = 0
    roctree%NRESIZE = 0
    roctree%NFREE   = 0

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_cpy3,T)(roctree, h_DataDest)

!<description>
    ! This subroutine copies the content of the octree to a handle.
    ! If the handle is not associated, then a new handle with correct
    ! size is allocated. Otherwise, the handle is reallocated if
    ! it does not provide enough memory.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree
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
      Isize = (/3, roctree%NVT/)
      call storage_new('otree_cpy3','p_DataDest',&
          Isize, T_STORAGE, h_DataDest, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_DataDest, Isize)
      if (Isize(2) < roctree%NVT) then
        call storage_realloc('otree_cpy3',&
            roctree%NVT, h_DataDest, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! Set pointers
    call storage_getbase(h_DataDest, p_DataDest)

    ! Call copy routine
    call otree_copy(roctree, p_DataDest)
#else
    call output_line('Octree does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'otree_cpy3')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_cpy4,T)(roctree, DataDest)

!<description>
    ! This subroutine copies the content of the octree to an array.
    ! If the array is too small, then this subroutines terminates.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree
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
    if (Isize(1) .ne. 3 .or. Isize(2) < roctree%NVT) then
      call output_line('Array too small!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_cpy4')
      call sys_halt()
    end if

    ! Copy data
    do i = 1, roctree%NVT
      DataDest(:,i) = roctree%p_Data(:,i)
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_cpy1,T)(h_DataSrc, roctree)

!<description>
    ! This subroutine copies the content of a handle to the octree.
!</description>

!<input>
    ! Handle to the coordinate vector
    integer, intent(in) :: h_DataSrc
!</input>

!<inputoutput>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

#ifdef T_STORAGE
    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_DataSrc

    ! Set pointer
    call storage_getbase(h_DataSrc, p_DataSrc)

    ! Call copy routine
    call otree_copy(p_DataSrc, roctree)
#else
    call output_line('OCtree does not support storage handles!',&
        OU_CLASS_ERROR,OU_MODE_STD,'otree_cpy1')
    call sys_halt()
#endif

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_cpy2,T)(DataSrc, roctree)

!<description>
    ! This subroutine copies the content of an array to the octree.
!</description>

!<input>
    ! Coordinate vector
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), intent(in) :: DataSrc
!</input>

!<inputoutput>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: i

    if (size(DataSrc,1) .ne. 3) then
      call output_line('First dimension of array must be 3!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_cpy2')
      call sys_halt()
    end if

    ! Adjust dimension of data array
    roctree%NVT = size(DataSrc,2)
    if (roctree%NNVT .lt. roctree%NVT) then
      call resizeNVT(roctree, roctree%NVT)
    end if

    ! Copy data array
    do i = 1, roctree%NVT
      roctree%p_Data(:,i) = DataSrc(:,i)
    end do

    ! Rebuild structure
    call otree_rebuild(roctree)

  end subroutine

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(otree_insert,T)(roctree, data, ivt, inode) result(iresult)

!<description>
    ! This function inserts a new coordinate item to the octree. The
    ! new position IVT is returned. The optional value INODE serves
    ! as starting node in the octree.  If there is no space left in
    ! this node, then it is subdivided into eight leaves and the
    ! insertion procedure continues recursively. If this optional
    ! value is not present, then the starting node in the octree is
    ! searched for internally.
!</description>

!<input>
    ! Coordinates of the new vertex
    FEAT2_PP_TTYPE(T_TYPE), dimension(3), intent(in) :: data

    ! OPTIONAL: Number of the node to which vertex should be inserted.
    ! If there is no space left, then the next free position will be used
    integer, intent(in), optional :: inode
!</input>

!<inputoutput>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>

!<output>
    ! Number of the inserted vertex
    integer, intent(out) :: ivt
!</output>

!<result>
    ! Result of the insertion:
    !   OTREE_FAILED
    !   OTREE_FOUND
    !   OTREE_INSERTED
    integer :: iresult
!</result>
!</function>

    ! local variables
    integer :: isize,jnode,jpos,jvt

    if (present(inode)) then
      jnode = inode
    else
      ! Search potential candidate for insertion
      iresult = otree_find(roctree, data, jnode, jpos, jvt)
      if (iresult .eq. OTREE_FOUND) then
        ivt = jvt
        return
      end if
    end if

    ! Check if there is enough space left in the nodal component of the octree
    if (roctree%NVT .eq. roctree%NNVT) then
      isize = max(roctree%NNVT+1, ceiling(roctree%dfactor*roctree%NNVT))
      call resizeNVT(roctree, isize)
    end if

    ! Update values
    roctree%NVT = roctree%NVT+1
    ivt = roctree%NVT
    roctree%p_Data(:,ivt) = data

    ! Insert entry recursively
    iresult = insert(ivt, jnode)

    ! Check success
    if (iresult .eq. OTREE_FAILED) then
      roctree%NVT = roctree%NVT-1
      roctree%p_Data(:,ivt) = 0.0_DP
      ivt = 0
    end if

  contains

    !**************************************************************
    ! Here, the recursive insertion routine follows

    recursive function insert(ivt, inode) result(iresult)

      integer, intent(in) :: ivt,inode
      integer :: iresult

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE) :: xmin,ymin,zmin,xmax,ymax,zmax,xmid,ymid,zmid
      integer :: i,jvt,jnode,nnode,isize,jresult


      ! Check status of current node
      if (roctree%p_Knode(OTREE_STATUS, inode) .eq. roctree%ndata) then

        ! Node is full

        ! Check if there are deleted nodes
        if (roctree%NFREE .ne. OTREE_EMPTY) then

          ! Reuse memory of the eight nodes which have been deleted lately
          nnode = roctree%NFREE

          ! Update pointer to next free nodes
          roctree%NFREE = -roctree%p_Knode(OTREE_PARENT, nnode)

          ! Decrease starting position of first node by one
          nnode = nnode-1

        else

          ! Otherwise, create new memory if required
          if (roctree%nnode+OTREE_MAX > roctree%nnnode) then
            isize = max(roctree%nnode+OTREE_MAX,&
                        ceiling(roctree%dfactor*roctree%NNNODE))
            call resizeNNODE(roctree, isize)
          end if

          ! New nodes are stored after all existing nodes
          nnode = roctree%NNODE

        end if

        ! Increase number of nodes
        roctree%NNODE = roctree%NNODE+OTREE_MAX

        ! Compute spatial coordinates of bounding boxes
        xmin = roctree%p_BdBox(OTREE_XMIN,inode)
        ymin = roctree%p_BdBox(OTREE_YMIN,inode)
        zmin = roctree%p_BdBox(OTREE_ZMIN,inode)
        xmax = roctree%p_BdBox(OTREE_XMAX,inode)
        ymax = roctree%p_BdBox(OTREE_YMAX,inode)
        zmax = roctree%p_BdBox(OTREE_ZMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)

        ! NWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWF) = OTREE_NWF
        roctree%p_BdBox(:,nnode+OTREE_NWF)            = (/xmin,ymid,zmin,xmid,ymax,zmid/)

        ! SWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWF) = OTREE_SWF
        roctree%p_BdBox(:,nnode+OTREE_SWF)            = (/xmin,ymin,zmin,xmid,ymid,zmid/)

        ! SEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEF) = OTREE_SEF
        roctree%p_BdBox(:,nnode+OTREE_SEF)            = (/xmid,ymin,zmin,xmax,ymid,zmid/)

        ! NEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEF) = OTREE_NEF
        roctree%p_BdBox(:,nnode+OTREE_NEF)            = (/xmid,ymid,zmin,xmax,ymax,zmid/)

        ! NWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWB) = OTREE_NWB
        roctree%p_BdBox(:,nnode+OTREE_NWB)            = (/xmin,ymid,zmid,xmid,ymax,zmax/)

        ! SWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWB) = OTREE_SWB
        roctree%p_BdBox(:,nnode+OTREE_SWB)            = (/xmin,ymin,zmid,xmid,ymid,zmax/)

        ! SEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEB) = OTREE_SEB
        roctree%p_BdBox(:,nnode+OTREE_SEB)            = (/xmid,ymin,zmid,xmax,ymid,zmax/)

        ! NEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEB) = OTREE_NEB
        roctree%p_BdBox(:,nnode+OTREE_NEB)            = (/xmid,ymid,zmid,xmax,ymax,zmax/)

        ! Add the data values from INODE to the eight new nodes
        do i = 1, roctree%ndata
          jvt = roctree%p_Knode(i, inode)
          jnode = nnode+otree_getDirection(roctree, roctree%p_Data(:,jvt), inode)
          jresult = insert(jvt, jnode)

          if (jresult .eq. OTREE_FAILED) then
            call output_line('Internal error in insertion!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'otree_insert')
            call sys_halt()
          end if
        end do

        ! Mark the current nodes as subdivided and set pointers to its eight children
        roctree%p_Knode(OTREE_STATUS,inode) =     OTREE_SUBDIV
        roctree%p_Knode(1:OTREE_MAX, inode) = - (/nnode+OTREE_NWF, nnode+OTREE_SWF,&
                                                  nnode+OTREE_SEF, nnode+OTREE_NEF,&
                                                  nnode+OTREE_NWB, nnode+OTREE_SWB,&
                                                  nnode+OTREE_SEB, nnode+OTREE_NEB/)

        ! Add the new entry to the next position recursively
        jnode = nnode+otree_getDirection(roctree, roctree%p_Data(:,ivt),inode)
        iresult = insert(ivt, jnode)

      elseif (roctree%p_Knode(OTREE_STATUS, inode) .ge. OTREE_EMPTY) then

        ! There is still some space in the node
        roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)+1
        roctree%p_Knode(roctree%p_Knode(OTREE_STATUS, inode), inode) = ivt

        ! Set success flag
        iresult = OTREE_INSERTED

      elseif (roctree%p_Knode(OTREE_STATUS, inode) .eq. OTREE_SUBDIV) then

        ! Proceed to correcponding sub-tree
        jnode = -roctree%p_Knode(otree_getDirection(roctree, data, inode), inode)
        iresult = insert(ivt, jnode)

      else

        ! Set failure flag
        iresult = OTREE_FAILED

      end if

    end function insert

  end function

  ! ***************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(otree_delete1,T)(roctree, data, ivt) result(iresult)

!<description>
    ! This function deletes an item from the octree.
    ! The value IVT returns the number of the item which is
    ! moved to the position of the deleted item. If the deleted
    ! item was the last one, then IVT=NVT is returned.
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    FEAT2_PP_TTYPE(T_TYPE), dimension(3), intent(in) :: data
!</input>

!<inputoutput>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that is deleted
    integer, intent(out) :: ivt
!</output>

!<result>
    ! Result of the deletion:
    !   OTREE_FAILED
    !   OTREE_DELETED
    integer :: iresult
!</result>
!</function>

    ! Delete item startin at root
    iresult = delete(1, ivt)

  contains

    !**************************************************************
    ! Here, the recursive deletion routine follows

    recursive function delete(inode, ivt) result(iresult)

      integer, intent(in) :: inode
      integer, intent(out) :: ivt
      integer :: iresult

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE), dimension(3) :: DataTmp
      integer, dimension(OTREE_MAX) :: Knode
      integer :: jvt,i,jnode,ipos,jpos,nemptyChildren


      ! Check status of current node
      select case(roctree%p_Knode(OTREE_STATUS, inode))

      case (OTREE_SUBDIV)   ! Node is subdivided

        ! Compute child INODE which to look recursively.
        jnode = -roctree%p_Knode(otree_getDirection(roctree, data, inode), inode)
        iresult = delete(jnode, ivt)

        ! Save values from current node
        Knode = roctree%p_Knode(1:OTREE_MAX, inode)

        ! Check if three or more children are empty
        nemptyChildren = 0

        do i = 1, OTREE_MAX
          jnode = -Knode(i)
          if (roctree%p_Knode(OTREE_STATUS, jnode) .eq. OTREE_EMPTY) then
            nemptyChildren = nemptyChildren+1
          elseif (roctree%p_Knode(OTREE_STATUS, jnode) .eq. OTREE_SUBDIV) then
            ! If the child is not a leaf, then do not compress this sub-tree
            return
          end if
        end do

        if (nemptyChildren .ge. OTREE_MAX-1) then

          ! Mark node as empty
          roctree%p_Knode(OTREE_STATUS, inode) = OTREE_EMPTY

          ! Copy data from non-empty child (if any) and mark nodes as deleted
          do i = 1, OTREE_MAX
            jnode = -Knode(i)
            if (roctree%p_Knode(OTREE_STATUS, jnode) .gt. OTREE_EMPTY) then
              ! Copy status of node
              roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, jnode)

              ! Copy data
              roctree%p_Knode(1:roctree%NDATA, inode) = roctree%p_Knode(1:roctree%NDATA, jnode)
            end if

            ! Mark node as deleted
            roctree%p_Knode(OTREE_STATUS, jnode) = OTREE_DEL

            ! Set pointer to next free position
            roctree%p_Knode(OTREE_PARENT, jnode) = -roctree%NFREE

          end do

          ! Update pointer to next free position
          roctree%NFREE = -Knode(1)

          ! Reduce number of nodes
          roctree%NNODE = roctree%NNODE-OTREE_MAX

        end if


      case (OTREE_EMPTY)   ! Node is empty so it cannot contain the item

        iresult = OTREE_FAILED
        ivt = 0


      case (OTREE_DEL)   ! Node is deleted -> serious error

        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'otree_delete1')
        call sys_halt()


      case default   ! Node is a non-empty leaf

        ! Search for (x,y,z) in current node
        do ipos = 1, roctree%p_Knode(OTREE_STATUS, inode)

          ! Get vertex number
          ivt = roctree%p_Knode(ipos, inode)

          if (maxval(abs(roctree%p_Data(:, ivt)-data)) .le. SYS_EPSREAL_DP) then

            ! Physically remove the item IVT from node INODE
            jpos = roctree%p_Knode(OTREE_STATUS, inode)

            roctree%p_Knode(ipos, inode) = roctree%p_Knode(jpos, inode)
            roctree%p_Knode(jpos, inode) = 0
            roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)-1

            ! If IVT is not last item then find item with largest number JVT
            if (ivt .ne. roctree%NVT) then
              DataTmp(:) = roctree%p_Data(:,roctree%NVT)
              if (otree_find(roctree, DataTmp(:),&
                             jnode, jpos, jvt) .eq. OTREE_FOUND) then

                ! Move last item JVT to position IVT
                roctree%p_Data(:, ivt) = roctree%p_Data(:, jvt)
                roctree%p_Knode(jpos, jnode) = ivt
              else
                call output_line('Internal error in deletion!',&
                                 OU_CLASS_ERROR,OU_MODE_STD,'otree_delete')
                call sys_halt()
              end if

              ! Set number of removed vertex
              ivt = roctree%NVT
            end if

            ! Decrease number of vertices
            roctree%NVT = roctree%NVT-1

            ! We have found the item IVT in node INODE
            iresult = OTREE_DELETED

            ! That is it
            return
          end if
        end do

        ! We have not found the item
        iresult = OTREE_FAILED
        ivt = 0

      end select

    end function delete

  end function

  ! ***************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(otree_delete2,T)(roctree, ivt, ivtReplace) result(iresult)

!<description>
    ! This function deletes vertex with number IVT from the octree.
!</description>

!<input>
    ! Number of the vertex to be deleted
    integer, intent(in) :: ivt
!</input>

!<inputoutput>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that replaces the deleted vertex
    integer, intent(out) :: ivtReplace
!</output>

!<result>
    ! Result of the deletion:
    !   OTREE_FAILED
    !   OTREE_DELETED
    integer :: iresult
!</result>
!</function>

    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(3) :: data

    if (ivt .le. roctree%NVT) then
      ! Get coordinates and invoke deletion routine
      data   = roctree%p_Data(:,ivt)
      iresult = otree_delete(roctree, data, ivtReplace)
    else
      iresult = OTREE_FAILED
    end if

  end function

  ! ***************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(otree_find,T)(roctree, data, inode, ipos, ivt) result(iresult)

!<description>
    ! This subroutine searches for given coordinates in the octree.
    ! The result of the search operation is returned by the value IRESULT.
    ! If the item was found than INODE is the number of the node, IPOS
    ! is the position of the item in the node and IVT is number of the
    ! item in the data array. Otherwise, INODE is the number of the leaf,
    ! where the item would be placed in case of insertion.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree

    ! Coordinates that should be searched for
    FEAT2_PP_TTYPE(T_TYPE), dimension(3), intent(in) :: data
!</input>

!<output>
    ! Number of the node in which the given coordinates are
    integer, intent(out) :: inode

    ! Position of the coordinates in the node
    integer, intent(out) :: ipos

    ! Number of the vertex the coordinates correspond to
    integer, intent(out) :: ivt
!</output>

!<result>
    ! Result of the searching:
    !   OTREE_NOT_FOUND
    !   OTREE_FOUND
    integer :: iresult
!</result>
!</function>

    ! Initialise
    inode = 1; ipos = 1; ivt = 1

    ! Search for item
    iresult = search(inode, ipos, ivt)

  contains

    !**************************************************************
    ! Here, the recursive searching routine follows

    recursive function search(inode, ipos, ivt) result(iresult)

      integer, intent(inout) :: inode,ipos,ivt
      integer :: iresult


      ! Check status of current node
      select case(roctree%p_Knode(OTREE_STATUS, inode))

      case (OTREE_SUBDIV)   ! Node is subdivided

        ! Compute child INODE which to look recursively.
        inode = -roctree%p_Knode(otree_getDirection(roctree, data, inode), inode)
        iresult =  search(inode, ipos, ivt)


      case (OTREE_EMPTY)   ! Node is empty so it cannot contain the item

        iresult = OTREE_NOT_FOUND


      case (OTREE_DEL)   ! Node is deleted -> serious error

        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'otree_find')
        call sys_halt()


      case default   ! Node is a non-empty leaf

        ! Search for (x,y,z) in current node
        do ipos = 1, roctree%p_Knode(OTREE_STATUS, inode)

          ! Get vertex number
          ivt = roctree%p_Knode(ipos, inode)

          if (maxval(abs(roctree%p_Data(:, ivt)-data)) .le. SYS_EPSREAL_DP) then

            ! We have found the item IVT in node INODE
            iresult = OTREE_FOUND

            ! That is it
            return
          end if
        end do

        ! We have not found the item
        iresult = OTREE_NOT_FOUND

      end select

    end function search

  end function

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(otree_getDirection,T)(roctree, data, inode) result(idirection)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Data.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree

    ! Coordinates
    FEAT2_PP_TTYPE(T_TYPE), dimension(3), intent(in) :: data

    ! Number of node
    integer, intent(in) :: inode
!</input>

!<result>
    ! Further search direction
    integer :: idirection
!</result>
!</function>

    ! local variables
    FEAT2_PP_TTYPE(T_TYPE) :: xmid,ymid,zmid

    ! Compute midpoint of current node
    xmid=0.5*(roctree%p_BdBox(OTREE_XMIN, inode)+&
              roctree%p_BdBox(OTREE_XMAX, inode))
    ymid=0.5*(roctree%p_BdBox(OTREE_YMIN, inode)+&
              roctree%p_BdBox(OTREE_YMAX, inode))
    zmid=0.5*(roctree%p_BdBox(OTREE_ZMIN, inode)+&
              roctree%p_BdBox(OTREE_ZMAX, inode))

    ! Do we have to look in the front or back?
    if (data(3) < zmid) then

      if (data(1) > xmid) then
        if (data(2) > ymid) then
          idirection = OTREE_NEF
        else
          idirection = OTREE_SEF
        end if
      else
        if (data(2) > ymid) then
          idirection = OTREE_NWF
        else
          idirection = OTREE_SWF
        end if
      end if

    else

      if (data(1) > xmid) then
        if (data(2) > ymid) then
          idirection = OTREE_NEB
        else
          idirection = OTREE_SEB
        end if
      else
        if (data(2) > ymid) then
          idirection = OTREE_NWB
        else
          idirection = OTREE_SWB
        end if
      end if

    end if

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_print,T)(roctree, cfilename)

!<description>
    ! This subroutine writes the content of the octree to a file
    ! which can be visualised by means of Matlab
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree

    ! filename of the output file
    character(LEN=*), intent(in) :: cfilename
!</input>
!</subroutine>

    ! local variables
    FEAT2_PP_TTYPE(T_TYPE) :: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: iunit

    iunit = sys_getFreeUnit()
    open(UNIT=iunit, FILE=trim(adjustl(cfilename)))
    xmin = roctree%p_BdBox(OTREE_XMIN, 1)
    ymin = roctree%p_BdBox(OTREE_YMIN, 1)
    zmin = roctree%p_BdBox(OTREE_ZMIN, 1)
    xmax = roctree%p_BdBox(OTREE_XMAX, 1)
    ymax = roctree%p_BdBox(OTREE_YMAX, 1)
    zmax = roctree%p_BdBox(OTREE_ZMAX, 1)
    call print(xmin, ymin, zmin, xmax, ymax, zmax, 1)
    close(UNIT=iunit)

  contains

    !**************************************************************
    ! Here, the recursive print routine follows

    recursive subroutine print(xmin, ymin, zmin, xmax, ymax, zmax, inode)
      FEAT2_PP_TTYPE(T_TYPE), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax
      integer, intent(in) :: inode
      FEAT2_PP_TTYPE(T_TYPE) :: xmid,ymid,zmid
      integer :: i,j

      write(UNIT=iunit,FMT=*) 'hex'
      write(UNIT=iunit,FMT=*) xmin, ymin, zmin, xmax, ymax, zmax, inode

      if (roctree%p_Knode(OTREE_STATUS, inode) .eq. OTREE_SUBDIV) then

        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)

        call print(xmin, ymid, zmin, xmid, ymax, zmid,&
                   -roctree%p_Knode(OTREE_NWF, inode))
        call print(xmin, ymin, zmin, xmid, ymid, zmid,&
                   -roctree%p_Knode(OTREE_SWF, inode))
        call print(xmid, ymin, zmin, xmax, ymid, zmid,&
                   -roctree%p_Knode(OTREE_SEF, inode))
        call print(xmid, ymid, zmin, xmax, ymax, zmid,&
                   -roctree%p_Knode(OTREE_NEF, inode))

        call print(xmin, ymid, zmid, xmid, ymax, zmax,&
                   -roctree%p_Knode(OTREE_NWB, inode))
        call print(xmin, ymin, zmid, xmid, ymid, zmax,&
                   -roctree%p_Knode(OTREE_SWB, inode))
        call print(xmid, ymin, zmid, xmax, ymid, zmax,&
                   -roctree%p_Knode(OTREE_SEB, inode))
        call print(xmid, ymid, zmid, xmax, ymax, zmax,&
                   -roctree%p_Knode(OTREE_NEB, inode))

      elseif (roctree%p_Knode(OTREE_STATUS, inode) > OTREE_EMPTY) then

        do i = 1, roctree%p_Knode(OTREE_STATUS, inode)
          j = roctree%p_Knode(i, inode)
          write(UNIT=iunit, FMT=*) 'node'
          write(UNIT=iunit, FMT=*) roctree%p_Data(1, j), roctree%p_Data(2, j),&
                                   roctree%p_Data(3, j), j
        end do

      end if

    end subroutine print

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_info,T)(roctree)

!<description>
    ! This subroutine outputs statistical info about the octree
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree
!</input>
!</subroutine>

    call output_line('Octree:')
    call output_line('-------')
    call output_line('NVT:     '//trim(sys_siL(roctree%NVT,15)))
    call output_line('NNVT:    '//trim(sys_siL(roctree%NNVT,15)))
    call output_line('NNODE:   '//trim(sys_siL(roctree%NNODE,15)))
    call output_line('NNNODE:  '//trim(sys_siL(roctree%NNNODE,15)))
    call output_line('NDATA:   '//trim(sys_siL(roctree%NDATA,15)))
    call output_line('NRESIZE: '//trim(sys_siL(roctree%NRESIZE,5)))
    call output_line('dfactor: '//trim(sys_sdL(roctree%dfactor,2)))
#ifdef T_STORAGE
    call output_line('h_Data: '//trim(sys_siL(roctree%h_Data,15)))
    call output_line('h_BdBox: '//trim(sys_siL(roctree%h_BdBox,15)))
#endif
    call output_line('h_Knode: '//trim(sys_siL(roctree%h_Knode,15)))
    call output_lbrk()
    call output_line('Current data memory usage: '//&
        trim(sys_sdL(100*roctree%NVT/real(roctree%NNVT,DP),2))//'%')
    call output_line('Current node memory usage: '//&
        trim(sys_sdL(100*roctree%NNODE/real(roctree%NNNODE,DP),2))//'%')

  end subroutine

  !************************************************************************

!<function>

  pure function FEAT2_PP_TEMPLATE_T(otree_getSize,T)(roctree) result(nvt)

!<description>
    ! This function returns the number of vertices stored in the octree
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree
!</input>

!<result>
    ! Number of vertices in octree
    integer :: nvt
!</result>
!</function>

    nvt = roctree%NVT

  end function

  !************************************************************************

!<function>

  function FEAT2_PP_TEMPLATE_T(otree_getBoundingBox,T)(roctree, inode) result(BdBox)

!<description>
    ! This function returns the bounding box of the specified node.
    ! If no node number is given, then the outer bounding box is returned.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree

    ! OPTIONAL: number of node for which bounding box should be returned
    integer, intent(in), optional :: inode
!</input>

!<result>
    ! bounding box
    FEAT2_PP_TTYPE(T_TYPE), dimension(6) :: BdBox
!</result>
!</function>

    if (present(inode)) then
      if (inode > roctree%NVT) then
        call output_line('Node number exceeds octree dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'otree_getBoundingBox')
        call sys_halt()
      end if
      BdBox = roctree%p_BdBox(OTREE_XMIN:OTREE_ZMAX, inode)
    else
      BdBox = roctree%p_BdBox(OTREE_XMIN:OTREE_ZMAX, 1)
    end if

  end function

  !************************************************************************

!<function>

  elemental function FEAT2_PP_TEMPLATE_T(otree_getX,T)(roctree, ivt) result(x)

!<description>
    ! This function returns the X-value at the given position.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree

    ! position in the octree
    integer, intent(in) :: ivt
!</input>

!<result>
    FEAT2_PP_TTYPE(T_TYPE) :: x
!</result>
!</function>

    x = roctree%p_Data(1,ivt)

  end function

  !************************************************************************

!<function>

  elemental function FEAT2_PP_TEMPLATE_T(otree_getY,T)(roctree, ivt) result(y)

!<description>
    ! This function returns the Y-value at the given position.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree

    ! position in the octree
    integer, intent(in) :: ivt
!</input>

!<result>
    FEAT2_PP_TTYPE(T_TYPE) :: y
!</result>
!</function>

    y = roctree%p_Data(2,ivt)

  end function

  !************************************************************************

!<function>

  elemental function FEAT2_PP_TEMPLATE_T(otree_getZ,T)(roctree, ivt) result(z)

!<description>
    ! This function returns the Z-value at the given position.
!</description>

!<input>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree

    ! position in the octree
    integer, intent(in) :: ivt
!</input>

!<result>
    FEAT2_PP_TTYPE(T_TYPE) :: z
!</result>
!</function>

    z = roctree%p_Data(3,ivt)

  end function

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_duplicate,T)(roctree, roctreeBackup)

!<description>
    ! This subroutine makes a copy of an octree in memory.
    ! It does not make sense to share some information between octrees,
    ! so each vectors is physically copied from the source octree
    ! to the destination octree.
!</description>

!<input>
    ! Source octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctree
!</input>

!<inputoutput>
    ! Destination octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctreeBackup
!</inputoutput>
!</subroutine>

   ! Release backup octree
    call otree_release(roctreeBackup)

    ! Copy all data
    roctreeBackup = roctree

    ! Reset handles
    roctreeBackup%h_Data  = ST_NOHANDLE
    roctreeBackup%h_BdBox = ST_NOHANDLE
    roctreeBackup%h_Knode = ST_NOHANDLE

    ! Copy storage blocks
    if (roctree%h_Knode .ne. ST_NOHANDLE) then
      call storage_copy(roctree%h_Knode, roctreeBackup%h_Knode)
      call storage_getbase(roctreeBackup%h_Knode,&
                           roctreeBackup%p_Knode)
    end if

#ifdef T_STORAGE
    if (roctree%h_Data .ne. ST_NOHANDLE) then
      call storage_copy(roctree%h_Data, roctreeBackup%h_Data)
      call storage_getbase(roctreeBackup%h_Data,&
                           roctreeBackup%p_Data)
    end if
#else
    if (associated(roctree%p_Data)) then
      allocate(roctreeBackup%h_Data(size(roctree%p_Data,1),&
                                    size(roctree%p_Data,2))
      roctreeBackup%p_Data = roctree%p_Data
    end if
#endif

#ifdef T_STORAGE
    if (roctree%h_BdBox .ne. ST_NOHANDLE) then
      call storage_copy(roctree%h_BdBox, roctreeBackup%h_BdBox)
      call storage_getbase(roctreeBackup%h_BdBox,&
                           roctreeBackup%p_BdBox)
    end if
#else
    if (associated(roctree%p_BdBox)) then
      allocate(roctreeBackup%h_Data(size(roctree%p_BdBox,1),&
                                    size(roctree%p_BdBox,2))
      roctreeBackup%p_BdBox = roctree%p_BdBox
    end if
#endif 

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_restore,T)(roctreeBackup, roctree)

!<description>
    ! This subroutine restores an octree from a previous backup.
!</description>

!<input>
    ! Backup of an octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(in) :: roctreeBackup
!</input>

!<inputoutput>
    ! Destination octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    ! Release octree
    call otree_release(roctree)

    ! Duplicate the backup
    call otree_duplicate(roctreeBackup, roctree)

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine FEAT2_PP_TEMPLATE_T(otree_rebuild,T)(roctree)

!<description>
    ! This subroutine rebuilds the structure of an octree
!</description>

!<inputoutput>
    ! Octree
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Inodes
    integer :: h_Inodes

    ! Initialization
    roctree%NNODE = 1
    roctree%NFREE = 0

    roctree%NNNODE = size(roctree%p_Knode,2)
    roctree%NNVT   = size(roctree%p_Data,2)

    ! Check if octree contains data
    if (roctree%NVT .eq. 0) then

      ! Initialise first node
      roctree%nnode = 1
      roctree%p_Knode(OTREE_STATUS,1) = OTREE_EMPTY
      roctree%p_Knode(OTREE_PARENT,1) = OTREE_EMPTY
      roctree%p_Knode(OTREE_PARPOS,1) = OTREE_EMPTY
      roctree%p_Knode(1:roctree%NDATA, 1) = 0

    else

      ! Create temporary storage
      h_Inodes = ST_NOHANDLE
      call storage_new('otree_rebuild', 'Inodes', roctree%NVT,&
                       ST_INT, h_Inodes, ST_NEWBLOCK_ORDERED)
      call storage_getbase_int(h_Inodes, p_Inodes)

      ! Create octree top-down
      call rebuild(1, 1, roctree%NVT)

      ! Release temporary storage
      call storage_free(h_Inodes)

    end if

  contains

    !**************************************************************
    ! Here, the recursive top-down rebuild routine follows

    recursive subroutine rebuild(inode, istart, iend)

      integer, intent(in) :: inode,istart, iend

      ! local variables
      FEAT2_PP_TTYPE(T_TYPE) :: xmin,ymin,zmin,xmax,ymax,zmax,xmid,ymid,zmid
      integer :: i,isize,jnode,nnode,imid1,imid2,imid3,imid4,imid5,imid6,imid7


      ! Check if istart > iend then the node is empty
      if (istart .gt. iend) return

      ! Check if current partition can be stored in a leave
      if (iend-istart+1 .le. roctree%NDATA) then

        ! Store vertex numbers
        do i = istart, iend
          roctree%p_Knode(i-istart+1,inode) = p_Inodes(i)
        end do

        ! Store number of vertices
        roctree%p_Knode(OTREE_STATUS, inode) = iend-istart+1

      else

        ! Otherwise, the current partition is subdivided into eight subpartitions

        if (roctree%nnode+OTREE_MAX > roctree%nnnode) then
          isize = max(roctree%nnode+OTREE_MAX,&
                      ceiling(roctree%dfactor*roctree%NNNODE))
          call resizeNNODE(roctree, isize)
        end if

        ! Store the number of nodes
        nnode = roctree%NNODE
        roctree%NNODE = nnode+OTREE_MAX

        ! Mark the current node as subdivided and set pointers to its eight children
        roctree%p_Knode(OTREE_STATUS,inode) =    OTREE_SUBDIV
        roctree%p_Knode(1:OTREE_MAX,inode)  = -(/nnode+OTREE_NWF, nnode+OTREE_SWF,&
                                                 nnode+OTREE_SEF, nnode+OTREE_NEF,&
                                                 nnode+OTREE_NWB, nnode+OTREE_SWB,&
                                                 nnode+OTREE_SEB, nnode+OTREE_NEB/)

        ! Determine coordinates for bounding boxes
        xmin = roctree%p_BdBox(OTREE_XMIN,inode)
        ymin = roctree%p_BdBox(OTREE_YMIN,inode)
        zmin = roctree%p_BdBox(OTREE_ZMIN,inode)
        xmax = roctree%p_BdBox(OTREE_XMAX,inode)
        ymax = roctree%p_BdBox(OTREE_YMAX,inode)
        zmax = roctree%p_BdBox(OTREE_ZMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)

        ! NWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWF) = OTREE_NWF
        roctree%p_BdBox(:,nnode+OTREE_NWF)            = (/xmin,ymid,zmin,xmid,ymax,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_NWF) = 0

        ! SWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWF) = OTREE_SWF
        roctree%p_BdBox(:,nnode+OTREE_SWF)            = (/xmin,ymin,zmin,xmid,ymid,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SWF) = 0

        ! SEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEF) = OTREE_SEF
        roctree%p_BdBox(:,nnode+OTREE_SEF)            = (/xmid,ymin,zmin,xmax,ymid,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SEF) = 0

        ! NEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEF) = OTREE_NEF
        roctree%p_BdBox(:,nnode+OTREE_NEF)            = (/xmid,ymid,zmin,xmax,ymax,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_NEF) = 0

        ! NWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWB) = OTREE_NWB
        roctree%p_BdBox(:,nnode+OTREE_NWB)            = (/xmin,ymid,zmid,xmid,ymax,zmax/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_NWB) = 0

        ! SWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWB) = OTREE_SWB
        roctree%p_BdBox(:,nnode+OTREE_SWB)            = (/xmin,ymin,zmid,xmid,ymid,zmax/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SWB) = 0

        ! SEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEB) = OTREE_SEB
        roctree%p_BdBox(:,nnode+OTREE_SEB)            = (/xmid,ymin,zmid,xmax,ymid,zmax/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SEB) = 0

        ! NEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEB) = OTREE_NEB
        roctree%p_BdBox(:,nnode+OTREE_NEB)            = (/xmid,ymid,zmid,xmax,ymax,zmax/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_NEB) = 0

        ! Swap nodes with respect to x-axis
        imid1 = partition(istart, iend, 1, xmid)

        ! Swap nodes with respect to y-axis
        imid2 = partition(istart, imid1, 2, ymid)
        imid3 = partition(imid1+1, iend, 2, ymid)

        ! Swap nodes with respect to z-axis
        imid4 = partition(istart,  imid2, 3, zmid)
        imid5 = partition(imid2+1, imid1, 3, zmid)
        imid6 = partition(imid1+1, imid3, 3, zmid)
        imid7 = partition(imid3+1, iend,  3, zmid)

        if (istart .lt. imid4) then
          jnode = -roctree%p_Knode(OTREE_SWF, inode)
          call rebuild(jnode, istart, imid4-1)
        end if

        if (imid4 .lt. imid2) then
          jnode = -roctree%p_Knode(OTREE_SWB, inode)
          call rebuild(jnode, imid4, imid2-1)
        end if

        if (imid2 .lt. imid5) then
          jnode = -roctree%p_Knode(OTREE_NWF, inode)
          call rebuild(jnode, imid2, imid5-1)
        end if

        if (imid5 .lt. imid1) then
          jnode = -roctree%p_Knode(OTREE_NWB, inode)
          call rebuild(jnode, imid5, imid1-1)
        end if

        if (imid1 .lt. imid6) then
          jnode = -roctree%p_Knode(OTREE_SEF, inode)
          call rebuild(jnode, imid1, imid6-1)
        end if

        if (imid6 .lt. imid3) then
          jnode = -roctree%p_Knode(OTREE_SEB, inode)
          call rebuild(jnode, imid6, imid3-1)
        end if

        if (imid3 .lt. imid7) then
          jnode = -roctree%p_Knode(OTREE_NEF, inode)
          call rebuild(jnode, imid3, imid7-1)
        end if

        if (imid7 .le. iend) then
          jnode = -roctree%p_Knode(OTREE_NEB, inode)
          call rebuild(jnode, imid7, iend)
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
        imid = merge(iend, iend+1, roctree%p_Data(idim, p_Inodes(iend)) .gt. dmid)

      else

        ! Otherwise, sort items in partition
        call quicksort(istart, iend, idim)

        ! Find location of first item belonging to "IS GREATER" partition
        do i = istart, iend
          if (roctree%p_Data(idim, p_Inodes(i)) .gt. dmid) then
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
      dpivot = roctree%p_Data(idim, p_Inodes(iend))

      do
        do while((i .lt. iend) .and. (roctree%p_Data(idim, p_Inodes(i)) .le. dpivot))
          i = i+1
        end do

        do while((j .gt. istart) .and. (roctree%p_Data(idim, p_Inodes(j)) .ge. dpivot))
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

!<subroutine>

  subroutine resizeNVT(roctree, nnvt)

!<description>
    ! This subroutine reallocates memory for an existing octree
!</description>

!<input>
    ! New number of vertices that should be stored in the octree
    integer, intent(in) :: nnvt
!</input>

!<inputoutput>
    ! Octree that should be resized
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)) :: roctree
!</inputoutput>
!</subroutine>

#ifndef T_STORAGE
    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_Data
#endif

#ifdef T_STORAGE
    call storage_realloc('resizeNVT', nnvt, roctree%h_Data,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase(roctree%h_Data, roctree%p_Data)
#else
    allocate(p_Ddata(size(roctree%p_Data,1),size(roctree%p_Data,2))
    p_Ddata = roctree%p_Data
    deallocate(roctree%p_Data)
    allocate(roctree%p_Data(size(p_Ddata,1),nnvt))
    roctree%p_Data(:,1:roctree%NNVT) = p_Data
    deallocate(p_Data)
#endif

    roctree%NNVT    = nnvt
    roctree%NRESIZE = roctree%NRESIZE+1

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine resizeNNODE(roctree, nnnode)

!<description>
    ! This subroutine reallocates memory for an existing octree
!</description>

!<input>
    ! New number of nodes that should be stored in the octree
    integer, intent(in) :: nnnode
!</input>

!<inputoutput>
    ! Octree that should be resized
    type(FEAT2_PP_TEMPLATE_T(t_octree,T)), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

#ifndef T_STORAGE
    ! local variables
    FEAT2_PP_TTYPE(T_TYPE), dimension(:,:), pointer :: p_BdBox
#endif

    call storage_realloc('resizeNNODE', nnnode, roctree%h_BdBox,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase(roctree%h_BdBox, roctree%p_BdBox)

#ifdef T_STORAGE
    call storage_realloc('resizeNNODE', nnnode, roctree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase(roctree%h_Knode, roctree%p_Knode)
#else
    allocate(p_Ddata(size(roctree%p_BdBox,1),size(roctree%p_BdBox,2))
    p_BdBox = roctree%p_BdBox
    deallocate(roctree%p_BdBox)
    allocate(roctree%p_BdBox(size(p_BdBox,1),nnnode))
    roctree%p_BdBox(:,1:roctree%NNNODE) = p_BdBox
    deallocate(p_BdBox)
#endif

    roctree%NNNODE  = nnnode
    roctree%NRESIZE = roctree%NRESIZE+1

  end subroutine

#endif

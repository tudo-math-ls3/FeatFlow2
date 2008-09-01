!##############################################################################
!# ****************************************************************************
!# <name> octree </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a (linear) octree.
!#
!# The following routines are available:
!#
!# 1.) otree_createOctree
!#     -> Create a new octree structure
!#
!# 2.) otree_releaseOctree
!#     -> Release an existing octree
!#
!# 3.) otree_resizeOctree
!#     -> Resize an existing octree
!#
!# 4.) otree_copyToOctree
!#     -> Copy data to an octree
!#
!# 5.) otree_copyFromOctree
!#     -> Copy data from the octree
!#
!# 6.) otree_insertIntoOctree
!#     -> Insert data into octree
!#
!# 7.) otree_deleteFromOctree
!#     -> Delete data from octree
!#
!# 8.) otree_searchInOctree
!#     -> Search data in octree
!#
!# 9.) otree_getDirection
!#     -> Get direction for the next node
!#
!# 10.) otree_printOctree
!#      -> Write octree to file
!#
!# 11.) otree_infoOctree
!#      -> Output info about octree
!#
!# 12.) otree_getsize
!#      -> Return number of vertices in octree
!#
!# 13.) otree_getBoundingBox
!#      -> Return the outer bounding box
!#
!# 14.) otree_getX
!#      -> Return the X-value at a given position
!#
!# 15.) otree_getY
!#      -> Return the Y-value at a given position
!#
!# 16.) otree_getZ
!#      -> Return the Z-value at a given position
!#
!# 17.) otree_duplicateOctree
!#      -> Create a duplicate / backup of an octree
!#
!# 18.) otree_restoreOctree
!#      -> Restore an octree from a previous backup
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

module octree
  use fsystem
  use storage
  use genoutput
  implicit none
  
  private
  public :: t_octree
  public :: otree_createOctree
  public :: otree_releaseOctree
  public :: otree_copyToOctree
  public :: otree_copyFromOctree
  public :: otree_insertIntoOctree
  public :: otree_deleteFromOctree
  public :: otree_searchInOctree
  public :: otree_getDirection
  public :: otree_printOctree
  public :: otree_infoOctree
  public :: otree_getsize
  public :: otree_getBoundingBox
  public :: otree_getX
  public :: otree_getY
  public :: otree_getZ
  public :: otree_duplicateOctree
  public :: otree_restoreOctree

!<constants>

!<constantblock description="KIND values for octree data">
  
  ! Kind value for indices in quadtree
  integer, parameter, public :: PREC_OTREEIDX = I32

!</constantblock>

!<constantblock description="Constants for octree structure">
  
  ! Position of the status information
  integer, parameter :: OTREE_STATUS = 11
  
  ! Position of the parent information
  integer, parameter :: OTREE_PARENT = 10

  ! Position of the "position" of the parent information
  integer, parameter :: OTREE_PARPOS = 9

  ! Number of free positions in cube
  integer, parameter :: OTREE_FREE   = 9

  ! Item in "North-East-Back" position
  integer, parameter :: OTREE_NEB    = 8

  ! Item in "South-East-Back" position
  integer, parameter :: OTREE_SEB    = 7

  ! Item in "South-West-Back" position
  integer, parameter :: OTREE_SWB     = 6

  ! Item in "North-West-Back" position
  integer, parameter :: OTREE_NWB    = 5

  ! Item in "North-East-Front" position
  integer, parameter :: OTREE_NEF    = 4

  ! Item in "South-East-Front" position
  integer, parameter :: OTREE_SEF    = 3

  ! Item in "South-West-Front" position
  integer, parameter :: OTREE_SWF    = 2

  ! Item in "North-West-Front" position
  integer, parameter :: OTREE_NWF    = 1

  ! Identifier: Cube is empty
  integer, parameter :: OTREE_EMPTY  = 0

  ! Identifier: Status is subdivided
  integer, parameter :: OTREE_SUBDIV = -1
  
  ! Maximum number of items for each quad
  integer, parameter :: OTREE_MAX    = 8

!</constantblock> 

!<constantblock description="Constants for octree bounding-box">

  ! Position of the x-minimal value
  integer, parameter :: OTREE_XMIN   =  1

  ! Position of the y-minimal value
  integer, parameter :: OTREE_YMIN   =  2

  ! Position of the z-minimal value
  integer, parameter :: OTREE_ZMIN   =  3

  ! Position of the x-maximal value
  integer, parameter :: OTREE_XMAX   =  4

  ! Position of the y-maximal value
  integer, parameter :: OTREE_YMAX   =  5

  ! Position of the z-maximal value
  integer, parameter :: OTREE_ZMAX   =  6

!</constantblock>   
  
!<constantblock description="Constants for octreetree operations">

  ! Item could not be found in the octree
  integer, parameter, public :: OTREE_NOT_FOUND = -1

  ! Item could be found in the octree
  integer, parameter, public :: OTREE_FOUND     =  0
!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! A linear octree implemented as array

  type t_octree

    ! Remark: The content of this derived data type is declared private.
    ! Hence, it cannot be accessed outside of this module. This allows us
    ! to use techniques such as the pointers which are used to increase
    ! the performace. These pointers cannot be modified externally so that
    ! data consistency is guaranteed.
    private

    ! Handle to data vector
    integer :: h_Ddata             = ST_NOHANDLE

    ! Handle to bounding box
    integer :: h_Dbbox             = ST_NOHANDLE

    ! Handle to octree structure
    !   KNODE(1:11,1:NNODE)
    !   KNODE( 11,INODE)    : < 0, the node has been subdivided
    !                         = 0, the node is empty
    !                         > 0, the number of points stored in the node
    !   KNODE( 10,INODE)    : > 0, the node the present node came from
    !   KNODE(  9,INODE)    : > 0, the position in the node the present 
    !                              node came from
    !   KNODE(1:8,INODE)    : for KNODE(11,INODE) > 0 : the points stored in the node
    !                         for KNODE(11,INODE) < 0 : the nodes into which the present 
    !                                                   node was subdivided
    integer :: h_Knode             = ST_NOHANDLE

    ! Data vectors
    ! NOTE: This array is introduced to increase performance. It should not be touched 
    ! by the user. To this end, all quantities of this derived type are PRIVATE.
    ! If the handle h_Ddata would be dereferences for each operation such as search, 
    ! delete, performance would be very poor.
    real(DP), dimension(:,:), pointer :: p_Ddata

    ! Coordinates of the bounding box
    ! NOTE: This array is introduced to increase performance (see above).
    real(DP), dimension(:,:), pointer :: p_Dbbox

    ! Octree structure
    ! NOTE: This array is introduced to increase performance (see above).
    integer(PREC_OTREEIDX), dimension(:,:), pointer :: p_Knode

    ! Number of vertices currently stored in the octree
    integer(PREC_OTREEIDX) :: NVT    = 0

    ! Total number of vertices that can be stored  in the octree
    integer(PREC_OTREEIDX) :: NNVT   = 0

    ! Number of nodes currently store in the octree
    integer(PREC_OTREEIDX) :: NNODE  = 0

    ! Total number of nodes that can be stored in the octree
    integer(PREC_OTREEIDX) :: NNNODE = 0

    ! Total number of resize operations
    integer :: NRESIZE               = 0

    ! Factor by which the octree is enlarged if new storage is allocate
    real(DP) :: dfactor              = 1.5_DP
  end type t_octree
  
!</typeblock>
!</types>
  
  interface otree_copyToOctree
    module procedure otree_copyToOctree_handle
    module procedure otree_copyToOctree_array
  end interface

  interface otree_copyFromOctree
    module procedure otree_copyFromOctree_handle
    module procedure otree_copyFromOctree_array
  end interface

  interface otree_deleteFromOctree
    module procedure otree_deleteFromOctree
    module procedure otree_deleteFromOctreeByNumber
  end interface

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine otree_createOctree(roctree, nnvt, nnnode,&
                                xmin, ymin, zmin, xmax, ymax, zmax, dfactor)
  
!<description>
    ! This subroutine creates a new octree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the octree
    integer(PREC_OTREEIDX), intent(IN) :: nnvt

    ! Total number of nodes that should be stored in the octree
    integer(PREC_OTREEIDX), intent(IN) :: nnnode

    ! Dimensions of the initial bounding box
    real(DP), intent(IN) :: xmin,ymin,zmin,xmax,ymax,zmax

    ! OPTIONAL: Factor by which the octree should be enlarged if
    ! new storage has to be allocated
    real(DP), intent(IN), optional :: dfactor
!</input>

!<output>
    ! Octree structure
    type(t_octree), intent(OUT) :: roctree
!</output>
!</subroutine>
    
    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1_DP) roctree%dfactor=dfactor
    end if
    
    ! Set values
    roctree%NNNODE  = nnnode
    roctree%NNVT    = nnvt
    roctree%NNODE   = 0
    roctree%NVT     = 0
    roctree%NRESIZE = 0
    
    ! Allocate memory and associate pointers
    call storage_new('otree_createOctree', 'p_Ddata', (/3,nnvt/),&
                     ST_DOUBLE, roctree%h_Ddata, ST_NEWBLOCK_ZERO)
    call storage_new('otree_createOctree', 'p_Dbbox', (/6,nnnode/),&
                     ST_DOUBLE, roctree%h_Dbbox, ST_NEWBLOCK_ZERO)
    call storage_new('otree_createOctree', 'p_Knode', (/11,nnnode/),&
                     ST_INT, roctree%h_Knode, ST_NEWBLOCK_ZERO)
    call storage_getbase_double2D(roctree%h_Ddata, roctree%p_Ddata)
    call storage_getbase_double2D(roctree%h_Dbbox, roctree%p_Dbbox)
    call storage_getbase_int2D(roctree%h_Knode,    roctree%p_Knode)

    ! Initialize first quad
    roctree%nnode                   = 1
    roctree%p_Knode(OTREE_STATUS,1) = OTREE_EMPTY
    roctree%p_Knode(OTREE_PARENT,1) = OTREE_EMPTY
    roctree%p_Knode(OTREE_FREE,1)   = 1
    
    ! Set co-ordinates of the bounding box
    roctree%p_Dbbox(OTREE_XMIN,1) = xmin
    roctree%p_Dbbox(OTREE_YMIN,1) = ymin
    roctree%p_Dbbox(OTREE_ZMIN,1) = zmin
    roctree%p_Dbbox(OTREE_XMAX,1) = xmax
    roctree%p_Dbbox(OTREE_YMAX,1) = ymax
    roctree%p_Dbbox(OTREE_ZMAX,1) = zmax
  end subroutine otree_createOctree

  !************************************************************************
  
!<subroutine>
  
  subroutine otree_releaseOctree(roctree)

!<description>
    ! This subroutine releases an existing octree
!</description>

!<inputoutput>
    ! quadtree
    type(t_octree), intent(INOUT) :: roctree
!</inputoutput>
!</subroutine>
    
    ! Free memory
    if (roctree%h_Ddata .ne. ST_NOHANDLE) call storage_free(roctree%h_Ddata)
    if (roctree%h_Dbbox .ne. ST_NOHANDLE) call storage_free(roctree%h_Dbbox)
    if (roctree%h_Knode .ne. ST_NOHANDLE) call storage_free(roctree%h_Knode)
    nullify(roctree%p_Knode, roctree%p_Dbbox, roctree%p_Ddata)

    ! Reset values
    roctree%NNNODE  = 0
    roctree%NNVT    = 0
    roctree%NNODE   = 0
    roctree%NVT     = 0
    roctree%NRESIZE = 0
  end subroutine otree_releaseOctree

  !************************************************************************

!<subroutine>
  
  subroutine otree_copyFromOctree_handle(roctree, h_Ddata)

!<description>
    ! This subroutine copies the content of the octree to a handle.
    ! If the handle is not associated, then a new handle with correct
    ! size is allocated. Otherwise, the handle is reallocated if 
    ! it does not provide enough memory.
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree
!</input>

!<inputoutput>
    ! Handle to the coordinate vector
    integer, intent(INOUT) :: h_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata,p_DdataTmp
    integer(PREC_OTREEIDX), dimension(2) :: Isize
    
    ! Check if handle is associated
    if (h_Ddata .eq. ST_NOHANDLE) then
      Isize = (/3,roctree%NVT/)
      call storage_new('otree_copyFromOctree_handle','p_Ddata',&
                       Isize, ST_DOUBLE, h_Ddata, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_Ddata, Isize)
      if (Isize(2) < roctree%NVT) then
        call storage_realloc('otree_copyFromOctree_handle',&
                              roctree%NVT, h_Ddata, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Set pointers
    call storage_getbase_double2D(h_Ddata, p_Ddata)
    call storage_getbase_double2D(roctree%h_Ddata, p_DdataTmp)

    ! Copy data
    call DCOPY(3*roctree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  end subroutine otree_copyFromOctree_handle

  !************************************************************************

!<subroutine>

  subroutine otree_copyFromOctree_array(roctree, p_Ddata)

!<description>
    ! This subroutine copies the content of the octree to an array.
    ! If the array is too small, then this subroutines terminates.
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree
!</input>

!<inputoutput>
    ! Coordinate vector
    real(DP), dimension(:,:), intent(INOUT) :: p_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DdataTmp
    integer(PREC_OTREEIDX), dimension(2) :: Isize

    ! Check size of array
    Isize = shape(p_Ddata)
    if (Isize(1) .ne. 3 .or. Isize(2) < roctree%NVT) then
      call output_line('Array too small!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyFromOctree_array')
      call sys_halt()
    end if

    ! Set pointers
    call storage_getbase_double2D(roctree%h_Ddata, p_DdataTmp)

    ! Copy data
    call DCOPY(3*roctree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  end subroutine otree_copyFromOctree_array

  !************************************************************************

!<subroutine>
  
  subroutine otree_copyToOctree_handle(h_Ddata, roctree)

!<description>
    ! This subroutine copies the content of a handle to the octree.
!</description>

!<input>
    ! Handle to the coordinate vector
    integer, intent(IN) :: h_Ddata
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(INOUT) :: roctree
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata
    integer(PREC_OTREEIDX) :: ivt,jvt,inode,ipos
    
    ! Set pointer and check its shape
    call storage_getbase_double2D(h_Ddata, p_Ddata)
    if (size(p_Ddata, 1) .ne. 3) then
      call output_line('First dimension of array must be 3!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyToOctree_handle')
      call sys_halt()
    end if
    
    do ivt = 1, size(p_Ddata, 2)
      if (otree_searchInOctree(roctree, p_Ddata(:,ivt),&
                               inode, ipos,jvt) .eq. OTREE_NOT_FOUND)&
          call otree_insertIntoOctree(roctree, ivt, p_Ddata(:,ivt), inode)
    end do
  end subroutine otree_copyToOctree_handle

  !************************************************************************

!<subroutine>
  
  subroutine otree_copyToOctree_array(p_Ddata, roctree)

!<description>
    ! This subroutine copies the content of an array to the octree.
!</description>

!<input>
    ! Coordinate vector
    real(DP), dimension(:,:), intent(IN) :: p_Ddata
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(INOUT) :: roctree
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(PREC_OTREEIDX) :: ivt,jvt,inode,ipos
    
    if (size(p_Ddata,1) .ne. 3) then
      call output_line('First dimension of array must be 3!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyToOctree_array')
      call sys_halt()
    end if

    do ivt = 1, size(p_Ddata,2)
      if (otree_searchInOctree(roctree, p_Ddata(:,ivt),&
                               inode, ipos, jvt) .eq. OTREE_NOT_FOUND)&
          call otree_insertIntoOctree(roctree, ivt, p_Ddata(:,ivt), inode)
    end do
  end subroutine otree_copyToOctree_array

  !************************************************************************
  
!<subroutine>
  
  subroutine otree_insertIntoOctree(roctree, ivt, Ddata, inode)

!<description>
    ! This subroutine inserts a new coordinate item into the octree
!</description>

!<input>
    ! Number of the inserted vertex
    integer(PREC_OTREEIDX), intent(IN) :: ivt

    ! Number of the node into which vertex is inserted
    integer(PREC_OTREEIDX), intent(IN) :: inode
    
    ! Coordinates of the new vertex
    real(DP), dimension(3), intent(IN) :: Ddata
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(INOUT)      :: roctree
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer(PREC_OTREEIDX) :: isize

    ! Check if there is enough space left in the vertex component of the octree
    if (roctree%NVT .eq. roctree%NNVT) then
      isize = max(roctree%NNVT+1, ceiling(roctree%dfactor*roctree%NNVT))
      call resizeNVT(roctree, isize)
    end if
    
    ! Update values and add new entry recursively
    roctree%NVT            = roctree%NVT+1
    roctree%p_Ddata(:,ivt) = Ddata
    call insert(ivt, inode)
    
  contains  
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    recursive subroutine insert(ivt, inode)
      integer(PREC_OTREEIDX), intent(IN) :: ivt,inode
      real(DP) :: xmin,ymin,zmin,xmax,ymax,zmax,xmid,ymid,zmid
      integer(PREC_OTREEIDX) :: i,jnode,knode,nnode,isize
      
      if (roctree%p_Knode(OTREE_STATUS, inode) .eq. OTREE_MAX) then
        
        if (roctree%nnode+OTREE_MAX > roctree%nnnode) then
          isize = max(roctree%nnode+OTREE_MAX, ceiling(roctree%dfactor*roctree%NNNODE))
          call resizeNNODE(roctree, isize)
        end if
        
        ! Node is full and needs to be refined into eight new nodes
        xmin = roctree%p_Dbbox(OTREE_XMIN,inode)
        ymin = roctree%p_Dbbox(OTREE_YMIN,inode)
        zmin = roctree%p_Dbbox(OTREE_ZMIN,inode)
        xmax = roctree%p_Dbbox(OTREE_XMAX,inode)
        ymax = roctree%p_Dbbox(OTREE_YMAX,inode)
        zmax = roctree%p_Dbbox(OTREE_ZMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)
        
        ! Store the number of nodes
        nnode         = roctree%NNODE
        roctree%NNODE = nnode+OTREE_MAX
        
        ! NWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWF) = OTREE_NWF
        roctree%p_Dbbox(:,nnode+OTREE_NWF)            = (/xmin,ymid,zmin,xmid,ymax,zmid/)
        
        ! SWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWF) = OTREE_SWF
        roctree%p_Dbbox(:,nnode+OTREE_SWF)            = (/xmin,ymin,zmin,xmid,ymid,zmid/)
        
        ! SEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEF) = OTREE_SEF
        roctree%p_Dbbox(:,nnode+OTREE_SEF)            = (/xmid,ymin,zmin,xmax,ymid,zmid/)
        
        ! NEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEF) = OTREE_NEF
        roctree%p_Dbbox(:,nnode+OTREE_NEF)            = (/xmid,ymid,zmin,xmax,ymax,zmid/)
        
        ! NWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWB) = OTREE_NWB
        roctree%p_Dbbox(:,nnode+OTREE_NWB)            = (/xmin,ymid,zmid,xmid,ymax,zmax/)
        
        ! SWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWB) = OTREE_SWB
        roctree%p_Dbbox(:,nnode+OTREE_SWB)            = (/xmin,ymin,zmid,xmid,ymid,zmax/)
        
        ! SEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEB) = OTREE_SEB
        roctree%p_Dbbox(:,nnode+OTREE_SEB)            = (/xmid,ymin,zmid,xmax,ymid,zmax/)
        
        ! NEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEB) = OTREE_NEB
        roctree%p_Dbbox(:,nnode+OTREE_NEB)            = (/xmid,ymid,zmid,xmax,ymax,zmax/)

        ! Add the eight values from INODE to the eight new nodes
        ! NNODE+1:NNODE+8 recursively
        do i = 1, OTREE_MAX
          knode = roctree%p_Knode(i,inode)
          jnode = nnode+otree_getDirection(roctree, roctree%p_Ddata(:,knode),inode)
          call insert(knode, jnode)
        end do
        
        ! Mark the current nodes as subdivided and set pointers to its eight children 
        roctree%p_Knode(OTREE_STATUS,inode) = OTREE_SUBDIV
        roctree%p_Knode(1:OTREE_MAX, inode) = (/nnode+OTREE_NWF, nnode+OTREE_SWF,&
                                                nnode+OTREE_SEF, nnode+OTREE_NEF,&
                                                nnode+OTREE_NWB, nnode+OTREE_SWB,&
                                                nnode+OTREE_SEB, nnode+OTREE_NEB/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+otree_getDirection(roctree, roctree%p_Ddata(:,ivt),inode)
        call insert(ivt, jnode)
        
      else
        
        ! Node is not full, so new items can be stored
        roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)+1
        roctree%p_Knode(roctree%p_Knode(OTREE_STATUS, inode), inode) = ivt
      end if
    end subroutine insert
  end subroutine otree_insertIntoOctree
  
  ! ***************************************************************************
  
!<function>
  
  function otree_deleteFromOctree(roctree, Ddata, ivt) result(f)

!<description>
    ! This function deletes an item from the octree
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    real(DP), dimension(3), intent(IN) :: Ddata
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(INOUT) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that is deleted
    integer(PREC_OTREEIDX), intent(OUT) :: ivt
!</output>

!<result>
    ! Result of the deletion: OTREE_NOT_FOUND / OTREE_FOUND
    integer :: f
!</result>
!</function>
    
    ! local variables
    integer :: inode,ipos,jpos,jvt
    real(DP), dimension(3) :: DdataTmp
    
    ! Search for the given coordinates
    f = otree_searchInOctree(roctree, Ddata, inode, ipos, ivt)
    
    ! What can we do from the searching
    if (f .eq. OTREE_FOUND) then
      
      ! Remove item IVT from node INODE
      do jpos = ipos+1, roctree%p_Knode(OTREE_STATUS, inode)
        roctree%p_Knode(jpos-1, inode) = roctree%p_Knode(jpos, inode)
      end do
      roctree%p_Knode(roctree%p_Knode(OTREE_STATUS, inode), inode) = 0
      roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)-1
      
      ! If IVT is not last item move last item NVT to position IVT
      if (ivt .ne. roctree%NVT) then
        DdataTmp(1:3) = roctree%p_Ddata(1:3, roctree%NVT)
        if (otree_searchInOctree(roctree, DdataTmp(:),&
                                 inode, ipos, jvt) .eq. OTREE_FOUND) then
          roctree%p_Ddata(:, ivt) = roctree%p_Ddata(:, roctree%NVT)
          roctree%p_Knode(ipos, inode) = ivt
        end if
        ivt = roctree%NVT
      end if
      roctree%NVT = roctree%NVT-1
    end if
  end function otree_deleteFromOctree

  ! ***************************************************************************

!<function>

  function otree_deleteFromOctreeByNumber(roctree, ivt, ivtReplace) result(f)

!<description>
    ! This function deletes vertex with number IVT from the octree
!</description>

!<input>
    ! Number of the vertex to be deleted
    integer(PREC_OTREEIDX), intent(IN) :: ivt
!</input>

!<inputoutput>
    ! octree
    type(t_octree), intent(INOUT) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that replaces the deleted vertex
    integer(PREC_OTREEIDX), intent(OUT) :: ivtReplace
!</output>

!<result>
    ! Result of the deletion: OTREE_NOT_FOUND / OTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    real(DP), dimension(3) :: Ddata
    
    if (ivt .le. roctree%NVT) then
      ! Get coordinates and invoke deletion routine
      Ddata = roctree%p_Ddata(:,ivt)
      f     = otree_deleteFromOctree(roctree, Ddata, ivtReplace)
    else
      call output_line('Invalid vertex number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'otree_deleteFromOctreeByNumber')
      call sys_halt()
    end if
  end function otree_deleteFromOctreeByNumber

  ! ***************************************************************************

!<function>
  
  function otree_searchInOctree(roctree, Ddata, inode, ipos, ivt) result(f)

!<description>
    ! This subroutine searches for given coordinates in the octree
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree
    
    ! Coordinates that should be searched for
    real(DP), dimension(3), intent(IN) :: Ddata
!</input>

!<output>
    ! Number of the node in which the given coordinates are
    integer(PREC_OTREEIDX), intent(OUT) :: inode

    ! Position of the coordinates in the node
    integer(PREC_OTREEIDX), intent(OUT) :: ipos

    ! Number of the vertex the coordinates correspond to
    integer(PREC_OTREEIDX), intent(OUT) :: ivt
!</output>

!<result>
    ! Result of the searching: OTREE_NOT_FOUND / OTREE_FOUND
    integer :: f
!</result>
!</function>
    
    ! Initialize
    inode = 1
    ipos  = 1
    ivt   = 1
    f     = search(inode, ipos, ivt)
    
  contains
    
    !**************************************************************
    ! Here, the recursive searching routine follows

    recursive function search(inode, ipos, ivt) result(f)
      integer(PREC_OTREEIDX), intent(INOUT) :: inode,ipos,ivt
      integer :: f
      
      f = OTREE_NOT_FOUND
      if (roctree%p_Knode(OTREE_STATUS, inode) .eq. OTREE_SUBDIV) then
        
        ! Node is subdivided. Compute child INODE which to look recursively.
        inode = roctree%p_Knode(otree_getDirection(roctree, Ddata, inode), inode)
        f     = search(inode,ipos, ivt)
        
      else
        
        ! Node is not subdivided. Search for (x,y,z) in current node
        do ipos = 1, roctree%p_Knode(OTREE_STATUS, inode)
          ivt = roctree%p_Knode(ipos, inode)
          if (sqrt(sum((roctree%p_Ddata(:, ivt)-Ddata)**2)) .le. SYS_EPSREAL) then
            f = OTREE_FOUND; return
          end if
        end do
        
      end if
    end function search
  end function otree_searchInOctree

  !************************************************************************
  
!<function>
  
  pure function otree_getDirection(roctree, Ddata, inode) result(d)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Ddata
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree
    
    ! Coordinates
    real(DP), dimension(3), intent(IN) :: Ddata

    ! Number of node
    integer(PREC_OTREEIDX), intent(IN) :: inode
!</input>

!<result>
    ! Further search direction
    integer :: d
!</result>
!</function>
    
    ! local variables
    real(DP) :: xmid,ymid,zmid

    ! Compute midpoint of current node
    xmid=0.5*(roctree%p_Dbbox(OTREE_XMIN, inode)+&
              roctree%p_Dbbox(OTREE_XMAX, inode))
    ymid=0.5*(roctree%p_Dbbox(OTREE_YMIN, inode)+&
              roctree%p_Dbbox(OTREE_YMAX, inode))
    zmid=0.5*(roctree%p_Dbbox(OTREE_ZMIN, inode)+&
              roctree%p_Dbbox(OTREE_ZMAX, inode))

    ! Do we have to look in the front or back?
    if (Ddata(3) < zmid) then
      
      if (Ddata(1) > xmid) then
        if (Ddata(2) > ymid) then
          d = OTREE_NEF; return
        else
          d = OTREE_SEF; return
        end if
      else
        if (Ddata(2) > ymid) then
          d = OTREE_NWF; return
        else
          d = OTREE_SWF; return
        end if
      end if

    else

      if (Ddata(1) > xmid) then
        if (Ddata(2) > ymid) then
          d = OTREE_NEB; return
        else
          d = OTREE_SEB; return
        end if
      else
        if (Ddata(2) > ymid) then
          d = OTREE_NWB; return
        else
          d = OTREE_SWB; return
        end if
      end if

    end if
  end function otree_getDirection

  !************************************************************************
  
!<subroutine>

  subroutine otree_printOctree(roctree, cfilename)

!<description>
    ! This subroutine writes the content of the octree to a file
    ! which can be visualized by means of Matlab
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree

    ! filename of the output file
    character(LEN=*), intent(IN) :: cfilename
!</input>
!</subroutine>
    
    ! local variables
    real(DP) :: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: iunit
    
    iunit = sys_getFreeUnit()
    open(UNIT=iunit, FILE=trim(adjustl(cfilename)))
    xmin = roctree%p_Dbbox(OTREE_XMIN, 1)
    ymin = roctree%p_Dbbox(OTREE_YMIN, 1)
    zmin = roctree%p_Dbbox(OTREE_ZMIN, 1)
    xmax = roctree%p_Dbbox(OTREE_XMAX, 1)
    ymax = roctree%p_Dbbox(OTREE_YMAX, 1)
    zmax = roctree%p_Dbbox(OTREE_ZMAX, 1)
    call print(xmin, ymin, zmin, xmax, ymax, zmax, 1)
    close(UNIT=iunit)
    
  contains
    
    !**************************************************************
    ! Here, the recursive print routine follows
    
    recursive subroutine print(xmin, ymin, zmin, xmax, ymax, zmax, inode)
      real(DP), intent(IN) :: xmin,ymin,zmin,xmax,ymax,zmax
      integer(PREC_OTREEIDX), intent(IN) :: inode
      real(DP)               :: xmid,ymid,zmid
      integer(PREC_OTREEIDX) :: i

      write(UNIT=iunit,FMT=*) 'hex'
      write(UNIT=iunit,FMT=10) xmin, ymin, zmin, xmax, ymax, zmax
      
      if (roctree%p_Knode(OTREE_STATUS, inode) .eq. OTREE_SUBDIV) then
        
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)
        
        call print(xmin, ymid, zmin, xmid, ymax, zmid,&
                   roctree%p_Knode(OTREE_NWF, inode))
        call print(xmin, ymin, zmin, xmid, ymid, zmid,&
                   roctree%p_Knode(OTREE_SWF, inode))
        call print(xmid, ymin, zmin, xmax, ymid, zmid,&
                   roctree%p_Knode(OTREE_SEF, inode))
        call print(xmid, ymid, zmin, xmax, ymax, zmid,&
                   roctree%p_Knode(OTREE_NEF, inode))

        call print(xmin, ymid, zmid, xmid, ymax, zmax,&
                   roctree%p_Knode(OTREE_NWB, inode))
        call print(xmin, ymin, zmid, xmid, ymid, zmax,&
                   roctree%p_Knode(OTREE_SWB, inode))
        call print(xmid, ymin, zmid, xmax, ymid, zmax,&
                   roctree%p_Knode(OTREE_SEB, inode))
        call print(xmid, ymid, zmid, xmax, ymax, zmax,&
                   roctree%p_Knode(OTREE_NEB, inode))

      elseif (roctree%p_Knode(OTREE_STATUS, inode) > OTREE_EMPTY) then
        
        do i = 1, roctree%p_Knode(OTREE_STATUS, inode)
          write(UNIT=iunit, FMT=*) 'node'
          write(UNIT=iunit, FMT=20) roctree%p_Ddata(:, roctree%p_Knode(i, inode))
        end do
        
      end if
      
10    format(6E15.6E3)
20    format(3E15.6E3)
    end subroutine print
  end subroutine otree_printOctree

  !************************************************************************

!<subroutine>

  subroutine otree_infoOctree(roctree)

!<description>
    ! This subroutine outputs statistical info about the octree
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree
!</input>
!</subroutine>

    call output_line('Octree:')
    call output_line('-------')
    call output_line('NVT:     '//trim(sys_siL(roctree%NVT,15)))
    call output_line('NNVT:    '//trim(sys_siL(roctree%NNVT,15)))
    call output_line('NNODE:   '//trim(sys_siL(roctree%NNODE,15)))
    call output_line('NNNODE:  '//trim(sys_siL(roctree%NNNODE,15)))
    call output_line('NRESIZE: '//trim(sys_siL(roctree%NRESIZE,5)))
    call output_line('dfactor: '//trim(sys_sdL(roctree%dfactor,2)))
    call output_line('h_Ddata: '//trim(sys_siL(roctree%h_Ddata,15)))
    call output_line('h_Dbbox: '//trim(sys_siL(roctree%h_Dbbox,15)))
    call output_line('h_Knode: '//trim(sys_siL(roctree%h_Knode,15)))
    call output_lbrk()
    write(*,*)
    call output_line('Current data memory usage: '//&
        trim(sys_sdL(100*roctree%NVT/real(roctree%NNVT,DP),2))//'%')
    call output_line('Current node memory usage: '//&
        trim(sys_sdL(100*roctree%NNODE/real(roctree%NNNODE,DP),2))//'%')
  end subroutine otree_infoOctree

  !************************************************************************

!<function>

  pure function otree_getSize(roctree) result(nvt)

!<description>
    ! This function returns the number of vertices stored in the octree
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree
!</input>

!<result>
    ! Number of vertices in octree
    integer(PREC_OTREEIDX) :: nvt
!</result>
!</function>

    nvt = roctree%NVT
  end function otree_getSize

  !************************************************************************

!<function>

  function otree_getBoundingBox(roctree, inode) result(bbox)
    
!<description>
    ! This function returns the bounding box of the specified node.
    ! If no node number is given, then the outer bounding box is returned.
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree

    ! OPTIONAL: number of node for which bounding box should be returned
    integer, intent(IN), optional :: inode
!</input>

!<result>
    ! bounding box
    real(DP), dimension(6) :: bbox
!</result>
!</function>
    
    if (present(inode)) then
      if (inode > roctree%NVT) then
        call output_line('Node number exceeds octree dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'otree_getBoundingBox')
        call sys_halt()
      end if
      bbox = roctree%p_Dbbox(OTREE_XMIN:OTREE_ZMAX, inode)
    else
      bbox = roctree%p_Dbbox(OTREE_XMIN:OTREE_ZMAX, 1)
    end if
  end function otree_getBoundingBox

  !************************************************************************

!<function>

  elemental function otree_getX(roctree, ivt) result(x)

!<description>
    ! This function returns the X-value at the given position.
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree

    ! position in the octree
    integer(PREC_OTREEIDX), intent(IN) :: ivt
!</input>

!<result>
    real(DP) :: x
!</result>
!</function>

    x = roctree%p_Ddata(1,ivt)
  end function otree_getX

  !************************************************************************

!<function>

  elemental function otree_getY(roctree, ivt) result(y)

!<description>
    ! This function returns the Y-value at the given position.
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree

    ! position in the octree
    integer(PREC_OTREEIDX), intent(IN) :: ivt
!</input>

!<result>
    real(DP) :: y
!</result>
!</function>

    y = roctree%p_Ddata(2,ivt)
  end function otree_getY

  !************************************************************************

!<function>

  elemental function otree_getZ(roctree, ivt) result(z)

!<description>
    ! This function returns the Z-value at the given position.
!</description>

!<input>
    ! Octree
    type(t_octree), intent(IN) :: roctree

    ! position in the octree
    integer(PREC_OTREEIDX), intent(IN) :: ivt
!</input>

!<result>
    real(DP) :: z
!</result>
!</function>

    z = roctree%p_Ddata(3,ivt)
  end function otree_getZ

  !************************************************************************

!<subroutine>

  subroutine otree_duplicateOctree(roctree, roctreeBackup)

!<description>
    ! This subroutine makes a copy of an octree in memory.
    ! It does not make sense to share some information between octrees,
    ! so each vectors is physically copied from the source octree
    ! to the destination octree.
!</description>

!<input>
    ! Source octree
    type(t_octree), intent(IN) :: roctree
!</input>

!<inputoutput>
    ! Destination octree
    type(t_octree), intent(INOUT) :: roctreeBackup
!</inputoutput>
!</subroutine>

    ! Release backup octree
    call otree_releaseOctree(roctreeBackup)

    ! Copy all data
    roctreeBackup = roctree

    ! Reset handles
    roctreeBackup%h_Ddata = ST_NOHANDLE
    roctreeBackup%h_Dbbox = ST_NOHANDLE
    roctreeBackup%h_Knode = ST_NOHANDLE

    ! Copy storage blocks
    if (roctree%h_Ddata .ne. ST_NOHANDLE) then
      call storage_copy(roctree%h_Ddata, roctreeBackup%h_Ddata)
      call storage_getbase_double2D(roctreeBackup%h_Ddata,&
                                    roctreeBackup%p_Ddata)
    end if

    if (roctree%h_Dbbox .ne. ST_NOHANDLE) then
      call storage_copy(roctree%h_Dbbox, roctreeBackup%h_Dbbox)
      call storage_getbase_double2D(roctreeBackup%h_Dbbox,&
                                    roctreeBackup%p_Dbbox)
    end if

    if (roctree%h_Knode .ne. ST_NOHANDLE) then
      call storage_copy(roctree%h_Knode, roctreeBackup%h_Knode)
      call storage_getbase_int2D(roctreeBackup%h_Knode,&
                                 roctreeBackup%p_Knode)
    end if
  end subroutine otree_duplicateOctree

  !************************************************************************

!<subroutine>

  subroutine otree_restoreOctree(roctreeBackup, roctree)

!<description>
    ! This subroutine restores an octree from a previous backup.
!</description>

!<input>
    ! Backup of an octree
    type(t_octree), intent(IN) :: roctreeBackup
!</input>

!<inputoutput>
    ! Destination octree
    type(t_octree), intent(INOUT) :: roctree
!</inputoutput>
!</subroutine>

    ! Release octree
    call otree_releaseOctree(roctree)

    ! Duplicate the backup
    call otree_duplicateOctree(roctreeBackup, roctree)
  end subroutine otree_restoreOctree

  !************************************************************************
  
!<subroutine>
  
  subroutine resizeNVT(roctree, nnvt)

!<description>
    ! This subroutine reallocates memory for an existing octree
!</description>

!<input>
    ! New number of vertices that should be stored in the octree
    integer(PREC_OTREEIDX), intent(IN) :: nnvt
!</input>

!<inputoutput>
    ! Octree that should be resized
    type(t_octree) :: roctree
!</inputoutput>
!</subroutine>

    call storage_realloc('resizeNVT', nnvt, roctree%h_Ddata,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_double2D(roctree%h_Ddata, roctree%p_Ddata)
    
    roctree%NNVT    = nnvt
    roctree%NRESIZE = roctree%NRESIZE+1
  end subroutine resizeNVT
  
  !************************************************************************
  
!<subroutine>
  
  subroutine resizeNNODE(roctree, nnnode)

!<description>
    ! This subroutine reallocates memory for an existing octree
!</description>

!<input>
    ! New number of nodes that should be stored in the octree
    integer(PREC_OTREEIDX), intent(IN) :: nnnode
!</input>

!<inputoutput>
    ! Octree that should be resized
    type(t_octree) :: roctree
!</inputoutput>
!</subroutine>

    call storage_realloc('resizeNNODE', nnnode, roctree%h_Dbbox,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_realloc('resizeNNODE', nnnode, roctree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_double2D(roctree%h_Dbbox, roctree%p_Dbbox)
    call storage_getbase_int2D(roctree%h_Knode,    roctree%p_Knode)
    
    roctree%NNNODE  = nnnode
    roctree%NRESIZE =roctree%NRESIZE+1
  end subroutine resizeNNODE
end module octree


!##############################################################################
!# ****************************************************************************
!# <name> octree </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module implements a (linear) octree. The implementation is
!# based on the description of octrees by
!#
!# R. Lohner, Applied CFD Techniques. An Introduction based on
!#            Finite Element Methods, Wiley, 2008
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
!#     -> Copy data to the octree
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
!# 19.) otree_rebuildOctree
!#      -> Rebuilds the structure of an octree
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
  public :: otree_rebuildOctree

!<constants>

!<constantblock description="Constants for octree structure">
  
  ! Maximum number of items for each node
  integer, parameter :: OTREE_MAX    = 8

  ! Item in "North-West-Front" position
  integer, parameter :: OTREE_NWF    = 1

  ! Item in "South-West-Front" position
  integer, parameter :: OTREE_SWF    = 2

  ! Item in "South-East-Front" position
  integer, parameter :: OTREE_SEF    = 3

  ! Item in "North-East-Front" position
  integer, parameter :: OTREE_NEF    = 4

  ! Item in "North-West-Back" position
  integer, parameter :: OTREE_NWB    = 5

  ! Item in "South-West-Back" position
  integer, parameter :: OTREE_SWB    = 6

  ! Item in "South-East-Back" position
  integer, parameter :: OTREE_SEB    = 7

  ! Item in "North-East-Back" position
  integer, parameter :: OTREE_NEB    = 8
  
  ! Position of the status information
  integer, parameter :: OTREE_STATUS = 0
  
  ! Position of the parent information
  integer, parameter :: OTREE_PARENT = -1

  ! Position of the "position" of the parent information
  integer, parameter :: OTREE_PARPOS = -2

  ! Identifier: Node is empty
  integer, parameter :: OTREE_EMPTY  =  0

  ! Identifier: Node is subdivided
  integer, parameter :: OTREE_SUBDIV = -1

  ! Identifier: Node is deleted
  integer, parameter :: OTREE_DEL    = -2

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
  
!<constantblock description="Constants for octree operations">

  ! Operation on octree failed
  integer, parameter, public :: OTREE_FAILED    = -2

  ! Item could not be found in the octree
  integer, parameter, public :: OTREE_NOT_FOUND = -1

  ! Item could be found in the octree
  integer, parameter, public :: OTREE_FOUND     =  0

  ! Item was inserted into the octree
  integer, parameter, public :: OTREE_INSERTED  =  1

  ! Item was deleted from the octree
  integer, parameter, public :: OTREE_DELETED   =  2
  
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
    integer :: h_Ddata = ST_NOHANDLE

    ! Handle to bounding box
    integer :: h_Dbbox = ST_NOHANDLE

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
    integer, dimension(:,:), pointer :: p_Knode

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
                                xmin, ymin, zmin, xmax, ymax, zmax, dfactor, ndata)
  
!<description>
    ! This subroutine creates a new octree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the octree
    integer, intent(in) :: nnvt

    ! Total number of nodes that should be stored in the octree
    integer, intent(in) :: nnnode

    ! Dimensions of the initial bounding box
    real(DP), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax

    ! OPTIONAL: Factor by which the octree should be enlarged if
    ! new storage has to be allocated
    real(DP), intent(in), optional :: dfactor

    ! OPTIONAL: Number of data items stored per node
    integer, optional :: ndata
!</input>

!<output>
    ! Octree structure
    type(t_octree), intent(out) :: roctree
!</output>
!</subroutine>
    
    ! local variables
    integer, dimension(2) :: Isize, Ilbound, Iubound

    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1_DP) roctree%dfactor=dfactor
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
    Isize = (/3, nnvt/)
    call storage_new('otree_createOctree', 'p_Ddata',&
                     Isize, ST_DOUBLE, roctree%h_Ddata, ST_NEWBLOCK_ZERO)

    Isize = (/6, nnnode/)
    call storage_new('otree_createOctree', 'p_Dbbox',&
                     Isize, ST_DOUBLE, roctree%h_Dbbox, ST_NEWBLOCK_ZERO)

    Ilbound = (/OTREE_PARPOS,1/); Iubound = (/roctree%NDATA,nnnode/)
    call storage_new('otree_createOctree', 'p_Knode', Ilbound,&
                     Iubound, ST_INT, roctree%h_Knode, ST_NEWBLOCK_ZERO)
    call storage_getbase_double2D(roctree%h_Ddata, roctree%p_Ddata)
    call storage_getbase_double2D(roctree%h_Dbbox, roctree%p_Dbbox)
    call storage_getbase_int2D(roctree%h_Knode,    roctree%p_Knode)

    ! Initialize first node
    roctree%nnode = 1
    roctree%p_Knode(OTREE_STATUS,1) = OTREE_EMPTY
    roctree%p_Knode(OTREE_PARENT,1) = OTREE_EMPTY
    roctree%p_Knode(OTREE_PARPOS,1) = 1
    roctree%p_Knode(1:roctree%ndata, 1) = 0
    
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
    ! Octree
    type(t_octree), intent(inout) :: roctree
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
    roctree%NFREE   = 0
    
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
    type(t_octree), intent(in) :: roctree
!</input>

!<inputoutput>
    ! Handle to the coordinate vector
    integer, intent(inout) :: h_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata
    integer, dimension(2) :: Isize
    
    ! Check if handle is associated
    if (h_Ddata .eq. ST_NOHANDLE) then
      Isize = (/3, roctree%NVT/)
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

    ! Call copy routine
    call otree_copyFromOctree(roctree, p_Ddata)

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
    type(t_octree), intent(in) :: roctree
!</input>

!<inputoutput>
    ! Coordinate vector
    real(DP), dimension(:,:), intent(inout) :: p_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize

    ! Check size of array
    Isize = shape(p_Ddata)
    if (Isize(1) .ne. 3 .or. Isize(2) < roctree%NVT) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyFromOctree_array')
      call sys_halt()
    end if

    ! Copy data
    call DCOPY(3*roctree%NVT, roctree%p_Ddata, 1, p_Ddata, 1)

  end subroutine otree_copyFromOctree_array

  !************************************************************************

!<subroutine>
  
  subroutine otree_copyToOctree_handle(h_Ddata, roctree)

!<description>
    ! This subroutine copies the content of a handle to the octree.
!</description>

!<input>
    ! Handle to the coordinate vector
    integer, intent(in) :: h_Ddata
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata
    
    ! Set pointer
    call storage_getbase_double2D(h_Ddata, p_Ddata)

    ! Call copy routine
    call otree_copyToOctree(p_Ddata, roctree)
    
  end subroutine otree_copyToOctree_handle

  !************************************************************************

!<subroutine>
  
  subroutine otree_copyToOctree_array(p_Ddata, roctree)

!<description>
    ! This subroutine copies the content of an array to the octree.
!</description>

!<input>
    ! Coordinate vector
    real(DP), dimension(:,:), intent(in) :: p_Ddata
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    
    if (size(p_Ddata,1) .ne. 3) then
      call output_line('First dimension of array must be 3!',&
                       OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyToOctree_array')
      call sys_halt()
    end if

    ! Adjust dimension of data array
    roctree%NVT = size(p_Ddata,2)
    if (roctree%NNVT .lt. roctree%NVT) then
      call resizeNVT(roctree, roctree%NVT)
    end if
    
    ! Copy data array
    call DCOPY(3*roctree%NVT, p_Ddata, 1, roctree%p_Ddata, 1)

    ! Rebuild structure
    call otree_rebuildOctree(roctree)

  end subroutine otree_copyToOctree_array

  !************************************************************************
  
!<function>
  
  function otree_insertIntoOctree(roctree, Ddata, ivt, inode) result(iresult)

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
    real(DP), dimension(3), intent(in) :: Ddata

    ! OPTIONAL: Number of the node to which vertex should be inserted.
    ! If there is no space left, then the next free position will be used
    integer, intent(in), optional :: inode
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>

!<output>
    ! Number of the inserted vertex
    integer, intent(out) :: ivt
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
      iresult = otree_searchInOctree(roctree, Ddata, jnode, jpos, jvt)
      if (iresult .eq. OTREE_FOUND) then
        ivt = jvt
        return
      end if
    end if

    ! Check if there is enough space left in the vertex component of the octree
    if (roctree%NVT .eq. roctree%NNVT) then
      isize = max(roctree%NNVT+1, ceiling(roctree%dfactor*roctree%NNVT))
      call resizeNVT(roctree, isize)
    end if
    
    ! Update values
    roctree%NVT = roctree%NVT+1
    ivt = roctree%NVT
    roctree%p_Ddata(:,ivt) = Ddata

    ! Insert entry recursively
    iresult = insert(ivt, jnode)

    ! Check success
    if (iresult .eq. OTREE_FAILED) then
      roctree%NVT = roctree%NVT-1
      roctree%p_Ddata(:,ivt) = 0.0_DP
      ivt = 0
    end if
    
  contains
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    recursive function insert(ivt, inode) result(iresult)

      integer, intent(in) :: ivt,inode
      integer :: iresult

      ! local variables
      real(DP) :: xmin,ymin,zmin,xmax,ymax,zmax,xmid,ymid,zmid
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

          ! Decrease starting position of first nodes by one
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
        roctree%NNODE = nnode+OTREE_MAX

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

        ! Add the data values from INODE to the eight new nodes
        do i = 1, roctree%ndata
          jvt = roctree%p_Knode(i, inode)
          jnode = nnode+otree_getDirection(roctree, roctree%p_Ddata(:,jvt), inode)
          jresult = insert(jvt, jnode)

          if (jresult .eq. OTREE_FAILED) then
            call output_line('Internal error in insertion!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'otree_insertIntoOctree')
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
        jnode = nnode+otree_getDirection(roctree, roctree%p_Ddata(:,ivt),inode)
        iresult = insert(ivt, jnode)
        
      elseif (roctree%p_Knode(OTREE_STATUS, inode) .ge. OTREE_EMPTY) then

        ! There is still some space in the node
        roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)+1
        roctree%p_Knode(roctree%p_Knode(OTREE_STATUS, inode), inode) = ivt

        ! Set success flag
        iresult = OTREE_INSERTED

      elseif (roctree%p_Knode(OTREE_STATUS, inode) .eq. OTREE_SUBDIV) then

        ! Proceed to correcponding sub-tree
        jnode = -roctree%p_Knode(otree_getDirection(roctree, Ddata, inode), inode)
        iresult = insert(ivt, jnode)

      else
        
        ! Set failure flag
        iresult = OTREE_FAILED

      end if

    end function insert

  end function otree_insertIntoOctree
  
  ! ***************************************************************************
  
!<function>
  
  function otree_deleteFromOctree(roctree, Ddata, ivt) result(iresult)

!<description>
    ! This function deletes an item from the octree.
    ! The value IVT returns the number of the item which is
    ! moved to the position of the deleted item. If the deleted
    ! item was the last one, then IVT=NVT is returned.
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    real(DP), dimension(3), intent(in) :: Ddata
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that is deleted
    integer, intent(out) :: ivt
!</output>

!<result>
    ! Result of the deletion:
    !   QTREE_FAILED
    !   QTREE_DELETED
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
      real(DP), dimension(3) :: DdataTmp
      integer, dimension(OTREE_MAX) :: Knode
      integer :: jvt,i,jnode,ipos,jpos,nemptyChildren
      
      
      ! Check status of current node
      select case(roctree%p_Knode(OTREE_STATUS, inode))

      case (OTREE_SUBDIV)   ! Node is subdivided
        
        ! Compute child INODE which to look recursively.
        jnode = -roctree%p_Knode(otree_getDirection(roctree, Ddata, inode), inode)
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
                         OU_CLASS_ERROR,OU_MODE_STD,'otree_deleteFromOctree')
        call sys_halt()


      case default   ! Node is a non-empty leaf

        ! Search for (x,y,z) in current node
        do ipos = 1, roctree%p_Knode(OTREE_STATUS, inode)
          
          ! Get vertex number
          ivt = roctree%p_Knode(ipos, inode)
          
          if (maxval(abs(roctree%p_Ddata(:, ivt)-Ddata)) .le. SYS_EPSREAL_DP) then
          
            ! Physically remove the item IVT from node INODE
            jpos = roctree%p_Knode(OTREE_STATUS, inode)
            
            roctree%p_Knode(ipos, inode) = roctree%p_Knode(jpos, inode)
            roctree%p_Knode(jpos, inode) = 0
            roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)-1
            
            ! If IVT is not last item then find item with largest number JVT
            if (ivt .ne. roctree%NVT) then
              DdataTmp(:) = roctree%p_Ddata(:,roctree%NVT)
              if (otree_searchInOctree(roctree, DdataTmp(:),&
                                       jnode, jpos, jvt) .eq. OTREE_FOUND) then
                
                ! Move last item JVT to position IVT
                roctree%p_Ddata(:, ivt) = roctree%p_Ddata(:, roctree%NVT)
                roctree%p_Knode(jpos, jnode) = ivt
              else
                call output_line('Internal error in deletion!',&
                                 OU_CLASS_ERROR,OU_MODE_STD,'otree_deleteFromOctree')
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

  end function otree_deleteFromOctree

  ! ***************************************************************************

!<function>

  function otree_deleteFromOctreeByNumber(roctree, ivt, ivtReplace) result(iresult)

!<description>
    ! This function deletes vertex with number IVT from the octree.
!</description>

!<input>
    ! Number of the vertex to be deleted
    integer, intent(in) :: ivt
!</input>

!<inputoutput>
    ! Octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that replaces the deleted vertex
    integer, intent(out) :: ivtReplace
!</output>

!<result>
    ! Result of the deletion:
    !   QTREE_FAILED
    !   QTREE_DELETED
    integer :: iresult
!</result>
!</function>

    ! local variables
    real(DP), dimension(3) :: Ddata
    
    if (ivt .le. roctree%NVT) then
      ! Get coordinates and invoke deletion routine
      Ddata   = roctree%p_Ddata(:,ivt)
      iresult = otree_deleteFromOctree(roctree, Ddata, ivtReplace)
    else
      iresult = OTREE_FAILED
    end if

  end function otree_deleteFromOctreeByNumber

  ! ***************************************************************************

!<function>
  
  function otree_searchInOctree(roctree, Ddata, inode, ipos, ivt) result(iresult)

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
    type(t_octree), intent(in) :: roctree
    
    ! Coordinates that should be searched for
    real(DP), dimension(3), intent(in) :: Ddata
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
    !   QTREE_NOT_FOUND
    !   QTREE_FOUND
    integer :: iresult
!</result>
!</function>
    
    ! Initialize
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
        inode = -roctree%p_Knode(otree_getDirection(roctree, Ddata, inode), inode)
        iresult =  search(inode, ipos, ivt)

        
      case (OTREE_EMPTY)   ! Node is empty so it cannot contain the item
        
        iresult = OTREE_NOT_FOUND

        
      case (OTREE_DEL)   ! Node is deleted -> serious error
        
        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'otree_searchInOctree')
        call sys_halt()


      case default   ! Node is a non-empty leaf
        
        ! Search for (x,y,z) in current node
        do ipos = 1, roctree%p_Knode(OTREE_STATUS, inode)
          
          ! Get vertex number
          ivt = roctree%p_Knode(ipos, inode)
          
          if (maxval(abs(roctree%p_Ddata(:, ivt)-Ddata)) .le. SYS_EPSREAL_DP) then
            
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

  end function otree_searchInOctree

  !************************************************************************
  
!<function>
  
  pure function otree_getDirection(roctree, Ddata, inode) result(idirection)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Ddata.
!</description>

!<input>
    ! Octree
    type(t_octree), intent(in) :: roctree
    
    ! Coordinates
    real(DP), dimension(3), intent(in) :: Ddata

    ! Number of node
    integer, intent(in) :: inode
!</input>

!<result>
    ! Further search direction
    integer :: idirection
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
          idirection = OTREE_NEF
        else
          idirection = OTREE_SEF
        end if
      else
        if (Ddata(2) > ymid) then
          idirection = OTREE_NWF
        else
          idirection = OTREE_SWF
        end if
      end if

    else

      if (Ddata(1) > xmid) then
        if (Ddata(2) > ymid) then
          idirection = OTREE_NEB
        else
          idirection = OTREE_SEB
        end if
      else
        if (Ddata(2) > ymid) then
          idirection = OTREE_NWB
        else
          idirection = OTREE_SWB
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
    type(t_octree), intent(in) :: roctree

    ! filename of the output file
    character(LEN=*), intent(in) :: cfilename
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
      real(DP), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax
      integer, intent(in) :: inode
      real(DP) :: xmid,ymid,zmid
      integer :: i,j

      write(UNIT=iunit,FMT=*) 'hex'
      write(UNIT=iunit,FMT=10) xmin, ymin, zmin, xmax, ymax, zmax, inode
      
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
          write(UNIT=iunit, FMT=20) roctree%p_Ddata(1, j), roctree%p_Ddata(2, j),&
                                    roctree%p_Ddata(3, j), j
        end do
        
      end if
      
10    format(6(E15.6E3,1X),I10)
20    format(3(E15.6E3,1X),I10)

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
    type(t_octree), intent(in) :: roctree
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
    call output_line('h_Ddata: '//trim(sys_siL(roctree%h_Ddata,15)))
    call output_line('h_Dbbox: '//trim(sys_siL(roctree%h_Dbbox,15)))
    call output_line('h_Knode: '//trim(sys_siL(roctree%h_Knode,15)))
    call output_lbrk()
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
    type(t_octree), intent(in) :: roctree
!</input>

!<result>
    ! Number of vertices in octree
    integer :: nvt
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
    type(t_octree), intent(in) :: roctree

    ! OPTIONAL: number of node for which bounding box should be returned
    integer, intent(in), optional :: inode
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
    type(t_octree), intent(in) :: roctree

    ! position in the octree
    integer, intent(in) :: ivt
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
    type(t_octree), intent(in) :: roctree

    ! position in the octree
    integer, intent(in) :: ivt
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
    type(t_octree), intent(in) :: roctree

    ! position in the octree
    integer, intent(in) :: ivt
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
    type(t_octree), intent(in) :: roctree
!</input>

!<inputoutput>
    ! Destination octree
    type(t_octree), intent(inout) :: roctreeBackup
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
    type(t_octree), intent(in) :: roctreeBackup
!</input>

!<inputoutput>
    ! Destination octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    ! Release octree
    call otree_releaseOctree(roctree)

    ! Duplicate the backup
    call otree_duplicateOctree(roctreeBackup, roctree)

  end subroutine otree_restoreOctree

  !************************************************************************

!<subroutine>

  subroutine otree_rebuildOctree(roctree)

!<description>
    ! This subroutine rebuilds the structure of an octree
!</description>

!<inputoutput>
    ! Octree
    type(t_octree), intent(inout) :: roctree
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Inodes
    integer, dimension(2) :: Isize
    integer :: h_Inodes

    ! Initialization
    roctree%NNODE = 1
    roctree%NFREE = 0

    call storage_getsize(roctree%h_Knode, Isize); roctree%NNNODE = Isize(2)
    call storage_getsize(roctree%h_Ddata, Isize); roctree%NNVT   = Isize(2)

    ! Check if octree contains data
    if (roctree%NVT .eq. 0) then

      ! Initialize first node
      roctree%nnode = 1
      roctree%p_Knode(OTREE_STATUS,1) = OTREE_EMPTY
      roctree%p_Knode(OTREE_PARENT,1) = OTREE_EMPTY
      roctree%p_Knode(OTREE_PARPOS,1) = OTREE_EMPTY
      roctree%p_Knode(1:roctree%NDATA, 1) = 0
      
    else
      
      ! Create temporary storage
      h_Inodes = ST_NOHANDLE
      call storage_new('otree_rebuildOctree', 'Inodes', roctree%NVT,&
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
      real(DP) :: xmin,ymin,zmin,xmax,ymax,zmax,xmid,ymid,zmid
      integer :: i,isize,jnode,nnode,ivt,imid1,imid2,imid3,imid4,imid5,imid6,imid7

      
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
        xmin = roctree%p_Dbbox(OTREE_XMIN,inode)
        ymin = roctree%p_Dbbox(OTREE_YMIN,inode)
        zmin = roctree%p_Dbbox(OTREE_ZMIN,inode)
        xmax = roctree%p_Dbbox(OTREE_XMAX,inode)
        ymax = roctree%p_Dbbox(OTREE_YMAX,inode)
        zmax = roctree%p_Dbbox(OTREE_ZMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)
        
        ! NWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWF) = OTREE_NWF
        roctree%p_Dbbox(:,nnode+OTREE_NWF)            = (/xmin,ymid,zmin,xmid,ymax,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_NWF) = 0
        
        ! SWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWF) = OTREE_SWF
        roctree%p_Dbbox(:,nnode+OTREE_SWF)            = (/xmin,ymin,zmin,xmid,ymid,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SWF) = 0
        
        ! SEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEF) = OTREE_SEF
        roctree%p_Dbbox(:,nnode+OTREE_SEF)            = (/xmid,ymin,zmin,xmax,ymid,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SEF) = 0
        
        ! NEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEF) = OTREE_NEF
        roctree%p_Dbbox(:,nnode+OTREE_NEF)            = (/xmid,ymid,zmin,xmax,ymax,zmid/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_NEF) = 0
        
        ! NWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWB) = OTREE_NWB
        roctree%p_Dbbox(:,nnode+OTREE_NWB)            = (/xmin,ymid,zmid,xmid,ymax,zmax/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_NWB) = 0
        
        ! SWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWB) = OTREE_SWB
        roctree%p_Dbbox(:,nnode+OTREE_SWB)            = (/xmin,ymin,zmid,xmid,ymid,zmax/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SWB) = 0
        
        ! SEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEB) = OTREE_SEB
        roctree%p_Dbbox(:,nnode+OTREE_SEB)            = (/xmid,ymin,zmid,xmax,ymid,zmax/)
        roctree%p_Knode(1:roctree%NDATA,nnode+OTREE_SEB) = 0
        
        ! NEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEB) = OTREE_NEB
        roctree%p_Dbbox(:,nnode+OTREE_NEB)            = (/xmid,ymid,zmid,xmax,ymax,zmax/)
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

      real(DP), intent(in) :: dmid
      integer, intent(in) :: istart, iend, idim
      integer :: imid

      ! local variables
      integer :: i

      
      if (istart .gt. iend) then
        
        ! If istart > iend then there are no items in this partition
        imid = iend+1
        
      else if (istart .eq. iend) then
        
        ! If istart = iend then there is one item in this partition
        imid = merge(iend, iend+1, roctree%p_Ddata(idim, p_Inodes(iend)) .gt. dmid)

      else
        
        ! Otherwise, sort items in partition
        call quicksort(istart, iend, idim)
        
        ! Find location of first item belonging to "IS GREATER" partition
        do i = istart, iend
          if (roctree%p_Ddata(idim, p_Inodes(i)) .gt. dmid) then
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
      real(DP), dimension(2) :: Daux
      real(DP) :: dpivot
      integer :: i, j, iaux
      
      ! Initialization
      i = istart
      j = iend-1
      dpivot = roctree%p_Ddata(idim, p_Inodes(iend))

      do
        do while((i .lt. iend) .and. (roctree%p_Ddata(idim, p_Inodes(i)) .le. dpivot))
          i = i+1
        end do
        
        do while((j .gt. istart) .and. (roctree%p_Ddata(idim, p_Inodes(j)) .ge. dpivot))
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

  end subroutine otree_rebuildOctree

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
    integer, intent(in) :: nnnode
!</input>

!<inputoutput>
    ! Octree that should be resized
    type(t_octree) :: roctree
!</inputoutput>
!</subroutine>

    call storage_realloc('resizeNNODE', nnnode, roctree%h_Dbbox,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_realloc('resizeNNODE', 1, nnnode, roctree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_double2D(roctree%h_Dbbox, roctree%p_Dbbox)
    call storage_getbase_int2D(roctree%h_Knode,    roctree%p_Knode)
    
    roctree%NNNODE  = nnnode
    roctree%NRESIZE =roctree%NRESIZE+1

  end subroutine resizeNNODE

end module octree


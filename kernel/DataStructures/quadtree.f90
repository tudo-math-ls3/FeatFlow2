!##############################################################################
!# ****************************************************************************
!# <name> quadtree </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a (linear) quadtree.
!#
!# The following routines are available:
!#
!# 1.) qtree_createQuadtree
!#     -> Create a new quadtree structure
!#
!# 2.) qtree_releaseQuadtree
!#     -> Release an existing quadtree
!#
!# 3.) qtree_resizeQuadtree
!#     -> Resize an existing quadtree
!#
!# 4.) qtree_copyToQuadtree
!#     -> Copy data to a quadtree
!#
!# 5.) qtree_copyFromQuadtree
!#     -> Copy data from the quadtree
!#
!# 6.) qtree_insertIntoQuadtree
!#     -> Insert data into quadtree
!#
!# 7.) qtree_deleteFromQuadtree
!#     -> Delete data from quadtree
!#
!# 8.) qtree_searchInQuadtree
!#     -> Search data in quadtree
!#
!# 9.) qtree_getDirection
!#     -> Get direction for the next node
!#
!# 10.) qtree_printQuadtree
!#      -> Write quadtree to file
!#
!# 11.) qtree_infoQuadtree
!#      -> Output info about quadtree
!#
!# 12.) qtree_getsize
!#      -> Return number of vertices in quadtree
!#
!# 13.) qtree_getBoundingBox
!#      -> Return the outer bounding box
!#
!# 14.) qtree_getX
!#      -> Return the X-value at a given position
!#
!# 15.) qtree_getY
!#      -> Return the Y-value at a given position
!#
!# 16.) qtree_duplicateQuadtree
!#      -> Create a duplicate / backup of a quadtree
!#
!# 17.) qtree_restoreQuadtree
!#      -> Restore a quadtree from a previous backup
!#
!# For the internal use the following routines are available:
!#
!# 1.) resizeNVT
!#     -> Resize the number of stored values
!#
!# 2.) resizeNNODE
!#     -> Resize the number of stored quads
!#
!# </purpose>
!##############################################################################

module quadtree
  use fsystem
  use storage
  use genoutput
  implicit none
  
  private
  public :: t_quadtree
  public :: qtree_createQuadtree
  public :: qtree_releaseQuadtree
  public :: qtree_copyToQuadtree
  public :: qtree_copyFromQuadtree
  public :: qtree_insertIntoQuadtree
  public :: qtree_deleteFromQuadtree
  public :: qtree_searchInQuadtree
  public :: qtree_getDirection
  public :: qtree_printQuadtree
  public :: qtree_infoQuadtree
  public :: qtree_getsize
  public :: qtree_getBoundingBox
  public :: qtree_getX
  public :: qtree_getY
  public :: qtree_duplicateQuadtree
  public :: qtree_restoreQuadtree

!<constants>

!<constantblock description="Constants for quadtree structure">
  
  ! Position of the status information
  integer, parameter :: QTREE_STATUS = 7
  
  ! Position of the parent information
  integer, parameter :: QTREE_PARENT = 6

  ! Position of the "position" of the parent information
  integer, parameter :: QTREE_PARPOS = 5

  ! Number of free positions in quad
  integer, parameter :: QTREE_FREE   = 5

  ! Item in "North-East" position
  integer, parameter :: QTREE_NE     = 4

  ! Item in "South-East" position
  integer, parameter :: QTREE_SE     = 3

  ! Item in "South-West" position
  integer, parameter :: QTREE_SW     = 2

  ! Item in "North-West" position
  integer, parameter :: QTREE_NW     = 1

  ! Identifier: Quad is empty
  integer, parameter :: QTREE_EMPTY  =  0

  ! Identifier: Status is subdivided
  integer, parameter :: QTREE_SUBDIV = -1
  
  ! Maximum number of items for each quad
  integer, parameter :: QTREE_MAX    =  4

!</constantblock> 

!<constantblock description="Constants for quadtree bounding-box">

  ! Position of the x-minimal value
  integer, parameter :: QTREE_XMIN   =  1

  ! Position of the y-minimal value
  integer, parameter :: QTREE_YMIN   =  2

  ! Position of the x-maximal value
  integer, parameter :: QTREE_XMAX   =  3

  ! Position of the y-maximal value
  integer, parameter :: QTREE_YMAX   =  4

!</constantblock>   
  
!<constantblock description="Constants for quadtree operations">

  ! Item could not be found in the quadtree
  integer, parameter, public :: QTREE_NOT_FOUND = -1

  ! Item could be found in the quadtree
  integer, parameter, public :: QTREE_FOUND     =  0
!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! A linear quadtree implemented as array

  type t_quadtree

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

    ! Handle to quadtree structure
    !   KNODE(1:7,1:NNNODE)
    !   KNODE(  7,INODE)    : < 0, the quad has been subdivided
    !                         = 0, the quad is empty
    !                         > 0, the number of points stored in the quad
    !   KNODE(  6,INODE)    : > 0, the quad the present quad came from
    !   KNODE(  5,INODE)    : > 0, the position in the quad the present 
    !                           quad came from
    !   KNODE(1:4,INODE)    : for KNODE(7,INODE) > 0 : the points stored in the quad
    !                         for KNODE(7,INODE) < 0 : the quads into
    !                                                  which the present quad was subdivided
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

    ! Quadtree structure
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:,:), pointer :: p_Knode

    ! Number of vertices currently stored in the quadtree
    integer :: NVT    = 0

    ! Total number of vertices that can be stored  in the quadtree
    integer :: NNVT   = 0

    ! Number of subdivisions currently store in the quadtree
    integer :: NNODE  = 0

    ! Total number of subdivision that can be stored in the quadtree
    integer :: NNNODE = 0

    ! Total number of resize operations
    integer :: NRESIZE            = 0

    ! Factor by which the quadtree is enlarged if new storage is allocate
    real(DP) :: dfactor           = 1.5_DP
  end type t_quadtree
  
!</typeblock>
!</types>
  
  interface qtree_copyToQuadtree
    module procedure qtree_copyToQuadtree_handle
    module procedure qtree_copyToQuadtree_array
  end interface

  interface qtree_copyFromQuadtree
    module procedure qtree_copyFromQuadtree_handle
    module procedure qtree_copyFromQuadtree_array
  end interface
   
  interface qtree_deleteFromQuadtree
    module procedure qtree_deleteFromQuadtree
    module procedure qtree_deleteFromQtreeByNumber
  end interface
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine qtree_createQuadtree(rquadtree, nnvt, nnnode,&
                                  xmin, ymin, xmax, ymax, dfactor)
  
!<description>
    ! This subroutine creates a new quadtree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the quadtree
    integer, intent(IN) :: nnvt

    ! Total number of subdivisions that should be stored in the quadtree
    integer, intent(IN) :: nnnode

    ! Dimensions of the initial bounding box
    real(DP), intent(IN) :: xmin,ymin,xmax,ymax

    ! OPTIONAL: Factor by which the quadtree should be enlarged if
    ! new storage has to be allocated
    real(DP), intent(IN), optional :: dfactor
!</input>

!<output>
    ! Quadtree structure
    type(t_quadtree), intent(OUT) :: rquadtree
!</output>
!</subroutine>
    
    integer(I32), dimension(2) :: Isize
    
    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1_DP) rquadtree%dfactor=dfactor
    end if
    
    ! Set values
    rquadtree%NNNODE = nnnode
    rquadtree%NNVT   = nnvt
    rquadtree%NNODE  = 0
    rquadtree%NVT    = 0
    rquadtree%NRESIZE= 0
    
    ! Allocate memory and associate pointers
    Isize = (/2,nnvt/)
    call storage_new('qtree_createQuadtree', 'p_Ddata',&
                     Isize, ST_DOUBLE, rquadtree%h_Ddata, ST_NEWBLOCK_ZERO)
    Isize = (/4,nnnode/)
    call storage_new('qtree_createQuadtree', 'p_Dbbox',&
                     Isize, ST_DOUBLE, rquadtree%h_Dbbox, ST_NEWBLOCK_ZERO)
    Isize = (/7,nnnode/)
    call storage_new('qtree_createQuadtree', 'p_Knode',&
                     Isize, ST_INT, rquadtree%h_Knode, ST_NEWBLOCK_ZERO)
    call storage_getbase_double2D(rquadtree%h_Ddata, rquadtree%p_Ddata)
    call storage_getbase_double2D(rquadtree%h_Dbbox, rquadtree%p_Dbbox)
    call storage_getbase_int2D(rquadtree%h_Knode,    rquadtree%p_Knode)

    ! Initialize first quad
    rquadtree%nnode = 1
    rquadtree%p_Knode(QTREE_STATUS,1) = QTREE_EMPTY
    rquadtree%p_Knode(QTREE_PARENT,1) = QTREE_EMPTY
    rquadtree%p_Knode(QTREE_FREE,1)   = 1
    
    ! Set co-ordinates of the bounding box
    rquadtree%p_Dbbox(QTREE_XMIN,1) = xmin
    rquadtree%p_Dbbox(QTREE_YMIN,1) = ymin
    rquadtree%p_Dbbox(QTREE_XMAX,1) = xmax
    rquadtree%p_Dbbox(QTREE_YMAX,1) = ymax
  end subroutine qtree_createQuadtree

  !************************************************************************
  
!<subroutine>
  
  subroutine qtree_releaseQuadtree(rquadtree)

!<description>
    ! This subroutine releases an existing quadtree
!</description>

!<inputoutput>
    ! quadtree
    type(t_quadtree), intent(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! Free memory
    if (rquadtree%h_Ddata .ne. ST_NOHANDLE) call storage_free(rquadtree%h_Ddata)
    if (rquadtree%h_Dbbox .ne. ST_NOHANDLE) call storage_free(rquadtree%h_Dbbox)
    if (rquadtree%h_Knode .ne. ST_NOHANDLE) call storage_free(rquadtree%h_Knode)
    nullify(rquadtree%p_Knode, rquadtree%p_Dbbox, rquadtree%p_Ddata)

    ! Reset values
    rquadtree%NNNODE  = 0
    rquadtree%NNVT    = 0
    rquadtree%NNODE   = 0
    rquadtree%NVT     = 0
    rquadtree%NRESIZE = 0
  end subroutine qtree_releaseQuadtree

  !************************************************************************

!<subroutine>
  
  subroutine qtree_copyFromQuadtree_handle(rquadtree, h_Ddata)

!<description>
    ! This subroutine copies the content of the quadtree to a handle.
    ! If the handle is not associated, then a new handle with correct
    ! size is allocated. Otherwise, the handle is reallocated if 
    ! it does not provide enough memory.
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree
!</input>

!<inputoutput>
    ! Handle to the coordinate vector
    integer, intent(INOUT) :: h_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata,p_DdataTmp
    integer, dimension(2) :: Isize
    
    ! Check if handle is associated
    if (h_Ddata .eq. ST_NOHANDLE) then
      Isize = (/2, rquadtree%NVT/)
      call storage_new('qtree_copyFromQuadtree_handle', 'p_Ddata',&
                       Isize, ST_DOUBLE, h_Ddata, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_Ddata, Isize)
      if (Isize(2) < rquadtree%NVT) then
        call storage_realloc('qtree_copyFromQuadtree_handle',&
                             rquadtree%NVT, h_Ddata, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Set pointers
    call storage_getbase_double2D(h_Ddata, p_Ddata)
    call storage_getbase_double2D(rquadtree%h_Ddata, p_DdataTmp)

    ! Copy data
    call DCOPY(2*rquadtree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  end subroutine qtree_copyFromQuadtree_handle

  !************************************************************************

!<subroutine>

  subroutine qtree_copyFromQuadtree_array(rquadtree, p_Ddata)

!<description>
    ! This subroutine copies the content of the quadtree to an array.
    ! If the array is too small, then this subroutines terminates.
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree
!</input>

!<inputoutput>
    ! Coordinate vector
    real(DP), dimension(:,:), intent(INOUT) :: p_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DdataTmp
    integer, dimension(2) :: Isize

    ! Check size of array
    Isize = shape(p_Ddata)
    if (Isize(1) .ne. 2 .or. Isize(2) < rquadtree%NVT) then
      call output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyFromQuadtree_array')
      call sys_halt()
    end if

    ! Set pointers
    call storage_getbase_double2D(rquadtree%h_Ddata, p_DdataTmp)

    ! Copy data
    call DCOPY(2*rquadtree%NVT, p_DdataTmp, 1, p_Ddata, 1)
    
  end subroutine qtree_copyFromQuadtree_array

  !************************************************************************

!<subroutine>
  
  subroutine qtree_copyToQuadtree_handle(h_Ddata, rquadtree)

!<description>
    ! This subroutine copies the content of a handle to the quadtree.
!</description>

!<input>
    ! Handle to the coordinate vector
    integer, intent(IN) :: h_Ddata
!</input>

!<inputoutput>
    ! quadtree
    type(t_quadtree), intent(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata
    integer :: ivt,jvt,inode,ipos
    
    ! Set pointer and check its shape
    call storage_getbase_double2D(h_Ddata, p_Ddata)
    if (size(p_Ddata,1) .ne. 2) then
      call output_line('First dimension of array must be 2!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyToQuadtree_handle')
      call sys_halt()
    end if
    
    do ivt = 1, size(p_Ddata, 2)
      if (qtree_searchInQuadtree(rquadtree, p_Ddata(:,ivt),&
                                 inode, ipos, jvt) .eq. QTREE_NOT_FOUND)&
          call qtree_insertIntoQuadtree(rquadtree, ivt, p_Ddata(:,ivt), inode)
    end do
  end subroutine qtree_copyToQuadtree_handle

  !************************************************************************

!<subroutine>
  
  subroutine qtree_copyToQuadtree_array(p_Ddata, rquadtree)

!<description>
    ! This subroutine copies the content of an array to the quadtree.
!</description>

!<input>
    ! Coordinate vector
    real(DP), dimension(:,:), intent(IN) :: p_Ddata
!</input>

!<inputoutput>
    ! quadtree
    type(t_quadtree), intent(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: ivt,jvt,inode,ipos
    
    if (size(p_Ddata,1) .ne. 2) then
      call output_line('First dimension of array must be 2!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyToQuadtree_handle')
      call sys_halt()
    end if

    do ivt = 1, size(p_Ddata, 2)
      if (qtree_searchInQuadtree(rquadtree, p_Ddata(:,ivt),&
                                 inode, ipos, jvt) .eq. QTREE_NOT_FOUND)&
          call qtree_insertIntoQuadtree(rquadtree, ivt, p_Ddata(:,ivt), inode)
    end do
  end subroutine qtree_copyToQuadtree_array

  !************************************************************************
  
!<subroutine>
  
  subroutine qtree_insertIntoQuadtree(rquadtree, ivt, Ddata, inode)

!<description>
    ! This subroutine inserts a new coordinate item to the quadtree
!</description>

!<input>
    ! Number of the inserted vertex
    integer, intent(IN) :: ivt

    ! Number of the quad to which vertex is inserted
    integer, intent(IN) :: inode
    
    ! Coordinates of the new vertex
    real(DP), dimension(2), intent(IN) :: Ddata
!</input>

!<inputoutput>
    ! quadtree
    type(t_quadtree), intent(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: isize
    
    ! Check if there is enough space left in the nodal component of the quadtree
    if (rquadtree%NVT .eq. rquadtree%NNVT) then
      isize = max(rquadtree%NNVT+1, ceiling(rquadtree%dfactor*rquadtree%NNVT))
      call resizeNVT(rquadtree, isize)
    end if
    
    ! Update values and add new entry recursively
    rquadtree%NVT            = rquadtree%NVT+1
    rquadtree%p_Ddata(:,ivt) = Ddata
    call insert(ivt, inode)
    
  contains  
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    recursive subroutine insert(ivt,inode)
      integer, intent(IN) :: ivt,inode
      real(DP) :: xmin,ymin,xmax,ymax,xmid,ymid
      integer :: i,jnode,nnode,knode,isize
            
      if (rquadtree%p_Knode(QTREE_STATUS,inode) .eq. QTREE_MAX) then
        
        if (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode) then
          isize = max(rquadtree%nnode+QTREE_MAX, ceiling(rquadtree%dfactor*rquadtree%NNNODE))
          call resizeNNODE(rquadtree, isize)
        end if
        
        ! Quad is full and needs to be refined into four new quads
        xmin = rquadtree%p_Dbbox(QTREE_XMIN,inode)
        ymin = rquadtree%p_Dbbox(QTREE_YMIN,inode)
        xmax = rquadtree%p_Dbbox(QTREE_XMAX,inode)
        ymax = rquadtree%p_Dbbox(QTREE_YMAX,inode)
        xmid = (xmin+xmax)/2._DP
        ymid = (ymin+ymax)/2._DP
        
        ! Store the number of quads
        nnode           = rquadtree%NNODE
        rquadtree%NNODE = nnode+QTREE_MAX
        
        ! NW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NW) = QTREE_NW
        rquadtree%p_Dbbox(:,nnode+QTREE_NW)            = (/xmin,ymid,xmid,ymax/)
        
        ! SW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_Dbbox(:,nnode+QTREE_SW)            = (/xmin,ymin,xmid,ymid/)
        
        ! SE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_Dbbox(:,nnode+QTREE_SE)            = (/xmid,ymin,xmax,ymid/)
        
        ! NE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_Dbbox(:,nnode+QTREE_NE)            = (/xmid,ymid,xmax,ymax/)
        
        ! Add the four values from INODE to the four new quads 
        ! NNODE+1:NNODE+4 recursively
        do i = 1, QTREE_MAX
          knode = rquadtree%p_Knode(i, inode)
          jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,knode), inode)
          call insert(knode, jnode)
        end do
        
        ! Mark the current quad as subdivided and set pointers to its four children 
        rquadtree%p_Knode(QTREE_STATUS,inode) = QTREE_SUBDIV
        rquadtree%p_Knode(1:QTREE_MAX,inode)  = (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                  nnode+QTREE_SE, nnode+QTREE_NE/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,ivt), inode)
        call insert(ivt, jnode)
        
      else
        
        ! Quad is not full, so new items can be stored
        rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)+1
        rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = ivt
      end if
    end subroutine insert
  end subroutine qtree_insertIntoQuadtree
  
  ! ***************************************************************************
  
!<function>
  
  function qtree_deleteFromQuadtree(rquadtree, Ddata, ivt) result(f)

!<description>
    ! This function deletes an item from the quadtree
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    real(DP), dimension(2), intent(IN) :: Ddata
!</input>

!<inputoutput>
    ! quadtree
    type(t_quadtree), intent(INOUT) :: rquadtree
!</inputoutput>

!<output>
    ! Number of the vertex that is deleted
    integer, intent(OUT) :: ivt
!</output>

!<result>
    ! Result of the deletion: QTREE_NOT_FOUND / QTREE_FOUND
    integer :: f
!</result>
!</function>
    
    ! local variables
    integer :: inode,ipos,jpos,jvt
    real(DP), dimension(2) :: DdataTmp
    
    ! Search for the given coordinates
    f = qtree_searchInQuadtree(rquadtree, Ddata, inode, ipos, ivt)
    
    ! What can we do from the searching
    if (f .eq. QTREE_FOUND) then
      
      ! Remove item IVT from node INODE
      do jpos = ipos+1, rquadtree%p_Knode(QTREE_STATUS, inode)
        rquadtree%p_Knode(jpos-1, inode) = rquadtree%p_Knode(jpos, inode)
      end do
      rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = 0
      rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)-1
      
      ! If IVT is not last item move last item to position IVT
      if (ivt .ne. rquadtree%NVT) then
        DdataTmp(1:2) = rquadtree%p_Ddata(1:2,rquadtree%NVT)
        if (qtree_searchInQuadtree(rquadtree, DdataTmp(:),&
                                   inode, ipos, jvt) .eq. QTREE_FOUND) then
          rquadtree%p_Ddata(:, ivt)      = rquadtree%p_Ddata(:, rquadtree%NVT)
          rquadtree%p_Knode(ipos, inode) = ivt
        end if

        ! Set number of removed vertex
        ivt = rquadtree%NVT
      end if
      
      ! Decrease number of vertices
      rquadtree%NVT = rquadtree%NVT-1
    end if
  end function qtree_deleteFromQuadtree

  ! ***************************************************************************

!<function>

  function qtree_deleteFromQtreeByNumber(rquadtree, ivt, ivtReplace) result(f)

!<description>
    ! This function deletes vertex with number IVT from the quadtree
!</description>

!<input>
    ! Number of the vertex to be deleted
    integer, intent(IN) :: ivt
!</input>

!<inputoutput>
    ! quadtree
    type(t_quadtree), intent(INOUT) :: rquadtree
!</inputoutput>

!<output>
    ! Number of the vertex that replaces the deleted vertex
    integer, intent(OUT) :: ivtReplace
!</output>

!<result>
    ! Result of the deletion: QTREE_NOT_FOUND / QTREE_FOUND
    integer :: f
!</result>
!</function>

    ! local variables
    real(DP), dimension(2) :: Ddata
    
    if (ivt .le. rquadtree%NVT) then
      ! Get coordinates and invoke deletion routine
      Ddata = rquadtree%p_Ddata(:,ivt)
      f     = qtree_deleteFromQuadtree(rquadtree, Ddata, ivtReplace)
    else
      call output_line('Invalid vertex number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_deleteFromQtreeByNumber')
      call sys_halt()
    end if
  end function qtree_deleteFromQtreeByNumber
 
  ! ***************************************************************************

!<function>
  
  function qtree_searchInQuadtree(rquadtree, Ddata, inode, ipos, ivt) result(f)

!<description>
    ! This subroutine searches for given coordinates in the quadtree
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree
    
    ! coordinates that should be searched
    real(DP), dimension(2), intent(IN) :: Ddata
!</input>

!<output>
    ! Number of the quad in which the given coordinates are
    integer, intent(OUT) :: inode

    ! Position of the coordinates in the quad
    integer, intent(OUT) :: ipos

    ! Number of the vertex the coordinates correspond to
    integer, intent(OUT) :: ivt
!</output>

!<result>
    ! Result of the searching: QTREE_NOT_FOUND / QTREE_FOUND
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
      integer, intent(INOUT) :: inode,ipos,ivt
      integer :: f
      
      f = QTREE_NOT_FOUND
      if (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. QTREE_SUBDIV) then
        
        ! Quad is subdivided. Compute child INODE which to look recursively.
        inode = rquadtree%p_Knode(qtree_getDirection(rquadtree, Ddata, inode), inode)
        f     = search(inode, ipos, ivt)
        
      else
        
        ! Quad is not subdivided. Search for (x,y) in current quad
        do ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          ivt = rquadtree%p_Knode(ipos, inode)
          if (maxval(abs(rquadtree%p_Ddata(:, ivt)-Ddata)) .le. SYS_EPSREAL) then
            f = QTREE_FOUND; return
          end if
        end do
        
      end if
    end function search
  end function qtree_searchInQuadtree

  !************************************************************************
  
!<function>
  
  pure function qtree_getDirection(rquadtree, Ddata, inode) result(d)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Ddata
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree
    
    ! Coordinates
    real(DP), dimension(2), intent(IN) :: Ddata

    ! Number of quad
    integer, intent(IN) :: inode
!</input>

!<result>
    ! Further search direction
    integer :: d
!</result>
!</function>
    
    ! local variables
    real(DP) :: xmid,ymid

    ! Compute midpoint of current quad
    xmid = 0.5*(rquadtree%p_Dbbox(QTREE_XMIN, inode)+&
                rquadtree%p_Dbbox(QTREE_XMAX, inode))
    ymid = 0.5*(rquadtree%p_Dbbox(QTREE_YMIN, inode)+&
                rquadtree%p_Dbbox(QTREE_YMAX, inode))
    
    if (Ddata(1) > xmid) then
      if (Ddata(2) > ymid) then
        d = QTREE_NE; return
      else
        d = QTREE_SE; return
      end if
    else
      if (Ddata(2) > ymid) then
        d = QTREE_NW; return
      else
        d = QTREE_SW; return
      end if
    end if
  end function qtree_getDirection

  !************************************************************************
  
!<subroutine>

  subroutine qtree_printQuadtree(rquadtree, cfilename)

!<description>
    ! This subroutine writes the content of the quadtree to a file
    ! which can be visualized by means of Matlab
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree

    ! filename of the output file
    character(LEN=*), intent(IN) :: cfilename
!</input>
!</subroutine>
    
    ! local variables
    real(DP) :: xmin,xmax,ymin,ymax
    integer :: iunit
    
    iunit = sys_getFreeUnit()
    open(UNIT=iunit, FILE=trim(adjustl(cfilename)))
    xmin = rquadtree%p_Dbbox(QTREE_XMIN, 1)
    ymin = rquadtree%p_Dbbox(QTREE_YMIN, 1)
    xmax = rquadtree%p_Dbbox(QTREE_XMAX, 1)
    ymax = rquadtree%p_Dbbox(QTREE_YMAX, 1)
    call print(xmin, ymin, xmax, ymax, 1)
    close(UNIT=iunit)
    
  contains
    
    !**************************************************************
    ! Here, the recursive print routine follows
    
    recursive subroutine print(xmin, ymin, xmax, ymax, inode)
      real(DP), intent(IN) :: xmin,ymin,xmax,ymax
      integer, intent(IN) :: inode
      real(DP) :: xmid,ymid
      integer :: i

      write(UNIT=iunit,FMT=*) 'rect'
      write(UNIT=iunit,FMT=10) xmin, ymin, xmax, ymax
      
      if (rquadtree%p_Knode(QTREE_STATUS,inode) .eq. QTREE_SUBDIV) then
        
        xmid = (xmin+xmax)/2._DP
        ymid = (ymin+ymax)/2._DP
        
        call print(xmin, ymid, xmid, ymax, rquadtree%p_Knode(QTREE_NW, inode))
        call print(xmin, ymin, xmid, ymid, rquadtree%p_Knode(QTREE_SW, inode))
        call print(xmid, ymin, xmax, ymid, rquadtree%p_Knode(QTREE_SE, inode))
        call print(xmid, ymid, xmax, ymax, rquadtree%p_Knode(QTREE_NE, inode))
        
      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) > QTREE_EMPTY) then
        
        do i = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          write(UNIT=iunit, FMT=*) 'node'
          write(UNIT=iunit, FMT=20) rquadtree%p_Ddata(:, rquadtree%p_Knode(i, inode))
        end do
        
      end if
      
10    format(4E15.6E3)
20    format(2E15.6E3)
    end subroutine print
  end subroutine qtree_printQuadtree

  !************************************************************************

!<subroutine>

  subroutine qtree_infoQuadtree(rquadtree)

!<description>
    ! This subroutine outputs statistical info about the quadtree
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree
!</input>
!</subroutine>

    call output_line('Quadtree:')
    call output_line('---------')
    call output_line('NVT:     '//trim(sys_siL(rquadtree%NVT,15)))
    call output_line('NNVT:    '//trim(sys_siL(rquadtree%NNVT,15)))
    call output_line('NNODE:   '//trim(sys_siL(rquadtree%NNODE,15)))
    call output_line('NNNODE:  '//trim(sys_siL(rquadtree%NNNODE,15)))
    call output_line('NRESIZE: '//trim(sys_siL(rquadtree%NRESIZE,5)))
    call output_line('dfactor: '//trim(sys_sdL(rquadtree%dfactor,2)))
    call output_line('h_Ddata: '//trim(sys_siL(rquadtree%h_Ddata,15)))
    call output_line('h_Dbbox: '//trim(sys_siL(rquadtree%h_Dbbox,15)))
    call output_line('h_Knode: '//trim(sys_siL(rquadtree%h_Knode,15)))
    call output_lbrk()
    write(*,*)
    call output_line('Current data memory usage: '//&
        trim(sys_sdL(100*rquadtree%NVT/real(rquadtree%NNVT,DP),2))//'%')
    call output_line('Current node memory usage: '//&
        trim(sys_sdL(100*rquadtree%NNODE/real(rquadtree%NNNODE,DP),2))//'%')
  end subroutine qtree_infoQuadtree

  !************************************************************************

!<function>

  pure function qtree_getSize(rquadtree) result(nvt)

!<description>
    ! This function returns the number of vertices stored in the quadtree
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree
!</input>

!<result>
    ! number of vertices in quadtree
    integer :: nvt
!</result>
!</function>

    nvt = rquadtree%NVT
  end function qtree_getSize

  !************************************************************************

!<function>

  function qtree_getBoundingBox(rquadtree, inode) result(bbox)
    
!<description>
    ! This function returns the bounding box of the specified quad.
    ! If no quad number is given, then the outer bounding box is
    ! returned.
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree

    ! OPTIONAL: number of quad for which bounding box should be
    ! returned
    integer, intent(IN), optional :: inode
!</input>

!<result>
    ! bounding box
    real(DP), dimension(4) :: bbox
!</result>
!</function>
    
    if (present(inode)) then
      if (inode > rquadtree%NVT) then
        call output_line('Node number exceeds quadtree dimension',&
            OU_CLASS_ERROR,OU_MODE_STD,'qtree_getBoundingBox')
        call sys_halt()
      end if
      bbox = rquadtree%p_Dbbox(QTREE_XMIN:QTREE_YMAX, inode)
    else
      bbox = rquadtree%p_Dbbox(QTREE_XMIN:QTREE_YMAX, 1)
    end if
  end function qtree_getBoundingBox

  !************************************************************************

!<function>

  elemental function qtree_getX(rquadtree, ivt) result(x)

!<description>
    ! This function returns the X-value at the given position.
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree

    ! position in the quadtree
    integer, intent(IN) :: ivt
!</input>

!<result>
    real(DP) :: x
!</result>
!</function>

    x = rquadtree%p_Ddata(1, ivt)
  end function qtree_getX

  !************************************************************************

!<function>

  elemental function qtree_getY(rquadtree, ivt) result(y)

!<description>
    ! This function returns the Y-value at the given position.
!</description>

!<input>
    ! quadtree
    type(t_quadtree), intent(IN) :: rquadtree

    ! position in the quadtree
    integer, intent(IN) :: ivt
!</input>

!<result>
    real(DP) :: y
!</result>
!</function>

    y = rquadtree%p_Ddata(2, ivt)
  end function qtree_getY

  !************************************************************************
  
!<subroutine>

  subroutine qtree_duplicateQuadtree(rquadtree, rquadtreeBackup)

!<description>
    ! This subroutine makes a copy of a quadtree in memory.
    ! It does not make sense to share some information between quadtrees,
    ! so each vectors is physically copied from the source quadtree
    ! to the destination quadtree.
!</description>

!<input>
    ! Source quadtree
    type(t_quadtree), intent(IN) :: rquadtree
!</input>

!<inputoutput>
    ! Destination quadtree
    type(t_quadtree), intent(INOUT) :: rquadtreeBackup
!</inputoutput>
!</subroutine>

    ! Release backup quadtree
    call qtree_releaseQuadtree(rquadtreeBackup)

    ! Copy all data
    rquadtreeBackup = rquadtree

    ! Reset handles
    rquadtreeBackup%h_Ddata = ST_NOHANDLE
    rquadtreeBackup%h_Dbbox = ST_NOHANDLE
    rquadtreeBackup%h_Knode = ST_NOHANDLE

    ! Copy storage blocks
    if (rquadtree%h_Ddata .ne. ST_NOHANDLE) then
      call storage_copy(rquadtree%h_Ddata, rquadtreeBackup%h_Ddata)
      call storage_getbase_double2D(rquadtreeBackup%h_Ddata,&
                                    rquadtreeBackup%p_Ddata)
    end if

    if (rquadtree%h_Dbbox .ne. ST_NOHANDLE) then
      call storage_copy(rquadtree%h_Dbbox, rquadtreeBackup%h_Dbbox)
      call storage_getbase_double2D(rquadtreeBackup%h_Dbbox,&
                                    rquadtreeBackup%p_Dbbox)
    end if

    if (rquadtree%h_Knode .ne. ST_NOHANDLE) then
      call storage_copy(rquadtree%h_Knode, rquadtreeBackup%h_Knode)
      call storage_getbase_int2D(rquadtreeBackup%h_Knode,&
                                 rquadtreeBackup%p_Knode)
    end if
  end subroutine qtree_duplicateQuadtree

  !************************************************************************

!<subroutine>

  subroutine qtree_restoreQuadtree(rquadtreeBackup, rquadtree)

!<description>
    ! This subroutine restores a quadtree from a previous backup.
!</description>

!<input>
    ! Backup of an quadtree
    type(t_quadtree), intent(IN) :: rquadtreeBackup
!</input>

!<inputoutput>
    ! Destination quadtree
    type(t_quadtree), intent(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>

    ! Release quadtree
    call qtree_releaseQuadtree(rquadtree)

    ! Duplicate the backup
    call qtree_duplicateQuadtree(rquadtreeBackup, rquadtree)
  end subroutine qtree_restoreQuadtree

  !************************************************************************

!<subroutine>
  
  subroutine resizeNVT(rquadtree, nnvt)

!<description>
    ! This subroutine reallocates memory for an existing quadtree
!</description>

!<input>
    ! New number of vertices that should be stored in the quadtree
    integer, intent(IN) :: nnvt
!</input>

!<inputoutput>
    ! quadtree that should be resized
    type(t_quadtree) :: rquadtree
!</inputoutput>
!</subroutine>

    call storage_realloc('resizeNVT', nnvt, rquadtree%h_Ddata,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_double2D(rquadtree%h_Ddata, rquadtree%p_Ddata)
    
    rquadtree%NNVT    = nnvt
    rquadtree%NRESIZE = rquadtree%NRESIZE+1
  end subroutine resizeNVT

  !************************************************************************
  
!<subroutine>
  
  subroutine resizeNNODE(rquadtree, nnnode)

!<description>
    ! This subroutine reallocates memory for an existing quadtree
!</description>

!<input>
    ! New number of quads that should be stored in the quadtree
    integer, intent(IN) :: nnnode
!</input>

!<inputoutput>
    ! quadtree that should be resized
    type(t_quadtree) :: rquadtree
!</inputoutput>
!</subroutine>

    call storage_realloc('resizeNNODE', nnnode, rquadtree%h_Dbbox,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_realloc('resizeNNODE', nnnode, rquadtree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_double2D(rquadtree%h_Dbbox, rquadtree%p_Dbbox)
    call storage_getbase_int2D(rquadtree%h_Knode,    rquadtree%p_Knode)
    
    rquadtree%NNNODE  = nnnode
    rquadtree%NRESIZE = rquadtree%NRESIZE+1
  end subroutine resizeNNODE
end module quadtree

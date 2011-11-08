!##############################################################################
!# ****************************************************************************
!# <name> quadtree </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module implements a (linear) quadtree. The implementation is
!# based on the description of quadtrees by
!#
!# R. Lohner, Applied CFD Techniques. An Introduction based on
!#            Finite Element Methods, Wiley, 2008
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
!#     -> Copy data to the quadtree
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
!# 18.) qtree_rebuildQuadtree
!#      -> Rebuilds the structure of a quadtree
!#
!# 19.) qtree_moveInQuadtree
!#      -> Reposition item in quadtree
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
  public :: qtree_rebuildQuadtree
  public :: qtree_moveInQuadtree

!<constants>

!<constantblock description="Constants for quadtree structure">
  
  ! Maximum number of items for each quad
  integer, parameter :: QTREE_MAX    = 4

  ! Item in "North-West" position
  integer, parameter :: QTREE_NW     = 1
  
  ! Item in "South-West" position
  integer, parameter :: QTREE_SW     = 2

  ! Item in "South-East" position
  integer, parameter :: QTREE_SE     = 3

  ! Item in "North-East" position
  integer, parameter :: QTREE_NE     = 4
  
  ! Position of the status information
  integer, parameter :: QTREE_STATUS = 0
  
  ! Position of the parent information
  integer, parameter :: QTREE_PARENT = -1

  ! Position of the "position" of the parent information
  integer, parameter :: QTREE_PARPOS = -2

  ! Identifier: Quad is empty
  integer, parameter :: QTREE_EMPTY  =  0

  ! Identifier: Quad is subdivided
  integer, parameter :: QTREE_SUBDIV = -1

  ! Identifier: Quad is deleted
  integer, parameter :: QTREE_DEL    = -2
  
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

  ! Operation on quadtree failed
  integer, parameter, public :: QTREE_FAILED    = -2

  ! Item could not be found in the quadtree
  integer, parameter, public :: QTREE_NOT_FOUND = -1

  ! Item could be found in the quadtree
  integer, parameter, public :: QTREE_FOUND     =  0

  ! Item was inserted into the quadtree
  integer, parameter, public :: QTREE_INSERTED  =  1

  ! Item was deleted from the quadtree
  integer, parameter, public :: QTREE_DELETED   =  2

  ! Item was moved in the quadtree
  integer, parameter, public :: QTREE_MOVED     =  3
  
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
    integer :: h_Ddata = ST_NOHANDLE

    ! Handle to bounding box
    integer :: h_Dbbox = ST_NOHANDLE

    ! Handle to quadtree structure
    !   KNODE(QTREE_STATUS,INODE) : < 0, the quad has been subdivided
    !                               = 0, the quad is empty
    !                               > 0, the number of points stored in the quad
    !   KNODE(QTREE_PARENT,INODE) : > 0, the quad the present quad came from
    !                               < 0, position of the next free quad which has been deleted
    !   KNODE(QTREE_PARPOS,INODE) : > 0, the position in the quad the present quad came from
    !   KNODE(1:NDATA,INODE)      : for KNODE(QTREE_STATUS,INODE) > 0 : the points stored in the quad
    !                               for KNODE(QTREE_STATUS,INODE) < 0 : the quads into which the present quad was subdivided
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

    ! Quadtree structure
    ! NOTE: This array is introduced to increase performance (see above).
    integer, dimension(:,:), pointer :: p_Knode

    ! Number of next free quad
    integer :: NFREE = 0

    ! Number of data items stored per quad
    ! Default value: QTREE_MAX
    integer :: NDATA = QTREE_MAX

    ! Number of vertices currently stored in the quadtree
    integer :: NVT = 0

    ! Total number of vertices that can be stored  in the quadtree
    integer :: NNVT = 0

    ! Number of subdivisions currently store in the quadtree
    integer :: NNODE = 0

    ! Total number of subdivision that can be stored in the quadtree
    integer :: NNNODE = 0

    ! Total number of resize operations
    integer :: NRESIZE = 0

    ! Factor by which the quadtree is enlarged if new storage is allocate
    real(DP) :: dfactor = 1.5_DP
    
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
    module procedure qtree_deleteFromQuadtreeDefault
    module procedure qtree_deleteFromQtreeByNumber
  end interface
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine qtree_createQuadtree(rquadtree, nnvt, nnnode,&
                                  xmin, ymin, xmax, ymax, dfactor, ndata)
  
!<description>
    ! This subroutine creates a new quadtree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the quadtree
    integer, intent(in) :: nnvt

    ! Total number of subdivisions that should be stored in the quadtree
    integer, intent(in) :: nnnode

    ! Dimensions of the initial bounding box
    real(DP), intent(in) :: xmin,ymin,xmax,ymax

    ! OPTIONAL: Factor by which the quadtree should be enlarged if
    ! new storage has to be allocated
    real(DP), intent(in), optional :: dfactor

    ! OPTIONAL: Number of data items stored per quad
    integer, optional :: ndata
!</input>

!<output>
    ! Quadtree structure
    type(t_quadtree), intent(out) :: rquadtree
!</output>
!</subroutine>
    
    ! local variables
    integer, dimension(2) :: Isize, Ilbound, Iubound
    
    ! Set factor
    if (present(dfactor)) then
      if (dfactor > 1_DP) rquadtree%dfactor=dfactor
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
    Isize = (/2, nnvt/)
    call storage_new('qtree_createQuadtree', 'p_Ddata',&
                     Isize, ST_DOUBLE, rquadtree%h_Ddata, ST_NEWBLOCK_ZERO)
    
    Isize = (/4, nnnode/)
    call storage_new('qtree_createQuadtree', 'p_Dbbox',&
                     Isize, ST_DOUBLE, rquadtree%h_Dbbox, ST_NEWBLOCK_ZERO)

    Ilbound = (/QTREE_PARPOS,1/); Iubound = (/rquadtree%NDATA,nnnode/)
    call storage_new('qtree_createQuadtree', 'p_Knode', Ilbound,&
                     Iubound, ST_INT, rquadtree%h_Knode, ST_NEWBLOCK_ZERO)
    call storage_getbase_double2D(rquadtree%h_Ddata, rquadtree%p_Ddata)
    call storage_getbase_double2D(rquadtree%h_Dbbox, rquadtree%p_Dbbox)
    call storage_getbase_int2D(rquadtree%h_Knode,    rquadtree%p_Knode)

    ! Initialize first quad
    rquadtree%nnode = 1
    rquadtree%p_Knode(QTREE_STATUS,1) = QTREE_EMPTY
    rquadtree%p_Knode(QTREE_PARENT,1) = QTREE_EMPTY
    rquadtree%p_Knode(QTREE_PARPOS,1) = QTREE_EMPTY
    rquadtree%p_Knode(1:rquadtree%ndata, 1) = 0
    
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
    ! This subroutine releases an existing quadtree.
!</description>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
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
    rquadtree%NFREE   = 0

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
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree
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

    ! Call copy routine
    call qtree_copyFromQuadtree(rquadtree, p_Ddata)
   
  end subroutine qtree_copyFromQuadtree_handle

  !************************************************************************

!<subroutine>

  subroutine qtree_copyFromQuadtree_array(rquadtree, p_Ddata)

!<description>
    ! This subroutine copies the content of the quadtree to an array.
    ! If the array is too small, then this subroutines terminates.
!</description>

!<input>
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree
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
    if (Isize(1) .ne. 2 .or. Isize(2) < rquadtree%NVT) then
      call output_line('Array too small!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyFromQuadtree_array')
      call sys_halt()
    end if

    ! Copy data
    call DCOPY(2*rquadtree%NVT, rquadtree%p_Ddata, 1, p_Ddata, 1)
    
  end subroutine qtree_copyFromQuadtree_array

  !************************************************************************

!<subroutine>
  
  subroutine qtree_copyToQuadtree_handle(h_Ddata, rquadtree)

!<description>
    ! This subroutine copies the content of a handle to the quadtree.
!</description>

!<input>
    ! Handle to the coordinate vector
    integer, intent(in) :: h_Ddata
!</input>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata
    
    ! Set pointer
    call storage_getbase_double2D(h_Ddata, p_Ddata)

    ! Call copy routine
    call qtree_copyToQuadtree(p_Ddata, rquadtree)
    
  end subroutine qtree_copyToQuadtree_handle

  !************************************************************************

!<subroutine>
  
  subroutine qtree_copyToQuadtree_array(p_Ddata, rquadtree)

!<description>
    ! This subroutine copies the content of an array to the quadtree.
!</description>

!<input>
    ! Coordinate vector
    real(DP), dimension(:,:), intent(in) :: p_Ddata
!</input>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>
    
    
    if (size(p_Ddata,1) .ne. 2) then
      call output_line('First dimension of array must be 2!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyToQuadtree_handle')
      call sys_halt()
    end if

    ! Adjust dimension of data array
    rquadtree%NVT = size(p_Ddata,2)
    if (rquadtree%NNVT .lt. rquadtree%NVT) then
      call resizeNVT(rquadtree, rquadtree%NVT)
    end if

    ! Copy data array
    call DCOPY(2*rquadtree%NVT, p_Ddata, 1, rquadtree%p_Ddata, 1)
 
    ! Rebuild structure
    call qtree_rebuildQuadtree(rquadtree)
      
  end subroutine qtree_copyToQuadtree_array

  !************************************************************************
  
!<function>
  
  function qtree_insertIntoQuadtree(rquadtree, Ddata, ivt, inode) result(iresult)

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
    real(DP), dimension(2), intent(in) :: Ddata

    ! OPTIONAL: Number of the quad to which vertex should be inserted.
    ! If there is no space left, then the next free position will be used
    integer, intent(in), optional :: inode
!</input>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
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
      iresult = qtree_searchInQuadtree(rquadtree, Ddata, jnode, jpos, jvt)
      if (iresult .eq. QTREE_FOUND) then
        ivt = jvt
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
    ivt = rquadtree%NVT
    rquadtree%p_Ddata(:,ivt) = Ddata
    
    ! Insert entry recursively
    iresult = insert(ivt, jnode)

    ! Check success
    if (iresult .eq. QTREE_FAILED) then
      rquadtree%NVT = rquadtree%NVT-1
      rquadtree%p_Ddata(:,ivt) = 0.0_DP
      ivt = 0
    end if

  contains
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    recursive function insert(ivt, inode) result(iresult)
      
      integer, intent(in) :: ivt,inode
      integer :: iresult

      ! local variables
      real(DP) :: xmin,ymin,xmax,ymax,xmid,ymid
      integer :: i,jvt,jnode,nnode,isize,jresult


      ! Check status of current quad
      if (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. rquadtree%ndata) then
      
        ! Quad is full

        ! Check if there are deleted quads
        if (rquadtree%NFREE .ne. QTREE_EMPTY) then

          ! Reuse memory of the four quads which have been deleted lately
          nnode = rquadtree%NFREE

          ! Update pointer to next free quads
          rquadtree%NFREE = -rquadtree%p_Knode(QTREE_PARENT, nnode)

          ! Decrease starting position of first quad by one
          nnode = nnode-1

        else

          ! Otherwise, create new memory if required
          if (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode) then
            isize = max(rquadtree%nnode+QTREE_MAX,&
                        ceiling(rquadtree%dfactor*rquadtree%NNNODE))
            call resizeNNODE(rquadtree, isize)
          end if

          ! New quads are stored after all existing quads
          nnode = rquadtree%NNODE
          
        end if
        
        ! Increase number of quads
        rquadtree%NNODE = rquadtree%NNODE+QTREE_MAX
        
        ! Compute spatial coordinates of bounding boxes
        xmin = rquadtree%p_Dbbox(QTREE_XMIN,inode)
        ymin = rquadtree%p_Dbbox(QTREE_YMIN,inode)
        xmax = rquadtree%p_Dbbox(QTREE_XMAX,inode)
        ymax = rquadtree%p_Dbbox(QTREE_YMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        
        ! NW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NW) = QTREE_NW
        rquadtree%p_Dbbox(:,nnode+QTREE_NW)            = (/xmin,ymid,xmid,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NW) = 0
        
        ! SW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_Dbbox(:,nnode+QTREE_SW)            = (/xmin,ymin,xmid,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SW) = 0
        
        ! SE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_Dbbox(:,nnode+QTREE_SE)            = (/xmid,ymin,xmax,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SE) = 0
        
        ! NE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_Dbbox(:,nnode+QTREE_NE)            = (/xmid,ymid,xmax,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NE) = 0
        
        ! Add the data values from INODE to the four new quads
        do i = 1, rquadtree%ndata
          jvt = rquadtree%p_Knode(i, inode)
          jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,jvt), inode)
          jresult = insert(jvt, jnode)
          
          if (jresult .eq. QTREE_FAILED) then
            call output_line('Internal error in insertion!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'qtree_insertIntoQuadtree')
            call sys_halt()
          end if
        end do
        
        ! Mark the current quad as subdivided and set pointers to its four children
        rquadtree%p_Knode(QTREE_STATUS, inode) =     QTREE_SUBDIV
        rquadtree%p_Knode(1:QTREE_MAX, inode)  = - (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                    nnode+QTREE_SE, nnode+QTREE_NE/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,ivt), inode)
        iresult = insert(ivt, inode)
        
      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .ge. QTREE_EMPTY) then
        
        ! There is still some space in the quad
        rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)+1
        rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = ivt
        
        ! Set success flag
        iresult = QTREE_INSERTED

      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. QTREE_SUBDIV) then

        ! Proceed to correcponding sub-tree
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, Ddata, inode), inode)
        iresult = insert(ivt, jnode)

      else

        ! Set failure flag
        iresult = QTREE_FAILED

      end if
      
    end function insert
    
  end function qtree_insertIntoQuadtree
  
  ! ***************************************************************************
  
!<function>
  
  function qtree_deleteFromQuadtreeDefault(rquadtree, Ddata, ivt) result(iresult)

!<description>
    ! This function deletes an item from the quadtree.
    ! The value IVT returns the number of the item which is
    ! moved to the position of the deleted item. If the deleted
    ! item was the last one, then IVT=NVT is returned.
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    real(DP), dimension(2), intent(in) :: Ddata
!</input>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
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
    
    
    ! Delete item starting at root
    iresult = delete(1, ivt)

  contains

    !**************************************************************
    ! Here, the recursive deletion routine follows

    recursive function delete(inode, ivt) result(iresult)
      
      integer, intent(in) :: inode
      integer, intent(inout) :: ivt
      integer :: iresult
      
      ! local variables
      real(DP), dimension(2) :: DdataTmp
      integer, dimension(QTREE_MAX) :: Knode
      integer :: i,jvt,jnode,ipos,jpos,nemptyChildren
      
      
      ! Check status of current quad
      select case(rquadtree%p_Knode(QTREE_STATUS, inode))

      case (QTREE_SUBDIV)   ! Quad is subdivided
        
        ! Compute child INODE which to look recursively.
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, Ddata, inode), inode)
        iresult = delete(jnode, ivt)
        
        ! Save values from current quad
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

          ! Mark quad as empty
          rquadtree%p_Knode(QTREE_STATUS, inode) = QTREE_EMPTY

          ! Copy data from non-empty child (if any) and mark nodes as deleted
          do i = 1, QTREE_MAX
            jnode = -Knode(i)
            if (rquadtree%p_Knode(QTREE_STATUS, jnode) .gt. QTREE_EMPTY) then
              ! Copy status of quad
              rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, jnode)
              
              ! Copy data
              rquadtree%p_Knode(1:rquadtree%NDATA, inode) = rquadtree%p_Knode(1:rquadtree%NDATA, jnode)
            end if
            
            ! Mark quad as deleted
            rquadtree%p_Knode(QTREE_STATUS, jnode) = QTREE_DEL

            ! Set pointer to next free position
            rquadtree%p_Knode(QTREE_PARENT, jnode) = -rquadtree%NFREE

          end do
          
          ! Update pointer to next free position
          rquadtree%NFREE = -Knode(1)
          
          ! Reduce number of quads
          rquadtree%NNODE = rquadtree%NNODE-QTREE_MAX

        end if
                
  
      case (QTREE_EMPTY)   ! Quad is empty so it cannot contain the item

        iresult = QTREE_FAILED
        ivt = 0


      case (QTREE_DEL)   ! Quad is deleted -> serious error

        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_deleteFromQuadtreeDefault')
        call sys_halt()
        

      case default   ! Quad is a non-empty leaf
        
        ! Search for (x,y) in current quad
        do ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          
          ! Get vertex number
          ivt = rquadtree%p_Knode(ipos, inode)

          if (maxval(abs(rquadtree%p_Ddata(:, ivt)-Ddata)) .le. SYS_EPSREAL_DP) then
           
            ! Physically remove the item IVT from node INODE
            jpos = rquadtree%p_Knode(QTREE_STATUS, inode)
            
            rquadtree%p_Knode(ipos, inode) = rquadtree%p_Knode(jpos, inode)
            rquadtree%p_Knode(jpos, inode) = 0
            rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)-1
            
            ! If IVT is not last item then find item with largest number JVT
            if (ivt .ne. rquadtree%NVT) then
              DdataTmp(:) = rquadtree%p_Ddata(:,rquadtree%NVT)
              if (qtree_searchInQuadtree(rquadtree, DdataTmp(:),&
                                         jnode, jpos, jvt) .eq. QTREE_FOUND) then

                ! Move last item JVT to position IVT
                rquadtree%p_Ddata(:, ivt) = rquadtree%p_Ddata(:, jvt)
                rquadtree%p_Knode(jpos, jnode) = ivt
              else
                call output_line('Internal error in deletion!',&
                                 OU_CLASS_ERROR,OU_MODE_STD,'qtree_deleteFromQuadtreeDefault')
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

  end function qtree_deleteFromQuadtreeDefault

  ! ***************************************************************************

!<function>

  function qtree_deleteFromQtreeByNumber(rquadtree, ivt, ivtReplace) result(iresult)

!<description>
    ! This function deletes vertex with number IVT from the quadtree.
!</description>

!<input>
    ! Number of the vertex to be deleted
    integer, intent(in) :: ivt
!</input>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
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
    real(DP), dimension(2) :: Ddata
    
    if (ivt .le. rquadtree%NVT) then
      ! Get coordinates and invoke deletion routine
      Ddata   = rquadtree%p_Ddata(:,ivt)
      iresult = qtree_deleteFromQuadtree(rquadtree, Ddata, ivtReplace)
    else
      iresult = QTREE_FAILED
    end if

  end function qtree_deleteFromQtreeByNumber
 
  ! ***************************************************************************

!<function>
  
  function qtree_searchInQuadtree(rquadtree, Ddata, inode, ipos, ivt) result(iresult)

!<description>
    ! This subroutine searches for given coordinates in the quadtree.
    ! The result of the search operation is returned by the value IRESULT.
    ! If the item was found than INODE is the number of the quad, IPOS
    ! is the position of the item in the quad and IVT is number of the
    ! item in the data array. Otherwise, INODE is the number of the leaf,
    ! where the item would be placed in case of insertion.
!</description>

!<input>
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree
    
    ! coordinates that should be searched
    real(DP), dimension(2), intent(in) :: Ddata
!</input>

!<output>
    ! Number of the quad in which the given coordinates are
    integer, intent(out) :: inode

    ! Position of the coordinates in the quad
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
      
      
      ! Check status of current quad
      select case(rquadtree%p_Knode(QTREE_STATUS, inode))

      case (QTREE_SUBDIV)   ! Quad is subdivided
        
        ! Compute child INODE which to look recursively.
        inode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, Ddata, inode), inode)
        iresult = search(inode, ipos, ivt)
        
        
      case (QTREE_EMPTY)   ! Quad is empty so it cannot contain the item
        
        iresult = QTREE_NOT_FOUND
        
        
      case (QTREE_DEL)   ! Quad is deleted -> serious error
        
        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_searchInQuadtree')
        call sys_halt()
      

      case default   ! Quad is a non-empty leaf

        ! Search for (x,y) in current quad
        do ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          
          ! Get vertex number
          ivt = rquadtree%p_Knode(ipos, inode)
          
          if (maxval(abs(rquadtree%p_Ddata(:, ivt)-Ddata)) .le. SYS_EPSREAL_DP) then

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

  end function qtree_searchInQuadtree

  !************************************************************************
  
!<function>
  
  pure function qtree_getDirection(rquadtree, Ddata, inode) result(idirection)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Ddata.
!</description>

!<input>
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree
    
    ! Coordinates
    real(DP), dimension(2), intent(in) :: Ddata

    ! Number of quad
    integer, intent(in) :: inode
!</input>

!<result>
    ! Further search direction
    integer :: idirection
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
        idirection = QTREE_NE
      else
        idirection = QTREE_SE
      end if
    else
      if (Ddata(2) > ymid) then
        idirection = QTREE_NW
      else
        idirection = QTREE_SW
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
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree

    ! filename of the output file
    character(LEN=*), intent(in) :: cfilename
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
      real(DP), intent(in) :: xmin,ymin,xmax,ymax
      integer, intent(in) :: inode
      real(DP) :: xmid,ymid
      integer :: i,j

      write(UNIT=iunit,FMT=*) 'rect'
      write(UNIT=iunit,FMT=10) xmin, ymin, xmax, ymax, inode
      
      if (rquadtree%p_Knode(QTREE_STATUS,inode) .eq. QTREE_SUBDIV) then
        
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        
        call print(xmin, ymid, xmid, ymax, -rquadtree%p_Knode(QTREE_NW, inode))
        call print(xmin, ymin, xmid, ymid, -rquadtree%p_Knode(QTREE_SW, inode))
        call print(xmid, ymin, xmax, ymid, -rquadtree%p_Knode(QTREE_SE, inode))
        call print(xmid, ymid, xmax, ymax, -rquadtree%p_Knode(QTREE_NE, inode))
        
      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) > QTREE_EMPTY) then
        
        do i = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          write(UNIT=iunit, FMT=*) 'node'
          j = rquadtree%p_Knode(i, inode)
          write(UNIT=iunit, FMT=20) rquadtree%p_Ddata(1, j), rquadtree%p_Ddata(2, j), j
        end do
        
      end if
      
10    format(4(E15.6E3,1X),I10)
20    format(2(E15.6E3,1X),I10)

    end subroutine print

  end subroutine qtree_printQuadtree

  !************************************************************************

!<subroutine>

  subroutine qtree_infoQuadtree(rquadtree)

!<description>
    ! This subroutine outputs statistical info about the quadtree
!</description>

!<input>
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree
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
    call output_line('h_Ddata: '//trim(sys_siL(rquadtree%h_Ddata,15)))
    call output_line('h_Dbbox: '//trim(sys_siL(rquadtree%h_Dbbox,15)))
    call output_line('h_Knode: '//trim(sys_siL(rquadtree%h_Knode,15)))
    call output_lbrk()
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
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree
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
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree

    ! OPTIONAL: number of quad for which bounding box should be
    ! returned
    integer, intent(in), optional :: inode
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
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree

    ! position in the quadtree
    integer, intent(in) :: ivt
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
    ! Quadtree
    type(t_quadtree), intent(in) :: rquadtree

    ! position in the quadtree
    integer, intent(in) :: ivt
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
    type(t_quadtree), intent(in) :: rquadtree
!</input>

!<inputoutput>
    ! Destination quadtree
    type(t_quadtree), intent(inout) :: rquadtreeBackup
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
    type(t_quadtree), intent(in) :: rquadtreeBackup
!</input>

!<inputoutput>
    ! Destination quadtree
    type(t_quadtree), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    ! Release quadtree
    call qtree_releaseQuadtree(rquadtree)

    ! Duplicate the backup
    call qtree_duplicateQuadtree(rquadtreeBackup, rquadtree)

  end subroutine qtree_restoreQuadtree

  !************************************************************************

!<subroutine>

  subroutine qtree_rebuildQuadtree(rquadtree)

!<description>
    ! This subroutine rebuilds the structure of a quadtree
!</description>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Inodes
    integer, dimension(2) :: Isize
    integer :: h_Inodes

    ! Initialization
    rquadtree%NNODE = 1
    rquadtree%NFREE = 0

    call storage_getsize(rquadtree%h_Knode, Isize); rquadtree%NNNODE = Isize(2)
    call storage_getsize(rquadtree%h_Ddata, Isize); rquadtree%NNVT   = Isize(2)

    ! Check if quadtree contains data
    if (rquadtree%NVT .eq. 0) then

      ! Initialize first quad
      rquadtree%nnode = 1
      rquadtree%p_Knode(QTREE_STATUS, 1) = QTREE_EMPTY
      rquadtree%p_Knode(QTREE_PARENT, 1) = QTREE_EMPTY
      rquadtree%p_Knode(QTREE_PARPOS, 1) = QTREE_EMPTY
      rquadtree%p_Knode(1:rquadtree%NDATA, 1) = 0
      
    else
      
      ! Create temporary storage
      h_Inodes = ST_NOHANDLE
      call storage_new('qtree_rebuildQuadtree', 'Inodes', rquadtree%NVT,&
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
      real(DP) :: xmin,ymin,xmax,ymax,xmid,ymid
      integer :: i,isize,jnode,nnode,ivt,imid1,imid2,imid3

      
      ! Check if istart > iend then the quad is empty
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

        ! Store the number of quads
        nnode = rquadtree%NNODE
        rquadtree%NNODE = nnode+QTREE_MAX

        ! Mark the current quad as subdivided and set pointers to its four children
        rquadtree%p_Knode(QTREE_STATUS,inode) =     QTREE_SUBDIV
        rquadtree%p_Knode(1:QTREE_MAX,inode)  = - (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                    nnode+QTREE_SE, nnode+QTREE_NE/)
        
        ! Determine coordinates for bounding boxes
        xmin = rquadtree%p_Dbbox(QTREE_XMIN,inode)
        ymin = rquadtree%p_Dbbox(QTREE_YMIN,inode)
        xmax = rquadtree%p_Dbbox(QTREE_XMAX,inode)
        ymax = rquadtree%p_Dbbox(QTREE_YMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)

        ! NW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NW) = QTREE_NW
        rquadtree%p_Dbbox(:,nnode+QTREE_NW)            = (/xmin,ymid,xmid,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NW) = 0
        
        ! SW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_Dbbox(:,nnode+QTREE_SW)            = (/xmin,ymin,xmid,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SW) = 0
        
        ! SE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_Dbbox(:,nnode+QTREE_SE)            = (/xmid,ymin,xmax,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SE) = 0
        
        ! NE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_Dbbox(:,nnode+QTREE_NE)            = (/xmid,ymid,xmax,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NE) = 0
        
        
        ! Partition nodes with respect to x-axis
        imid1 = partition(istart, iend, 1, xmid)
          
        ! Partition nodes with respect to y-axis
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
        imid = merge(iend, iend+1, rquadtree%p_Ddata(idim, p_Inodes(iend)) .gt. dmid)

      else
        
        ! Otherwise, sort items in partition
        call quicksort(istart, iend, idim)
        
        ! Find location of first item belonging to "IS GREATER" partition
        do i = istart, iend
          if (rquadtree%p_Ddata(idim, p_Inodes(i)) .gt. dmid) then
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
      dpivot = rquadtree%p_Ddata(idim, p_Inodes(iend))

      do
        do while((i .lt. iend) .and. (rquadtree%p_Ddata(idim, p_Inodes(i)) .le. dpivot))
          i = i+1
        end do
        
        do while((j .gt. istart) .and. (rquadtree%p_Ddata(idim, p_Inodes(j)) .ge. dpivot))
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

  end subroutine qtree_rebuildQuadtree

  !************************************************************************
  
!<function>
  
  function qtree_moveInQuadtree(rquadtree, DdataOld, DdataNew) result(iresult)

!<description>
    ! This function modifies the coordinates of an item in the quadtree.
    ! First, the item with coordinates DDATA is searched. If the new
    ! coordinates belong to the same node (e.g. they are comprised in the
    ! nodes bounding box), then the coordinates are simply updated.
    ! Otherwise, the vertex is deleted from the quadtree and re-inserted
    ! with the new coordinates, whereby the vertex number does not change.
!</description>

!<input>
    ! Old coordinates of the vertex
    real(DP), dimension(2), intent(in) :: DdataOld

    ! New coordinates of the vertex
    real(DP), dimension(2), intent(in) :: DdataNew
!</input>

!<inputoutput>
    ! Quadtree
    type(t_quadtree), intent(inout) :: rquadtree
!</inputoutput>

!<result>
    ! Result of the insertion:
    !   QTREE_FAILED
    !   QTREE_MOVED
    integer :: iresult
!</result>
!</function>

    ! local variables
    integer :: ivt,jvt,jnode,jpos

    ! Move item starting at root
    iresult = move(1, ivt)

    if (iresult .eq. QTREE_DELETED) then
      ! Search potential candidate for insertion
      iresult = qtree_searchInQuadtree(rquadtree, DdataNew, jnode, jpos, jvt)
      if (iresult .eq. QTREE_FOUND) then
        call output_line('Duplicate entry in quadtree!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_moveInQuadtree')
        call sys_halt()
      end if
      
      ! Update values
      rquadtree%p_Ddata(:,ivt) = DdataNew

      ! Insert entry recursively
      iresult = insert(ivt, jnode)
    end if

  contains

    !**************************************************************
    ! Here, the recursive move routine follows

    recursive function move(inode, ivt) result(iresult)
      
      integer, intent(in) :: inode
      integer, intent(inout) :: ivt
      integer :: iresult
      
      ! local variables
      real(DP), dimension(2) :: DdataTmp
      integer, dimension(QTREE_MAX) :: Knode
      integer :: i,jvt,jnode,ipos,jpos,nemptyChildren
      
      
      ! Check status of current quad
      select case(rquadtree%p_Knode(QTREE_STATUS, inode))

      case (QTREE_SUBDIV)   ! Quad is subdivided
        
        ! Compute child INODE which to look recursively.
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, DdataOld, inode), inode)
        iresult = move(jnode, ivt)
        
        ! Save values from current quad
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

          ! Mark quad as empty
          rquadtree%p_Knode(QTREE_STATUS, inode) = QTREE_EMPTY

          ! Copy data from non-empty child (if any) and mark nodes as deleted
          do i = 1, QTREE_MAX
            jnode = -Knode(i)
            if (rquadtree%p_Knode(QTREE_STATUS, jnode) .gt. QTREE_EMPTY) then
              ! Copy status of quad
              rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, jnode)
              
              ! Copy data
              rquadtree%p_Knode(1:rquadtree%NDATA, inode) = rquadtree%p_Knode(1:rquadtree%NDATA, jnode)
            end if
            
            ! Mark quad as deleted
            rquadtree%p_Knode(QTREE_STATUS, jnode) = QTREE_DEL

            ! Set pointer to next free position
            rquadtree%p_Knode(QTREE_PARENT, jnode) = -rquadtree%NFREE

          end do
          
          ! Update pointer to next free position
          rquadtree%NFREE = -Knode(1)
          
          ! Reduce number of quads
          rquadtree%NNODE = rquadtree%NNODE-QTREE_MAX

        end if
                
  
      case (QTREE_EMPTY)   ! Quad is empty so it cannot contain the item

        iresult = QTREE_FAILED
        ivt = 0


      case (QTREE_DEL)   ! Quad is deleted -> serious error

        call output_line('Internal error in deletion!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'qtree_moveInQuadtree')
        call sys_halt()
        

      case default   ! Quad is a non-empty leaf
        
        ! Search for (x,y) in current quad
        do ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          
          ! Get vertex number
          ivt = rquadtree%p_Knode(ipos, inode)

          if (maxval(abs(rquadtree%p_Ddata(:, ivt)-DdataOld)) .le. SYS_EPSREAL_DP) then
           
            ! Check if the new coordinates can be stored in the same node
            if ((rquadtree%p_Dbbox(QTREE_XMIN,inode) .le. DdataNew(1)) .and.&
                (rquadtree%p_Dbbox(QTREE_XMAX,inode) .ge. DdataNew(1)) .and.&
                (rquadtree%p_Dbbox(QTREE_YMIN,inode) .le. DdataNew(2)) .and.&
                (rquadtree%p_Dbbox(QTREE_YMAX,inode) .ge. DdataNew(2))) then

              ! Just update coordinates
              rquadtree%p_Ddata(:,ivt) = DdataNew

              ! We have updated the item IVT in node INODE
              iresult = QTREE_MOVED
              
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
      real(DP) :: xmin,ymin,xmax,ymax,xmid,ymid
      integer :: i,jvt,jnode,nnode,isize,jresult


      ! Check status of current quad
      if (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. rquadtree%ndata) then
      
        ! Quad is full

        ! Check if there are deleted quads
        if (rquadtree%NFREE .ne. QTREE_EMPTY) then

          ! Reuse memory of the four quads which have been deleted lately
          nnode = rquadtree%NFREE

          ! Update pointer to next free quads
          rquadtree%NFREE = -rquadtree%p_Knode(QTREE_PARENT, nnode)

          ! Decrease starting position of first quad by one
          nnode = nnode-1

        else

          ! Otherwise, create new memory if required
          if (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode) then
            isize = max(rquadtree%nnode+QTREE_MAX,&
                        ceiling(rquadtree%dfactor*rquadtree%NNNODE))
            call resizeNNODE(rquadtree, isize)
          end if

          ! New quads are stored after all existing quads
          nnode = rquadtree%NNODE
          
        end if
        
        ! Increase number of quads
        rquadtree%NNODE = rquadtree%NNODE+QTREE_MAX
        
        ! Compute spatial coordinates of bounding boxes
        xmin = rquadtree%p_Dbbox(QTREE_XMIN,inode)
        ymin = rquadtree%p_Dbbox(QTREE_YMIN,inode)
        xmax = rquadtree%p_Dbbox(QTREE_XMAX,inode)
        ymax = rquadtree%p_Dbbox(QTREE_YMAX,inode)
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        
        ! NW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NW) = QTREE_NW
        rquadtree%p_Dbbox(:,nnode+QTREE_NW)            = (/xmin,ymid,xmid,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NW) = 0
        
        ! SW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_Dbbox(:,nnode+QTREE_SW)            = (/xmin,ymin,xmid,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SW) = 0
        
        ! SE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_Dbbox(:,nnode+QTREE_SE)            = (/xmid,ymin,xmax,ymid/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_SE) = 0
        
        ! NE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_Dbbox(:,nnode+QTREE_NE)            = (/xmid,ymid,xmax,ymax/)
        rquadtree%p_Knode(1:rquadtree%NDATA, nnode+QTREE_NE) = 0
        
        ! Add the data values from INODE to the four new quads
        do i = 1, rquadtree%ndata
          jvt = rquadtree%p_Knode(i, inode)
          jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,jvt), inode)
          jresult = insert(jvt, jnode)
          
          if (jresult .eq. QTREE_FAILED) then
            call output_line('Internal error in insertion!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'qtree_moveInQuadtree')
            call sys_halt()
          end if
        end do
        
        ! Mark the current quad as subdivided and set pointers to its four children
        rquadtree%p_Knode(QTREE_STATUS, inode) =     QTREE_SUBDIV
        rquadtree%p_Knode(1:QTREE_MAX, inode)  = - (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                    nnode+QTREE_SE, nnode+QTREE_NE/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,ivt), inode)
        iresult = insert(ivt, inode)
        
      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .ge. QTREE_EMPTY) then
        
        ! There is still some space in the quad
        rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)+1
        rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = ivt
        
        ! Set success flag
        iresult = QTREE_INSERTED

      elseif (rquadtree%p_Knode(QTREE_STATUS, inode) .eq. QTREE_SUBDIV) then

        ! Proceed to correcponding sub-tree
        jnode = -rquadtree%p_Knode(qtree_getDirection(rquadtree, DdataNew, inode), inode)
        iresult = insert(ivt, jnode)

      else

        ! Set failure flag
        iresult = QTREE_FAILED

      end if
      
    end function insert

  end function qtree_moveInQuadtree
  
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
    type(t_quadtree), intent(inout) :: rquadtree
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
    integer, intent(in) :: nnnode
!</input>

!<inputoutput>
    ! Quadtree that should be resized
    type(t_quadtree), intent(inout) :: rquadtree
!</inputoutput>
!</subroutine>

    call storage_realloc('resizeNNODE', nnnode, rquadtree%h_Dbbox,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_realloc('resizeNNODE', 1, nnnode, rquadtree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .true.)
    call storage_getbase_double2D(rquadtree%h_Dbbox, rquadtree%p_Dbbox)
    call storage_getbase_int2D(rquadtree%h_Knode,    rquadtree%p_Knode)
    
    rquadtree%NNNODE  = nnnode
    rquadtree%NRESIZE = rquadtree%NRESIZE+1

  end subroutine resizeNNODE

end module quadtree

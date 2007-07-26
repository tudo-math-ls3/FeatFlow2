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
!# 9.) otree_printOctree
!#     -> Write octree to file
!#
!# 10.) otree_infoOctree
!#      -> Output info about octree
!#
!# 11.) otree_getsize
!#      -> Return number of vertices in octree
!#
!# 12.) otree_getBoundingBox
!#      -> Return the outer bounding box
!#
!# 13.) otree_getX
!#      -> Return the X-value at a given position
!#
!# 14.) otree_getY
!#      -> Return the Y-value at a given position
!#
!# 15.) otree_getZ
!#      -> Return the Z-value at a given position
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

MODULE octree
  USE fsystem
  USE storage
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_octree
  PUBLIC :: otree_createOctree
  PUBLIC :: otree_releaseOctree
  PUBLIC :: otree_copyToOctree
  PUBLIC :: otree_copyFromOctree
  PUBLIC :: otree_insertIntoOctree
  PUBLIC :: otree_deleteFromOctree
  PUBLIC :: otree_searchInOctree
  PUBLIC :: otree_printOctree
  PUBLIC :: otree_infoOctree
  PUBLIC :: otree_getsize
  PUBLIC :: otree_getBoundingBox
  PUBLIC :: otree_getX
  PUBLIC :: otree_getY
  PUBLIC :: otree_getZ

!<constants>

!<constantblock description="KIND values for octree data">
  
  ! Kind value for indices in quadtree
  INTEGER, PARAMETER, PUBLIC :: PREC_OTREEIDX = I32

!</constantblock>

!<constantblock description="Constants for octree structure">
  
  ! Position of the status information
  INTEGER, PARAMETER :: OTREE_STATUS = 11
  
  ! Position of the parent information
  INTEGER, PARAMETER :: OTREE_PARENT = 10

  ! Position of the "position" of the parent information
  INTEGER, PARAMETER :: OTREE_PARPOS = 9

  ! Number of free positions in cube
  INTEGER, PARAMETER :: OTREE_FREE   = 9

  ! Item in "North-East-Back" position
  INTEGER, PARAMETER :: OTREE_NEB    = 8

  ! Item in "South-East-Back" position
  INTEGER, PARAMETER :: OTREE_SEB    = 7

  ! Item in "South-West-Back" position
  INTEGER, PARAMETER :: OTREE_SWB     = 6

  ! Item in "North-West-Back" position
  INTEGER, PARAMETER :: OTREE_NWB    = 5

  ! Item in "North-East-Front" position
  INTEGER, PARAMETER :: OTREE_NEF    = 4

  ! Item in "South-East-Front" position
  INTEGER, PARAMETER :: OTREE_SEF    = 3

  ! Item in "South-West-Front" position
  INTEGER, PARAMETER :: OTREE_SWF    = 2

  ! Item in "North-West-Front" position
  INTEGER, PARAMETER :: OTREE_NWF    = 1

  ! Identifier: Cube is empty
  INTEGER, PARAMETER :: OTREE_EMPTY  = 0

  ! Identifier: Status is subdivided
  INTEGER, PARAMETER :: OTREE_SUBDIV = -1
  
  ! Maximum number of items for each quad
  INTEGER, PARAMETER :: OTREE_MAX    = 8

!</constantblock> 

!<constantblock description="Constants for octree bounding-box">

  ! Position of the x-minimal value
  INTEGER, PARAMETER :: OTREE_XMIN   =  1

  ! Position of the y-minimal value
  INTEGER, PARAMETER :: OTREE_YMIN   =  2

  ! Position of the z-minimal value
  INTEGER, PARAMETER :: OTREE_ZMIN   =  3

  ! Position of the x-maximal value
  INTEGER, PARAMETER :: OTREE_XMAX   =  4

  ! Position of the y-maximal value
  INTEGER, PARAMETER :: OTREE_YMAX   =  5

  ! Position of the z-maximal value
  INTEGER, PARAMETER :: OTREE_ZMAX   =  6

!</constantblock>   
  
!<constantblock description="Constants for octreetree operations">

  ! Item could not be found in the octree
  INTEGER, PARAMETER, PUBLIC :: OTREE_NOT_FOUND = -1

  ! Item could be found in the octree
  INTEGER, PARAMETER, PUBLIC :: OTREE_FOUND     =  0
!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! A linear octree implemented as array

  TYPE t_octree

    ! Remark: The content of this derived data type is declared private.
    ! Hence, it cannot be accessed outside of this module. This allows us
    ! to use techniques such as the pointers which are used to increase
    ! the performace. These pointers cannot be modified externally so that
    ! data consistency is guaranteed.
    PRIVATE

    ! Handle to data vector
    INTEGER :: h_Ddata             = ST_NOHANDLE

    ! Handle to bounding box
    INTEGER :: h_Dbbox             = ST_NOHANDLE

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
    INTEGER :: h_Knode             = ST_NOHANDLE

    ! Data vectors
    ! NOTE: This array is introduced to increase performance. It should not be touched 
    ! by the user. To this end, all quantities of this derived type are PRIVATE.
    ! If the handle h_Ddata would be dereferences for each operation such as search, 
    ! delete, performance would be very poor.
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata

    ! Coordinates of the bounding box
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dbbox

    ! Octree structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_OTREEIDX), DIMENSION(:,:), POINTER :: p_Knode

    ! Number of vertices currently stored in the octree
    INTEGER(PREC_OTREEIDX) :: NVT    = 0

    ! Total number of vertices that can be stored  in the octree
    INTEGER(PREC_OTREEIDX) :: NNVT   = 0

    ! Number of nodes currently store in the octree
    INTEGER(PREC_OTREEIDX) :: NNODE  = 0

    ! Total number of nodes that can be stored in the octree
    INTEGER(PREC_OTREEIDX) :: NNNODE = 0

    ! Total number of resize operations
    INTEGER :: NRESIZE               = 0

    ! Factor by which the octree is enlarged if new storage is allocate
    REAL(DP) :: dfactor              = 1.5_DP
  END TYPE t_octree
  
!</typeblock>
!</types>

  INTERFACE otree_createOctree
    MODULE PROCEDURE t_octree_create
  END INTERFACE
  
  INTERFACE otree_releaseOctree
    MODULE PROCEDURE t_octree_release
  END INTERFACE
  
  INTERFACE otree_copyToOctree
    MODULE PROCEDURE t_octree_copyto_handle
    MODULE PROCEDURE t_octree_copyto_array
  END INTERFACE

  INTERFACE otree_copyFromOctree
    MODULE PROCEDURE t_octree_copyfrom_handle
    MODULE PROCEDURE t_octree_copyfrom_array
  END INTERFACE
  
  INTERFACE otree_insertIntoOctree
    MODULE PROCEDURE t_octree_insert
  END INTERFACE
  INTERFACE insert   ! for internal use
    MODULE PROCEDURE t_octree_insert
  END INTERFACE

  INTERFACE otree_deleteFromOctree
    MODULE PROCEDURE t_octree_delete
  END INTERFACE
  INTERFACE delete   ! for internal use
    MODULE PROCEDURE t_octree_delete
  END INTERFACE
  
  INTERFACE otree_searchInOctree
    MODULE PROCEDURE t_octree_search
  END INTERFACE
  INTERFACE search   ! for internal use
    MODULE PROCEDURE t_octree_search
  END INTERFACE
  
  INTERFACE direction   ! for internal use
    MODULE PROCEDURE t_octree_direction
  END INTERFACE

  INTERFACE otree_printOctree
    MODULE PROCEDURE t_octree_print
  END INTERFACE

  INTERFACE otree_infoOctree
    MODULE PROCEDURE t_octree_info
  END INTERFACE
  
  INTERFACE otree_getsize
    MODULE PROCEDURE t_octree_getsize
  END INTERFACE

  INTERFACE otree_getBoundingBox
    MODULE PROCEDURE t_octree_getboundingbox
  END INTERFACE

  INTERFACE otree_getX
    MODULE PROCEDURE t_octree_getxvalue
  END INTERFACE

  INTERFACE otree_getY
    MODULE PROCEDURE t_octree_getyvalue
  END INTERFACE

  INTERFACE otree_getZ
    MODULE PROCEDURE t_octree_getzvalue
  END INTERFACE

  INTERFACE resizeNVT
    MODULE PROCEDURE t_octree_resize_nvt
  END INTERFACE

  INTERFACE resizeNNODE
    MODULE PROCEDURE t_octree_resize_nnode
  END INTERFACE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_octree_create(roctree,nnvt,nnnode,xmin,ymin,zmin,xmax,ymax,zmax,dfactor)
  
!<description>
    ! This subroutine creates a new octree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the octree
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: nnvt

    ! Total number of nodes that should be stored in the octree
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: nnnode

    ! Dimensions of the initial bounding box
    REAL(DP), INTENT(IN) :: xmin,ymin,zmin,xmax,ymax,zmax

    ! OPTIONAL: Factor by which the octree should be enlarged if
    ! new storage has to be allocated
    REAL(DP), INTENT(IN), OPTIONAL :: dfactor
!</input>

!<output>
    ! Octree structure
    TYPE(t_octree), INTENT(OUT) :: roctree
!</output>
!</subroutine>
    
    ! Set factor
    IF (PRESENT(dfactor)) THEN
      IF (dfactor > 1_DP) roctree%dfactor=dfactor
    END IF
    
    ! Set values
    roctree%NNNODE  = nnnode
    roctree%NNVT    = nnvt
    roctree%NNODE   = 0
    roctree%NVT     = 0
    roctree%NRESIZE = 0
    
    ! Allocate memory and associate pointers
    CALL storage_new('t_octree_create', 'p_Ddata', (/3,nnvt/),&
        ST_DOUBLE, roctree%h_Ddata, ST_NEWBLOCK_ZERO)
    CALL storage_new('t_octree_create', 'p_Dbbox', (/6,nnnode/),&
        ST_DOUBLE, roctree%h_Dbbox, ST_NEWBLOCK_ZERO)
    CALL storage_new('t_octree_create', 'p_Knode', (/11,nnnode/),&
        ST_INT, roctree%h_Knode, ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double2D(roctree%h_Ddata,roctree%p_Ddata)
    CALL storage_getbase_double2D(roctree%h_Dbbox,roctree%p_Dbbox)
    CALL storage_getbase_int2D(roctree%h_Knode,   roctree%p_Knode)

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
  END SUBROUTINE t_octree_create

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_octree_release(roctree)

!<description>
    ! This subroutine releases an existing octree
!</description>

!<inputoutput>
    ! quadtree
    TYPE(t_octree), INTENT(INOUT) :: roctree
!</inputoutput>
!</subroutine>
    
    ! Free memory
    IF (roctree%h_Ddata /= ST_NOHANDLE) CALL storage_free(roctree%h_Ddata)
    IF (roctree%h_Dbbox /= ST_NOHANDLE) CALL storage_free(roctree%h_Dbbox)
    IF (roctree%h_Knode /= ST_NOHANDLE) CALL storage_free(roctree%h_Knode)
    NULLIFY(roctree%p_Knode,roctree%p_Dbbox,roctree%p_Ddata)

    ! Reset values
    roctree%NNNODE = 0
    roctree%NNVT   = 0
    roctree%NNODE  = 0
    roctree%NVT    = 0
    roctree%NRESIZE= 0
  END SUBROUTINE t_octree_release
  
  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_octree_resize_nvt(roctree,nnvt)

!<description>
    ! This subroutine reallocates memory for an existing octree
!</description>

!<input>
    ! New number of vertices that should be stored in the octree
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: nnvt
!</input>

!<inputoutput>
    ! Octree that should be resized
    TYPE(t_octree) :: roctree
!</inputoutput>
!</subroutine>

    CALL storage_realloc('t_octree_resize_nvt', nnvt,roctree%h_Ddata,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_double2D(roctree%h_Ddata,roctree%p_Ddata)
    
    roctree%NNVT   = nnvt
    roctree%NRESIZE=roctree%NRESIZE+1
  END SUBROUTINE t_octree_resize_nvt

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_octree_resize_nnode(roctree,nnnode)

!<description>
    ! This subroutine reallocates memory for an existing octree
!</description>

!<input>
    ! New number of nodes that should be stored in the octree
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: nnnode
!</input>

!<inputoutput>
    ! Octree that should be resized
    TYPE(t_octree) :: roctree
!</inputoutput>
!</subroutine>

    CALL storage_realloc('t_octree_resize_nnode',nnnode,roctree%h_Dbbox,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_realloc('t_octree_resize_nnode',nnnode,roctree%h_Knode,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_double2D(roctree%h_Dbbox,roctree%p_Dbbox)
    CALL storage_getbase_int2D(roctree%h_Knode,   roctree%p_Knode)
    
    roctree%NNNODE = nnnode
    roctree%NRESIZE=roctree%NRESIZE+1
  END SUBROUTINE t_octree_resize_nnode

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_octree_copyfrom_handle(roctree,h_Ddata)

!<description>
    ! This subroutine copies the content of the octree to a handle.
    ! If the handle is not associated, then a new handle with correct
    ! size is allocated. Otherwise, the handle is reallocated if 
    ! it does not provide enough memory.
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree
!</input>

!<inputoutput>
    ! Handle to the coordinate vector
    INTEGER, INTENT(INOUT) :: h_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata,p_DdataTmp
    INTEGER(PREC_OTREEIDX), DIMENSION(2) :: Isize
    
    ! Check if handle is associated
    IF (h_Ddata == ST_NOHANDLE) THEN
      Isize = (/3,roctree%NVT/)
      CALL storage_new('t_octree_copyfrom_handle','p_Ddata',&
          Isize, ST_DOUBLE, h_Ddata, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Ddata,Isize)
      IF (Isize(2) < roctree%NVT) THEN
        CALL storage_realloc('t_octree_copyfrom_handle',&
            roctree%NVT, h_Ddata, ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Set pointers
    CALL storage_getbase_double2D(h_Ddata,p_Ddata)
    CALL storage_getbase_double2D(roctree%h_Ddata,p_DdataTmp)

    ! Copy data
    CALL DCOPY(3*roctree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  END SUBROUTINE t_octree_copyfrom_handle

  !************************************************************************

!<subroutine>

  SUBROUTINE t_octree_copyfrom_array(roctree,p_Ddata)

!<description>
    ! This subroutine copies the content of the octree to an array.
    ! If the array is too small, then this subroutines terminates.
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree
!</input>

!<inputoutput>
    ! Coordinate vector
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DdataTmp
    INTEGER(PREC_OTREEIDX), DIMENSION(2) :: Isize

    ! Check size of array
    Isize = SHAPE(p_Ddata)
    IF (Isize(1) /= 3 .OR. Isize(2) < roctree%NVT) THEN
      PRINT *, "t_octree_copyfrom_array: Array too small!"
      CALL sys_halt()
    END IF

    ! Set pointers
    CALL storage_getbase_double2D(roctree%h_Ddata,p_DdataTmp)

    ! Copy data
    CALL DCOPY(3*roctree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  END SUBROUTINE t_octree_copyfrom_array

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_octree_copyto_handle(h_Ddata,roctree)

!<description>
    ! This subroutine copies the content of a handle to the octree.
!</description>

!<input>
    ! Handle to the coordinate vector
    INTEGER, INTENT(IN) :: h_Ddata
!</input>

!<inputoutput>
    ! Octree
    TYPE(t_octree), INTENT(INOUT) :: roctree
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata
    INTEGER(PREC_OTREEIDX) :: ivt,jvt,inode,ipos
    
    ! Set pointer and check its shape
    CALL storage_getbase_double2D(h_Ddata,p_Ddata)
    IF (SIZE(p_Ddata,1) /= 3) THEN
      PRINT *, "t_octree_copyto_handle: First dimension of array must be 3!"
      CALL sys_halt()
    END IF
    
    DO ivt=1,SIZE(p_Ddata,2)
      IF (search(roctree,p_Ddata(:,ivt),inode,ipos,jvt) == OTREE_NOT_FOUND)&
          CALL insert(roctree,ivt,p_Ddata(:,ivt),inode)
    END DO
  END SUBROUTINE t_octree_copyto_handle

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_octree_copyto_array(p_Ddata,roctree)

!<description>
    ! This subroutine copies the content of an array to the octree.
!</description>

!<input>
    ! Coordinate vector
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: p_Ddata
!</input>

!<inputoutput>
    ! Octree
    TYPE(t_octree), INTENT(INOUT) :: roctree
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_OTREEIDX) :: ivt,jvt,inode,ipos
    
    IF (SIZE(p_Ddata,1) /= 3) THEN
      PRINT *, "t_octree_copyto_handle: First dimension of array must be 3!"
      CALL sys_halt()
    END IF

    DO ivt=1,SIZE(p_Ddata,2)
      IF (search(roctree,p_Ddata(:,ivt),inode,ipos,jvt) == OTREE_NOT_FOUND)&
          CALL insert(roctree,ivt,p_Ddata(:,ivt),inode)
    END DO
  END SUBROUTINE t_octree_copyto_array

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_octree_insert(roctree,ivt,Ddata,inode)

!<description>
    ! This subroutine inserts a new coordinate item into the octree
!</description>

!<input>
    ! Number of the inserted vertex
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: ivt

    ! Number of the node into which vertex is inserted
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: inode
    
    ! Coordinates of the new vertex
    REAL(DP), DIMENSION(3), INTENT(IN) :: Ddata
!</input>

!<inputoutput>
    ! Octree
    TYPE(t_octree), INTENT(INOUT)      :: roctree
!</inputoutput>
!</subroutine>
    
    ! Check if there is enough space left in the vertex component of the octree
    IF (roctree%NVT == roctree%NNVT)&
        CALL resizeNVT(roctree,CEILING(roctree%dfactor*roctree%NNVT))
    
    ! Update values and add new entry recursively
    roctree%NVT            = roctree%NVT+1
    roctree%p_Ddata(:,ivt) = Ddata
    CALL insert(ivt,inode)
    
  CONTAINS  
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    RECURSIVE SUBROUTINE insert(ivt,inode)
      INTEGER(PREC_OTREEIDX), INTENT(IN) :: ivt,inode
      REAL(DP) :: xmin,ymin,zmin,xmax,ymax,zmax,xmid,ymid,zmid
      INTEGER(PREC_OTREEIDX) :: i,jnode,nnode
      
      IF (roctree%p_Knode(OTREE_STATUS,inode) == OTREE_MAX) THEN
        
        IF (roctree%nnode+OTREE_MAX > roctree%nnnode)&
            CALL resizeNNODE(roctree,CEILING(roctree%dfactor*roctree%NNNODE))
        
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
        roctree%p_Dbbox(:,nnode+OTREE_NWF) = (/xmin,ymid,zmin,xmid,ymax,zmid/)
        
        ! SWF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWF) = OTREE_SWF
        roctree%p_Dbbox(:,nnode+OTREE_SWF) = (/xmin,ymin,zmin,xmid,ymid,zmid/)
        
        ! SEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEF) = OTREE_SEF
        roctree%p_Dbbox(:,nnode+OTREE_SEF) = (/xmid,ymin,zmin,xmax,ymid,zmid/)
        
        ! NEF-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEF) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEF) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEF) = OTREE_NEF
        roctree%p_Dbbox(:,nnode+OTREE_NEF) = (/xmid,ymid,zmin,xmax,ymax,zmid/)
        
        ! NWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NWB) = OTREE_NWB
        roctree%p_Dbbox(:,nnode+OTREE_NWB) = (/xmin,ymid,zmid,xmid,ymax,zmax/)
        
        ! SWB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SWB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SWB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SWB) = OTREE_SWB
        roctree%p_Dbbox(:,nnode+OTREE_SWB) = (/xmin,ymin,zmid,xmid,ymid,zmax/)
        
        ! SEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_SEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_SEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_SEB) = OTREE_SEB
        roctree%p_Dbbox(:,nnode+OTREE_SEB) = (/xmid,ymin,zmid,xmax,ymid,zmax/)
        
        ! NEB-node
        roctree%p_Knode(OTREE_STATUS,nnode+OTREE_NEB) = OTREE_EMPTY
        roctree%p_Knode(OTREE_PARENT,nnode+OTREE_NEB) = inode
        roctree%p_Knode(OTREE_PARPOS,nnode+OTREE_NEB) = OTREE_NEB
        roctree%p_Dbbox(:,nnode+OTREE_NEB) = (/xmid,ymid,zmid,xmax,ymax,zmax/)

        ! Add the eight values from INODE to the eight new nodes
        ! NNODE+1:NNODE+8 recursively
        DO i=1,OTREE_MAX
          jnode=nnode+direction(roctree,roctree%p_Ddata(:,roctree%p_Knode(i,inode)),inode)
          CALL insert(roctree%p_Knode(i,inode),jnode)
        END DO
        
        ! Mark the current nodes as subdivided and set pointers to its eight children 
        roctree%p_Knode(OTREE_STATUS,inode) = OTREE_SUBDIV
        roctree%p_Knode(1:OTREE_MAX, inode) = (/nnode+OTREE_NWF,nnode+OTREE_SWF,&
            nnode+OTREE_SEF,nnode+OTREE_NEF,nnode+OTREE_NWB,nnode+OTREE_SWB,&
            nnode+OTREE_SEB,nnode+OTREE_NEB/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+direction(roctree,roctree%p_Ddata(:,ivt),inode)
        CALL insert(ivt,jnode)
        
      ELSE
        
        ! Node is not full, so new items can be stored
        roctree%p_Knode(OTREE_STATUS,inode) = roctree%p_Knode(OTREE_STATUS,inode)+1
        roctree%p_Knode(roctree%p_Knode(OTREE_STATUS,inode),inode) = ivt
      END IF
    END SUBROUTINE insert
  END SUBROUTINE t_octree_insert
  
  ! ***************************************************************************
  
!<function>
  
  FUNCTION t_octree_delete(roctree,Ddata,ivt) RESULT(f)

!<description>
    ! This function deletes an item from the octree
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    REAL(DP), DIMENSION(3), INTENT(IN) :: Ddata
!</input>

!<inputoutput>
    ! Octree
    TYPE(t_octree), INTENT(INOUT) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that is deleted
    INTEGER(PREC_OTREEIDX), INTENT(OUT) :: ivt
!</output>

!<result>
    ! Result of the deletion: OTREE_NOT_FOUND / OTREE_FOUND
    INTEGER :: f
!</result>
!</function>
    
    ! local variables
    INTEGER :: inode,ipos,jpos,jvt
    
    ! Search for the given coordinates
    f=search(roctree,Ddata,inode,ipos,ivt)
    
    ! What can we do from the searching
    IF (f == OTREE_FOUND) THEN
      
      ! Remove item IVT from node INODE
      DO jpos=ipos+1,roctree%p_Knode(OTREE_STATUS,inode)
        roctree%p_Knode(jpos-1,inode) = roctree%p_Knode(jpos,inode)
      END DO
      roctree%p_Knode(roctree%p_Knode(OTREE_STATUS,inode),inode) = 0
      roctree%p_Knode(OTREE_STATUS,inode) = roctree%p_Knode(OTREE_STATUS,inode)-1
      
      ! If IVT is not last item move last item NVT to position IVT
      IF (ivt /= roctree%NVT) THEN
        IF (search(roctree,roctree%p_Ddata(1:3,roctree%NVT),inode,ipos,jvt) == OTREE_FOUND) THEN
          roctree%p_Ddata(:,ivt) = roctree%p_Ddata(:,roctree%NVT)
          roctree%p_Knode(ipos,inode) = ivt
        END IF
        ivt=roctree%NVT
      END IF
      roctree%NVT = roctree%NVT-1
    END IF
  END FUNCTION t_octree_delete
  
  ! ***************************************************************************

!<function>
  
  FUNCTION t_octree_search(roctree,Ddata,inode,ipos,ivt) RESULT(f)

!<description>
    ! This subroutine searches for given coordinates in the octree
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree
    
    ! Coordinates that should be searched for
    REAL(DP), DIMENSION(3), INTENT(IN) :: Ddata
!</input>

!<output>
    ! Number of the node in which the given coordinates are
    INTEGER(PREC_OTREEIDX), INTENT(OUT) :: inode

    ! Position of the coordinates in the node
    INTEGER(PREC_OTREEIDX), INTENT(OUT) :: ipos

    ! Number of the vertex the coordinates correspond to
    INTEGER(PREC_OTREEIDX), INTENT(OUT) :: ivt
!</output>

!<result>
    ! Result of the searching: OTREE_NOT_FOUND / OTREE_FOUND
    INTEGER :: f
!</result>
!</function>
    
    ! Initialize
    inode=1; ipos=1; ivt=1; f=search(inode,ipos,ivt)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive searching routine follows

    RECURSIVE FUNCTION search(inode,ipos,ivt) RESULT(f)
      INTEGER(PREC_OTREEIDX), INTENT(INOUT) :: inode,ipos,ivt
      INTEGER :: f
      
      f=OTREE_NOT_FOUND
      IF (roctree%p_Knode(OTREE_STATUS,inode) == OTREE_SUBDIV) THEN
        
        ! Node is subdivided. Compute child INODE which to look recursively.
        inode = roctree%p_Knode(direction(roctree,Ddata,inode),inode)
        f=search(inode,ipos,ivt)
        
      ELSE
        
        ! Node is not subdivided. Search for (x,y,z) in current node
        DO ipos=1,roctree%p_Knode(OTREE_STATUS,inode)
          ivt = roctree%p_Knode(ipos,inode)
          IF (SQRT(SUM((roctree%p_Ddata(:,ivt)-Ddata)**2)) <= SYS_EPSREAL) THEN
            f=OTREE_FOUND
            RETURN
          END IF
        END DO
        
      END IF
    END FUNCTION search
  END FUNCTION t_octree_search

  !************************************************************************
  
!<function>
  
  PURE FUNCTION t_octree_direction(roctree,Ddata,inode) RESULT(d)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Ddata
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree
    
    ! Coordinates
    REAL(DP), DIMENSION(3), INTENT(IN) :: Ddata

    ! Number of node
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: inode
!</input>

!<result>
    ! Further search direction
    INTEGER :: d
!</result>
!</function>
    
    ! local variables
    REAL(DP) :: xmid,ymid,zmid

    ! Compute midpoint of current node
    xmid=0.5*(roctree%p_Dbbox(OTREE_XMIN,inode)+&
              roctree%p_Dbbox(OTREE_XMAX,inode))
    ymid=0.5*(roctree%p_Dbbox(OTREE_YMIN,inode)+&
              roctree%p_Dbbox(OTREE_YMAX,inode))
    zmid=0.5*(roctree%p_Dbbox(OTREE_ZMIN,inode)+&
              roctree%p_Dbbox(OTREE_ZMAX,inode))

    ! Do we have to look in the front or back?
    IF (Ddata(3) < zmid) THEN
      
      IF (Ddata(1) > xmid) THEN
        IF (Ddata(2) > ymid) THEN
          d=OTREE_NEF; RETURN
        ELSE
          d=OTREE_SEF; RETURN
        END IF
      ELSE
        IF (Ddata(2) > ymid) THEN
          d=OTREE_NWF; RETURN
        ELSE
          d=OTREE_SWF; RETURN
        END IF
      END IF

    ELSE

      IF (Ddata(1) > xmid) THEN
        IF (Ddata(2) > ymid) THEN
          d=OTREE_NEB; RETURN
        ELSE
          d=OTREE_SEB; RETURN
        END IF
      ELSE
        IF (Ddata(2) > ymid) THEN
          d=OTREE_NWB; RETURN
        ELSE
          d=OTREE_SWB; RETURN
        END IF
      END IF

    END IF
  END FUNCTION t_octree_direction

  !************************************************************************
  
!<subroutine>

  SUBROUTINE t_octree_print(roctree,cfilename)

!<description>
    ! This subroutine writes the content of the octree to a file
    ! which can be visualized by means of Matlab
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree

    ! filename of the output file
    CHARACTER(LEN=*), INTENT(IN) :: cfilename
!</input>
!</subroutine>
    
    ! local variables
    REAL(DP) :: xmin,xmax,ymin,ymax,zmin,zmax
    INTEGER :: iunit
    
    iunit=sys_getFreeUnit()
    OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(cfilename)))
    xmin = roctree%p_Dbbox(OTREE_XMIN,1)
    ymin = roctree%p_Dbbox(OTREE_YMIN,1)
    zmin = roctree%p_Dbbox(OTREE_ZMIN,1)
    xmax = roctree%p_Dbbox(OTREE_XMAX,1)
    ymax = roctree%p_Dbbox(OTREE_YMAX,1)
    zmax = roctree%p_Dbbox(OTREE_ZMAX,1)
    CALL print(xmin,ymin,zmin,xmax,ymax,zmax,1)
    CLOSE(UNIT=iunit)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive print routine follows
    
    RECURSIVE SUBROUTINE print(xmin,ymin,zmin,xmax,ymax,zmax,inode)
      REAL(DP), INTENT(IN) :: xmin,ymin,zmin,xmax,ymax,zmax
      INTEGER(PREC_OTREEIDX), INTENT(IN) :: inode
      REAL(DP)               :: xmid,ymid,zmid
      INTEGER(PREC_OTREEIDX) :: i

      WRITE(UNIT=iunit,FMT=*) 'hex'
      WRITE(UNIT=iunit,FMT=10) xmin,ymin,zmin,xmax,ymax,zmax
      
      IF (roctree%p_Knode(OTREE_STATUS,inode) == OTREE_SUBDIV) THEN
        
        xmid=0.5*(xmin+xmax)
        ymid=0.5*(ymin+ymax)
        zmid=0.5*(zmin+zmax)
        
        CALL print(xmin,ymid,zmin,xmid,ymax,zmid,roctree%p_Knode(OTREE_NWF,inode))
        CALL print(xmin,ymin,zmin,xmid,ymid,zmid,roctree%p_Knode(OTREE_SWF,inode))
        CALL print(xmid,ymin,zmin,xmax,ymid,zmid,roctree%p_Knode(OTREE_SEF,inode))
        CALL print(xmid,ymid,zmin,xmax,ymax,zmid,roctree%p_Knode(OTREE_NEF,inode))

        CALL print(xmin,ymid,zmid,xmid,ymax,zmax,roctree%p_Knode(OTREE_NWB,inode))
        CALL print(xmin,ymin,zmid,xmid,ymid,zmax,roctree%p_Knode(OTREE_SWB,inode))
        CALL print(xmid,ymin,zmid,xmax,ymid,zmax,roctree%p_Knode(OTREE_SEB,inode))
        CALL print(xmid,ymid,zmid,xmax,ymax,zmax,roctree%p_Knode(OTREE_NEB,inode))

      ELSEIF (roctree%p_Knode(OTREE_STATUS,inode) > OTREE_EMPTY) THEN
        
        DO i=1,roctree%p_Knode(OTREE_STATUS,inode)
          WRITE(UNIT=iunit,FMT=*) 'node'
          WRITE(UNIT=iunit,FMT=20) roctree%p_Ddata(:,roctree%p_Knode(i,inode))
        END DO
        
      END IF
      
10    FORMAT(6E15.6E3)
20    FORMAT(3E15.6E3)
    END SUBROUTINE print
  END SUBROUTINE t_octree_print

  !************************************************************************

!<subroutine>

  SUBROUTINE t_octree_info(roctree)

!<description>
    ! This subroutine outputs statistical info about the octree
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree
!</input>
!</subroutine>

    WRITE(*,FMT=*) ' Octree:'
    WRITE(*,FMT=*) ' ======='
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Ddata =',roctree%h_Ddata
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Dbbox =',roctree%h_Dbbox
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Knode =',roctree%h_Knode
    WRITE(*,*)
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8)') '  NVT     =',roctree%NVT,'NNVT     =',roctree%NNVT
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8)') '  NNODE   =',roctree%NNODE,'NNNODE   =',roctree%NNNODE
    WRITE(*,FMT='(1X,A,1X,I8,A,4X,F5.1,A,4X,F5.1,A)') '  NRESIZE =',roctree%NRESIZE, "   NODES    =", &
        100*roctree%NNODE/REAL(roctree%NNNODE,DP),'%  FILLING  =',100*roctree%NVT/REAL(roctree%NNVT,DP),'%'
    WRITE(*,*)
  END SUBROUTINE t_octree_info

  !************************************************************************

!<function>

  PURE FUNCTION t_octree_getsize(roctree) RESULT(nvt)

!<description>
    ! This function returns the number of vertices stored in the octree
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree
!</input>

!<result>
    ! Number of vertices in octree
    INTEGER(PREC_OTREEIDX) :: nvt
!</result>
!</function>

    nvt=roctree%NVT
  END FUNCTION t_octree_getsize

  !************************************************************************

!<function>

  FUNCTION t_octree_getboundingbox(roctree,inode) RESULT(bbox)
    
!<description>
    ! This function returns the bounding box of the specified node.
    ! If no node number is given, then the outer bounding box is returned.
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree

    ! OPTIONAL: number of node for which bounding box should be returned
    INTEGER, INTENT(IN), OPTIONAL :: inode
!</input>

!<result>
    ! bounding box
    REAL(DP), DIMENSION(6) :: bbox
!</result>
!</function>
    
    IF (PRESENT(inode)) THEN
      IF (inode > roctree%NVT) THEN
        PRINT *, "t_octree_getboundingbox: node number exceeds octree dimension"
        CALL sys_halt()
      END IF
      bbox=roctree%p_Dbbox(OTREE_XMIN:OTREE_ZMAX,inode)
    ELSE
      bbox=roctree%p_Dbbox(OTREE_XMIN:OTREE_ZMAX,1)
    END IF
  END FUNCTION t_octree_getboundingbox

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION t_octree_getXvalue(roctree,ivt) RESULT(x)

!<description>
    ! This function returns the X-value at the given position.
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree

    ! position in the octree
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: ivt
!</input>

!<result>
    REAL(DP) :: x
!</result>
!</function>

    x=roctree%p_Ddata(1,ivt)
  END FUNCTION t_octree_getXvalue

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION t_octree_getYvalue(roctree,ivt) RESULT(y)

!<description>
    ! This function returns the Y-value at the given position.
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree

    ! position in the octree
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: ivt
!</input>

!<result>
    REAL(DP) :: y
!</result>
!</function>

    y=roctree%p_Ddata(2,ivt)
  END FUNCTION t_octree_getYvalue

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION t_octree_getZvalue(roctree,ivt) RESULT(z)

!<description>
    ! This function returns the Z-value at the given position.
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree

    ! position in the octree
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: ivt
!</input>

!<result>
    REAL(DP) :: z
!</result>
!</function>

    z=roctree%p_Ddata(3,ivt)
  END FUNCTION t_octree_getZvalue
END MODULE octree

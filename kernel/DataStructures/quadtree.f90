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
!# 9.) qtree_printQuadtree
!#     -> Write quadtree to file
!#
!# 10.) qtree_infoQuadtree
!#      -> Output info about quadtree
!#
!# 11.) qtree_getsize
!#      -> Return number of vertices in quadtree
!#
!# 12.) qtree_getBoundingBox
!#      -> Return the outer bounding box
!#
!# 13.) qtree_getX
!#      -> Return the X-value at a given position
!#
!# 14.) qtree_getY
!#      -> Return the Y-value at a given position
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

MODULE quadtree
  USE fsystem
  USE storage
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: t_quadtree
  PUBLIC :: qtree_createQuadtree
  PUBLIC :: qtree_releaseQuadtree
  PUBLIC :: qtree_copyToQuadtree
  PUBLIC :: qtree_copyFromQuadtree
  PUBLIC :: qtree_insertIntoQuadtree
  PUBLIC :: qtree_deleteFromQuadtree
  PUBLIC :: qtree_searchInQuadtree
  PUBLIC :: qtree_printQuadtree
  PUBLIC :: qtree_infoQuadtree
  PUBLIC :: qtree_getsize
  PUBLIC :: qtree_getBoundingBox
  PUBLIC :: qtree_getX
  PUBLIC :: qtree_getY

!<constants>

!<constantblock description="KIND values for quadtree data">
  
  ! Kind value for indices in quadtree
  INTEGER, PARAMETER, PUBLIC :: PREC_QTREEIDX = I32

!</constantblock>

!<constantblock description="Constants for quadtree structure">
  
  ! Position of the status information
  INTEGER, PARAMETER :: QTREE_STATUS = 7
  
  ! Position of the parent information
  INTEGER, PARAMETER :: QTREE_PARENT = 6

  ! Position of the "position" of the parent information
  INTEGER, PARAMETER :: QTREE_PARPOS = 5

  ! Number of free positions in quad
  INTEGER, PARAMETER :: QTREE_FREE   = 5

  ! Item in "North-East" position
  INTEGER, PARAMETER :: QTREE_NE     = 4

  ! Item in "South-East" position
  INTEGER, PARAMETER :: QTREE_SE     = 3

  ! Item in "South-West" position
  INTEGER, PARAMETER :: QTREE_SW     = 2

  ! Item in "North-West" position
  INTEGER, PARAMETER :: QTREE_NW     = 1

  ! Identifier: Quad is empty
  INTEGER, PARAMETER :: QTREE_EMPTY  =  0

  ! Identifier: Status is subdivided
  INTEGER, PARAMETER :: QTREE_SUBDIV = -1
  
  ! Maximum number of items for each quad
  INTEGER, PARAMETER :: QTREE_MAX    =  4

!</constantblock> 

!<constantblock description="Constants for quadtree bounding-box">

  ! Position of the x-minimal value
  INTEGER, PARAMETER :: QTREE_XMIN   =  1

  ! Position of the y-minimal value
  INTEGER, PARAMETER :: QTREE_YMIN   =  2

  ! Position of the x-maximal value
  INTEGER, PARAMETER :: QTREE_XMAX   =  3

  ! Position of the y-maximal value
  INTEGER, PARAMETER :: QTREE_YMAX   =  4

!</constantblock>   
  
!<constantblock description="Constants for quadtree operations">

  ! Item could not be found in the quadtree
  INTEGER, PARAMETER, PUBLIC :: QTREE_NOT_FOUND = -1

  ! Item could be found in the quadtree
  INTEGER, PARAMETER, PUBLIC :: QTREE_FOUND     =  0
!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! A linear quadtree implemented as array

  TYPE t_quadtree

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

    ! Quadtree structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_QTREEIDX), DIMENSION(:,:), POINTER :: p_Knode

    ! Number of vertices currently stored in the quadtree
    INTEGER(PREC_QTREEIDX) :: NVT    = 0

    ! Total number of vertices that can be stored  in the quadtree
    INTEGER(PREC_QTREEIDX) :: NNVT   = 0

    ! Number of subdivisions currently store in the quadtree
    INTEGER(PREC_QTREEIDX) :: NNODE  = 0

    ! Total number of subdivision that can be stored in the quadtree
    INTEGER(PREC_QTREEIDX) :: NNNODE = 0

    ! Total number of resize operations
    INTEGER :: NRESIZE            = 0

    ! Factor by which the quadtree is enlarged if new storage is allocate
    REAL(DP) :: dfactor           = 1.5_DP
  END TYPE t_quadtree
  
!</typeblock>
!</types>

  INTERFACE qtree_createQuadtree
    MODULE PROCEDURE t_quadtree_create
  END INTERFACE
  
  INTERFACE qtree_releaseQuadtree
    MODULE PROCEDURE t_quadtree_release
  END INTERFACE
  
  INTERFACE qtree_copyToQuadtree
    MODULE PROCEDURE t_quadtree_copyto_handle
    MODULE PROCEDURE t_quadtree_copyto_array
  END INTERFACE

  INTERFACE qtree_copyFromQuadtree
    MODULE PROCEDURE t_quadtree_copyfrom_handle
    MODULE PROCEDURE t_quadtree_copyfrom_array
  END INTERFACE
  
  INTERFACE qtree_insertIntoQuadtree
    MODULE PROCEDURE t_quadtree_insert
  END INTERFACE
  INTERFACE insert   ! for internal use
    MODULE PROCEDURE t_quadtree_insert
  END INTERFACE

  INTERFACE qtree_deleteFromQuadtree
    MODULE PROCEDURE t_quadtree_delete
  END INTERFACE
  INTERFACE delete   ! for internal use
    MODULE PROCEDURE t_quadtree_delete
  END INTERFACE
  
  INTERFACE qtree_searchInQuadtree
    MODULE PROCEDURE t_quadtree_search
  END INTERFACE
  INTERFACE search   ! for internal use
    MODULE PROCEDURE t_quadtree_search
  END INTERFACE
  
  INTERFACE direction   ! for internal use
    MODULE PROCEDURE t_quadtree_direction
  END INTERFACE

  INTERFACE qtree_printQuadtree
    MODULE PROCEDURE t_quadtree_print
  END INTERFACE

  INTERFACE qtree_infoQuadtree
    MODULE PROCEDURE t_quadtree_info
  END INTERFACE
  
  INTERFACE qtree_getsize
    MODULE PROCEDURE t_quadtree_getsize
  END INTERFACE

  INTERFACE qtree_getBoundingBox
    MODULE PROCEDURE t_quadtree_getboundingbox
  END INTERFACE

  INTERFACE qtree_getX
    MODULE PROCEDURE t_quadtree_getxvalue
  END INTERFACE

  INTERFACE qtree_getY
    MODULE PROCEDURE t_quadtree_getyvalue
  END INTERFACE

  INTERFACE resizeNVT
    MODULE PROCEDURE t_quadtree_resize_nvt
  END INTERFACE

  INTERFACE resizeNNODE
    MODULE PROCEDURE t_quadtree_resize_nnode
  END INTERFACE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_quadtree_create(rquadtree,nnvt,nnnode,xmin,ymin,xmax,ymax,dfactor)
  
!<description>
    ! This subroutine creates a new quadtree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the quadtree
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: nnvt

    ! Total number of subdivisions that should be stored in the quadtree
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: nnnode

    ! Dimensions of the initial bounding box
    REAL(DP), INTENT(IN) :: xmin,ymin,xmax,ymax

    ! OPTIONAL: Factor by which the quadtree should be enlarged if
    ! new storage has to be allocated
    REAL(DP), INTENT(IN), OPTIONAL :: dfactor
!</input>

!<output>
    ! Quadtree structure
    TYPE(t_quadtree), INTENT(OUT) :: rquadtree
!</output>
!</subroutine>
    
    ! Set factor
    IF (PRESENT(dfactor)) THEN
      IF (dfactor > 1_DP) rquadtree%dfactor=dfactor
    END IF
    
    ! Set values
    rquadtree%NNNODE = nnnode
    rquadtree%NNVT   = nnvt
    rquadtree%NNODE  = 0
    rquadtree%NVT    = 0
    rquadtree%NRESIZE= 0
    
    ! Allocate memory and associate pointers
    CALL storage_new('t_quadtree_create','p_Ddata',&
        (/2,nnvt/),ST_DOUBLE,rquadtree%h_Ddata,ST_NEWBLOCK_ZERO)
    CALL storage_new('t_quadtree_create','p_Dbbox',&
        (/4,nnnode/),ST_DOUBLE,rquadtree%h_Dbbox,ST_NEWBLOCK_ZERO)
    CALL storage_new('t_quadtree_create','p_Knode',&
        (/7,nnnode/),ST_INT,rquadtree%h_Knode,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double2D(rquadtree%h_Ddata,rquadtree%p_Ddata)
    CALL storage_getbase_double2D(rquadtree%h_Dbbox,rquadtree%p_Dbbox)
    CALL storage_getbase_int2D(rquadtree%h_Knode,   rquadtree%p_Knode)

    ! Initialize first quad
    rquadtree%nnode=1
    rquadtree%p_Knode(QTREE_STATUS,1)=QTREE_EMPTY
    rquadtree%p_Knode(QTREE_PARENT,1)=QTREE_EMPTY
    rquadtree%p_Knode(QTREE_FREE,1)=1
    
    ! Set co-ordinates of the bounding box
    rquadtree%p_Dbbox(QTREE_XMIN,1) = xmin
    rquadtree%p_Dbbox(QTREE_YMIN,1) = ymin
    rquadtree%p_Dbbox(QTREE_XMAX,1) = xmax
    rquadtree%p_Dbbox(QTREE_YMAX,1) = ymax
  END SUBROUTINE t_quadtree_create

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_quadtree_release(rquadtree)

!<description>
    ! This subroutine releases an existing quadtree
!</description>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! Free memory
    IF (rquadtree%h_Ddata /= ST_NOHANDLE) CALL storage_free(rquadtree%h_Ddata)
    IF (rquadtree%h_Dbbox /= ST_NOHANDLE) CALL storage_free(rquadtree%h_Dbbox)
    IF (rquadtree%h_Knode /= ST_NOHANDLE) CALL storage_free(rquadtree%h_Knode)
    NULLIFY(rquadtree%p_Knode,rquadtree%p_Dbbox,rquadtree%p_Ddata)

    ! Reset values
    rquadtree%NNNODE = 0
    rquadtree%NNVT   = 0
    rquadtree%NNODE  = 0
    rquadtree%NVT    = 0
    rquadtree%NRESIZE= 0
  END SUBROUTINE t_quadtree_release
  
  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_quadtree_resize_nvt(rquadtree,nnvt)

!<description>
    ! This subroutine reallocates memory for an existing quadtree
!</description>

!<input>
    ! New number of vertices that should be stored in the quadtree
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: nnvt
!</input>

!<inputoutput>
    ! quadtree that should be resized
    TYPE(t_quadtree) :: rquadtree
!</inputoutput>
!</subroutine>

    CALL storage_realloc('t_quadtree_resize_nvt',&
        nnvt,rquadtree%h_Ddata,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_double2D(rquadtree%h_Ddata,rquadtree%p_Ddata)
    
    rquadtree%NNVT   = nnvt
    rquadtree%NRESIZE=rquadtree%NRESIZE+1
  END SUBROUTINE t_quadtree_resize_nvt

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_quadtree_resize_nnode(rquadtree,nnnode)

!<description>
    ! This subroutine reallocates memory for an existing quadtree
!</description>

!<input>
    ! New number of quads that should be stored in the quadtree
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: nnnode
!</input>

!<inputoutput>
    ! quadtree that should be resized
    TYPE(t_quadtree) :: rquadtree
!</inputoutput>
!</subroutine>

    CALL storage_realloc('t_quadtree_resize_nnode',&
        nnnode,rquadtree%h_Dbbox,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_realloc('t_quadtree_resize_nnode',&
        nnnode,rquadtree%h_Knode,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_double2D(rquadtree%h_Dbbox,rquadtree%p_Dbbox)
    CALL storage_getbase_int2D(rquadtree%h_Knode,   rquadtree%p_Knode)
    
    rquadtree%NNNODE = nnnode
    rquadtree%NRESIZE=rquadtree%NRESIZE+1
  END SUBROUTINE t_quadtree_resize_nnode

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_quadtree_copyfrom_handle(rquadtree,h_Ddata)

!<description>
    ! This subroutine copies the content of the quadtree to a handle.
    ! If the handle is not associated, then a new handle with correct
    ! size is allocated. Otherwise, the handle is reallocated if 
    ! it does not provide enough memory.
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
!</input>

!<inputoutput>
    ! Handle to the coordinate vector
    INTEGER, INTENT(INOUT) :: h_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata,p_DdataTmp
    INTEGER(PREC_QTREEIDX), DIMENSION(2) :: Isize
    
    ! Check if handle is associated
    IF (h_Ddata == ST_NOHANDLE) THEN
      Isize = (/2,rquadtree%NVT/)
      CALL storage_new('t_quadtree_copyfrom_handle','p_Ddata',&
          Isize,ST_DOUBLE,h_Ddata,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Ddata,Isize)
      IF (Isize(2) < rquadtree%NVT) THEN
        CALL storage_realloc('t_quadtree_copyfrom_handle',&
            rquadtree%NVT,h_Ddata,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF
    
    ! Set pointers
    CALL storage_getbase_double2D(h_Ddata,p_Ddata)
    CALL storage_getbase_double2D(rquadtree%h_Ddata,p_DdataTmp)

    ! Copy data
    CALL DCOPY(2*rquadtree%NVT,p_DdataTmp,1,p_Ddata,1)
  END SUBROUTINE t_quadtree_copyfrom_handle

  !************************************************************************

!<subroutine>

  SUBROUTINE t_quadtree_copyfrom_array(rquadtree,p_Ddata)

!<description>
    ! This subroutine copies the content of the quadtree to an array.
    ! If the array is too small, then this subroutines terminates.
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
!</input>

!<inputoutput>
    ! Coordinate vector
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p_Ddata
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DdataTmp
    INTEGER(PREC_QTREEIDX), DIMENSION(2) :: Isize

    ! Check size of array
    Isize = SHAPE(p_Ddata)
    IF (Isize(1) /= 2 .OR. Isize(2) < rquadtree%NVT) THEN
      PRINT *, "t_quadtree_copyfrom_array: Array too small!"
      CALL sys_halt()
    END IF

    ! Set pointers
    CALL storage_getbase_double2D(rquadtree%h_Ddata,p_DdataTmp)

    ! Copy data
    CALL DCOPY(2*rquadtree%NVT,p_DdataTmp,1,p_Ddata,1)
    
  END SUBROUTINE t_quadtree_copyfrom_array

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_quadtree_copyto_handle(h_Ddata,rquadtree)

!<description>
    ! This subroutine copies the content of a handle to the quadtree.
!</description>

!<input>
    ! Handle to the coordinate vector
    INTEGER, INTENT(IN) :: h_Ddata
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata
    INTEGER(PREC_QTREEIDX) :: ivt,jvt,inode,ipos
    
    ! Set pointer and check its shape
    CALL storage_getbase_double2D(h_Ddata,p_Ddata)
    IF (SIZE(p_Ddata,1) /= 2) THEN
      PRINT *, "t_quadtree_copyto_handle: First dimension of array must be 2!"
      CALL sys_halt()
    END IF
    
    DO ivt=1,SIZE(p_Ddata,2)
      IF (search(rquadtree,p_Ddata(:,ivt),inode,ipos,jvt) == QTREE_NOT_FOUND)&
          CALL insert(rquadtree,ivt,p_Ddata(:,ivt),inode)
    END DO
  END SUBROUTINE t_quadtree_copyto_handle

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_quadtree_copyto_array(p_Ddata,rquadtree)

!<description>
    ! This subroutine copies the content of an array to the quadtree.
!</description>

!<input>
    ! Coordinate vector
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: p_Ddata
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_QTREEIDX) :: ivt,jvt,inode,ipos
    
    IF (SIZE(p_Ddata,1) /= 2) THEN
      PRINT *, "t_quadtree_copyto_handle: First dimension of array must be 2!"
      CALL sys_halt()
    END IF

    DO ivt=1,SIZE(p_Ddata,2)
      IF (search(rquadtree,p_Ddata(:,ivt),inode,ipos,jvt) == QTREE_NOT_FOUND)&
          CALL insert(rquadtree,ivt,p_Ddata(:,ivt),inode)
    END DO
  END SUBROUTINE t_quadtree_copyto_array

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_quadtree_insert(rquadtree,ivt,Ddata,inode)

!<description>
    ! This subroutine inserts a new coordinate item to the quadtree
!</description>

!<input>
    ! Number of the inserted vertex
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: ivt

    ! Number of the quad to which vertex is inserted
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: inode
    
    ! Coordinates of the new vertex
    REAL(DP), DIMENSION(2), INTENT(IN) :: Ddata
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! Check if there is enough space left in the nodal component of the quadtree
    IF (rquadtree%NVT == rquadtree%NNVT)&
        CALL resizeNVT(rquadtree,CEILING(rquadtree%dfactor*rquadtree%NNVT))
    
    ! Update values and add new entry recursively
    rquadtree%NVT            = rquadtree%NVT+1
    rquadtree%p_Ddata(:,ivt) = Ddata
    CALL insert(ivt,inode)
    
  CONTAINS  
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    RECURSIVE SUBROUTINE insert(ivt,inode)
      INTEGER(PREC_QTREEIDX), INTENT(IN) :: ivt,inode
      REAL(DP) :: xmin,ymin,xmax,ymax,xmid,ymid
      INTEGER(PREC_QTREEIDX) :: i,jnode,nnode
      
      IF (rquadtree%p_Knode(QTREE_STATUS,inode) == QTREE_MAX) THEN
        
        IF (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode)&
            CALL resizeNNODE(rquadtree,CEILING(rquadtree%dfactor*rquadtree%NNNODE))
        
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
        rquadtree%p_Dbbox(:,nnode+QTREE_NW)       = (/xmin,ymid,xmid,ymax/)
        
        ! SW-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SW) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SW) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SW) = QTREE_SW
        rquadtree%p_Dbbox(:,nnode+QTREE_SW)       = (/xmin,ymin,xmid,ymid/)
        
        ! SE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_SE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_SE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_SE) = QTREE_SE
        rquadtree%p_Dbbox(:,nnode+QTREE_SE)       = (/xmid,ymin,xmax,ymid/)
        
        ! NE-quad
        rquadtree%p_Knode(QTREE_STATUS,nnode+QTREE_NE) = QTREE_EMPTY
        rquadtree%p_Knode(QTREE_PARENT,nnode+QTREE_NE) = inode
        rquadtree%p_Knode(QTREE_PARPOS,nnode+QTREE_NE) = QTREE_NE
        rquadtree%p_Dbbox(:,nnode+QTREE_NE)       = (/xmid,ymid,xmax,ymax/)
        
        ! Add the four values from INODE to the four new quads 
        ! NNODE+1:NNODE+4 recursively
        DO i=1,QTREE_MAX
          jnode=nnode+direction(rquadtree,rquadtree%p_Ddata(:,rquadtree%p_Knode(i,inode)),inode)
          CALL insert(rquadtree%p_Knode(i,inode),jnode)
        END DO
        
        ! Mark the current quad as subdivided and set pointers to its four children 
        rquadtree%p_Knode(QTREE_STATUS,inode) = QTREE_SUBDIV
        rquadtree%p_Knode(1:4,inode)     = (/nnode+QTREE_NW,nnode+QTREE_SW,nnode+QTREE_SE,nnode+QTREE_NE/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+direction(rquadtree,rquadtree%p_Ddata(:,ivt),inode)
        CALL insert(ivt,jnode)
        
      ELSE
        
        ! Quad is not full, so new items can be stored
        rquadtree%p_Knode(QTREE_STATUS,inode) = rquadtree%p_Knode(QTREE_STATUS,inode)+1
        rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS,inode),inode) = ivt
      END IF
    END SUBROUTINE insert
  END SUBROUTINE t_quadtree_insert
  
  ! ***************************************************************************
  
!<function>
  
  FUNCTION t_quadtree_delete(rquadtree,Ddata,ivt) RESULT(f)

!<description>
    ! This function deletes an item from the quadtree
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    REAL(DP), DIMENSION(2), INTENT(IN) :: Ddata
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>

!<output>
    ! Number of the vertex that is deleted
    INTEGER(PREC_QTREEIDX), INTENT(OUT) :: ivt
!</output>

!<result>
    ! Result of the deletion: QTREE_NOT_FOUND / QTREE_FOUND
    INTEGER :: f
!</result>
!</function>
    
    ! local variables
    INTEGER :: inode,ipos,jpos,jvt
    
    ! Search for the given coordinates
    f=search(rquadtree,Ddata,inode,ipos,ivt)
    
    ! What can we do from the searching
    IF (f == QTREE_FOUND) THEN
      
      ! Remove item IVT from quad INODE
      DO jpos=ipos+1,rquadtree%p_Knode(QTREE_STATUS,inode)
        rquadtree%p_Knode(jpos-1,inode) = rquadtree%p_Knode(jpos,inode)
      END DO
      rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS,inode),inode) = 0
      rquadtree%p_Knode(QTREE_STATUS,inode) = rquadtree%p_Knode(QTREE_STATUS,inode)-1
      
      ! If IVT is not last item move last item NVT to position IVT
      IF (ivt /= rquadtree%NVT) THEN
        IF (search(rquadtree,rquadtree%p_Ddata(1:2,rquadtree%NVT),inode,ipos,jvt) == QTREE_FOUND) THEN
          rquadtree%p_Ddata(:,ivt) = rquadtree%p_Ddata(:,rquadtree%NVT)
          rquadtree%p_Knode(ipos,inode) = ivt
        END IF
        ivt=rquadtree%NVT
      END IF
      rquadtree%NVT = rquadtree%NVT-1
    END IF
  END FUNCTION t_quadtree_delete
  
  ! ***************************************************************************

!<function>
  
  FUNCTION t_quadtree_search(rquadtree,Ddata,inode,ipos,ivt) RESULT(f)

!<description>
    ! This subroutine searches for given coordinates in the quadtree
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
    
    ! coordinates that should be searched
    REAL(DP), DIMENSION(2), INTENT(IN) :: Ddata
!</input>

!<output>
    ! Number of the quad in which the given coordinates are
    INTEGER(PREC_QTREEIDX), INTENT(OUT) :: inode

    ! Position of the coordinates in the quad
    INTEGER(PREC_QTREEIDX), INTENT(OUT) :: ipos

    ! Number of the vertex the coordinates correspond to
    INTEGER(PREC_QTREEIDX), INTENT(OUT) :: ivt
!</output>

!<result>
    ! Result of the searching: QTREE_NOT_FOUND / QTREE_FOUND
    INTEGER :: f
!</result>
!</function>
    
    ! Initialize
    inode=1; ipos=1; ivt=1; f=search(inode,ipos,ivt)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive searching routine follows

    RECURSIVE FUNCTION search(inode,ipos,ivt) RESULT(f)
      INTEGER(PREC_QTREEIDX), INTENT(INOUT) :: inode,ipos,ivt
      INTEGER :: f
      
      f=QTREE_NOT_FOUND
      IF (rquadtree%p_Knode(QTREE_STATUS,inode) == QTREE_SUBDIV) THEN
        
        ! Quad is subdivided. Compute child INODE which to look recursively.
        inode = rquadtree%p_Knode(direction(rquadtree,Ddata,inode),inode)
        f=search(inode,ipos,ivt)
        
      ELSE
        
        ! Quad is not subdivided. Search for (x,y) in current quad
        DO ipos=1,rquadtree%p_Knode(QTREE_STATUS,inode)
          ivt = rquadtree%p_Knode(ipos,inode)
          IF (MAXVAL(ABS(rquadtree%p_Ddata(:,ivt)-Ddata)) <= SYS_EPSREAL) THEN
            f=QTREE_FOUND; RETURN
          END IF
        END DO
        
      END IF
    END FUNCTION search
  END FUNCTION t_quadtree_search

  !************************************************************************
  
!<function>
  
  PURE FUNCTION t_quadtree_direction(rquadtree,Ddata,inode) RESULT(d)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. Ddata
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
    
    ! Coordinates
    REAL(DP), DIMENSION(2), INTENT(IN) :: Ddata

    ! Number of quad
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: inode
!</input>

!<result>
    ! Further search direction
    INTEGER :: d
!</result>
!</function>
    
    ! local variables
    REAL(DP) :: xmid,ymid

    ! Compute midpoint of current quad
    xmid=(rquadtree%p_Dbbox(QTREE_XMIN,inode)+rquadtree%p_Dbbox(QTREE_XMAX,inode))/2._DP
    ymid=(rquadtree%p_Dbbox(QTREE_YMIN,inode)+rquadtree%p_Dbbox(QTREE_YMAX,inode))/2._DP
    
    IF (Ddata(1) > xmid) THEN
      IF (Ddata(2) > ymid) THEN
        d=QTREE_NE; RETURN
      ELSE
        d=QTREE_SE; RETURN
      END IF
    ELSE
      IF (Ddata(2) > ymid) THEN
        d=QTREE_NW; RETURN
      ELSE
        d=QTREE_SW; RETURN
      END IF
    END IF
  END FUNCTION t_quadtree_direction

  !************************************************************************
  
!<subroutine>

  SUBROUTINE t_quadtree_print(rquadtree,cfilename)

!<description>
    ! This subroutine writes the content of the quadtree to a file
    ! which can be visualized by means of Matlab
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree

    ! filename of the output file
    CHARACTER(LEN=*), INTENT(IN) :: cfilename
!</input>
!</subroutine>
    
    ! local variables
    REAL(DP) :: xmin,xmax,ymin,ymax
    INTEGER :: iunit
    
    iunit=sys_getFreeUnit()
    OPEN(UNIT=iunit,FILE=TRIM(ADJUSTL(cfilename)))
    xmin = rquadtree%p_Dbbox(QTREE_XMIN,1)
    ymin = rquadtree%p_Dbbox(QTREE_YMIN,1)
    xmax = rquadtree%p_Dbbox(QTREE_XMAX,1)
    ymax = rquadtree%p_Dbbox(QTREE_YMAX,1)
    CALL print(xmin,ymin,xmax,ymax,1)
    CLOSE(UNIT=iunit)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive print routine follows
    
    RECURSIVE SUBROUTINE print(xmin,ymin,xmax,ymax,inode)
      REAL(DP), INTENT(IN) :: xmin,ymin,xmax,ymax
      INTEGER(PREC_QTREEIDX), INTENT(IN) :: inode
      REAL(DP) :: xmid,ymid
      INTEGER(PREC_QTREEIDX) :: i

      WRITE(UNIT=iunit,FMT=*) 'rect'
      WRITE(UNIT=iunit,FMT=10) xmin,ymin,xmax,ymax
      
      IF (rquadtree%p_Knode(QTREE_STATUS,inode) == QTREE_SUBDIV) THEN
        
        xmid=(xmin+xmax)/2._DP
        ymid=(ymin+ymax)/2._DP
        
        CALL print(xmin,ymid,xmid,ymax,rquadtree%p_Knode(QTREE_NW,inode))
        CALL print(xmin,ymin,xmid,ymid,rquadtree%p_Knode(QTREE_SW,inode))
        CALL print(xmid,ymin,xmax,ymid,rquadtree%p_Knode(QTREE_SE,inode))
        CALL print(xmid,ymid,xmax,ymax,rquadtree%p_Knode(QTREE_NE,inode))
        
      ELSEIF (rquadtree%p_Knode(QTREE_STATUS,inode) > QTREE_EMPTY) THEN
        
        DO i=1,rquadtree%p_Knode(QTREE_STATUS,inode)
          WRITE(UNIT=iunit,FMT=*) 'node'
          WRITE(UNIT=iunit,FMT=20) rquadtree%p_Ddata(:,rquadtree%p_Knode(i,inode))
        END DO
        
      END IF
      
10    FORMAT(4E15.6E3)
20    FORMAT(2E15.6E3)
    END SUBROUTINE print
  END SUBROUTINE t_quadtree_print

  !************************************************************************

!<subroutine>

  SUBROUTINE t_quadtree_info(rquadtree)

!<description>
    ! This subroutine outputs statistical info about the quadtree
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
!</input>
!</subroutine>

    WRITE(*,FMT=*) ' Quadtree:'
    WRITE(*,FMT=*) ' ========='
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Ddata =',rquadtree%h_Ddata
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Dbbox =',rquadtree%h_Dbbox
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Knode =',rquadtree%h_Knode
    WRITE(*,*)
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8)') '  NVT     =',rquadtree%NVT,'NNVT     =',rquadtree%NNVT
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8)') '  NNODE   =',rquadtree%NNODE,'NNNODE   =',rquadtree%NNNODE
    WRITE(*,FMT='(1X,A,1X,I8,A,4X,F5.1,A,4X,F5.1,A)') '  NRESIZE =',rquadtree%NRESIZE, "   QUADS    =", &
        100*rquadtree%NNODE/REAL(rquadtree%NNNODE,DP),'%  FILLING  =',100*rquadtree%NVT/REAL(rquadtree%NNVT,DP),'%'
    WRITE(*,*)
  END SUBROUTINE t_quadtree_info

  !************************************************************************

!<function>

  PURE FUNCTION t_quadtree_getsize(rquadtree) RESULT(nvt)

!<description>
    ! This function returns the number of vertices stored in the quadtree
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
!</input>

!<result>
    ! number of vertices in quadtree
    INTEGER(PREC_QTREEIDX) :: nvt
!</result>
!</function>

    nvt=rquadtree%NVT
  END FUNCTION t_quadtree_getsize

  !************************************************************************

!<function>

  FUNCTION t_quadtree_getboundingbox(rquadtree,inode) RESULT(bbox)
    
!<description>
    ! This function returns the bounding box of the specified quad.
    ! If no quad number is given, then the outer bounding box is
    ! returned.
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree

    ! OPTIONAL: number of quad for which bounding box should be
    ! returned
    INTEGER, INTENT(IN), OPTIONAL :: inode
!</input>

!<result>
    ! bounding box
    REAL(DP), DIMENSION(4) :: bbox
!</result>
!</function>
    
    IF (PRESENT(inode)) THEN
      IF (inode > rquadtree%NVT) THEN
        PRINT *, "t_quadtree_getboundingbox: quad number exceeds quadtree dimension"
        CALL sys_halt()
      END IF
      bbox=rquadtree%p_Dbbox(QTREE_XMIN:QTREE_YMAX,inode)
    ELSE
      bbox=rquadtree%p_Dbbox(QTREE_XMIN:QTREE_YMAX,1)
    END IF
  END FUNCTION t_quadtree_getboundingbox

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION t_quadtree_getXvalue(rquadtree,ivt) RESULT(x)

!<description>
    ! This function returns the X-value at the given position.
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree

    ! position in the quadtree
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: ivt
!</input>

!<result>
    REAL(DP) :: x
!</result>
!</function>

    x=rquadtree%p_Ddata(1,ivt)
  END FUNCTION t_quadtree_getXvalue

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION t_quadtree_getYvalue(rquadtree,ivt) RESULT(y)

!<description>
    ! This function returns the Y-value at the given position.
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree

    ! position in the quadtree
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: ivt
!</input>

!<result>
    REAL(DP) :: y
!</result>
!</function>

    y=rquadtree%p_Ddata(2,ivt)
  END FUNCTION t_quadtree_getYvalue
END MODULE quadtree

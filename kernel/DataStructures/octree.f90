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

MODULE octree
  USE fsystem
  USE storage
  USE genoutput
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
  PUBLIC :: otree_getDirection
  PUBLIC :: otree_printOctree
  PUBLIC :: otree_infoOctree
  PUBLIC :: otree_getsize
  PUBLIC :: otree_getBoundingBox
  PUBLIC :: otree_getX
  PUBLIC :: otree_getY
  PUBLIC :: otree_getZ
  PUBLIC :: otree_duplicateOctree
  PUBLIC :: otree_restoreOctree

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
  
  INTERFACE otree_copyToOctree
    MODULE PROCEDURE otree_copyToOctree_handle
    MODULE PROCEDURE otree_copyToOctree_array
  END INTERFACE

  INTERFACE otree_copyFromOctree
    MODULE PROCEDURE otree_copyFromOctree_handle
    MODULE PROCEDURE otree_copyFromOctree_array
  END INTERFACE

  INTERFACE otree_deleteFromOctree
    MODULE PROCEDURE otree_deleteFromOctree
    MODULE PROCEDURE otree_deleteFromOctreeByNumber
  END INTERFACE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE otree_createOctree(roctree, nnvt, nnnode,&
                                xmin, ymin, zmin, xmax, ymax, zmax, dfactor)
  
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
    CALL storage_new('otree_createOctree', 'p_Ddata', (/3,nnvt/),&
                     ST_DOUBLE, roctree%h_Ddata, ST_NEWBLOCK_ZERO)
    CALL storage_new('otree_createOctree', 'p_Dbbox', (/6,nnnode/),&
                     ST_DOUBLE, roctree%h_Dbbox, ST_NEWBLOCK_ZERO)
    CALL storage_new('otree_createOctree', 'p_Knode', (/11,nnnode/),&
                     ST_INT, roctree%h_Knode, ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double2D(roctree%h_Ddata, roctree%p_Ddata)
    CALL storage_getbase_double2D(roctree%h_Dbbox, roctree%p_Dbbox)
    CALL storage_getbase_int2D(roctree%h_Knode,    roctree%p_Knode)

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
  END SUBROUTINE otree_createOctree

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE otree_releaseOctree(roctree)

!<description>
    ! This subroutine releases an existing octree
!</description>

!<inputoutput>
    ! quadtree
    TYPE(t_octree), INTENT(INOUT) :: roctree
!</inputoutput>
!</subroutine>
    
    ! Free memory
    IF (roctree%h_Ddata .NE. ST_NOHANDLE) CALL storage_free(roctree%h_Ddata)
    IF (roctree%h_Dbbox .NE. ST_NOHANDLE) CALL storage_free(roctree%h_Dbbox)
    IF (roctree%h_Knode .NE. ST_NOHANDLE) CALL storage_free(roctree%h_Knode)
    NULLIFY(roctree%p_Knode, roctree%p_Dbbox, roctree%p_Ddata)

    ! Reset values
    roctree%NNNODE  = 0
    roctree%NNVT    = 0
    roctree%NNODE   = 0
    roctree%NVT     = 0
    roctree%NRESIZE = 0
  END SUBROUTINE otree_releaseOctree

  !************************************************************************

!<subroutine>
  
  SUBROUTINE otree_copyFromOctree_handle(roctree, h_Ddata)

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
    IF (h_Ddata .EQ. ST_NOHANDLE) THEN
      Isize = (/3,roctree%NVT/)
      CALL storage_new('otree_copyFromOctree_handle','p_Ddata',&
                       Isize, ST_DOUBLE, h_Ddata, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Ddata, Isize)
      IF (Isize(2) < roctree%NVT) THEN
        CALL storage_realloc('otree_copyFromOctree_handle',&
                              roctree%NVT, h_Ddata, ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Set pointers
    CALL storage_getbase_double2D(h_Ddata, p_Ddata)
    CALL storage_getbase_double2D(roctree%h_Ddata, p_DdataTmp)

    ! Copy data
    CALL DCOPY(3*roctree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  END SUBROUTINE otree_copyFromOctree_handle

  !************************************************************************

!<subroutine>

  SUBROUTINE otree_copyFromOctree_array(roctree, p_Ddata)

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
    IF (Isize(1) .NE. 3 .OR. Isize(2) < roctree%NVT) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyFromOctree_array')
      CALL sys_halt()
    END IF

    ! Set pointers
    CALL storage_getbase_double2D(roctree%h_Ddata, p_DdataTmp)

    ! Copy data
    CALL DCOPY(3*roctree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  END SUBROUTINE otree_copyFromOctree_array

  !************************************************************************

!<subroutine>
  
  SUBROUTINE otree_copyToOctree_handle(h_Ddata, roctree)

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
    CALL storage_getbase_double2D(h_Ddata, p_Ddata)
    IF (SIZE(p_Ddata, 1) .NE. 3) THEN
      CALL output_line('First dimension of array must be 3!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyToOctree_handle')
      CALL sys_halt()
    END IF
    
    DO ivt = 1, SIZE(p_Ddata, 2)
      IF (otree_searchInOctree(roctree, p_Ddata(:,ivt),&
                               inode, ipos,jvt) .EQ. OTREE_NOT_FOUND)&
          CALL otree_insertIntoOctree(roctree, ivt, p_Ddata(:,ivt), inode)
    END DO
  END SUBROUTINE otree_copyToOctree_handle

  !************************************************************************

!<subroutine>
  
  SUBROUTINE otree_copyToOctree_array(p_Ddata, roctree)

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
    
    IF (SIZE(p_Ddata,1) .NE. 3) THEN
      CALL output_line('First dimension of array must be 3!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'otree_copyToOctree_array')
      CALL sys_halt()
    END IF

    DO ivt = 1, SIZE(p_Ddata,2)
      IF (otree_searchInOctree(roctree, p_Ddata(:,ivt),&
                               inode, ipos, jvt) .EQ. OTREE_NOT_FOUND)&
          CALL otree_insertIntoOctree(roctree, ivt, p_Ddata(:,ivt), inode)
    END DO
  END SUBROUTINE otree_copyToOctree_array

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE otree_insertIntoOctree(roctree, ivt, Ddata, inode)

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
    
    ! local variables
    INTEGER(PREC_OTREEIDX) :: isize

    ! Check if there is enough space left in the vertex component of the octree
    IF (roctree%NVT .EQ. roctree%NNVT) THEN
      isize = MAX(roctree%NNVT+1, CEILING(roctree%dfactor*roctree%NNVT))
      CALL resizeNVT(roctree, isize)
    END IF
    
    ! Update values and add new entry recursively
    roctree%NVT            = roctree%NVT+1
    roctree%p_Ddata(:,ivt) = Ddata
    CALL insert(ivt, inode)
    
  CONTAINS  
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    RECURSIVE SUBROUTINE insert(ivt, inode)
      INTEGER(PREC_OTREEIDX), INTENT(IN) :: ivt,inode
      REAL(DP) :: xmin,ymin,zmin,xmax,ymax,zmax,xmid,ymid,zmid
      INTEGER(PREC_OTREEIDX) :: i,jnode,knode,nnode,isize
      
      IF (roctree%p_Knode(OTREE_STATUS, inode) .EQ. OTREE_MAX) THEN
        
        IF (roctree%nnode+OTREE_MAX > roctree%nnnode) THEN
          isize = MAX(roctree%nnode+OTREE_MAX, CEILING(roctree%dfactor*roctree%NNNODE))
          CALL resizeNNODE(roctree, isize)
        END IF
        
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
        DO i = 1, OTREE_MAX
          knode = roctree%p_Knode(i,inode)
          jnode = nnode+otree_getDirection(roctree, roctree%p_Ddata(:,knode),inode)
          CALL insert(knode, jnode)
        END DO
        
        ! Mark the current nodes as subdivided and set pointers to its eight children 
        roctree%p_Knode(OTREE_STATUS,inode) = OTREE_SUBDIV
        roctree%p_Knode(1:OTREE_MAX, inode) = (/nnode+OTREE_NWF, nnode+OTREE_SWF,&
                                                nnode+OTREE_SEF, nnode+OTREE_NEF,&
                                                nnode+OTREE_NWB, nnode+OTREE_SWB,&
                                                nnode+OTREE_SEB, nnode+OTREE_NEB/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+otree_getDirection(roctree, roctree%p_Ddata(:,ivt),inode)
        CALL insert(ivt, jnode)
        
      ELSE
        
        ! Node is not full, so new items can be stored
        roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)+1
        roctree%p_Knode(roctree%p_Knode(OTREE_STATUS, inode), inode) = ivt
      END IF
    END SUBROUTINE insert
  END SUBROUTINE otree_insertIntoOctree
  
  ! ***************************************************************************
  
!<function>
  
  FUNCTION otree_deleteFromOctree(roctree, Ddata, ivt) RESULT(f)

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
    REAL(DP), DIMENSION(3) :: DdataTmp
    
    ! Search for the given coordinates
    f = otree_searchInOctree(roctree, Ddata, inode, ipos, ivt)
    
    ! What can we do from the searching
    IF (f .EQ. OTREE_FOUND) THEN
      
      ! Remove item IVT from node INODE
      DO jpos = ipos+1, roctree%p_Knode(OTREE_STATUS, inode)
        roctree%p_Knode(jpos-1, inode) = roctree%p_Knode(jpos, inode)
      END DO
      roctree%p_Knode(roctree%p_Knode(OTREE_STATUS, inode), inode) = 0
      roctree%p_Knode(OTREE_STATUS, inode) = roctree%p_Knode(OTREE_STATUS, inode)-1
      
      ! If IVT is not last item move last item NVT to position IVT
      IF (ivt .NE. roctree%NVT) THEN
        DdataTmp(1:3) = roctree%p_Ddata(1:3, roctree%NVT)
        IF (otree_searchInOctree(roctree, DdataTmp(:),&
                                 inode, ipos, jvt) .EQ. OTREE_FOUND) THEN
          roctree%p_Ddata(:, ivt) = roctree%p_Ddata(:, roctree%NVT)
          roctree%p_Knode(ipos, inode) = ivt
        END IF
        ivt = roctree%NVT
      END IF
      roctree%NVT = roctree%NVT-1
    END IF
  END FUNCTION otree_deleteFromOctree

  ! ***************************************************************************

!<function>

  FUNCTION otree_deleteFromOctreeByNumber(roctree, ivt, ivtReplace) RESULT(f)

!<description>
    ! This function deletes vertex with number IVT from the octree
!</description>

!<input>
    ! Number of the vertex to be deleted
    INTEGER(PREC_OTREEIDX), INTENT(IN) :: ivt
!</input>

!<inputoutput>
    ! octree
    TYPE(t_octree), INTENT(INOUT) :: roctree
!</inputoutput>

!<output>
    ! Number of the vertex that replaces the deleted vertex
    INTEGER(PREC_OTREEIDX), INTENT(OUT) :: ivtReplace
!</output>

!<result>
    ! Result of the deletion: OTREE_NOT_FOUND / OTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    REAL(DP), DIMENSION(3) :: Ddata
    
    IF (ivt .LE. roctree%NVT) THEN
      ! Get coordinates and invoke deletion routine
      Ddata = roctree%p_Ddata(:,ivt)
      f     = otree_deleteFromOctree(roctree, Ddata, ivtReplace)
    ELSE
      CALL output_line('Invalid vertex number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'otree_deleteFromOctreeByNumber')
      CALL sys_halt()
    END IF
  END FUNCTION otree_deleteFromOctreeByNumber

  ! ***************************************************************************

!<function>
  
  FUNCTION otree_searchInOctree(roctree, Ddata, inode, ipos, ivt) RESULT(f)

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
    inode = 1
    ipos  = 1
    ivt   = 1
    f     = search(inode, ipos, ivt)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive searching routine follows

    RECURSIVE FUNCTION search(inode, ipos, ivt) RESULT(f)
      INTEGER(PREC_OTREEIDX), INTENT(INOUT) :: inode,ipos,ivt
      INTEGER :: f
      
      f = OTREE_NOT_FOUND
      IF (roctree%p_Knode(OTREE_STATUS, inode) .EQ. OTREE_SUBDIV) THEN
        
        ! Node is subdivided. Compute child INODE which to look recursively.
        inode = roctree%p_Knode(otree_getDirection(roctree, Ddata, inode), inode)
        f     = search(inode,ipos, ivt)
        
      ELSE
        
        ! Node is not subdivided. Search for (x,y,z) in current node
        DO ipos = 1, roctree%p_Knode(OTREE_STATUS, inode)
          ivt = roctree%p_Knode(ipos, inode)
          IF (SQRT(SUM((roctree%p_Ddata(:, ivt)-Ddata)**2)) .LE. SYS_EPSREAL) THEN
            f = OTREE_FOUND; RETURN
          END IF
        END DO
        
      END IF
    END FUNCTION search
  END FUNCTION otree_searchInOctree

  !************************************************************************
  
!<function>
  
  PURE FUNCTION otree_getDirection(roctree, Ddata, inode) RESULT(d)

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
    xmid=0.5*(roctree%p_Dbbox(OTREE_XMIN, inode)+&
              roctree%p_Dbbox(OTREE_XMAX, inode))
    ymid=0.5*(roctree%p_Dbbox(OTREE_YMIN, inode)+&
              roctree%p_Dbbox(OTREE_YMAX, inode))
    zmid=0.5*(roctree%p_Dbbox(OTREE_ZMIN, inode)+&
              roctree%p_Dbbox(OTREE_ZMAX, inode))

    ! Do we have to look in the front or back?
    IF (Ddata(3) < zmid) THEN
      
      IF (Ddata(1) > xmid) THEN
        IF (Ddata(2) > ymid) THEN
          d = OTREE_NEF; RETURN
        ELSE
          d = OTREE_SEF; RETURN
        END IF
      ELSE
        IF (Ddata(2) > ymid) THEN
          d = OTREE_NWF; RETURN
        ELSE
          d = OTREE_SWF; RETURN
        END IF
      END IF

    ELSE

      IF (Ddata(1) > xmid) THEN
        IF (Ddata(2) > ymid) THEN
          d = OTREE_NEB; RETURN
        ELSE
          d = OTREE_SEB; RETURN
        END IF
      ELSE
        IF (Ddata(2) > ymid) THEN
          d = OTREE_NWB; RETURN
        ELSE
          d = OTREE_SWB; RETURN
        END IF
      END IF

    END IF
  END FUNCTION otree_getDirection

  !************************************************************************
  
!<subroutine>

  SUBROUTINE otree_printOctree(roctree, cfilename)

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
    
    iunit = sys_getFreeUnit()
    OPEN(UNIT=iunit, FILE=TRIM(ADJUSTL(cfilename)))
    xmin = roctree%p_Dbbox(OTREE_XMIN, 1)
    ymin = roctree%p_Dbbox(OTREE_YMIN, 1)
    zmin = roctree%p_Dbbox(OTREE_ZMIN, 1)
    xmax = roctree%p_Dbbox(OTREE_XMAX, 1)
    ymax = roctree%p_Dbbox(OTREE_YMAX, 1)
    zmax = roctree%p_Dbbox(OTREE_ZMAX, 1)
    CALL print(xmin, ymin, zmin, xmax, ymax, zmax, 1)
    CLOSE(UNIT=iunit)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive print routine follows
    
    RECURSIVE SUBROUTINE print(xmin, ymin, zmin, xmax, ymax, zmax, inode)
      REAL(DP), INTENT(IN) :: xmin,ymin,zmin,xmax,ymax,zmax
      INTEGER(PREC_OTREEIDX), INTENT(IN) :: inode
      REAL(DP)               :: xmid,ymid,zmid
      INTEGER(PREC_OTREEIDX) :: i

      WRITE(UNIT=iunit,FMT=*) 'hex'
      WRITE(UNIT=iunit,FMT=10) xmin, ymin, zmin, xmax, ymax, zmax
      
      IF (roctree%p_Knode(OTREE_STATUS, inode) .EQ. OTREE_SUBDIV) THEN
        
        xmid = 0.5*(xmin+xmax)
        ymid = 0.5*(ymin+ymax)
        zmid = 0.5*(zmin+zmax)
        
        CALL print(xmin, ymid, zmin, xmid, ymax, zmid,&
                   roctree%p_Knode(OTREE_NWF, inode))
        CALL print(xmin, ymin, zmin, xmid, ymid, zmid,&
                   roctree%p_Knode(OTREE_SWF, inode))
        CALL print(xmid, ymin, zmin, xmax, ymid, zmid,&
                   roctree%p_Knode(OTREE_SEF, inode))
        CALL print(xmid, ymid, zmin, xmax, ymax, zmid,&
                   roctree%p_Knode(OTREE_NEF, inode))

        CALL print(xmin, ymid, zmid, xmid, ymax, zmax,&
                   roctree%p_Knode(OTREE_NWB, inode))
        CALL print(xmin, ymin, zmid, xmid, ymid, zmax,&
                   roctree%p_Knode(OTREE_SWB, inode))
        CALL print(xmid, ymin, zmid, xmax, ymid, zmax,&
                   roctree%p_Knode(OTREE_SEB, inode))
        CALL print(xmid, ymid, zmid, xmax, ymax, zmax,&
                   roctree%p_Knode(OTREE_NEB, inode))

      ELSEIF (roctree%p_Knode(OTREE_STATUS, inode) > OTREE_EMPTY) THEN
        
        DO i = 1, roctree%p_Knode(OTREE_STATUS, inode)
          WRITE(UNIT=iunit, FMT=*) 'node'
          WRITE(UNIT=iunit, FMT=20) roctree%p_Ddata(:, roctree%p_Knode(i, inode))
        END DO
        
      END IF
      
10    FORMAT(6E15.6E3)
20    FORMAT(3E15.6E3)
    END SUBROUTINE print
  END SUBROUTINE otree_printOctree

  !************************************************************************

!<subroutine>

  SUBROUTINE otree_infoOctree(roctree)

!<description>
    ! This subroutine outputs statistical info about the octree
!</description>

!<input>
    ! Octree
    TYPE(t_octree), INTENT(IN) :: roctree
!</input>
!</subroutine>

    CALL output_line('Octree:')
    CALL output_line('-------')
    CALL output_line('NVT:     '//TRIM(sys_siL(roctree%NVT,15)))
    CALL output_line('NNVT:    '//TRIM(sys_siL(roctree%NNVT,15)))
    CALL output_line('NNODE:   '//TRIM(sys_siL(roctree%NNODE,15)))
    CALL output_line('NNNODE:  '//TRIM(sys_siL(roctree%NNNODE,15)))
    CALL output_line('NRESIZE: '//TRIM(sys_siL(roctree%NRESIZE,5)))
    CALL output_line('dfactor: '//TRIM(sys_sdL(roctree%dfactor,2)))
    CALL output_line('h_Ddata: '//TRIM(sys_siL(roctree%h_Ddata,15)))
    CALL output_line('h_Dbbox: '//TRIM(sys_siL(roctree%h_Dbbox,15)))
    CALL output_line('h_Knode: '//TRIM(sys_siL(roctree%h_Knode,15)))
    CALL output_lbrk()
    WRITE(*,*)
    CALL output_line('Current data memory usage: '//&
        TRIM(sys_sdL(100*roctree%NVT/REAL(roctree%NNVT,DP),2))//'%')
    CALL output_line('Current node memory usage: '//&
        TRIM(sys_sdL(100*roctree%NNODE/REAL(roctree%NNNODE,DP),2))//'%')
  END SUBROUTINE otree_infoOctree

  !************************************************************************

!<function>

  PURE FUNCTION otree_getSize(roctree) RESULT(nvt)

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

    nvt = roctree%NVT
  END FUNCTION otree_getSize

  !************************************************************************

!<function>

  FUNCTION otree_getBoundingBox(roctree, inode) RESULT(bbox)
    
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
        CALL output_line('Node number exceeds octree dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'otree_getBoundingBox')
        CALL sys_halt()
      END IF
      bbox = roctree%p_Dbbox(OTREE_XMIN:OTREE_ZMAX, inode)
    ELSE
      bbox = roctree%p_Dbbox(OTREE_XMIN:OTREE_ZMAX, 1)
    END IF
  END FUNCTION otree_getBoundingBox

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION otree_getX(roctree, ivt) RESULT(x)

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

    x = roctree%p_Ddata(1,ivt)
  END FUNCTION otree_getX

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION otree_getY(roctree, ivt) RESULT(y)

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

    y = roctree%p_Ddata(2,ivt)
  END FUNCTION otree_getY

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION otree_getZ(roctree, ivt) RESULT(z)

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

    z = roctree%p_Ddata(3,ivt)
  END FUNCTION otree_getZ

  !************************************************************************

!<subroutine>

  SUBROUTINE otree_duplicateOctree(roctree, roctreeBackup)

!<description>
    ! This subroutine makes a copy of an octree in memory.
    ! It does not make sense to share some information between octrees,
    ! so each vectors is physically copied from the source octree
    ! to the destination octree.
!</description>

!<input>
    ! Source octree
    TYPE(t_octree), INTENT(IN) :: roctree
!</input>

!<inputoutput>
    ! Destination octree
    TYPE(t_octree), INTENT(INOUT) :: roctreeBackup
!</inputoutput>
!</subroutine>

    ! Release backup octree
    CALL otree_releaseOctree(roctreeBackup)

    ! Copy all data
    roctreeBackup = roctree

    ! Reset handles
    roctreeBackup%h_Ddata = ST_NOHANDLE
    roctreeBackup%h_Dbbox = ST_NOHANDLE
    roctreeBackup%h_Knode = ST_NOHANDLE

    ! Copy storage blocks
    IF (roctree%h_Ddata .NE. ST_NOHANDLE) THEN
      CALL storage_copy(roctree%h_Ddata, roctreeBackup%h_Ddata)
      CALL storage_getbase_double2D(roctreeBackup%h_Ddata,&
                                    roctreeBackup%p_Ddata)
    END IF

    IF (roctree%h_Dbbox .NE. ST_NOHANDLE) THEN
      CALL storage_copy(roctree%h_Dbbox, roctreeBackup%h_Dbbox)
      CALL storage_getbase_double2D(roctreeBackup%h_Dbbox,&
                                    roctreeBackup%p_Dbbox)
    END IF

    IF (roctree%h_Knode .NE. ST_NOHANDLE) THEN
      CALL storage_copy(roctree%h_Knode, roctreeBackup%h_Knode)
      CALL storage_getbase_int2D(roctreeBackup%h_Knode,&
                                 roctreeBackup%p_Knode)
    END IF
  END SUBROUTINE otree_duplicateOctree

  !************************************************************************

!<subroutine>

  SUBROUTINE otree_restoreOctree(roctreeBackup, roctree)

!<description>
    ! This subroutine restores an octree from a previous backup.
!</description>

!<input>
    ! Backup of an octree
    TYPE(t_octree), INTENT(IN) :: roctreeBackup
!</input>

!<inputoutput>
    ! Destination octree
    TYPE(t_octree), INTENT(INOUT) :: roctree
!</inputoutput>
!</subroutine>

    ! Release octree
    CALL otree_releaseOctree(roctree)

    ! Duplicate the backup
    CALL otree_duplicateOctree(roctreeBackup, roctree)
  END SUBROUTINE otree_restoreOctree

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE resizeNVT(roctree, nnvt)

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

    CALL storage_realloc('resizeNVT', nnvt, roctree%h_Ddata,&
                         ST_NEWBLOCK_ZERO, .TRUE.)
    CALL storage_getbase_double2D(roctree%h_Ddata, roctree%p_Ddata)
    
    roctree%NNVT    = nnvt
    roctree%NRESIZE = roctree%NRESIZE+1
  END SUBROUTINE resizeNVT
  
  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE resizeNNODE(roctree, nnnode)

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

    CALL storage_realloc('resizeNNODE', nnnode, roctree%h_Dbbox,&
                         ST_NEWBLOCK_ZERO, .TRUE.)
    CALL storage_realloc('resizeNNODE', nnnode, roctree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .TRUE.)
    CALL storage_getbase_double2D(roctree%h_Dbbox, roctree%p_Dbbox)
    CALL storage_getbase_int2D(roctree%h_Knode,    roctree%p_Knode)
    
    roctree%NNNODE  = nnnode
    roctree%NRESIZE =roctree%NRESIZE+1
  END SUBROUTINE resizeNNODE
END MODULE octree


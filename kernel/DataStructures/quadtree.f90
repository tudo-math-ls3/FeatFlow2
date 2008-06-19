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
!#      -> Write quadtree to file
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
!# 15.) qtree_duplicateQuadtree
!#      -> Create a duplicate / backup of a quadtree
!#
!# 16.) qtree_restoreQuadtree
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

MODULE quadtree
  USE fsystem
  USE storage
  USE genoutput
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
  PUBLIC :: qtree_duplicateQuadtree
  PUBLIC :: qtree_restoreQuadtree

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
  
  INTERFACE qtree_copyToQuadtree
    MODULE PROCEDURE qtree_copyToQuadtree_handle
    MODULE PROCEDURE qtree_copyToQuadtree_array
  END INTERFACE

  INTERFACE qtree_copyFromQuadtree
    MODULE PROCEDURE qtree_copyFromQuadtree_handle
    MODULE PROCEDURE qtree_copyFromQuadtree_array
  END INTERFACE
   
  INTERFACE qtree_deleteFromQuadtree
    MODULE PROCEDURE qtree_deleteFromQuadtree
    MODULE PROCEDURE qtree_deleteFromQtreeByNumber
  END INTERFACE
  
CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE qtree_createQuadtree(rquadtree, nnvt, nnnode,&
                                  xmin, ymin, xmax, ymax, dfactor)
  
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
    
    INTEGER(I32), DIMENSION(2) :: Isize
    
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
    Isize = (/2,nnvt/)
    CALL storage_new('qtree_createQuadtree', 'p_Ddata',&
                     Isize, ST_DOUBLE, rquadtree%h_Ddata, ST_NEWBLOCK_ZERO)
    Isize = (/4,nnnode/)
    CALL storage_new('qtree_createQuadtree', 'p_Dbbox',&
                     Isize, ST_DOUBLE, rquadtree%h_Dbbox, ST_NEWBLOCK_ZERO)
    Isize = (/7,nnnode/)
    CALL storage_new('qtree_createQuadtree', 'p_Knode',&
                     Isize, ST_INT, rquadtree%h_Knode, ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double2D(rquadtree%h_Ddata, rquadtree%p_Ddata)
    CALL storage_getbase_double2D(rquadtree%h_Dbbox, rquadtree%p_Dbbox)
    CALL storage_getbase_int2D(rquadtree%h_Knode,    rquadtree%p_Knode)

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
  END SUBROUTINE qtree_createQuadtree

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE qtree_releaseQuadtree(rquadtree)

!<description>
    ! This subroutine releases an existing quadtree
!</description>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! Free memory
    IF (rquadtree%h_Ddata .NE. ST_NOHANDLE) CALL storage_free(rquadtree%h_Ddata)
    IF (rquadtree%h_Dbbox .NE. ST_NOHANDLE) CALL storage_free(rquadtree%h_Dbbox)
    IF (rquadtree%h_Knode .NE. ST_NOHANDLE) CALL storage_free(rquadtree%h_Knode)
    NULLIFY(rquadtree%p_Knode, rquadtree%p_Dbbox, rquadtree%p_Ddata)

    ! Reset values
    rquadtree%NNNODE  = 0
    rquadtree%NNVT    = 0
    rquadtree%NNODE   = 0
    rquadtree%NVT     = 0
    rquadtree%NRESIZE = 0
  END SUBROUTINE qtree_releaseQuadtree

  !************************************************************************

!<subroutine>
  
  SUBROUTINE qtree_copyFromQuadtree_handle(rquadtree, h_Ddata)

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
    IF (h_Ddata .EQ. ST_NOHANDLE) THEN
      Isize = (/2, rquadtree%NVT/)
      CALL storage_new('qtree_copyFromQuadtree_handle', 'p_Ddata',&
                       Isize, ST_DOUBLE, h_Ddata, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Ddata, Isize)
      IF (Isize(2) < rquadtree%NVT) THEN
        CALL storage_realloc('qtree_copyFromQuadtree_handle',&
                             rquadtree%NVT, h_Ddata, ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Set pointers
    CALL storage_getbase_double2D(h_Ddata, p_Ddata)
    CALL storage_getbase_double2D(rquadtree%h_Ddata, p_DdataTmp)

    ! Copy data
    CALL DCOPY(2*rquadtree%NVT, p_DdataTmp, 1, p_Ddata, 1)
  END SUBROUTINE qtree_copyFromQuadtree_handle

  !************************************************************************

!<subroutine>

  SUBROUTINE qtree_copyFromQuadtree_array(rquadtree, p_Ddata)

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
    IF (Isize(1) .NE. 2 .OR. Isize(2) < rquadtree%NVT) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyFromQuadtree_array')
      CALL sys_halt()
    END IF

    ! Set pointers
    CALL storage_getbase_double2D(rquadtree%h_Ddata, p_DdataTmp)

    ! Copy data
    CALL DCOPY(2*rquadtree%NVT, p_DdataTmp, 1, p_Ddata, 1)
    
  END SUBROUTINE qtree_copyFromQuadtree_array

  !************************************************************************

!<subroutine>
  
  SUBROUTINE qtree_copyToQuadtree_handle(h_Ddata, rquadtree)

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
    CALL storage_getbase_double2D(h_Ddata, p_Ddata)
    IF (SIZE(p_Ddata,1) .NE. 2) THEN
      CALL output_line('First dimension of array must be 2!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyToQuadtree_handle')
      CALL sys_halt()
    END IF
    
    DO ivt = 1, SIZE(p_Ddata, 2)
      IF (qtree_searchInQuadtree(rquadtree, p_Ddata(:,ivt),&
                                 inode, ipos, jvt) .EQ. QTREE_NOT_FOUND)&
          CALL qtree_insertIntoQuadtree(rquadtree, ivt, p_Ddata(:,ivt), inode)
    END DO
  END SUBROUTINE qtree_copyToQuadtree_handle

  !************************************************************************

!<subroutine>
  
  SUBROUTINE qtree_copyToQuadtree_array(p_Ddata, rquadtree)

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
    
    IF (SIZE(p_Ddata,1) .NE. 2) THEN
      CALL output_line('First dimension of array must be 2!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_copyToQuadtree_handle')
      CALL sys_halt()
    END IF

    DO ivt = 1, SIZE(p_Ddata, 2)
      IF (qtree_searchInQuadtree(rquadtree, p_Ddata(:,ivt),&
                                 inode, ipos, jvt) .EQ. QTREE_NOT_FOUND)&
          CALL qtree_insertIntoQuadtree(rquadtree, ivt, p_Ddata(:,ivt), inode)
    END DO
  END SUBROUTINE qtree_copyToQuadtree_array

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE qtree_insertIntoQuadtree(rquadtree, ivt, Ddata, inode)

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

    ! local variables
    INTEGER(PREC_QTREEIDX) :: isize
    
    ! Check if there is enough space left in the nodal component of the quadtree
    IF (rquadtree%NVT .EQ. rquadtree%NNVT) THEN
      isize = MAX(rquadtree%NNVT+1, CEILING(rquadtree%dfactor*rquadtree%NNVT))
      CALL resizeNVT(rquadtree, isize)
    END IF
    
    ! Update values and add new entry recursively
    rquadtree%NVT            = rquadtree%NVT+1
    rquadtree%p_Ddata(:,ivt) = Ddata
    CALL insert(ivt, inode)
    
  CONTAINS  
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    RECURSIVE SUBROUTINE insert(ivt,inode)
      INTEGER(PREC_QTREEIDX), INTENT(IN) :: ivt,inode
      REAL(DP)                           :: xmin,ymin,xmax,ymax,xmid,ymid
      INTEGER(PREC_QTREEIDX)             :: i,jnode,nnode,knode,isize
            
      IF (rquadtree%p_Knode(QTREE_STATUS,inode) .EQ. QTREE_MAX) THEN
        
        IF (rquadtree%nnode+QTREE_MAX > rquadtree%nnnode) THEN
          isize = MAX(rquadtree%nnode+QTREE_MAX, CEILING(rquadtree%dfactor*rquadtree%NNNODE))
          CALL resizeNNODE(rquadtree, isize)
        END IF
        
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
        DO i=1,QTREE_MAX
          knode = rquadtree%p_Knode(i,inode)
          jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,knode), inode)
          CALL insert(rquadtree%p_Knode(i, inode), jnode)
        END DO
        
        ! Mark the current quad as subdivided and set pointers to its four children 
        rquadtree%p_Knode(QTREE_STATUS,inode) = QTREE_SUBDIV
        rquadtree%p_Knode(1:4,inode)          = (/nnode+QTREE_NW, nnode+QTREE_SW,&
                                                  nnode+QTREE_SE, nnode+QTREE_NE/)
        
        ! Add the new entry to the next position recursively
        jnode = nnode+qtree_getDirection(rquadtree, rquadtree%p_Ddata(:,ivt), inode)
        CALL insert(ivt, jnode)
        
      ELSE
        
        ! Quad is not full, so new items can be stored
        rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)+1
        rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = ivt
      END IF
    END SUBROUTINE insert
  END SUBROUTINE qtree_insertIntoQuadtree
  
  ! ***************************************************************************
  
!<function>
  
  FUNCTION qtree_deleteFromQuadtree(rquadtree, Ddata, ivt) RESULT(f)

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
    REAL(DP), DIMENSION(2) :: DdataTmp
    
    ! Search for the given coordinates
    f = qtree_searchInQuadtree(rquadtree, Ddata, inode, ipos, ivt)
    
    ! What can we do from the searching
    IF (f .EQ. QTREE_FOUND) THEN
      
      ! Remove item IVT from node INODE
      DO jpos = ipos+1, rquadtree%p_Knode(QTREE_STATUS, inode)
        rquadtree%p_Knode(jpos-1, inode) = rquadtree%p_Knode(jpos, inode)
      END DO
      rquadtree%p_Knode(rquadtree%p_Knode(QTREE_STATUS, inode), inode) = 0
      rquadtree%p_Knode(QTREE_STATUS, inode) = rquadtree%p_Knode(QTREE_STATUS, inode)-1
      
      ! If IVT is not last item move last item to position IVT
      IF (ivt .NE. rquadtree%NVT) THEN
        DdataTmp(1:2) = rquadtree%p_Ddata(1:2,rquadtree%NVT)
        IF (qtree_searchInQuadtree(rquadtree, DdataTmp(:),&
                                   inode, ipos, jvt) .EQ. QTREE_FOUND) THEN
          rquadtree%p_Ddata(:, ivt)      = rquadtree%p_Ddata(:, rquadtree%NVT)
          rquadtree%p_Knode(ipos, inode) = ivt
        END IF

        ! Set number of removed vertex
        ivt = rquadtree%NVT
      END IF
      
      ! Decrease number of vertices
      rquadtree%NVT = rquadtree%NVT-1
    END IF
  END FUNCTION qtree_deleteFromQuadtree

  ! ***************************************************************************

!<function>

  FUNCTION qtree_deleteFromQtreeByNumber(rquadtree, ivt, ivtReplace) RESULT(f)

!<description>
    ! This function deletes vertex with number IVT from the quadtree
!</description>

!<input>
    ! Number of the vertex to be deleted
    INTEGER(PREC_QTREEIDX), INTENT(IN) :: ivt
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>

!<output>
    ! Number of the vertex that replaces the deleted vertex
    INTEGER(PREC_QTREEIDX), INTENT(OUT) :: ivtReplace
!</output>

!<result>
    ! Result of the deletion: QTREE_NOT_FOUND / QTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    REAL(DP), DIMENSION(2) :: Ddata
    
    IF (ivt .LE. rquadtree%NVT) THEN
      ! Get coordinates and invoke deletion routine
      Ddata = rquadtree%p_Ddata(:,ivt)
      f     = qtree_deleteFromQuadtree(rquadtree, Ddata, ivtReplace)
    ELSE
      CALL output_line('Invalid vertex number!',&
          OU_CLASS_ERROR,OU_MODE_STD,'qtree_deleteFromQtreeByNumber')
      CALL sys_halt()
    END IF
  END FUNCTION qtree_deleteFromQtreeByNumber
 
  ! ***************************************************************************

!<function>
  
  FUNCTION qtree_searchInQuadtree(rquadtree, Ddata, inode, ipos, ivt) RESULT(f)

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
    inode = 1
    ipos  = 1 
    ivt   = 1
    f     = search(inode, ipos, ivt)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive searching routine follows

    RECURSIVE FUNCTION search(inode, ipos, ivt) RESULT(f)
      INTEGER(PREC_QTREEIDX), INTENT(INOUT) :: inode,ipos,ivt
      INTEGER :: f
      
      f = QTREE_NOT_FOUND
      IF (rquadtree%p_Knode(QTREE_STATUS, inode) .EQ. QTREE_SUBDIV) THEN
        
        ! Quad is subdivided. Compute child INODE which to look recursively.
        inode = rquadtree%p_Knode(qtree_getDirection(rquadtree, Ddata, inode), inode)
        f     = search(inode, ipos, ivt)
        
      ELSE
        
        ! Quad is not subdivided. Search for (x,y) in current quad
        DO ipos = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          ivt = rquadtree%p_Knode(ipos, inode)
          IF (MAXVAL(ABS(rquadtree%p_Ddata(:, ivt)-Ddata)) .LE. SYS_EPSREAL) THEN
            f = QTREE_FOUND; RETURN
          END IF
        END DO
        
      END IF
    END FUNCTION search
  END FUNCTION qtree_searchInQuadtree

  !************************************************************************
  
!<function>
  
  PURE FUNCTION qtree_getDirection(rquadtree, Ddata, inode) RESULT(d)

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
    xmid = (rquadtree%p_Dbbox(QTREE_XMIN, inode)+&
            rquadtree%p_Dbbox(QTREE_XMAX, inode))/2._DP
    ymid = (rquadtree%p_Dbbox(QTREE_YMIN, inode)+&
            rquadtree%p_Dbbox(QTREE_YMAX, inode))/2._DP
    
    IF (Ddata(1) > xmid) THEN
      IF (Ddata(2) > ymid) THEN
        d = QTREE_NE; RETURN
      ELSE
        d = QTREE_SE; RETURN
      END IF
    ELSE
      IF (Ddata(2) > ymid) THEN
        d = QTREE_NW; RETURN
      ELSE
        d = QTREE_SW; RETURN
      END IF
    END IF
  END FUNCTION qtree_getDirection

  !************************************************************************
  
!<subroutine>

  SUBROUTINE qtree_printQuadtree(rquadtree, cfilename)

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
    OPEN(UNIT=iunit, FILE=TRIM(ADJUSTL(cfilename)))
    xmin = rquadtree%p_Dbbox(QTREE_XMIN, 1)
    ymin = rquadtree%p_Dbbox(QTREE_YMIN, 1)
    xmax = rquadtree%p_Dbbox(QTREE_XMAX, 1)
    ymax = rquadtree%p_Dbbox(QTREE_YMAX, 1)
    CALL print(xmin, ymin, xmax, ymax, 1)
    CLOSE(UNIT=iunit)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive print routine follows
    
    RECURSIVE SUBROUTINE print(xmin, ymin, xmax, ymax, inode)
      REAL(DP), INTENT(IN) :: xmin,ymin,xmax,ymax
      INTEGER(PREC_QTREEIDX), INTENT(IN) :: inode
      REAL(DP) :: xmid,ymid
      INTEGER(PREC_QTREEIDX) :: i

      WRITE(UNIT=iunit,FMT=*) 'rect'
      WRITE(UNIT=iunit,FMT=10) xmin,ymin,xmax,ymax
      
      IF (rquadtree%p_Knode(QTREE_STATUS,inode) .EQ. QTREE_SUBDIV) THEN
        
        xmid=(xmin+xmax)/2._DP
        ymid=(ymin+ymax)/2._DP
        
        CALL print(xmin, ymid, xmid, ymax, rquadtree%p_Knode(QTREE_NW, inode))
        CALL print(xmin, ymin, xmid, ymid, rquadtree%p_Knode(QTREE_SW, inode))
        CALL print(xmid, ymin, xmax, ymid, rquadtree%p_Knode(QTREE_SE, inode))
        CALL print(xmid, ymid, xmax, ymax, rquadtree%p_Knode(QTREE_NE, inode))
        
      ELSEIF (rquadtree%p_Knode(QTREE_STATUS, inode) > QTREE_EMPTY) THEN
        
        DO i = 1, rquadtree%p_Knode(QTREE_STATUS, inode)
          WRITE(UNIT=iunit, FMT=*) 'node'
          WRITE(UNIT=iunit, FMT=20) rquadtree%p_Ddata(:, rquadtree%p_Knode(i, inode))
        END DO
        
      END IF
      
10    FORMAT(4E15.6E3)
20    FORMAT(2E15.6E3)
    END SUBROUTINE print
  END SUBROUTINE qtree_printQuadtree

  !************************************************************************

!<subroutine>

  SUBROUTINE qtree_infoQuadtree(rquadtree)

!<description>
    ! This subroutine outputs statistical info about the quadtree
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
!</input>
!</subroutine>

    CALL output_line('Quadtree:')
    CALL output_line('---------')
    CALL output_line('NVT:     '//TRIM(sys_siL(rquadtree%NVT,15)))
    CALL output_line('NNVT:    '//TRIM(sys_siL(rquadtree%NNVT,15)))
    CALL output_line('NNODE:   '//TRIM(sys_siL(rquadtree%NNODE,15)))
    CALL output_line('NNNODE:  '//TRIM(sys_siL(rquadtree%NNNODE,15)))
    CALL output_line('NRESIZE: '//TRIM(sys_siL(rquadtree%NRESIZE,5)))
    CALL output_line('dfactor: '//TRIM(sys_sdL(rquadtree%dfactor,2)))
    CALL output_line('h_Ddata: '//TRIM(sys_siL(rquadtree%h_Ddata,15)))
    CALL output_line('h_Dbbox: '//TRIM(sys_siL(rquadtree%h_Dbbox,15)))
    CALL output_line('h_Knode: '//TRIM(sys_siL(rquadtree%h_Knode,15)))
    CALL output_lbrk()
    WRITE(*,*)
    CALL output_line('Current data memory usage: '//&
        TRIM(sys_sdL(100*rquadtree%NVT/REAL(rquadtree%NNVT,DP),2))//'%')
    CALL output_line('Current ndoe memory usage: '//&
        TRIM(sys_sdL(100*rquadtree%NNODE/REAL(rquadtree%NNNODE,DP),2))//'%')
  END SUBROUTINE qtree_infoQuadtree

  !************************************************************************

!<function>

  PURE FUNCTION qtree_getSize(rquadtree) RESULT(nvt)

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

    nvt = rquadtree%NVT
  END FUNCTION qtree_getSize

  !************************************************************************

!<function>

  FUNCTION qtree_getBoundingBox(rquadtree, inode) RESULT(bbox)
    
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
        CALL output_line('Node number exceeds quadtree dimension',&
            OU_CLASS_ERROR,OU_MODE_STD,'qtree_getBoundingBox')
        CALL sys_halt()
      END IF
      bbox = rquadtree%p_Dbbox(QTREE_XMIN:QTREE_YMAX, inode)
    ELSE
      bbox = rquadtree%p_Dbbox(QTREE_XMIN:QTREE_YMAX, 1)
    END IF
  END FUNCTION qtree_getBoundingBox

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION qtree_getX(rquadtree, ivt) RESULT(x)

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

    x = rquadtree%p_Ddata(1, ivt)
  END FUNCTION qtree_getX

  !************************************************************************

!<function>

  ELEMENTAL FUNCTION qtree_getY(rquadtree, ivt) RESULT(y)

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

    y = rquadtree%p_Ddata(2, ivt)
  END FUNCTION qtree_getY

  !************************************************************************
  
!<subroutine>

  SUBROUTINE qtree_duplicateQuadtree(rquadtree, rquadtreeBackup)

!<description>
    ! This subroutine makes a copy of a quadtree in memory.
    ! It does not make sense to share some information between quadtrees,
    ! so each vectors is physically copied from the source quadtree
    ! to the destination quadtree.
!</description>

!<input>
    ! Source quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
!</input>

!<inputoutput>
    ! Destination quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtreeBackup
!</inputoutput>
!</subroutine>

    ! Release backup quadtree
    CALL qtree_releaseQuadtree(rquadtreeBackup)

    ! Copy all data
    rquadtreeBackup = rquadtree

    ! Reset handles
    rquadtreeBackup%h_Ddata = ST_NOHANDLE
    rquadtreeBackup%h_Dbbox = ST_NOHANDLE
    rquadtreeBackup%h_Knode = ST_NOHANDLE

    ! Copy storage blocks
    IF (rquadtree%h_Ddata .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rquadtree%h_Ddata, rquadtreeBackup%h_Ddata)
      CALL storage_getbase_double2D(rquadtreeBackup%h_Ddata,&
                                    rquadtreeBackup%p_Ddata)
    END IF

    IF (rquadtree%h_Dbbox .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rquadtree%h_Dbbox, rquadtreeBackup%h_Dbbox)
      CALL storage_getbase_double2D(rquadtreeBackup%h_Dbbox,&
                                    rquadtreeBackup%p_Dbbox)
    END IF

    IF (rquadtree%h_Knode .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rquadtree%h_Knode, rquadtreeBackup%h_Knode)
      CALL storage_getbase_int2D(rquadtreeBackup%h_Knode,&
                                 rquadtreeBackup%p_Knode)
    END IF
  END SUBROUTINE qtree_duplicateQuadtree

  !************************************************************************

!<subroutine>

  SUBROUTINE qtree_restoreQuadtree(rquadtreeBackup, rquadtree)

!<description>
    ! This subroutine restores a quadtree from a previous backup.
!</description>

!<input>
    ! Backup of an quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtreeBackup
!</input>

!<inputoutput>
    ! Destination quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>

    ! Release quadtree
    CALL qtree_releaseQuadtree(rquadtree)

    ! Duplicate the backup
    CALL qtree_duplicateQuadtree(rquadtreeBackup, rquadtree)
  END SUBROUTINE qtree_restoreQuadtree

  !************************************************************************

!<subroutine>
  
  SUBROUTINE resizeNVT(rquadtree, nnvt)

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

    CALL storage_realloc('resizeNVT', nnvt, rquadtree%h_Ddata,&
                         ST_NEWBLOCK_ZERO, .TRUE.)
    CALL storage_getbase_double2D(rquadtree%h_Ddata, rquadtree%p_Ddata)
    
    rquadtree%NNVT    = nnvt
    rquadtree%NRESIZE = rquadtree%NRESIZE+1
  END SUBROUTINE resizeNVT

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE resizeNNODE(rquadtree, nnnode)

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

    CALL storage_realloc('resizeNNODE', nnnode, rquadtree%h_Dbbox,&
                         ST_NEWBLOCK_ZERO, .TRUE.)
    CALL storage_realloc('resizeNNODE', nnnode, rquadtree%h_Knode,&
                         ST_NEWBLOCK_ZERO, .TRUE.)
    CALL storage_getbase_double2D(rquadtree%h_Dbbox, rquadtree%p_Dbbox)
    CALL storage_getbase_int2D(rquadtree%h_Knode,    rquadtree%p_Knode)
    
    rquadtree%NNNODE  = nnnode
    rquadtree%NRESIZE = rquadtree%NRESIZE+1
  END SUBROUTINE resizeNNODE
END MODULE quadtree

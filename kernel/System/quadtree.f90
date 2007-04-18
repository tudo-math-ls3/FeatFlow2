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
!# 1.) createQuadtree
!#     -> Create a new quadtree structure
!#
!# 2.) releaseQuadtree
!#     -> Release an existing quadtree
!#
!# 3.) resizeQuadtree
!#     -> Resize an existing quadtree
!#
!# 4.) copyQuadtree
!#     -> Copy data to/from a quadtree
!#
!# 5.) insertIntoQuadtree
!#     -> Insert data into quadtree
!#
!# 6.) deleteFromQuadtree
!#     -> Delete data from quadtree
!#
!# 7.) searchInQuadtree
!#     -> Search data in quadtree
!#
!# 8.) printQuadtree
!#     -> Write quadtree to file
!#
!# 9.) infoQuadtree
!#     -> Output info about quadtree
!#
!# 10.) getBoundingBox
!#      -> Return the outer bounding box
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
  PUBLIC :: qtree_copyQuadtree
  PUBLIC :: qtree_insertIntoQuadtree
  PUBLIC :: qtree_deleteFromQuadtree
  PUBLIC :: qtree_searchInQuadtree
  PUBLIC :: qtree_printQuadtree
  PUBLIC :: qtree_infoQuadtree
  PUBLIC :: qtree_getBoundingBox

!<constants>

!<constantblock description="KIND values for quadtree data">
  
  ! kind value for indices in quadtree
  INTEGER, PARAMETER, PUBLIC :: PREC_QTIDX = I32

!</constantblock>

!<constantblock description="Constants for quadtree structure">
  
  ! Position of the status information
  INTEGER, PARAMETER :: QSTATUS = 7
  
  ! Position of the parent information
  INTEGER, PARAMETER :: QPARENT = 6

  ! Position of the "position" of the parent information
  INTEGER, PARAMETER :: QPARPOS = 5

  ! Number of free positions in quad
  INTEGER, PARAMETER :: QFREE   = 5

  ! Item in "North-East" position
  INTEGER, PARAMETER :: QNE     = 4

  ! Item in "South-East" position
  INTEGER, PARAMETER :: QSE     = 3

  ! Item in "South-West" position
  INTEGER, PARAMETER :: QSW     = 2

  ! Item in "North-West" position
  INTEGER, PARAMETER :: QNW     = 1

  ! Identifier: Quad is empty
  INTEGER, PARAMETER :: QEMPTY  =  0

  ! Identifier: Status is subdivided
  INTEGER, PARAMETER :: QSUBDIV = -1
  
  ! Maximum number of items for each quad
  INTEGER, PARAMETER :: QMAX    =  4

!</constantblock> 

!<constantblock description="Constants for quadtree bounding-box">

  ! Position of the x-minimal value
  INTEGER, PARAMETER :: QXMIN   =  1

  ! Position of the y-minimal value
  INTEGER, PARAMETER :: QYMIN   =  2

  ! Position of the x-maximal value
  INTEGER, PARAMETER :: QXMAX   =  3

  ! Position of the y-maximal value
  INTEGER, PARAMETER :: QYMAX   =  4

!</constantblock>   
  
!<constantblock description="Constants for quadtree operations">

  ! Item could not be found in the quadtree
  INTEGER, PARAMETER, PUBLIC :: QNOT_FOUND = -1

  ! Item could be found in the quadtree
  INTEGER, PARAMETER, PUBLIC :: QFOUND     =  0
!</constantblock>
!</constants>

!<types>
!<typeblock>

  ! A linear quadtree implemented as array

  TYPE t_quadtree
    PRIVATE

    ! Handle to data vector
    INTEGER :: h_Ddata             = ST_NOHANDLE

    ! Handle to bounding box
    INTEGER :: h_Dbbox             = ST_NOHANDLE

    ! Handle to quadtree structure
    !   KQUAD(1:7,1:NNQUAD)
    !   KQUAD(  7,IQUAD)    : < 0, the quad has been subdivided
    !                         = 0, the quad is empty
    !                         > 0, the number of points stored in the quad
    !   KQUAD(  6,IQUAD)    : > 0, the quad the present quad came from
    !   KQUAD(  5,IQUAD)    : > 0, the position in the quad the present 
    !                           quad came from
    !   KQUAD(1:4,IQUAD)    : for KQUAD(7,IQUAD) > 0 : the points stored in the quad
    !                         for KQUAD(7,IQUAD) < 0 : the quads into
    !                                                  which the present quad was subdivided
    INTEGER :: h_Kquad             = ST_NOHANDLE

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
    INTEGER(PREC_QTIDX), DIMENSION(:,:), POINTER :: p_Kquad

    ! Number of vertices currently stored in the quadtree
    INTEGER(PREC_QTIDX) :: NVT    = 0

    ! Total number of vertices that can be stored  in the quadtree
    INTEGER(PREC_QTIDX) :: NNVT   = 0

    ! Number of subdivisions currently store in the quadtree
    INTEGER(PREC_QTIDX) :: NQUAD  = 0

    ! Total number of subdivision that can be stored in the quadtree
    INTEGER(PREC_QTIDX) :: NNQUAD = 0

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
  
  INTERFACE qtree_copyQuadtree
    MODULE PROCEDURE t_quadtree_copyto
    MODULE PROCEDURE t_quadtree_copyfrom
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

  INTERFACE qtree_getBoundingBox
    MODULE PROCEDURE t_quadtree_getboundingbox
  END INTERFACE

  INTERFACE resizeNVT
    MODULE PROCEDURE t_quadtree_resize_nvt
  END INTERFACE

  INTERFACE resizeNQUAD
    MODULE PROCEDURE t_quadtree_resize_nquad
  END INTERFACE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_quadtree_create(rquadtree,nnvt,nnquad,xmin,ymin,xmax,ymax,dfactor)
  
!<description>
    ! This subroutine creates a new quadtree.
!</description>

!<input>
    ! Total number of vertices that should be stored in the quadtree
    INTEGER(PREC_QTIDX), INTENT(IN) :: nnvt

    ! Total number of subdivisions that should be stored in the quadtree
    INTEGER(PREC_QTIDX), INTENT(IN) :: nnquad

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
    rquadtree%NNQUAD = nnquad
    rquadtree%NNVT   = nnvt
    rquadtree%NQUAD  = 0
    rquadtree%NVT    = 0
    rquadtree%NRESIZE= 0
    
    ! Allocate memory and associate pointers
    CALL storage_new('t_quadtree_init','p_Ddata',(/2,nnvt/),ST_DOUBLE,rquadtree%h_Ddata,ST_NEWBLOCK_ZERO)
    CALL storage_new('t_quadtree_init','p_Dbbox',(/4,nnquad/),ST_DOUBLE,rquadtree%h_Dbbox,ST_NEWBLOCK_ZERO)
    CALL storage_new('t_quadtree_init','p_Kquad',(/7,nnquad/),ST_INT,rquadtree%h_Kquad,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double2D(rquadtree%h_Ddata,rquadtree%p_Ddata)
    CALL storage_getbase_double2D(rquadtree%h_Dbbox, rquadtree%p_Dbbox)
    CALL storage_getbase_int2D(rquadtree%h_Kquad,rquadtree%p_Kquad)

    ! Initialize first quad
    rquadtree%nquad=1
    rquadtree%p_Kquad(QSTATUS,1)=QEMPTY
    rquadtree%p_Kquad(QPARENT,1)=QEMPTY
    rquadtree%p_Kquad(QFREE,1)=1
    
    ! Set co-ordinates of the bounding box
    rquadtree%p_Dbbox(QXMIN,1) = xmin
    rquadtree%p_Dbbox(QYMIN,1) = ymin
    rquadtree%p_Dbbox(QXMAX,1) = xmax
    rquadtree%p_Dbbox(QYMAX,1) = ymax
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
    IF (rquadtree%h_Dbbox  /= ST_NOHANDLE) CALL storage_free(rquadtree%h_Dbbox)
    IF (rquadtree%h_Kquad  /= ST_NOHANDLE) CALL storage_free(rquadtree%h_Kquad)
    NULLIFY(rquadtree%p_Kquad,rquadtree%p_Dbbox,rquadtree%p_Ddata)

    ! Reset values
    rquadtree%NNQUAD = 0
    rquadtree%NNVT   = 0
    rquadtree%NQUAD  = 0
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
    INTEGER(PREC_QTIDX), INTENT(IN) :: nnvt
!</input>

!<inputoutput>
    ! quadtree that should be resized
    TYPE(t_quadtree) :: rquadtree
!</inputoutput>
!</subroutine>

    CALL storage_realloc('t_quadtree_resize',nnvt,rquadtree%h_Ddata,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_double2D(rquadtree%h_Ddata,rquadtree%p_Ddata)
    
    rquadtree%NNVT   = nnvt
    rquadtree%NRESIZE=rquadtree%NRESIZE+1
  END SUBROUTINE t_quadtree_resize_nvt

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_quadtree_resize_nquad(rquadtree,nnquad)

!<description>
    ! This subroutine reallocates memory for an existing quadtree
!</description>

!<input>
    ! New number of quads that should be stored in the quadtree
    INTEGER(PREC_QTIDX), INTENT(IN) :: nnquad
!</input>

!<inputoutput>
    ! quadtree that should be resized
    TYPE(t_quadtree) :: rquadtree
!</inputoutput>
!</subroutine>

    CALL storage_realloc('t_quadtree_resize',nnquad,rquadtree%h_Dbbox,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_realloc('t_quadtree_resize',nnquad,rquadtree%h_Kquad,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_double2D(rquadtree%h_Dbbox,rquadtree%p_Dbbox)
    CALL storage_getbase_int2D(rquadtree%h_Kquad,rquadtree%p_Kquad)
    
    rquadtree%NNQUAD = nnquad
    rquadtree%NRESIZE=rquadtree%NRESIZE+1
  END SUBROUTINE t_quadtree_resize_nquad

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_quadtree_copyfrom(rquadtree,h_Ddata)

!<description>
    ! This subroutine copies the content of the quadtree to a handle
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
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dcorvg,p_DcorvgTmp
    
    ! Transform the content of the quadtree to h_Ddata
    IF (h_Ddata /= ST_NOHANDLE) THEN
      CALL storage_realloc('t_quadtree_copy',rquadtree%NVT,h_Ddata,ST_NEWBLOCK_NOINIT,.TRUE.)
    ELSE
      CALL storage_new('t_quadtree_copy','p_Ddata',(/2,rquadtree%NVT/),ST_DOUBLE,h_Ddata,ST_NEWBLOCK_NOINIT)
    END IF
    CALL storage_getbase_double2D(h_Ddata,p_Dcorvg)
    CALL storage_getbase_double2D(rquadtree%h_Ddata,p_DcorvgTmp)
    CALL DCOPY(2*rquadtree%NVT,p_DcorvgTmp,1,p_Dcorvg,1)
  END SUBROUTINE t_quadtree_copyfrom

  !************************************************************************

!<subroutine>
  
  SUBROUTINE t_quadtree_copyto(h_DdataSrc,rquadtree)

!<description>
    ! This subroutine copies the content of a handle
    ! to the quadtree and vice versa
!</description>

!<input>
    ! Handle to the coordinate vector
    INTEGER, INTENT(IN) :: h_DdataSrc
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata
    INTEGER(PREC_QTIDX) :: ivt,jvt,iquad,ipos

    ! Transform the content of h_Ddata to the quadtree
    CALL storage_getbase_double2D(h_DdataSrc,p_Ddata)
    DO ivt=1,SIZE(p_Ddata,2)
      IF (search(rquadtree,p_Ddata(:,ivt),iquad,ipos,jvt) ==&
          QNOT_FOUND) CALL insert(rquadtree,ivt,p_Ddata(:,ivt),iquad)
    END DO
  END SUBROUTINE t_quadtree_copyto

  !************************************************************************
  
!<subroutine>
  
  SUBROUTINE t_quadtree_insert(rquadtree,ivt,dcorvg,iquad)

!<description>
    ! This subroutine inserts a new coordinate item to the quadtree
!</description>

!<input>
    ! Number of the inserted vertex
    INTEGER(PREC_QTIDX), INTENT(IN) :: ivt

    ! Number of the quad to which vertex is inserted
    INTEGER(PREC_QTIDX), INTENT(IN) :: iquad
    
    ! Coordinates of the new vertex
    REAL(DP), DIMENSION(2), INTENT(IN) :: dcorvg
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>
!</subroutine>
    
    ! Check if there is enough space left in the nodal component of the quadtree
    IF (rquadtree%NVT == rquadtree%NNVT) CALL resizeNVT(rquadtree,CEILING(rquadtree%dfactor*rquadtree%NNVT))
    
    ! Update values and add new entry recursively
    rquadtree%NVT           = rquadtree%NVT+1
    rquadtree%p_Ddata(:,ivt) = dcorvg
    CALL insert(ivt,iquad)
    
  CONTAINS  
    
    !**************************************************************
    ! Here, the recursive insertion routine follows
    
    RECURSIVE SUBROUTINE insert(ivt,iquad)
      INTEGER(PREC_QTIDX), INTENT(IN) :: ivt,iquad
      REAL(DP) :: xmin,ymin,xmax,ymax,xmid,ymid
      INTEGER(PREC_QTIDX) :: i,jquad,nquad
      
      IF (rquadtree%p_Kquad(QSTATUS,iquad) == QMAX) THEN
        
        IF (rquadtree%nquad+4 > rquadtree%nnquad) CALL resizeNQUAD(rquadtree,CEILING(rquadtree%dfactor*rquadtree%NNQUAD))
        
        ! Quad is full and needs to be refined into four new quads
        xmin = rquadtree%p_Dbbox(QXMIN,iquad)
        ymin = rquadtree%p_Dbbox(QYMIN,iquad)
        xmax = rquadtree%p_Dbbox(QXMAX,iquad)
        ymax = rquadtree%p_Dbbox(QYMAX,iquad)
        xmid = (xmin+xmax)/2._DP
        ymid = (ymin+ymax)/2._DP
        
        ! Store the number of quads
        nquad           = rquadtree%NQUAD
        rquadtree%NQUAD = nquad+4
        
        ! NW-quad
        rquadtree%p_Kquad(QSTATUS,nquad+QNW) = QEMPTY
        rquadtree%p_Kquad(QPARENT,nquad+QNW) = iquad
        rquadtree%p_Kquad(QPARPOS,nquad+QNW) = QNW
        rquadtree%p_Dbbox(:,nquad+QNW)       = (/xmin,ymid,xmid,ymax/)
        
        ! SW-quad
        rquadtree%p_Kquad(QSTATUS,nquad+QSW) = QEMPTY
        rquadtree%p_Kquad(QPARENT,nquad+QSW) = iquad
        rquadtree%p_Kquad(QPARPOS,nquad+QSW) = QSW
        rquadtree%p_Dbbox(:,nquad+QSW)       = (/xmin,ymin,xmid,ymid/)
        
        ! SE-quad
        rquadtree%p_Kquad(QSTATUS,nquad+QSE) = QEMPTY
        rquadtree%p_Kquad(QPARENT,nquad+QSE) = iquad
        rquadtree%p_Kquad(QPARPOS,nquad+QSE) = QSE
        rquadtree%p_Dbbox(:,nquad+QSE)       = (/xmid,ymin,xmax,ymid/)
        
        ! NE-quad
        rquadtree%p_Kquad(QSTATUS,nquad+QNE) = QEMPTY
        rquadtree%p_Kquad(QPARENT,nquad+QNE) = iquad
        rquadtree%p_Kquad(QPARPOS,nquad+QNE) = QNE
        rquadtree%p_Dbbox(:,nquad+QNE)       = (/xmid,ymid,xmax,ymax/)
        
        ! Add the four values from IQUAD to the four new quads 
        ! NQUAD+1:NQUAD+4 recursively
        DO i=1,QMAX
          jquad=nquad+direction(rquadtree,rquadtree%p_Ddata(:,rquadtree%p_Kquad(i,iquad)),iquad)
          CALL insert(rquadtree%p_Kquad(i,iquad),jquad)
        END DO
        
        ! Mark the current quad as subdivided and set pointers to its four children 
        rquadtree%p_Kquad(QSTATUS,iquad) = QSUBDIV
        rquadtree%p_Kquad(1:4,iquad)     = (/nquad+QNW,nquad+QSW,nquad+QSE,nquad+QNE/)
        
        ! Add the new entry to the next position recursively
        jquad = nquad+direction(rquadtree,rquadtree%p_Ddata(:,ivt),iquad)
        CALL insert(ivt,jquad)
        
      ELSE
        
        ! Quad is not full, so new items can be stored
        rquadtree%p_Kquad(QSTATUS,iquad) = rquadtree%p_Kquad(QSTATUS,iquad)+1
        rquadtree%p_Kquad(rquadtree%p_Kquad(QSTATUS,iquad),iquad) = ivt
      END IF
    END SUBROUTINE insert
  END SUBROUTINE t_quadtree_insert
  
  ! ***************************************************************************
  
!<function>
  
  FUNCTION t_quadtree_delete(rquadtree,dcorvg,ivt) RESULT(f)

!<description>
    ! This function deletes an item from the quadtree
!</description>

!<input>
    ! Coordinates of the vertex that should be deleted
    REAL(DP), DIMENSION(2), INTENT(IN) :: dcorvg
!</input>

!<inputoutput>
    ! quadtree
    TYPE(t_quadtree), INTENT(INOUT) :: rquadtree
!</inputoutput>

!<output>
    ! Number of the vertex that is deleted
    INTEGER(PREC_QTIDX), INTENT(OUT) :: ivt
!</output>

!<result>
    ! Result of the deletion: QNOT_FOUND / QFOUND
    INTEGER :: f
!</result>
!</function>
    
    ! local variables
    INTEGER :: iquad,ipos,jpos,jvt
    
    ! Search for the given coordinates
    f=search(rquadtree,dcorvg,iquad,ipos,ivt)
    
    ! What can we do from the searching
    IF (f == QFOUND) THEN
      
      ! Remove item IVT from quad IQUAD
      DO jpos=ipos+1,rquadtree%p_Kquad(QSTATUS,iquad)
        rquadtree%p_Kquad(jpos-1,iquad) = rquadtree%p_Kquad(jpos,iquad)
      END DO
      rquadtree%p_Kquad(rquadtree%p_Kquad(QSTATUS,iquad),iquad) = 0
      rquadtree%p_Kquad(QSTATUS,iquad) = rquadtree%p_Kquad(QSTATUS,iquad)-1
      
      ! If IVT is not last item move last item NVT to position IVT
      IF (ivt /= rquadtree%NVT) THEN
        IF (search(rquadtree,rquadtree%p_Ddata(1:2,rquadtree%NVT),iquad,ipos,jvt) == QFOUND) THEN
          rquadtree%p_Ddata(:,ivt) = rquadtree%p_Ddata(:,rquadtree%NVT)
          rquadtree%p_Kquad(ipos,iquad) = ivt
        END IF
        ivt=rquadtree%NVT
      END IF
      rquadtree%NVT = rquadtree%NVT-1
    END IF
  END FUNCTION t_quadtree_delete
  
  ! ***************************************************************************

!<function>
  
  FUNCTION t_quadtree_search(rquadtree,dcorvg,iquad,ipos,ivt) RESULT(f)

!<description>
    ! This subroutine searches for given coordinates in the quadtree
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
    
    ! coordinates that should be searched
    REAL(DP), DIMENSION(2), INTENT(IN) :: dcorvg
!</input>

!<output>
    ! Number of the quad in which the given coordinates are
    INTEGER(PREC_QTIDX), INTENT(OUT) :: iquad

    ! Position of the coordinates in the quad
    INTEGER(PREC_QTIDX), INTENT(OUT) :: ipos

    ! Number of the vertex the coordinates correspond to
    INTEGER(PREC_QTIDX), INTENT(OUT) :: ivt
!</output>

!<result>
    ! Result of the searching: QNOT_FOUND / QFOUND
    INTEGER :: f
!</result>
!</function>
    
    ! Initialize
    iquad=1; ipos=1; ivt=1; f=search(iquad,ipos,ivt)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive searching routine follows

    RECURSIVE FUNCTION search(iquad,ipos,ivt) RESULT(f)
      INTEGER(PREC_QTIDX), INTENT(INOUT) :: iquad,ipos,ivt
      INTEGER :: f
      
      f=QNOT_FOUND
      IF (rquadtree%p_Kquad(QSTATUS,iquad) == QSUBDIV) THEN
        
        ! Quad is subdivided. Compute child IQUAD which to look recursively.
        iquad = rquadtree%p_Kquad(direction(rquadtree,dcorvg,iquad),iquad)
        f=search(iquad,ipos,ivt)
        
      ELSE
        
        ! Quad is not subdivided. Search for (x,y) in current quad
        DO ipos=1,rquadtree%p_Kquad(QSTATUS,iquad)
          ivt = rquadtree%p_Kquad(ipos,iquad)
          IF (MAXVAL(ABS(rquadtree%p_Ddata(:,ivt)-dcorvg)) <= SYS_EPSREAL) THEN
            f=QFOUND; RETURN
          END IF
        END DO
        
      END IF
    END FUNCTION search
  END FUNCTION t_quadtree_search

  !************************************************************************
  
!<function>
  
  PURE FUNCTION t_quadtree_direction(rquadtree,dcorvg,iquad) RESULT(d)

!<description>
    ! This subroutine determines the direction to preceed w.r.t. dcorvg
!</description>

!<input>
    ! quadtree
    TYPE(t_quadtree), INTENT(IN) :: rquadtree
    
    ! Coordinates
    REAL(DP), DIMENSION(2), INTENT(IN) :: dcorvg

    ! Number of quad
    INTEGER(PREC_QTIDX), INTENT(IN) :: iquad
!</input>

!<result>
    ! Further search direction
    INTEGER :: d
!</result>
!</function>
    
    ! local variables
    REAL(DP) :: xmid,ymid

    ! Compute midpoint of current quad
    xmid=(rquadtree%p_Dbbox(QXMIN,iquad)+rquadtree%p_Dbbox(QXMAX,iquad))/2._DP
    ymid=(rquadtree%p_Dbbox(QYMIN,iquad)+rquadtree%p_Dbbox(QYMAX,iquad))/2._DP
    
    IF (dcorvg(1) > xmid) THEN
      IF (dcorvg(2) > ymid) THEN
        d=QNE; RETURN
      ELSE
        d=QSE; RETURN
      END IF
    ELSE
      IF (dcorvg(2) > ymid) THEN
        d=QNW; RETURN
      ELSE
        d=QSW; RETURN
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
    xmin = rquadtree%p_Dbbox(QXMIN,1)
    ymin = rquadtree%p_Dbbox(QYMIN,1)
    xmax = rquadtree%p_Dbbox(QXMAX,1)
    ymax = rquadtree%p_Dbbox(QYMAX,1)
    CALL print(xmin,ymin,xmax,ymax,1)
    CLOSE(UNIT=iunit)
    
  CONTAINS
    
    !**************************************************************
    ! Here, the recursive print routine follows
    
    RECURSIVE SUBROUTINE print(xmin,ymin,xmax,ymax,iquad)
      REAL(DP), INTENT(IN) :: xmin,ymin,xmax,ymax
      INTEGER(PREC_QTIDX), INTENT(IN) :: iquad
      REAL(DP) :: xmid,ymid
      INTEGER(PREC_QTIDX) :: i

      WRITE(UNIT=iunit,FMT=*) 'rect'
      WRITE(UNIT=iunit,FMT=10) xmin,ymin,xmax,ymax
      
      IF (rquadtree%p_Kquad(QSTATUS,iquad) == QSUBDIV) THEN
        
        xmid=(xmin+xmax)/2._DP
        ymid=(ymin+ymax)/2._DP
        
        CALL print(xmin,ymid,xmid,ymax,rquadtree%p_Kquad(QNW,iquad))
        CALL print(xmin,ymin,xmid,ymid,rquadtree%p_Kquad(QSW,iquad))
        CALL print(xmid,ymin,xmax,ymid,rquadtree%p_Kquad(QSE,iquad))
        CALL print(xmid,ymid,xmax,ymax,rquadtree%p_Kquad(QNE,iquad))
        
      ELSEIF (rquadtree%p_Kquad(QSTATUS,iquad) > QEMPTY) THEN
        
        DO i=1,rquadtree%p_Kquad(QSTATUS,iquad)
          WRITE(UNIT=iunit,FMT=*) 'node'
          WRITE(UNIT=iunit,FMT=20) rquadtree%p_Ddata(:,rquadtree%p_Kquad(i,iquad))
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
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Kquad =',rquadtree%h_Kquad
    WRITE(*,*)
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8)') '  NVT     =',rquadtree%NVT,'NNVT     =',rquadtree%NNVT
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8)') '  NQUAD   =',rquadtree%NQUAD,'NNQUAD   =',rquadtree%NNQUAD
    WRITE(*,FMT='(1X,A,1X,I8,A,4X,F5.1,A,4X,F5.1,A)') '  NRESIZE =',rquadtree%NRESIZE, "   QUADS    =", &
        100*rquadtree%NQUAD/REAL(rquadtree%NNQUAD,DP),'%  FILLING  =',100*rquadtree%NVT/REAL(rquadtree%NNVT,DP),'%'
    WRITE(*,*)
  END SUBROUTINE t_quadtree_info

  !************************************************************************

!<function>

  FUNCTION t_quadtree_getboundingbox(rquadtree,iquad) RESULT(bbox)
    
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
    INTEGER, INTENT(IN), OPTIONAL :: iquad
!</input>

!<result>
    ! bounding box
    REAL(DP), DIMENSION(4) :: bbox
!</result>
!</function>
    
    IF (PRESENT(iquad)) THEN
      IF (iquad > rquadtree%NVT) THEN
        PRINT *, "(EE) t_quadtree_getboundingbox: quad number exceeds quadtree dimension"
        STOP
      END IF
      bbox=rquadtree%p_Dbbox(QXMIN:QYMAX,iquad)
    ELSE
      bbox=rquadtree%p_Dbbox(QXMIN:QYMAX,1)
    END IF
  END FUNCTION t_quadtree_getboundingbox
END MODULE quadtree

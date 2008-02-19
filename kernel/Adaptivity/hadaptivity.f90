!##############################################################################
!# ****************************************************************************
!# <name> hadaptivity </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform h-adaptivity,
!# namely, grid refinement and grid coarsening.
!# In order to apply h-adaptivity, the initial mesh must be conforming, 
!# that is, it is allowed to contain a mixture of triangular and quadrilateral
!# elements without 'hanging' nodes. After one grid adaptivity step, the
!# resulting mesh is also conforming without 'hanging' nodes.
!#
!# This module makes extensive use of dynamic data structures such as quadtrees,
!# binary search trees and so-called arraylists which need to be generated from 
!# the static mesh structure. After grid adaptivity, the dynamic data needs to 
!# be reconverted to the static structures required for the simulation.
!#
!# Note that this module is implemented as a 'black-box' tool for grid adaptivity.
!# One of the building blocks is the t_hadapt data structure which provides all
!# required information. 
!#
!# The following routines are available:
!#
!#  1.) hadapt_initFromParameterlist
!#      -> Initialize adaptivity structure from parameterlist
!#
!#  2.) hadapt_initFromTriangulation
!#      -> Initialize adaptivity structure from triangulation structure
!#
!#  3.) hadapt_generateRawMesh
!#      ->  Generate the raw mesh from the adaptivity structure
!#
!#  4.) hadapt_releaseAdaptation
!#      -> Release all internal adaptation structures
!#
!#  5.) hadapt_duplicateAdaptation
!#      -> Create a duplicate / backup of an adaptivity structure.
!#
!#  6.) hadapt_restoreAdaptation
!#      -> Restores an adaptivity structure previously backed up with 
!#         hadapt_duplicateAdapation
!#
!#  7.) hadapt_setVertexCoords2D
!#      -> Set the coordinates of vertices to the adaptivity structure in 2D
!#
!#  8.) hadapt_getVertexCoords2D
!#      -> Get the coordinates of vertices from the adaptivity structure in 2D
!#
!#  9.) hadapt_setVertexCoords3D
!#      -> Set the coordinates of vertices to the adaptivity structure in 3D
!#
!# 10.) hadapt_getVertexCoords3D
!#      -> Get the coordinates of vertices from the adaptivity structure in 3D
!#
!# 11.) hadapt_setVerticesAtElement
!#      -> Set the "vertices-at-element" structure to the adaptivity structure
!#
!# 12.) hadapt_getVerticesAtElement
!#      -> Get the "vertices-at-element" structure from the adaptivity structure
!#
!# 13.) hadapt_setNeighboursAtElement
!#      -> Set the "neighbours-at-element" structure to the adaptivity structure
!#
!# 14.) hadapt_getNeighboursAtElement
!#      -> Get the "neighbours-at-element" structure from the adaptivity structure
!#
!# 15.) hadapt_setNelOfType
!#      -> Set the "InelOfType" array to the adaptivity structure
!#
!# 16.) hadapt_getNelOfType
!#      -> Get the "InelOfType" array from the adaptivity structure
!#
!# 17.) hadapt_setBoundary
!#      -> Set the boundary structure to the adaptivity structure
!#
!# 18.) hadapt_getBoundary
!#      -> Get the boundary structure form the adaptivity structure
!#
!# 19.) hadapt_setNodalProperty
!#      -> Set the "nodal property" list to the adaptivity structure
!#
!# 20.) hadapt_getNodalProperty
!#      -> Get the "nodal property" list from the adaptivity structure
!#
!# 21.) hadapt_genElementsAtVertex
!#      -> Generate the "elements-at-vertex" data structure
!#
!# 22.) hadapt_performAdaptation
!#      -> perform one step of grid adaptation
!#
!# 23.) hadapt_infoStatistics
!#      -> output information about the adaptivity structure
!#
!# 24.) hadapt_writeGridSVG
!#      -> write the adapted grid to file in SVG format
!#
!# 25.) hadapt_writeGridGMV
!#      -> write the adapted grid to file in GMV format
!#
!# 26.) hadapt_checkConsistency
!#      -> check the internal consistency of dynamic data structures
!#
!# 27.) hadapt_refreshAdaptation
!#      -> Refresh pointers of adaptation structure
!#
!# The following internal routines are available:
!#
!# 1.) redgreen_refine
!#     -> perform red-green refinement for marked elements
!#
!# 2.) redgreen_coarsen
!#     -> perform red-green coarsening for marked elements
!#
!#
!#  FAQ - Some explainations
!# --------------------------
!# 1.) So, how does red-green grid refinement work in detail?
!#
!#     In general, a conforming mesh in two spatial dimensions consists of 
!#     triangles and/or quadrilaterals which have common vertices, that is,
!#     there is no vertex located somewhere in the edge of one element.
!#     If you want to refine your mesh, a natural choice is to subdivide
!#     (a) one triangle into four similar triangles by placing new vertices
!#         at all midpoints of the three edges and connecting these new nodes
!#     (b) one quadrlateral into four similar quadrilaterals by placing
!#         four new vertices at the edge midpoints and one vertex at the
!#         center and connect all new nodes to the centered vertex.
!#     These two strategies are known as red-refinement.
!#
!#     In order to remove the hanging nodes in adjacent elements you have
!#     to refine these elements, too. Here, all all possibilities:
!#     (c) Refine one triangle into two triangles by placing one new node
!#         at the midpoint of one edge.
!#     (d) Refine one quadrilateral into two quadrilaterals by inserting two
!#         new vertices at the midpoints of opposite edges and connect them.
!#     (e) Refine one quadrilateral into three triangles by inserting one new
!#         vertex at the midpoint of one edge and connect it to the two nodes
!#         of the opposite edge.
!#     (f) Refine one quadrilateral into four triangles by inserting two nodes
!#         at the midpoint of adjacent edges, connect these vertices with 
!#         oneanother and with the opposite corner vertex.
!#     These four refinement strategies are known as green refinement.
!#
!#     Keep in mind, that there may be two hanging nodes for triangles or
!#     three hanging nodes for quadrilaterals. If these elements were refined 
!#     without introducing new vertices, we would perform blue refinement.
!#     However, this strategy should be avoided and an element marked for
!#     blue refinement should be subdivided into four similar elements
!#     introducing one additional vertex at the edge midpoint.
!#
!#     Red refinement does not deteriorate the grid quality but green does.
!#     If one would proceed like this, then green-refined element would tend
!#     to produce very sharp angles leading to nasty "needle" elements.
!#     In order to refine a green element, the original element from which
!#     the green element was created has to be restored. Then red refinement
!#     is performed for this element. 
!#
!#
!# 2.) What are the different numbers for the element marker MARK_xxx?
!#
!#     In two space dimensions, the maximum number of vertices per element and
!#     consequently the number of edges surrounding an element is two. Consider
!#     a 5-bitfield [bit4,bit3,bit2,bit1,bit0]. The bit0 is used to distinguish
!#     between triangles (bit0=0) and quadrilaterals (bit1=1). Each of bitX is
!#     associated with the X-th edge. If you want to mark an edge for refinement,
!#     then you have to set its bit to one. If you want to unmark it, then set 
!#     bitX=0. Finally, the integer associated to this 5-bitfiled is the element
!#     marker. Here, is the table of bitfields and markers
!#
!#     bit4  bit3  bit2  bit1  bit0   integer   description
!#      0     0     0     0     0        0       triangle, not marked
!#      0     0     0     0     1        1       quadrilateral, not marked
!#      0     0     0     1     0        2       triangle, marked for 1-tria : 2-tria
!#                                               refinement along the first edge
!#      0     0     0     1     1        3       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the first edge
!#      0     0     1     0     0        4       triangle, marked for 1-tria : 2-tria
!#                                               refinement along the second edge
!#      0     0     1     0     1        5       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the second edge
!#      0     0     1     1     0        6       triangle, marked for 1-tria : 3-tria
!#                                               refinement along first/second edge
!#      0     0     1     1     1        7       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along first/second edge
!#      0     1     0     0     0        8       triangle, marked for 1-tria : 2-tria
!#                                               refinement along the third edge
!#      0     1     0     0     1        9       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the third edge
!#      0     1     0     1     0       10       triangle, marked for 1-tria : 3-tria
!#                                               refinement along first/third edge
!#      0     1     0     1     1       11       quadrilateral, makred for 1-quad :
!#                                               2-quad along first/third edge
!#      0     1     1     0     0       12       triangle, marked for 1-tria : 3-tria
!#                                               refinement along second/third edge
!#      0     1     1     0     1       13       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along second/third edge
!#      0     1     1     1     0       14       triangle marked for 1:4 red refinement
!#      0     1     1     1     1       15       quadrilateral, marked for blue refinement
!#      1     0     0     0     0       16       n.a.
!#      1     0     0     0     1       17       quadtrilateral, marked for 1-quad :
!#                                               3-tria refinement along the fourth edge
!#      1     0     0     1     0       18       n.a
!#      1     0     0     1     1       19       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along first/fourth edge
!#      1     0     1     0     0       20       n.a.
!#      1     0     1     0     1       21       quadrilateral, makred for 1-quad :
!#                                               2-quad along second/fourth edge
!#      1     0     1     1     0       22       n.a.
!#      1     0     1     1     1       23       quadrilateral, marked for blue refinement
!#      1     1     0     0     0       24       n.a.
!#      1     1     0     0     1       25       quadrilateral, marked for 1-quad :
!#                                               4-tria refinement along third/fourth edge
!#      1     1     0     1     0       26       n.a.
!#      1     1     0     1     1       27       quadrilateral, marked for blue refinement
!#      1     1     1     0     0       28       n.a.
!#      1     1     1     0     1       29       quadrilateral, marked for blue refinement
!#      1     1     1     1     0       30       n.a.
!#      1     1     1     1     1       31       quadrilateral marked for 1:4 red refinement
!#
!#
!# 3.) Ok, and what does the state of an element mean?
!#
!#     The state of an element is computed from the age of its vertices.
!#     In 2D, an element can have four vertices at most
!#
!#     Typically, a more or less complicated
!#     data structure is required to maintain the father-son relationship.
!#     In this implementaiton, we use a new strategy based on element states
!#     which allows for an efficient computation of father-son relations.
!#     For each vertex, we save its creation age. All vertices of the initial
!#     grid have age 0 and cannot be removed. If a new vertex is introduced
!#     its age is computed as 1+MAX("Age of surrounding vertices"), that is,
!#     1+MAX("Age of left neighbor","Age of right neighbour") for edge midpoints
!#     and 1+MAX("Age of four corner vertices") for centered node.
!#
!# </purpose>
!##############################################################################

MODULE hadaptivity

  USE arraylist
  USE binarytree
  USE collection
  USE hadaptaux3d
  USE io
  USE linearsystemscalar
  USE list
  USE octree
  USE paramlist
  USE quadtree
  USE sort
  USE storage
  USE triangulation
  USE fsystem

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_hadapt
  PUBLIC :: hadapt_initFromParameterlist
  PUBLIC :: hadapt_initFromTriangulation
  PUBLIC :: hadapt_generateRawMesh
  PUBLIC :: hadapt_releaseAdaptation
  PUBLIC :: hadapt_duplicateAdaptation
  PUBLIC :: hadapt_restoreAdaptation
  PUBLIC :: hadapt_setVertexCoords2D
  PUBLIC :: hadapt_getVertexCoords2D
  PUBLIC :: hadapt_setVerticesAtElement
  PUBLIC :: hadapt_getVerticesAtElement
  PUBLIC :: hadapt_setNeighboursAtElement
  PUBLIC :: hadapt_getNeighboursAtElement
  PUBLIC :: hadapt_setNelOfType
  PUBLIC :: hadapt_getNelOfType
  PUBLIC :: hadapt_setBoundary
  PUBLIC :: hadapt_getBoundary
  PUBLIC :: hadapt_setNodalProperty
  PUBLIC :: hadapt_getNodalProperty
  PUBLIC :: hadapt_genElementsAtVertex
  PUBLIC :: hadapt_performAdaptation
  PUBLIC :: hadapt_infoStatistics
  PUBLIC :: hadapt_writeGridSVG
  PUBLIC :: hadapt_writeGridGMV
  PUBLIC :: hadapt_checkConsistency
  PUBLIC :: hadapt_refreshAdaptation

CONTAINS
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_initFromParameterlist(rhadapt, rparlist, ssection)

!<description>
    ! This subroutine initializes the adaptivity structure
    ! with the values supplied by the parameter list
!</description>

!<input>
    ! parameter list
    TYPE(t_parlist), INTENT(IN)  :: rparlist

    ! name of the section
    CHARACTER(LEN=*), INTENT(IN) :: ssection
!</input>

!<output>
    ! adaptivity structure
    TYPE(t_hadapt), INTENT(OUT)  :: rhadapt
!</output>
!</subroutine>

    ! Get mandatory parameters from list
    CALL parlst_getvalue_int   (rparlist,&
        ssection, "nsubdividemax", rhadapt%nsubdividemax)
    CALL parlst_getvalue_int   (rparlist,&
        ssection, "iadaptationStrategy", rhadapt%iadaptationStrategy)
    CALL parlst_getvalue_double(rparlist,&
        ssection, "drefinementTolerance", rhadapt%drefinementTolerance)
    CALL parlst_getvalue_double(rparlist,&
        ssection, "dcoarseningTolerance", rhadapt%dcoarseningTolerance)

    ! Initialize data
    rhadapt%iSpec = HADAPT_HAS_PARAMETERS
  END SUBROUTINE hadapt_initFromParameterlist

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_initFromTriangulation(rhadapt,rtriangulation)

!<description>
    ! This subroutine initializes all required components of the adaptativit 
    ! structure from the triangulation structure rtriangulation.
!</description>

!<input>
    ! Triangulation structure
    TYPE(t_triangulation), INTENT(IN) :: rtriangulation
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)     :: rhadapt
!</inputoutput>
!</subroutine>

    ! Initialize duplication flag
    rhadapt%iduplicationFlag=0

    ! Set dimension
    rhadapt%ndim=rtriangulation%ndim

    ! Set coordinates
    SELECT CASE(rhadapt%ndim)
    CASE(NDIM2D)
      CALL hadapt_setVertexCoords2D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE(NDIM3D)
      CALL hadapt_setVertexCoords3D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE DEFAULT
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_initFromTriangulation')
      CALL sys_halt()
    END SELECT

    ! Set nodal properties
    CALL hadapt_setNodalProperty(rhadapt,rtriangulation%h_InodalProperty)

    ! Set element numbers
    CALL hadapt_setNelOfType(rhadapt,rtriangulation%InelOfType)
    
    ! Set vertices at element
    CALL hadapt_setVerticesAtElement(rhadapt,&
        rtriangulation%h_IverticesAtElement,rtriangulation%NEL)

    ! Set elements adjacent to element
    CALL hadapt_setNeighboursAtElement(rhadapt,&
        rtriangulation%h_IneighboursAtElement)
    
    ! Set boundary
    CALL hadapt_setBoundary(rhadapt,rtriangulation%h_IboundaryCpIdx,&
        rtriangulation%h_IverticesAtBoundary,&
        rtriangulation%h_DvertexParameterValue,&
        rtriangulation%NBCT,rtriangulation%NVBD)

    ! Generate "elements-meeting-at-vertex" structure
    CALL hadapt_genElementsAtVertex(rhadapt)
    
    ! Create generation array and initialize all nodes with "age" 0
    CALL storage_new('hadapt_initFromTriangulation','p_IvertexAge',&
        rhadapt%NVT,ST_INT,rhadapt%h_IvertexAge,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(rhadapt%h_IvertexAge,rhadapt%p_IvertexAge)
  END SUBROUTINE hadapt_initFromTriangulation

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_generateRawMesh(rhadapt,rtriangulation)

!<description>
    ! This subroutine generates a raw mesh from the adaptivity structure
!</description>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)        :: rhadapt
    
    ! Triangulation structure
    TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
!</subroutine>

    ! Get dimensions
    rtriangulation%ndim=rhadapt%ndim

    ! Get coordinates
    SELECT CASE(rhadapt%ndim)
    CASE(NDIM2D)
      CALL hadapt_getVertexCoords2D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE(NDIM3D)
      CALL hadapt_getVertexCoords3D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE DEFAULT
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_generateRawMesh')
      CALL sys_halt()
    END SELECT

    ! Get number of elements
    CALL hadapt_getNelOfType(rhadapt,rtriangulation%InelOfType)

    ! Get vertices at element list
    CALL hadapt_getVerticesAtElement(rhadapt,&
        rtriangulation%h_IverticesAtElement,rtriangulation%NEL)

    ! Get element neighbours
    CALL hadapt_getNeighboursAtElement(rhadapt,rtriangulation%h_IneighboursAtElement)

    ! Get nodal property list
    CALL hadapt_getNodalProperty(rhadapt,rtriangulation%h_InodalProperty)

    ! Get boundary
    CALL hadapt_getBoundary(rhadapt,rtriangulation%h_IboundaryCpIdx,&
        rtriangulation%h_IverticesAtBoundary,rtriangulation%h_DvertexParameterValue,&
        rtriangulation%NVBD,rtriangulation%NBCT)
  END SUBROUTINE hadapt_generateRawMesh
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_releaseAdaptation(rhadapt)

!<description>
    ! This subroutine releases all internal structures of the
    ! adaptivity data structure rhadapt.
!</description>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32) :: idupflag
    INTEGER :: ibct

    idupflag = rhadapt%iduplicationFlag

    ! Check if quadtree exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).EQ.HADAPT_HAS_COORDS) THEN
      SELECT CASE(rhadapt%ndim)
      CASE(NDIM2D)
        CALL qtree_releaseQuadtree(rhadapt%rVertexCoordinates2D)
        
      CASE(NDIM3D)
        CALL otree_releaseOctree(rhadapt%rVertexCoordinates3D)

      CASE DEFAULT
        CALL output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_releaseAdaptation')
        CALL sys_halt()
      END SELECT
    END IF   

    ! Check if boundary structure exists
    IF (ASSOCIATED(rhadapt%rBoundary)) THEN
      DO ibct=1,SIZE(rhadapt%rBoundary,1)
        CALL btree_releaseTree(rhadapt%rBoundary(ibct))
      END DO
      DEALLOCATE(rhadapt%rBoundary)
      NULLIFY(rhadapt%rBoundary)
    END IF

    ! Release elements-meeting-at-vertex arraylist
    CALL arrlst_releaseArraylist(rhadapt%relementsAtVertex)

    ! Release storage which is no longer in use
    CALL checkAndRelease(idupflag, HADAPT_SHARE_IMARKER,&
        rhadapt%h_Imarker)
    CALL checkAndRelease(idupflag, HADAPT_SHARE_IVERTEXAGE,&
        rhadapt%h_IvertexAge)
    CALL checkAndRelease(idupflag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
        rhadapt%h_ImidneighboursAtElement)
    
    ! Nullify "performance-pointers"
    NULLIFY(rhadapt%p_IvertexAge)
    NULLIFY(rhadapt%p_InodalProperty)
    NULLIFY(rhadapt%p_IverticesAtElement)
    NULLIFY(rhadapt%p_IneighboursAtElement)
    NULLIFY(rhadapt%p_ImidneighboursAtElement)
    
    ! Clear parameters
    rhadapt%iadaptationStrategy  = HADAPT_NOADAPTATION
    rhadapt%drefinementTolerance = 0._DP
    rhadapt%dcoarseningTolerance = 0._DP

    ! Clear data
    rhadapt%iSpec            = HADAPT_UNDEFINED
    rhadapt%nRefinementSteps = 0
    rhadapt%nCoarseningSteps = 0
    rhadapt%nSmoothingSteps  = 0
    rhadapt%nGreenElements   = 0
    rhadapt%ndim             = 0
    rhadapt%NVT              = 0
    rhadapt%NVT0             = 0
    rhadapt%increaseNVT      = 0
    rhadapt%NVBD             = 0
    rhadapt%NVBD0            = 0
    rhadapt%NBCT             = 0
    rhadapt%NEL              = 0
    rhadapt%NEL0             = 0
    rhadapt%NELMAX           = 0
    rhadapt%InelOfType       = 0
    rhadapt%InelOfType0      = 0

  CONTAINS
    
    SUBROUTINE checkAndRelease (idupFlag,ibitfield,ihandle)
      INTEGER(I32), INTENT(IN) :: ibitfield
      INTEGER(I32), INTENT(IN) :: idupFlag
      INTEGER, INTENT(INOUT) :: ihandle
      
      IF (IAND(idupFlag,ibitfield) .NE. ibitfield) THEN
        IF (ihandle .NE. ST_NOHANDLE) CALL storage_free(ihandle)
      ELSE
        ihandle = ST_NOHANDLE
      END IF
      
    END SUBROUTINE checkAndRelease
  END SUBROUTINE hadapt_releaseAdaptation

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_duplicateAdaptation(rhadapt,rhadaptBackup,&
      iduplicationFlag,bupdate)

!<description>
    ! This subroutine makes a copy of an adaptivity structure in memory. The 
    ! variable iduplicationFlag decides on which arrays are copied in memory
    ! and which are not.
    !
    ! By setting the corresponding bit in iduplicationFlag to 0, the array is
    ! duplicated in memory, and any change to the new array will not harm
    ! the original one.
    ! By setting a flag HADAPT_SHARE_xxxx in iduplicationFlag, the corresponding 
    ! array is not duplicated. The handle of the original structure is simply put
    ! into the new structure to make the information accessable. Then this array
    ! is shared between two adaptivity structures!
!</description>

!<input>
    ! The source structure which provides all information
    TYPE(t_hadapt), INTENT(IN) :: rhadapt

    ! Bitfield that decides which handles are a copy of another structure, thus 
    ! which arrays are shared between the new and the old structure. 
    ! The bitfield contains a combination of HADAPT_SHARE_xxxx constants. Every 
    ! information whose flag is set in iduplicationFlag is shared between rhadapt
    ! and rhadaptBackup.
    ! Therefore e.g., iduplicationFlag=0 copies all arrays from rhadapt in memory,
    ! while HADAPT_SHARE_ALL will copy nothing, but will share everything between 
    ! rhadapt and rhadaptBackup.
    INTEGER(I32), INTENT(IN)   :: iduplicationFlag
  
    ! OPTIONAL. Defines how to create the backup.
    ! = .FALSE.: Treat rhadaptBackup as empty destination structure. If necessary, 
    !    information in rhadaptBackup is released. rhadaptBackup is rebuild 
    !    according to rhadapt and iduplicationFlag. 
    !    This is the standard setting if bupdate is not specified.
    ! = .TRUE. : Treat rhadaptBackup as existing copy of rhadapt which has to be 
    !    updated. Recover all arrays of rhadaptBackup by those of rhadapt.
    !    (I.e. those arrays which were duplicated by a previous
    !    call to iduplicationFlag with IUPD=0.)
    !    iduplicationFlag can used to specify which data do copy
    !    from rhadapt to rhadaptBackup. It is OR'ed with rhadaptBackup%iduplicationFlag 
    !    to get the actual duplication flag. This leads to the following interpretation
    !    of iduplicationFlag:
    !     =0:  Copy all data that was copied previously. This is the usual setting.
    !    <>0:  Copy all arrays where the corresponding flag in iduplicationFlag is not
    !          set and which exist as duplicates. Arrays corresponding to flags which
    !          are set or where the handles in rhadapt and rhadaptBackup coincide 
    !          are not touched.
    LOGICAL, INTENT(IN), OPTIONAL     :: bupdate
!</input>

!<inputoutput>
    ! The destination structure which receives the information
    TYPE(t_hadapt), INTENT(INOUT) :: rhadaptBackup
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(I32) :: idupFlag
    LOGICAL :: bupd

    bupd = .FALSE.
    IF (PRESENT(bupdate)) bupd = bupdate
    
    IF (.NOT. bupd) THEN
      ! Release any old data.
      CALL hadapt_releaseAdaptation(rhadaptBackup)

      rhadaptBackup%iSpec                = rhadapt%iSpec
      rhadaptBackup%iadaptationStrategy  = rhadapt%iadaptationStrategy
      rhadaptBackup%NSUBDIVIDEMAX        = rhadapt%NSUBDIVIDEMAX
      rhadaptBackup%nRefinementSteps     = rhadapt%nRefinementSteps
      rhadaptBackup%nCoarseningSteps     = rhadapt%nCoarseningSteps
      rhadaptBackup%nSmoothingSteps      = rhadapt%nSmoothingSteps
      rhadaptBackup%drefinementTolerance = rhadapt%dRefinementTolerance
      rhadaptBackup%dCoarseningTolerance = rhadapt%dCoarseningTolerance
      rhadaptBackup%ndim                 = rhadapt%ndim
      rhadaptBackup%NVT0                 = rhadapt%NVT0
      rhadaptBackup%NVT                  = rhadapt%NVT
      rhadaptBackup%increaseNVT          = rhadapt%increaseNVT
      rhadaptBackup%NVBD0                = rhadapt%NVBD0
      rhadaptBackup%NVBD                 = rhadapt%NVBD
      rhadaptBackup%NBCT                 = rhadapt%NBCT
      rhadaptBackup%NEL0                 = rhadapt%NEL0
      rhadaptBackup%NEL                  = rhadapt%NEL
      rhadaptBackup%NELMAX               = rhadapt%NELMAX
      rhadaptBackup%nGreenElements       = rhadapt%nGreenElements
      rhadaptBackup%InelOfType0          = rhadapt%InelOfType0
      rhadaptBackup%InelOfType           = rhadapt%InelOfType

      ! Decide on iduplicationFlag which arrays to copy
      rhadaptBackup%iduplicationFlag = iduplicationFlag
      idupFlag = iduplicationFlag

    ELSE

      ! Create a bitfiled what to copy by ORing iduplicationFlag with what
      ! we have in rhadaptBackup. That way, only arrays that exist as real
      ! duplicates are copied from rhadapt to rhadaptBackup.

      idupFlag = IOR(iduplicationFlag,rhadaptBackup%iduplicationFlag)

    END IF

    ! Call checkAndCopy for all components. This will either copy the handle
    ! or allocate new memory and copy the content of the component.

    ! Bit   0: Imarker
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IMARKER,&
        rhadapt%h_Imarker,&
        rhadaptBackup%h_Imarker)

    ! Bit   1: IvertexAge
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IVERTEXAGE,&
        rhadapt%h_IvertexAge,&
        rhadaptBackup%h_IvertexAge)
    CALL storage_getbase_int(rhadaptBackup%h_IvertexAge,&
        rhadaptBackup%p_IvertexAge)

    ! Bit   2: InodalProperty
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_INODALPROPERTY,&
        rhadapt%h_InodalProperty,&
        rhadaptBackup%h_InodalProperty)
    CALL storage_getbase_int(rhadaptBackup%h_InodalProperty,&
        rhadaptBackup%p_InodalProperty)

    ! Bit   3: IverticesAtElement
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IVERTICESATELEMENT,&
        rhadapt%h_IverticesAtElement,&
        rhadaptBackup%h_IverticesAtElement)
    CALL storage_getbase_int2D(rhadaptBackup%h_IverticesAtElement,&
        rhadaptBackup%p_IverticesAtElement)

    ! Bit   4: IneighboursAtElement
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_INEIGHATELEMENT,&
        rhadapt%h_IneighboursAtElement,&
        rhadaptBackup%h_IneighboursAtElement)
    CALL storage_getbase_int2D(rhadaptBackup%h_IneighboursAtElement,&
        rhadaptBackup%p_IneighboursAtElement)

    ! Bit   5: ImidneighboursAtElement
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
        rhadapt%h_ImidneighboursAtElement,&
        rhadaptBackup%h_ImidneighboursAtElement)
    CALL storage_getbase_int2D(rhadaptBackup%h_ImidneighboursAtElement,&
        rhadaptBackup%p_ImidneighboursAtElement)

    ! Bit   6: rVertexCoordinates
    IF (IAND(idupFlag, HADAPT_SHARE_RVERTEXCOORDINATES) .NE.&
        HADAPT_SHARE_RVERTEXCOORDINATES) THEN
      
      SELECT CASE(rhadaptBackup%ndim)
      CASE(NDIM2D)
        CALL qtree_duplicateQuadtree(rhadapt%rVertexCoordinates2D,&
            rhadaptBackup%rVertexCoordinates2D)
      CASE(NDIM3D)
        CALL otree_duplicateOctree(rhadapt%rVertexCoordinates3D,&
            rhadaptBackup%rVertexCoordinates3D)
      CASE DEFAULT
        CALL output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_duplicateAdaptation')
        CALL sys_halt()
      END SELECT
    END IF
        
    ! Bit   7: rBoundary
    
    
    ! Bit   8: rElementsAtVertex
    IF (IAND(idupFlag, HADAPT_SHARE_RELEMENTSATVERTEX) .NE.&
        HADAPT_SHARE_RELEMENTSATVERTEX) THEN
      CALL arrlst_duplicateArrayList(rhadapt%rElementsAtVertex,&
          rhadaptBackup%rElementsAtVertex)
    END IF

  CONTAINS
    
    SUBROUTINE checkAndCopy (idupFlag,ibitfield,isourcehandle,idesthandle)
      
      ! Checks if idupFlag has all bits ibitfield set.
      ! If yes, idesthandle is set to isourcehandle.
      ! Otherwise, the memory behind isourcehandle is duplicated in memory
      ! and idesthandle receives the handle to the new memory block.
      
      INTEGER(I32), INTENT(IN) :: ibitfield
      INTEGER(I32), INTENT(IN) :: idupFlag
      INTEGER, INTENT(IN) :: isourcehandle
      INTEGER, INTENT(INOUT) :: idesthandle
      
      IF (IAND(idupFlag,ibitfield) .NE. ibitfield) THEN
        IF (isourcehandle .NE. ST_NOHANDLE) THEN
          CALL storage_copy(isourcehandle,idesthandle)
        END IF
      ELSE
        idesthandle = isourcehandle
      END IF
      
    END SUBROUTINE checkAndCopy
  END SUBROUTINE hadapt_duplicateAdaptation

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE hadapt_restoreAdaptation(rhadaptBackup,rhadapt)

!<description>
    ! This subroutine restares data of an adaptivity structure. All information
    ! arrays which are not shared between rhadaptBackup and another adaptivity
    ! structure is copied (back) to rhadapt.
!</description>

!<input>
    ! Backup of an adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadaptBackup
!</input>

!<inputoutput>
    ! Destination adaptivity structure
    ! All components where a duplicates exist in rhadaptBackup are copied
    ! to rhadapt, overwriting the old information arrays.
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(I32) :: idupFlag
  
    idupFlag = rhadapt%iduplicationFlag

    rhadapt%iSpec                = rhadaptBackup%iSpec
    rhadapt%iadaptationStrategy  = rhadaptBackup%iadaptationStrategy
    rhadapt%NSUBDIVIDEMAX        = rhadaptBackup%NSUBDIVIDEMAX
    rhadapt%nRefinementSteps     = rhadaptBackup%nRefinementSteps
    rhadapt%nCoarseningSteps     = rhadaptBackup%nCoarseningSteps
    rhadapt%nSmoothingSteps      = rhadaptBackup%nSmoothingSteps
    rhadapt%drefinementTolerance = rhadaptBackup%dRefinementTolerance
    rhadapt%dCoarseningTolerance = rhadaptBackup%dCoarseningTolerance
    rhadapt%ndim                 = rhadaptBackup%ndim
    rhadapt%NVT0                 = rhadaptBackup%NVT0
    rhadapt%NVT                  = rhadaptBackup%NVT
    rhadapt%increaseNVT          = rhadaptBackup%increaseNVT
    rhadapt%NVBD0                = rhadaptBackup%NVBD0
    rhadapt%NVBD                 = rhadaptBackup%NVBD
    rhadapt%NBCT                 = rhadaptBackup%NBCT
    rhadapt%NEL0                 = rhadaptBackup%NEL0
    rhadapt%NEL                  = rhadaptBackup%NEL
    rhadapt%NELMAX               = rhadaptBackup%NELMAX
    rhadapt%nGreenElements       = rhadaptBackup%nGreenElements
    rhadapt%InelOfType0          = rhadaptBackup%InelOfType0
    rhadapt%InelOfType           = rhadaptBackup%InelOfType

    ! Call checkAndCopy for all components. This will either copy the handle
    ! or allocate new memory and copy the content of the component.
    
    ! Bit   0: Imarker
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IMARKER,&
        rhadapt%h_Imarker,&
        rhadaptBackup%h_Imarker)
   
    ! Bit   1: IvertexAge
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IVERTEXAGE,&
        rhadapt%h_IvertexAge,&
        rhadaptBackup%h_IvertexAge)
    CALL storage_getbase_int(rhadapt%h_IvertexAge,&
        rhadapt%p_IvertexAge)

    ! Bit   2: InodalProperty
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_INODALPROPERTY,&
        rhadapt%h_InodalProperty,&
        rhadaptBackup%h_InodalProperty)
    CALL storage_getbase_int(rhadapt%h_InodalProperty,&
        rhadapt%p_InodalProperty)

    ! Bit   3: IverticesAtElement
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IVERTICESATELEMENT,&
        rhadapt%h_IverticesAtElement,&
        rhadaptBackup%h_IverticesAtElement)
    CALL storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
        rhadapt%p_IverticesAtElement)

    ! Bit   4: IneighboursAtElement
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_INEIGHATELEMENT,&
        rhadapt%h_IneighboursAtElement,&
        rhadaptBackup%h_IneighboursAtElement)
    CALL storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
        rhadapt%p_IneighboursAtElement)

    ! Bit   5: ImidneighboursAtElement
    CALL checkAndCopy(idupFlag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
        rhadapt%h_ImidneighboursAtElement,&
        rhadaptBackup%h_ImidneighboursAtElement)
    CALL storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
        rhadapt%p_ImidneighboursAtElement)

    ! Bit   6: rVertexCoordinates
    IF (IAND(idupFlag, HADAPT_SHARE_RVERTEXCOORDINATES) .NE.&
        HADAPT_SHARE_RVERTEXCOORDINATES) THEN
      
      SELECT CASE(rhadaptBackup%ndim)
      CASE(NDIM2D)
        CALL qtree_restoreQuadtree(rhadaptBackup%rVertexCoordinates2D,&
            rhadapt%rVertexCoordinates2D)
      CASE(NDIM3D)
        CALL otree_restoreOctree(rhadaptBackup%rVertexCoordinates3D,&
            rhadapt%rVertexCoordinates3D)
      CASE DEFAULT
        CALL output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_restoreAdaptation')
        CALL sys_halt()
      END SELECT
    END IF

    ! Bit   7: rBoundary

    ! Bit   8: rElementsAtVertex
    IF (IAND(idupFlag, HADAPT_SHARE_RELEMENTSATVERTEX) .NE.&
        HADAPT_SHARE_RELEMENTSATVERTEX) THEN
      CALL arrlst_restoreArrayList(rhadaptBackup%rElementsAtVertex,&
          rhadapt%rElementsAtVertex)
    END IF
  
  CONTAINS
    
    SUBROUTINE checkAndCopy (idupFlag,ibitfield,idesthandle,isourcehandle)
      
      ! Checks if idupFlag has all bits ibitfield set.
      ! If not, the memory behind isourcehandle is copied to idesthandle
      ! overwriting all previous information.
      
      INTEGER(I32), INTENT(IN) :: ibitfield
      INTEGER(I32), INTENT(IN) :: idupFlag
      INTEGER, INTENT(IN) :: isourcehandle
      INTEGER, INTENT(INOUT) :: idesthandle
      
      IF (IAND(idupFlag,ibitfield) .NE. ibitfield) THEN
        IF (isourcehandle .NE. ST_NOHANDLE) THEN
          CALL storage_copy(isourcehandle,idesthandle)
        END IF
      END IF
      
    END SUBROUTINE checkAndCopy
  END SUBROUTINE hadapt_restoreAdaptation

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_setVertexCoords2D(rhadapt,h_DvertexCoords,nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the two-dimensional array 
    ! p_DvertexCoords and stores them in the quadtree.
!</description>

!<input>
    ! Handle to the vertex coordinates
    INTEGER, INTENT(IN)                 :: h_DvertexCoords

    ! Number of vertices
    INTEGER(PREC_VERTEXIDX), INTENT(IN) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)    :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_QTREEIDX)            :: nnode
    REAL(DP)                          :: xmin,xmax,ymin,ymax

    ! Check if handle is not empty
    IF (h_DvertexCoords .EQ. ST_NOHANDLE) THEN
      CALL output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords2D')
      CALL sys_halt()
    END IF

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).EQ.HADAPT_HAS_COORDS) THEN
      CALL qtree_releaseQuadtree(rhadapt%rVertexCoordinates2D)
    END IF
    
    ! Set pointer
    CALL storage_getbase_double2D(h_DvertexCoords,p_DvertexCoords,nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin=MINVAL(p_DvertexCoords(1,:))
    xmax=MAXVAL(p_DvertexCoords(1,:))
    ymin=MINVAL(p_DvertexCoords(2,:))
    ymax=MAXVAL(p_DvertexCoords(2,:))
    
    ! Estimate number of initial quadrilaterals
    nnode=INT(0.5_DP*nvt)

    ! Create quadtree for vertices
    CALL qtree_createQuadtree(rhadapt%rVertexCoordinates2D,nvt,&
        nnode,xmin,ymin,xmax,ymax)
    
    ! Copy vertex coordinates to quadtree
    CALL qtree_copyToQuadtree(p_DvertexCoords,rhadapt%rVertexCoordinates2D)

    ! Set specifier for quadtree
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim=NDIM2D
    rhadapt%NVT =nvt
  END SUBROUTINE hadapt_setVertexCoords2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getVertexCoords2D(rhadapt,h_DvertexCoords,nvt,ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the quadtree
    ! and stores them in the two-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will 
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    INTEGER, INTENT(INOUT)        :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    INTEGER(PREC_VERTEXIDX), INTENT(OUT), OPTIONAL :: nvt

    ! OPTIONAL: number of spatial dimensions
    INTEGER, INTENT(OUT), OPTIONAL :: ndim
!</output>
!</subroutine>

    ! Check if coordinates exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).NE.HADAPT_HAS_COORDS) THEN
      CALL output_line('Quadtree does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords2D')
      CALL sys_halt()
    END IF

    ! Check if coordinates are given in 2D
    IF (rhadapt%ndim .NE. NDIM2D) THEN
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords2D')
      CALL sys_halt()
    END IF

    ! Copy quadtree to handle h_DvertexCoords
    CALL qtree_copyFromQuadtree(rhadapt%rVertexCoordinates2D,h_DvertexCoords)

    ! Set dimension
    IF (PRESENT(ndim)) ndim=rhadapt%ndim
    IF (PRESENT(nvt))  nvt=rhadapt%NVT
  END SUBROUTINE hadapt_getVertexCoords2D
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_setVertexCoords3D(rhadapt,h_DvertexCoords,nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the three-dimensional array 
    ! p_DvertexCoords and stores them in the octree.
!</description>

!<input>
    ! Handle to the vertex coordinates
    INTEGER, INTENT(IN)                 :: h_DvertexCoords

    ! Number of vertices
    INTEGER(PREC_VERTEXIDX), INTENT(IN) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)    :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_QTREEIDX)            :: nnode
    REAL(DP)                          :: xmin,xmax,ymin,ymax,zmin,zmax

    ! Check if handle is not empty
    IF (h_DvertexCoords .EQ. ST_NOHANDLE) THEN
      CALL output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords3D')
      CALL sys_halt()
    END IF

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).EQ.HADAPT_HAS_COORDS) THEN
      CALL otree_releaseOctree(rhadapt%rVertexCoordinates3D)
    END IF
    
    ! Set pointer
    CALL storage_getbase_double2D(h_DvertexCoords,p_DvertexCoords,nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin=MINVAL(p_DvertexCoords(1,:))
    xmax=MAXVAL(p_DvertexCoords(1,:))
    ymin=MINVAL(p_DvertexCoords(2,:))
    ymax=MAXVAL(p_DvertexCoords(2,:))
    zmin=MINVAL(p_DvertexCoords(3,:))
    zmax=MAXVAL(p_DvertexCoords(3,:))
    
    ! Estimate number of initial quadrilaterals
    nnode=INT(0.5_DP*nvt)

    ! Create octree for vertices
    CALL otree_createOctree(rhadapt%rVertexCoordinates3D,nvt,&
        nnode,xmin,ymin,zmin,xmax,ymax,zmax)
    
    ! Copy vertex coordinates to octree
    CALL otree_copyToOctree(p_DvertexCoords,rhadapt%rVertexCoordinates3D)

    ! Set specifier for quadtree
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim=NDIM3D
    rhadapt%NVT =nvt
  END SUBROUTINE hadapt_setVertexCoords3D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getVertexCoords3D(rhadapt,h_DvertexCoords,nvt,ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the octree
    ! and stores them in the three-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will 
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    INTEGER, INTENT(INOUT)     :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    INTEGER(PREC_VERTEXIDX), INTENT(OUT), OPTIONAL :: nvt

    ! OPTIONAL: number of spatial dimensions
    INTEGER, INTENT(OUT), OPTIONAL :: ndim
!</output>
!</subroutine>

    ! Check if coordinates exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).NE.HADAPT_HAS_COORDS) THEN
      CALL output_line('Octree does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords3D')
      CALL sys_halt()
    END IF

    ! Check if coordinates are given in 2D
    IF (rhadapt%ndim .NE. NDIM3D) THEN
      CALL output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords3D')
      CALL sys_halt()
    END IF

    ! Copy octree to handle h_DvertexCoords
    CALL otree_copyFromOctree(rhadapt%rVertexCoordinates3D,h_DvertexCoords)

    ! Set dimension
    IF (PRESENT(ndim)) ndim=rhadapt%ndim
    IF (PRESENT(nvt))  nvt=rhadapt%NVT
  END SUBROUTINE hadapt_getVertexCoords3D

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE hadapt_setVerticesAtElement(rhadapt,h_IverticesAtElement,nel)

!<description>
    ! This routine assigns the handle to the "vertices-at-element" array
    ! to the adaptivity structure and sets up internal links.
!</description>

!<input>
    ! Handle to p_IverticesAtElement
    INTEGER, INTENT(IN) :: h_IverticesAtElement

    ! Total number of elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: nel
!</input>

!<inputoutput>
    ! adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    IF (h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVerticesAtElement')
      CALL sys_halt()
    END IF

    rhadapt%h_IverticesAtElement=h_IverticesAtElement
    CALL storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
        rhadapt%p_IverticesAtElement)
    
    ! Set specifier for IverticesAtElement
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_VERTATELEM)

    ! Set dimensions
    rhadapt%NEL =nel
  END SUBROUTINE hadapt_setVerticesAtElement

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getVerticesAtElement(rhadapt,h_IverticesAtElement,nel)

!<description>
    ! This routine assigns the "vertices-at-element" array from the
    ! adaptivity structure to the given handle.
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_IverticesAtElement
    INTEGER, INTENT(INOUT) :: h_IverticesAtElement
!</inputoutput>

!<output>
    ! OPTIONAL: number of elements
    INTEGER(PREC_ELEMENTIDX), INTENT(OUT), OPTIONAL :: nel
!</output>
!</subroutine>

    ! Check if "vertices-at-element" array exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_VERTATELEM).NE.HADAPT_HAS_VERTATELEM) THEN
      CALL output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVerticesAtElement')
      CALL sys_halt()
    END IF

    ! Check if handle needs to be freed first
    IF (h_IverticesAtElement /= ST_NOHANDLE .AND.&
        h_IverticesAtElement /= rhadapt%h_IverticesAtElement)&
        CALL storage_free(h_IverticesAtElement)

    ! Assign handle
    h_IverticesAtElement=rhadapt%h_IverticesAtElement

    ! Set dimensions
    IF(PRESENT(nel))   nel=rhadapt%NEL
  END SUBROUTINE hadapt_getVerticesAtElement
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE hadapt_setNeighboursAtElement(rhadapt,h_IneighboursAtElement)

!<description>
    ! This routine assigns the handle to the "neighbours-at-element" array
    ! to the adaptivity structure and sets up internal links.
!</description>

!<input>
    ! Handle to p_IneighboursAtElement
    INTEGER, INTENT(IN) :: h_IneighboursAtElement
!</input>

!<inputoutput>
    ! adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    IF (h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setNeighboursAtElement')
      CALL sys_halt()
    END IF

    rhadapt%h_IneighboursAtElement=h_IneighboursAtElement
    CALL storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
        rhadapt%p_IneighboursAtElement)
    
    ! Create structure of "mid-adjacent" neighbours
    CALL storage_copy(rhadapt%h_IneighboursAtElement,&
        rhadapt%h_ImidneighboursAtElement)
    CALL storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
        rhadapt%p_ImidneighboursAtElement)
    
    ! Set specifier for IverticesAtElement
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_NEIGHATELEM)
  END SUBROUTINE hadapt_setNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getNeighboursAtElement(rhadapt,&
      h_IneighboursAtElement,nel)

!<description>
    ! This routine assigns the "neighbours-at-element" array from the
    ! adaptivity structure to the given handle
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_IneighboursAtElement
    INTEGER, INTENT(INOUT) :: h_IneighboursAtElement
!</inputoutput>

!<output>
    ! OPTIONAL: number of elements
    INTEGER(PREC_ELEMENTIDX), INTENT(OUT), OPTIONAL :: nel
!</output>
!</subroutine>

    ! Check if "neighbours-at-element" array exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_NEIGHATELEM).NE.HADAPT_HAS_NEIGHATELEM) THEN
      CALL output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNeighboursAtElement')
      CALL sys_halt()
    END IF

    ! Check if handle needs to be freed first
    IF (h_IneighboursAtElement /= ST_NOHANDLE .AND.&
        h_IneighboursAtElement /= rhadapt%h_IneighboursAtElement)&
        CALL storage_free(h_IneighboursAtElement)

    ! Assign handle
    h_IneighboursAtElement=rhadapt%h_IneighboursAtElement

    ! Set dimensions
    IF(PRESENT(nel))   nel=rhadapt%NEL
  END SUBROUTINE hadapt_getNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_setNelOfType(rhadapt,InelOfType)

!<description>
    ! This subroutine sets the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</description>

!<input>
    ! Number of elements with a defined number of vertices per element.
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE), INTENT(IN) :: InelOfType
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    rhadapt%InelOfType=InelOfType

    ! Set specifier for InelOfType
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_NELOFTYPE)
  END SUBROUTINE hadapt_setNelOfType

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getNelOfType(rhadapt,InelOfType)

!<description>
    ! This subroutine returns the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>

!<output>
    ! Number of elements with a defined number of vertices per element.
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE), INTENT(OUT) :: InelOfType
!</output>
!</subroutine>

    ! Check if "InelOfType" array exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_NELOFTYPE).NE.HADAPT_HAS_NELOFTYPE) THEN
      CALL output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNelOfType')
      CALL sys_halt()
    END IF
    
    InelOfType=rhadapt%InelOfType
  END SUBROUTINE hadapt_getNelOfType

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_setBoundary(rhadapt,h_IboundaryCpIdx,&
      h_IverticesAtBoundary,h_DvertexParameterValue,nbct,nvbd)

!<description>
    ! This subroutine sets the boundary structure to the
    ! adaptivity structure and initializes internal data
!</description>

!<input>
    ! Handle to p_IboundaryCpIdx
    INTEGER, INTENT(IN) :: h_IboundaryCpIdx

    ! Handle to p_IverticesAtBoundary
    INTEGER, INTENT(IN) :: h_IverticesAtBoundary

    ! Handle to p_DvertexParameterValue
    INTEGER, INTENT(IN) :: h_DvertexParameterValue

    ! Number of boundary components
    INTEGER, INTENT(IN) :: nbct

    ! Number of vertices at the boundary 
    INTEGER, INTENT(IN) :: nvbd
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexParameterValue2D
    REAL(DP), DIMENSION(:), POINTER   :: p_DvertexParameterValue
    INTEGER, DIMENSION(:), POINTER    :: p_IboundaryCpIdx
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER    :: p_IverticesAtBoundary
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtBoundary
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER :: ioff,lvbd,ivbdStart,ivbdEnd,ibct,h_IneighboursAtBoundary
    
    ! Check if handle are not empty
    IF ((h_IboundaryCpIdx .EQ. ST_NOHANDLE) .OR.&
        (h_IverticesAtBoundary .EQ. ST_NOHANDLE) .OR.&
        (h_DvertexParameterValue .EQ. ST_NOHANDLE)) THEN
      CALL output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setBoundary')
      CALL sys_halt()
    END IF

    ! Check if boundary structure exists and remove it
    IF (ASSOCIATED(rhadapt%rBoundary)) THEN
      DO ibct=1,SIZE(rhadapt%rBoundary,1)
        CALL btree_releaseTree(rhadapt%rBoundary(ibct))
      END DO
      DEALLOCATE(rhadapt%rBoundary)
    END IF

    ! Set pointers
    CALL storage_getbase_double(h_DvertexParameterValue,p_DvertexParameterValue,nvbd)
    CALL convert_pointer(nvbd,p_DvertexParameterValue,p_DvertexParameterValue2D)
    CALL storage_getbase_int(h_IboundaryCpIdx,p_IboundaryCpIdx,nbct+1)
    CALL storage_getbase_int(h_IverticesAtBoundary,p_IverticesAtBoundary,nvbd)
    
    ! Allocate array of search trees
    ALLOCATE(rhadapt%rBoundary(nbct))
    
    ! Create auxiliary array
    Isize=(/2,nvbd/)
    CALL storage_new('hadapt_setBoundary','IData',Isize,ST_INT,&
        h_IneighboursAtBoundary,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int2D(h_IneighboursAtBoundary,p_IneighboursAtBoundary)

    ! Initialization
    ioff=0

    DO ibct=1,nbct

      ! Create a separate search tree for each boundary component
      lvbd=p_IboundaryCpIdx(ibct+1)-p_IboundaryCpIdx(ibct)
      CALL btree_createTree(rhadapt%rBoundary(ibct),lvbd,ST_INT,1,0,2)

      ! Set subdimensions
      ivbdStart=ioff+1
      ivbdEnd  =ioff+lvbd
      ioff     =ioff+lvbd

      ! Fill auxiliary working array p_IneighboursAtBoundary
      ! For each item ivbd we have:
      !   p_IneighboursAtBoundary(BdrPrev,ivbd) -> previous neighbour of ivbd
      !   p_IneighboursAtBoundary(BdrNext,ivbd) -> following neighbour of ivbd
      p_IneighboursAtBoundary(BdrPrev,ivbdStart)          =p_IverticesAtBoundary(ivbdEnd)
      p_IneighboursAtBoundary(BdrPrev,ivbdStart+1:ivbdEnd)=p_IverticesAtBoundary(ivbdStart:ivbdEnd-1)
      p_IneighboursAtBoundary(BdrNext,ivbdStart:ivbdEnd-1)=p_IverticesAtBoundary(ivbdStart+1:ivbdEnd)
      p_IneighboursAtBoundary(BdrPrev,ivbdEnd)            =p_IverticesAtBoundary(ivbdStart)

      ! Fill search tree
      CALL btree_copyToTree(p_IverticesAtBoundary(ivbdStart:ivbdEnd),rhadapt%rBoundary(ibct),&
          p_DData=p_DvertexParameterValue2D(:,ivbdStart:ivbdEnd),&
          p_IData=p_IneighboursAtBoundary(:,ivbdStart:ivbdEnd))
    END DO
    
    ! Set dimensions
    rhadapt%NBCT =nbct
    rhadapt%NVBD =nvbd
    rhadapt%NVBD0=nvbd

    ! Free auxiliary storage
    CALL storage_free(h_IneighboursAtBoundary)

    ! Set specifier for boundary
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_BOUNDARY)

  CONTAINS
    
    ! ****************************************
    ! The convertion routine 1D -> 2D
    
    SUBROUTINE convert_pointer(isize,ptr_1d,ptr_2d)
      INTEGER(I32), INTENT(IN)                 :: isize
      REAL(DP), DIMENSION(1,isize), INTENT(IN), TARGET :: ptr_1d
      REAL(DP), DIMENSION(:,:), POINTER                :: ptr_2d
      
      ptr_2d => ptr_1d
    END SUBROUTINE convert_pointer
  END SUBROUTINE hadapt_setBoundary
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getBoundary(rhadapt,h_IboundaryCpIdx,&
      h_IverticesAtBoundary,h_DvertexParameterValue,nvbd,nbct)

!<description>
    ! This subroutine extracts the boundary data from the adaptivity structure
    ! and generates the the arrays p_IboundaryCpIdx, p_IverticesAtBoundary and 
    ! p_DvertexParameterValue. If the optional parameter nvbd is given
    ! then the number of boundary vertices is assigned.
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Handle to p_IboundaryCpIdx
    INTEGER, INTENT(INOUT) :: h_IboundaryCpIdx
    
    ! Handle to p_IverticesAtBoundary
    INTEGER, INTENT(INOUT) :: h_IverticesAtBoundary

    ! Handle to p_DvertexParameterValue
    INTEGER, INTENT(INOUT) :: h_DvertexParameterValue

    ! OPTIONAL: number of vertices at the boundary
    INTEGER, INTENT(OUT), OPTIONAL :: nvbd

    ! OPTIONAL: number of boundary components
    INTEGER, INTENT(OUT), OPTIONAL :: nbct
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER, DIMENSION(:), POINTER  :: p_IboundaryCpIdx
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER   :: p_IverticesAtBoundary
    INTEGER(I32) :: isize
    INTEGER :: ioff,lvbd,ivbdStart,ivbdEnd,ibct

    ! Check if boundary data exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_BOUNDARY).NE.HADAPT_HAS_BOUNDARY) THEN
      CALL output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getBoundary')
      CALL sys_halt()
    END IF

    ! Due to the fact, that the boundary data may be stored in multiple 
    ! binary search trees (one for each component) the static arrays 
    ! are pre-allocated and filled step-by-step from the dynamic data
    ! in the boundary search tree(s).

    ! Check if the arrays p_IboundaryCpIdx, p_IverticesAtBoundary and 
    ! p_DvertexParameterValue have correct dimension
    IF (h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL storage_new('hadapt_getBoundary','p_IboundaryCpIdx',&
          rhadapt%NBCT+1,ST_INT,h_IboundaryCpIdx,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_IboundaryCpIdx,isize)
      IF (isize < rhadapt%NBCT+1) THEN
        CALL storage_realloc('hadapt_getBoundary',rhadapt%NBCT+1,&
            h_IboundaryCpIdx,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF

    IF (h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL storage_new('hadapt_getBoundary','p_IverticesAtBoundary',&
          rhadapt%NVBD,ST_INT,h_IverticesAtBoundary,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_IverticesAtBoundary,isize)
      IF (isize < rhadapt%NVBD) THEN
        CALL storage_realloc('hadapt_getBoundary',rhadapt%NVBD,&
            h_IverticesAtBoundary,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF

    IF (h_DvertexParameterValue .EQ. ST_NOHANDLE) THEN
      CALL storage_new('hadapt_getBoundary','p_DvertexParameterValue',&
          rhadapt%NVBD,ST_DOUBLE,h_DvertexParameterValue,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_DvertexParameterValue,isize)
      IF (isize < rhadapt%NVBD) THEN
        CALL storage_realloc('hadapt_getBoundary',rhadapt%NVBD,&
            h_DvertexParameterValue,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF

    ! Set pointers
    CALL storage_getbase_int(h_IboundaryCpIdx,p_IboundaryCpIdx,rhadapt%NBCT+1)
    CALL storage_getbase_int(h_IverticesAtBoundary,p_IverticesAtBoundary,rhadapt%NVBD)
    CALL storage_getbase_double(h_DvertexParameterValue,p_DvertexParameterValue,rhadapt%NVBD)

    ! Initialization
    p_IboundaryCpIdx(1)=1
    ioff=0

    DO ibct=1,rhadapt%NBCT
      
      ! Set subdimensions
      lvbd     =rhadapt%rBoundary(ibct)%NA
      ivbdStart=ioff+1
      ivbdEnd  =ioff+lvbd
      ioff     =ioff+lvbd

      ! Set index for next boundary component
      p_IboundaryCpIdx(ibct+1)=ioff+1

      ! Get static boundary data from tree
      CALL btree_copyFromTreeKey(rhadapt%rBoundary(ibct),p_IverticesAtBoundary(ivbdStart:ivbdEnd))
      CALL btree_copyFromTreeDble(rhadapt%rBoundary(ibct),p_DvertexParameterValue(ivbdStart:ivbdEnd),1)

      ! Sort array w.r.t. increasing parameter values
      CALL sort_dp(p_DvertexParameterValue(ivbdStart:ivbdEnd),SORT_QUICK,&
          p_IverticesAtBoundary(ivbdStart:ivbdEnd))
    END DO

    ! Set dimension
    IF(PRESENT(nvbd)) nvbd=rhadapt%NVBD
    IF(PRESENT(nbct)) nbct=rhadapt%NBCT
  END SUBROUTINE hadapt_getBoundary

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_setNodalProperty(rhadapt,h_InodalProperty)

!<description>
    ! This subroutine sets the nodal property list to the adaptivity structure
!</description>

!<input>
    ! Handle to p_InodalProperty
    INTEGER, INTENT(IN) :: h_InodalProperty
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    IF (h_InodalProperty .EQ. ST_NOHANDLE) THEN
      CALL output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setNodalProperty')
      CALL sys_halt()
    END IF

    rhadapt%h_InodalProperty=h_InodalProperty
    CALL storage_getbase_int(rhadapt%h_InodalProperty,rhadapt%p_InodalProperty)
    
    ! Set specifier for InodalProperty
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_NODALPROP)

  END SUBROUTINE hadapt_setNodalProperty

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getNodalProperty(rhadapt,h_InodalProperty)

!<description>
    ! This subroutine assignes the "nodal property" array from the
    ! adaptivity structure to the given handle
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_InodalProperty
    INTEGER, INTENT(INOUT) :: h_InodalProperty
!</inputoutput>
!</subroutine>

    ! Check if "nodal property" list exits
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_NODALPROP).NE.HADAPT_HAS_NODALPROP) THEN
      CALL output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNodalProperty')
      CALL sys_halt()
    END IF

    ! Check if handle needs to be freed first
    IF (h_InodalProperty /= ST_NOHANDLE .AND.&
        h_InodalProperty /= rhadapt%h_InodalProperty)&
        CALL storage_free(h_InodalProperty)

    ! Assign handle
    h_InodalProperty=rhadapt%h_InodalProperty
  END SUBROUTINE hadapt_getNodalProperty

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_genElementsAtVertex(rhadapt)

!<description>
    ! This subroutine generates a list of arrays for the "elements-at-vertex"
    ! structure. This is an internal structure of the adaptivity structure
    ! which is used, e.g., in the removel of vertices.
!</description>

!<inputoutput>
    ! Adaptivity structur
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: ive

    ! Check if "vertices-at-element" list exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_VERTATELEM).NE.HADAPT_HAS_VERTATELEM) THEN
      CALL output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_genElementsAtVertex')
      CALL sys_halt()
    END IF

    ! Create arraylist
    CALL arrlst_createArraylist(rhadapt%relementsAtVertex,&
        2*rhadapt%NVT,8*rhadapt%NEL,ST_INT,ARRAYLIST_UNORDERED)

    ! Fill arraylist
    DO iel=1,rhadapt%NEL
      DO ive=1,hadapt_getNVE(rhadapt,iel)
        CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,&
            rhadapt%p_IverticesAtElement(ive,iel),iel,ipos)
      END DO
    END DO
    
    ! Set specifier for relementsAtVertex
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_ELEMATVERTEX)
  END SUBROUTINE hadapt_genElementsAtVertex

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_performAdaptation(rhadapt, rindicator, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine performs the complete adaptation process.
    ! First, the internal data structures are generated and
    ! elements are marked for refinement and coarsening based on
    ! the indicator vector and the tolerances for refinement and
    ! coarsening, respecitvely.
    ! Next, the list of marked elements is updated so as to obtain
    ! a globally conforming triangulation.
    ! Finally, mesh refinement and coarsening is performed
!</description>

!<input>
    ! Indicator vector for refinement
    TYPE(t_vectorScalar), INTENT(IN)            :: rindicator

    ! callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL                                    :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(1) :: Ielements
    INTEGER(PREC_VERTEXIDX), DIMENSION(1) :: Ivertices
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_VERTEXIDX)    :: nvt
    INTEGER(PREC_ELEMENTIDX)   :: nel

    ! Check if dynamic data structures are available
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA) THEN
      CALL output_line('Dynamic data structures are not generated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
      CALL sys_halt()
    END IF

    ! Initialize initial dimensions
    CALL storage_getsize(rhadapt%h_IverticesAtElement, Isize)
    rhadapt%NELMAX = Isize(2)
    CALL storage_getsize(rhadapt%h_IneighboursAtElement, Isize)
    rhadapt%NELMAX = MIN(rhadapt%NELMAX,Isize(2))

    rhadapt%InelOfType0 = rhadapt%InelOfType
    rhadapt%NVT0        = rhadapt%NVT
    rhadapt%NEL0        = rhadapt%NEL
    rhadapt%NVBD0       = rhadapt%NVBD
    rhadapt%increaseNVT = 0
    
    ! What kind of grid refinement should be performed
    SELECT CASE(rhadapt%iadaptationStrategy)
    CASE (HADAPT_NOADAPTATION)   ! No grid refinement

      
    CASE (HADAPT_REDGREEN)       ! Red-green grid refinement

      ! Which spatial dimensions are we?
      SELECT CASE(rhadapt%ndim)
        
      CASE (NDIM2D)
        ! Mark elements for refinement based on indicator function
        CALL mark_refinement2D(rhadapt, rindicator)
        
        ! Mark additional elements to restore conformity
        CALL redgreen_mark_refinement2D(rhadapt, rcollection, fcb_hadaptCallback)
        
        ! Mark element for recoarsening based on indicator function
        CALL redgreen_mark_coarsening2D(rhadapt, rindicator)
        
      CASE DEFAULT
        CALL output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
        CALL sys_halt()
      END SELECT

      ! Compute new dimensions
      nvt = rhadapt%NVT+rhadapt%increaseNVT
      nel = NumberOfElements(rhadapt)

      ! Adjust vertex age array and nodal property array
      CALL storage_realloc('hadapt_performAdaptation', nvt,&
          rhadapt%h_IvertexAge, ST_NEWBLOCK_NOINIT, .TRUE.)
      CALL storage_realloc('hadapt_performAdaptation', nvt,&
          rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .TRUE.)
     
      ! Adjust elemental arrays
      CALL storage_realloc('hadapt_performAdaptation', nel,&
          rhadapt%h_IverticesAtElement, ST_NEWBLOCK_NOINIT, .TRUE.)
      CALL storage_realloc('hadapt_performAdaptation', nel,&
          rhadapt%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT, .TRUE.)
      CALL storage_realloc('hadapt_performAdaptation', nel,&
          rhadapt%h_ImidneighboursAtElement, ST_NEWBLOCK_NOINIT, .TRUE.)

      ! Reset pointers
      CALL storage_getbase_int(rhadapt%h_IvertexAge,&
                               rhadapt%p_IvertexAge)
      CALL storage_getbase_int(rhadapt%h_InodalProperty,&
                               rhadapt%p_InodalProperty)
      CALL storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                                 rhadapt%p_IverticesAtElement)
      CALL storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                                 rhadapt%p_IneighboursAtElement)
      CALL storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                                 rhadapt%p_ImidneighboursAtElement)

      ! Adjust dimension of solution vector
      IF (PRESENT(fcb_hadaptCallback) .AND. PRESENT(rcollection)) THEN
        Ivertices = (/nvt/)
        Ielements = (/0/)
        CALL fcb_hadaptCallback(rcollection, HADAPT_OPR_ADJUSTVERTEXDIM,&
                                Ivertices, Ielements)
      END IF

      ! Perform refinement
      CALL redgreen_refine(rhadapt, rcollection, fcb_hadaptCallback)

      ! Perform coarsening
      CALL redgreen_coarsen(rhadapt, rcollection, fcb_hadaptCallback)

      ! Adjust nodal property array
      CALL storage_realloc('hadapt_performAdaptation', rhadapt%NVT,&
          rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .TRUE.)
            
    CASE DEFAULT
      CALL output_line('Unsupported refinement strategy!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, some auxiliary routines follow
    
    !*************************************************************
    ! This function computes the number of elements after refinement

    FUNCTION NumberOfElements(rhadapt) RESULT(nel)
      
      TYPE(t_hadapt), INTENT(IN) :: rhadapt
      INTEGER(PREC_ELEMENTIDX)      :: nel
      
      ! local variables
      INTEGER, DIMENSION(:), POINTER :: p_Imarker
      INTEGER(PREC_ELEMENTIDX)       :: iel
      
      CALL storage_getbase_int(rhadapt%h_Imarker, p_Imarker)
      
      ! Initialize number of elements by current number
      nel = rhadapt%NEL0
      
      ! Loop over all elements and check marker
      DO iel = 1, rhadapt%NEL0
        
        SELECT CASE(p_Imarker(iel))
        CASE(MARK_REF_TRIA2TRIA_1,MARK_REF_TRIA2TRIA_2,MARK_REF_TRIA2TRIA_3,&
             MARK_REF_QUAD2QUAD_13,MARK_REF_QUAD2QUAD_24,MARK_CRS_2QUAD3TRIA)
          ! Interestingly enought, the coarsening of 2 quadrilaterals into three
          ! triangles which reduces the number of vertices by one also increases
          ! the number of elements by one which has to be taken into account here.
          nel = nel+1

        CASE(MARK_REF_QUAD3TRIA_1,MARK_REF_QUAD3TRIA_2,MARK_REF_QUAD3TRIA_3,&
             MARK_REF_QUAD3TRIA_4,MARK_REF_TRIA3TRIA_12,MARK_REF_TRIA3TRIA_23,&
             MARK_REF_TRIA3TRIA_13)
          nel = nel+2
          
        CASE(MARK_REF_TRIA4TRIA,MARK_REF_QUAD4QUAD,MARK_REF_QUAD4TRIA_12,&
             MARK_REF_QUAD4TRIA_23,MARK_REF_QUAD4TRIA_34,MARK_REF_QUAD4TRIA_14)
          nel = nel+3
        END SELECT
      END DO
    END FUNCTION NumberOfElements
  END SUBROUTINE hadapt_performAdaptation

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE hadapt_infoStatistics(rhadapt)

!<description>
    ! This subroutine outputs statistical info about the adaptivity data structure
!</description>

!<input>
    ! adaptivity data structure
    TYPE(t_hadapt), INTENT(IN) :: rhadapt
!</input>
!</subroutine>

    ! local variables
    INTEGER :: ibct
    
    CALL output_line('Adaptivity statistics:')
    CALL output_line('----------------------')
    CALL output_line('Total number of grid refinement steps:       '//&
        TRIM(sys_siL(rhadapt%nRefinementSteps,3)))
    CALL output_line('Total number of grid coarsening steps:       '//&
        TRIM(sys_siL(rhadapt%nCoarseningSteps,3)))
    CALL output_line('Total number of grid smoothing  steps:       '//&
        TRIM(sys_siL(rhadapt%nSmoothingSteps,3)))
    CALL output_line('Strategy for grid refinement and coarsening: '//&
        TRIM(sys_siL(rhadapt%iadaptationStrategy,3)))
    CALL output_line('Total/initial number of elements:            '//&
        TRIM(sys_siL(rhadapt%NEL,15))//"("//TRIM(sys_siL(rhadapt%NEL0,15))//")")
    CALL output_line('Total/initial number of vertices:            '//&
        TRIM(sys_siL(rhadapt%NVT,15))//"("//TRIM(sys_siL(rhadapt%NVT0,15))//")")
    CALL output_line('Total/initial number of vertices at boundary:'//&
        TRIM(sys_siL(rhadapt%NVBD,15))//"("//TRIM(sys_siL(rhadapt%NVBD0,15))//")")
    CALL output_lbrk

    CALL output_line('Handles:')
    CALL output_line('--------')
    CALL output_line('h_Imarker:                '//&
        TRIM(sys_siL(rhadapt%h_Imarker,15)))
    CALL output_line('h_IvertexAge:             '//&
        TRIM(sys_siL(rhadapt%h_IvertexAge,15)))
    CALL output_line('h_InodalProperty:         '//&
        TRIM(sys_siL(rhadapt%h_InodalProperty,15)))
    CALL output_line('h_IverticesAtElement:     '//&
        TRIM(sys_siL(rhadapt%h_IverticesAtElement,15)))
    CALL output_line('h_IneighboursAtelement:   '//&
        TRIM(sys_siL(rhadapt%h_IneighboursAtElement,15)))
    CALL output_line('h_ImidneighboursAtelement:'//&
        TRIM(sys_siL(rhadapt%h_ImidneighboursAtElement,15)))
    CALL output_lbrk

    CALL output_line('Coordinates:')
    CALL output_line('------------')
    SELECT CASE(rhadapt%ndim)
    CASE(NDIM2D)
      CALL qtree_infoQuadtree(rhadapt%rVertexCoordinates2D)
    CASE(NDIM3D)
      CALL otree_infoOctree(rhadapt%rVertexCoordinates3D)
    END SELECT
    CALL output_lbrk
    
    CALL output_line('Boundary:')
    CALL output_line('---------')
    DO ibct=1,rhadapt%NBCT
      CALL output_line('Boundary component '//TRIM(sys_siL(ibct,3))//":")
      CALL output_line('------------------------')
      CALL btree_infoTree(rhadapt%rBoundary(ibct))
    END DO
    CALL output_lbrk

    CALL output_line('Element at vertex structure:')
    CALL output_line('----------------------------')
    CALL arrlst_infoArraylist(rhadapt%relementsAtVertex)
  END SUBROUTINE hadapt_infoStatistics

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_writeGridSVG(rhadapt,coutputFile,width,height)

!<description>
    ! This subroutine outputs the current state of the adapted grid stored
    ! in the dynamic data structure to a given file in SVG format.
!</description>

!<input>
    ! Output file name w/o suffix .svg
    CHARACTER(LEN=*), INTENT(IN)     :: coutputFile

    ! OPTIONAL: Width of the generated file
    INTEGER, INTENT(IN), OPTIONAL :: width

    ! OPTIONAL: Heigh of the generated file
    INTEGER, INTENT(IN), OPTIONAL :: height
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local parameters
    INTEGER, PARAMETER :: defaultwidth  = 1280
    INTEGER, PARAMETER :: defaultheight = 1024
    INTEGER, PARAMETER :: xoffset       = 100
    INTEGER, PARAMETER :: yoffset       = 100
    INTEGER, PARAMETER :: font_size     = 18
    INTEGER, PARAMETER :: colsep        = 10
    INTEGER, PARAMETER :: linesep       = 32
    CHARACTER(LEN=*), PARAMETER :: font_family = "Arial"
    
    ! local variables
    REAL(DP), DIMENSION(2*NDIM2D) :: bbox
    REAL(DP) :: x0,y0,xdim,ydim,dscale
    
    INTEGER, DIMENSION(:), POINTER :: p_Imarker
    INTEGER(PREC_VERTEXIDX)  :: ivt
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_ARRAYLISTIDX):: ipos
    INTEGER :: xsize,ysize,iunit,ive,nve
    INTEGER, SAVE :: iout=0
    
    ! Check if dynamic data structures generated
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA) THEN
      CALL output_line('Dynamic data structures are not generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridSVG')
      CALL sys_halt()
    END IF
    
    ! Increment the sample number
    iout=iout+1
    
    ! Open output file for writing
    CALL io_openFileForWriting(TRIM(ADJUSTL(coutputFile))//'.'//&
        TRIM(sys_siL(iout,5))//'.svg',iunit,SYS_REPLACE,bformatted=.TRUE.)
    
    ! Write prolog, XML-declaration, document type declaration)
    WRITE(iunit,FMT='(A)') '<?xml version="1.0" encoding="utf-8" standalone="yes"?>'
    WRITE(iunit,FMT='(A)') '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN"'
    WRITE(iunit,FMT='(A)') ' "http://www.w3.org/Graphics/SVG/1.0/DTD/svg11.dtd">'
    WRITE(iunit,FMT='(A)')
    WRITE(iunit,FMT='(A)') '<!-- Created with Featflow2 (http://www.featflow.de/) -->'
    WRITE(iunit,FMT='(A)')
    WRITE(iunit,FMT='(A)') '<svg version="1.0"'
    WRITE(iunit,FMT='(A)') ' xmlns="http://www.w3.org/2000/svg"'
    WRITE(iunit,FMT='(A)') ' xmlns:xlink="http://www.w3.org/1999/xlink"'
    
    IF (PRESENT(width)) THEN
      xsize=width
    ELSE
      xsize=defaultwidth
    END IF

    IF (PRESENT(height)) THEN
      ysize=height
    ELSE
      ysize=defaultheight
    END IF

    ! Determine bounding box
    bbox  =qtree_getBoundingBox(rhadapt%rVertexCoordinates2D)
    xdim  =bbox(3)-bbox(1); x0=bbox(1)
    ydim  =bbox(4)-bbox(2); y0=bbox(2)
    dscale=MIN(xsize/xdim,ysize/ydim)

    ! Set height and width of image
    WRITE(iunit,FMT='(A)') ' width="100%" height="100%" xml:space="preserve"'
    WRITE(iunit,FMT='(A)') ' viewBox="'//&
        TRIM(sys_siL(-xoffset,10))//' '//&
        TRIM(sys_siL(-yoffset,10))//' '//&
        TRIM(sys_siL(2*xoffset+xsize,10))//' '//&
        TRIM(sys_siL(2*yoffset+ysize,10))//'">'

    ! Write embedded java-scripts to file
    WRITE(iunit,FMT='(A)') '<defs>'
    WRITE(iunit,FMT='(A)') '<script type="text/javascript">'
    WRITE(iunit,FMT='(A)') '<![CDATA['

    WRITE(iunit,FMT='(A)') '  function ShowVertexInfo(evt,ivt,age,elems) {'
    WRITE(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    WRITE(iunit,FMT='(A)') '    var x=evt.clientX; var y=evt.clientY;'
    WRITE(iunit,FMT='(A)') '    var m=svgdoc.getScreenCTM();'
    WRITE(iunit,FMT='(A)') '    var p=svgdoc.createSVGPoint();'
    WRITE(iunit,FMT='(A)') '    p.x=evt.clientX; p.y=evt.clientY;'
    WRITE(iunit,FMT='(A)') '    p=p.matrixTransform(m.inverse());'
    WRITE(iunit,FMT='(A)') '    p.x=p.x+15; p.y=p.y+15;'

    WRITE(iunit,FMT='(A)') '    var vinfo=svgdoc.getElementById("vinfo");'
    WRITE(iunit,FMT='(A)') '    var vinfo_box=svgdoc.getElementById("vinfo_box");'
    WRITE(iunit,FMT='(A)') '    var vinfo_text1=svgdoc.getElementById("vinfo_text1");'
    WRITE(iunit,FMT='(A)') '    var vinfo_text2=svgdoc.getElementById("vinfo_text2");'
    WRITE(iunit,FMT='(A)') '    var vinfo_text3=svgdoc.getElementById("vinfo_text3");'
    WRITE(iunit,FMT='(A)') '    var vinfo_text4=svgdoc.getElementById("vinfo_text4");'
    WRITE(iunit,FMT='(A)') '    var vinfo_text5=svgdoc.getElementById("vinfo_text5");'
    WRITE(iunit,FMT='(A)') '    var vinfo_text6=svgdoc.getElementById("vinfo_text6");'
    
    WRITE(iunit,FMT='(A)') '    vinfo_box.setAttribute("x",p.x);'
    WRITE(iunit,FMT='(A)') '    vinfo_box.setAttribute("y",p.y);'
    WRITE(iunit,FMT='(A)') '    vinfo_text1.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text1.setAttribute("y",p.y+'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text2.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text2.setAttribute("y",p.y+2*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text3.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text3.setAttribute("y",p.y+3*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text4.setAttribute("x",p.x+vinfo_text1.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    vinfo_text4.setAttribute("y",p.y+'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text4.firstChild.nodeValue = ivt;'
    WRITE(iunit,FMT='(A)') '    vinfo_text5.setAttribute("x",p.x+vinfo_text2.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    vinfo_text5.setAttribute("y",p.y+2*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text5.firstChild.nodeValue = age;'
    WRITE(iunit,FMT='(A)') '    vinfo_text6.setAttribute("x",p.x+vinfo_text3.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    vinfo_text6.setAttribute("y",p.y+3*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    vinfo_text6.firstChild.nodeValue = elems;'

    WRITE(iunit,FMT='(A)') '    var textlen = vinfo_text1.getComputedTextLength()+&
        &vinfo_text4.getComputedTextLength();'
    WRITE(iunit,FMT='(A)') '    textlen = Math.max(textlen,vinfo_text2.getComputedTextLength()&
        &+vinfo_text5.getComputedTextLength());'
    WRITE(iunit,FMT='(A)') '    textlen = Math.max(textlen,vinfo_text3.getComputedTextLength()&
        &+vinfo_text6.getComputedTextLength());'
    WRITE(iunit,FMT='(A)') '    vinfo_box.setAttribute("width",textlen+30);'

    WRITE(iunit,FMT='(A)') '    vinfo.setAttribute("style","visibility: visible");'
    WRITE(iunit,FMT='(A)') '  }'
    
    WRITE(iunit,FMT='(A)') '  function HideVertexInfo() {'
    WRITE(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    WRITE(iunit,FMT='(A)') '    var vinfo=svgdoc.getElementById("vinfo");'
    WRITE(iunit,FMT='(A)') '    vinfo.setAttribute("style","visibility: hidden");'
    WRITE(iunit,FMT='(A)') '  }'

    WRITE(iunit,FMT='(A)') '  function ShowElementInfo(evt,iel,state,marker,kadj,kmidadj,kvert) {'
    WRITE(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    WRITE(iunit,FMT='(A)') '    var x=evt.clientX; var y=evt.clientY;'
    WRITE(iunit,FMT='(A)') '    var m=svgdoc.getScreenCTM();'
    WRITE(iunit,FMT='(A)') '    var p=svgdoc.createSVGPoint();'
    WRITE(iunit,FMT='(A)') '    p.x=evt.clientX; p.y=evt.clientY;'
    WRITE(iunit,FMT='(A)') '    p=p.matrixTransform(m.inverse());'
    WRITE(iunit,FMT='(A)') '    p.x=p.x+15; p.y=p.y+15;'

    WRITE(iunit,FMT='(A)') '    var einfo=svgdoc.getElementById("einfo");'
    WRITE(iunit,FMT='(A)') '    var einfo_box=svgdoc.getElementById("einfo_box");'
    WRITE(iunit,FMT='(A)') '    var einfo_text1=svgdoc.getElementById("einfo_text1");'
    WRITE(iunit,FMT='(A)') '    var einfo_text2=svgdoc.getElementById("einfo_text2");'
    WRITE(iunit,FMT='(A)') '    var einfo_text3=svgdoc.getElementById("einfo_text3");'
    WRITE(iunit,FMT='(A)') '    var einfo_text4=svgdoc.getElementById("einfo_text4");'
    WRITE(iunit,FMT='(A)') '    var einfo_text5=svgdoc.getElementById("einfo_text5");'
    WRITE(iunit,FMT='(A)') '    var einfo_text6=svgdoc.getElementById("einfo_text6");'
    WRITE(iunit,FMT='(A)') '    var einfo_text7=svgdoc.getElementById("einfo_text7");'
    WRITE(iunit,FMT='(A)') '    var einfo_text8=svgdoc.getElementById("einfo_text8");'
    WRITE(iunit,FMT='(A)') '    var einfo_text9=svgdoc.getElementById("einfo_text9");'
    WRITE(iunit,FMT='(A)') '    var einfo_text10=svgdoc.getElementById("einfo_text10");'
    WRITE(iunit,FMT='(A)') '    var einfo_text11=svgdoc.getElementById("einfo_text11");'
    WRITE(iunit,FMT='(A)') '    var einfo_text12=svgdoc.getElementById("einfo_text12");'
    
    WRITE(iunit,FMT='(A)') '    einfo_box.setAttribute("x",p.x);'
    WRITE(iunit,FMT='(A)') '    einfo_box.setAttribute("y",p.y);'
    WRITE(iunit,FMT='(A)') '    einfo_text1.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text1.setAttribute("y",p.y+'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text2.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text2.setAttribute("y",p.y+2*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text3.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text3.setAttribute("y",p.y+3*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text4.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text4.setAttribute("y",p.y+4*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text5.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text5.setAttribute("y",p.y+5*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text6.setAttribute("x",p.x+'//TRIM(sys_siL(colsep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text6.setAttribute("y",p.y+6*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text7.setAttribute("x",p.x+einfo_text1.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    einfo_text7.setAttribute("y",p.y+'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text7.firstChild.nodeValue = iel;'
    WRITE(iunit,FMT='(A)') '    einfo_text8.setAttribute("x",p.x+einfo_text2.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    einfo_text8.setAttribute("y",p.y+2*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text8.firstChild.nodeValue = state;'
    WRITE(iunit,FMT='(A)') '    einfo_text9.setAttribute("x",p.x+einfo_text3.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    einfo_text9.setAttribute("y",p.y+3*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text9.firstChild.nodeValue = marker;'
    WRITE(iunit,FMT='(A)') '    einfo_text10.setAttribute("x",p.x+einfo_text4.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    einfo_text10.setAttribute("y",p.y+4*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text10.firstChild.nodeValue = kadj;'
    WRITE(iunit,FMT='(A)') '    einfo_text11.setAttribute("x",p.x+einfo_text5.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    einfo_text11.setAttribute("y",p.y+5*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text11.firstChild.nodeValue = kmidadj;'
    WRITE(iunit,FMT='(A)') '    einfo_text12.setAttribute("x",p.x+einfo_text6.getComputedTextLength()+20);'
    WRITE(iunit,FMT='(A)') '    einfo_text12.setAttribute("y",p.y+6*'//TRIM(sys_siL(linesep,10))//');'
    WRITE(iunit,FMT='(A)') '    einfo_text12.firstChild.nodeValue = kvert;'

    WRITE(iunit,FMT='(A)') '    var textlen = einfo_text1.getComputedTextLength()+&
        &einfo_text7.getComputedTextLength();'
    WRITE(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text2.getComputedTextLength()+&
        &einfo_text8.getComputedTextLength());'
    WRITE(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text3.getComputedTextLength()+&
        &einfo_text9.getComputedTextLength());'
    WRITE(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text4.getComputedTextLength()+&
        &einfo_text10.getComputedTextLength());'
    WRITE(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text5.getComputedTextLength()+&
        &einfo_text11.getComputedTextLength());'
    WRITE(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text6.getComputedTextLength()+&
        &einfo_text12.getComputedTextLength());'
    WRITE(iunit,FMT='(A)') '    einfo_box.setAttribute("width",textlen+30);'

    WRITE(iunit,FMT='(A)') '    einfo.setAttribute("style","visibility: visible");'
    WRITE(iunit,FMT='(A)') '  }'

    WRITE(iunit,FMT='(A)') '  function HideElementInfo() {'
    WRITE(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    WRITE(iunit,FMT='(A)') '    var einfo=svgdoc.getElementById("einfo");'
    WRITE(iunit,FMT='(A)') '    einfo.setAttribute("style","visibility: hidden");'
    WRITE(iunit,FMT='(A)') '  }'
    WRITE(iunit,FMT='(A)') ']]>'
    WRITE(iunit,FMT='(A)') '</script>'
    WRITE(iunit,FMT='(A)') '</defs>'
       
    !---------------------------------------------------------------------------
    ! Output all elements
    !---------------------------------------------------------------------------
    IF (IAND(rhadapt%iSpec,HADAPT_MARKEDREFINE) .EQ.HADAPT_MARKEDREFINE .OR.&
        IAND(rhadapt%iSpec,HADAPT_MARKEDCOARSEN).EQ.HADAPT_MARKEDCOARSEN) THEN

      ! Set pointer to marker
      CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)

      ! Output elements and color those which are marked for refinement
      DO iel=1,rhadapt%NEL
        
        ! Get number of vertices per elements
        nve=hadapt_getNVE(rhadapt,iel)
        
        IF (iel .LE. rhadapt%NEL0) THEN
          
          ! What kind of element are we?
          SELECT CASE(p_Imarker(iel))
            
          ! Element is neither marked for refinement nor coarsening
          CASE(MARK_ASIS_TRIA,MARK_ASIS_QUAD)
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="white" stroke="black" stroke-width="1"'
            
          ! Element is marked for green refinement
          CASE(MARK_REF_TRIA2TRIA_1,MARK_REF_TRIA2TRIA_2,MARK_REF_TRIA2TRIA_3,&
               MARK_REF_QUAD3TRIA_1,MARK_REF_QUAD3TRIA_2,MARK_REF_QUAD3TRIA_3,&
               MARK_REF_QUAD3TRIA_4,MARK_REF_QUAD4TRIA_12,MARK_REF_QUAD4TRIA_23,&
               MARK_REF_QUAD4TRIA_34,MARK_REF_QUAD4TRIA_14,MARK_REF_QUAD2QUAD_13,&
               MARK_REF_QUAD2QUAD_24)
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="green" stroke="black" stroke-width="1"'

          ! Element is marked for blue refinement
          CASE(MARK_REF_TRIA3TRIA_12,MARK_REF_TRIA3TRIA_23,MARK_REF_TRIA3TRIA_13,&
               MARK_REF_QUADBLUE_412,MARK_REF_QUADBLUE_234,MARK_REF_QUADBLUE_123,&
               MARK_REF_QUADBLUE_341)
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="blue" stroke="black" stroke-width="1"'

          ! Element is marked for red refinement
          CASE(MARK_REF_TRIA4TRIA,MARK_REF_QUAD4QUAD)
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="red" stroke="black" stroke-width="1"'

          ! Element is marked for generic coarsening
          CASE(MARK_CRS_GENERIC)
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="lightgray" stroke="black" stroke-width="1"'
            
          ! Element is marked for specific coarsening
          CASE(MARK_CRS_4QUAD4TRIA:MARK_CRS_4TRIA1TRIA)
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="yellow" stroke="black" stroke-width="1"'

          ! Unknown element marker
          CASE DEFAULT
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="hotpink" stroke="black" stroke-width="1"'
          END SELECT
          
          ! Write data which is common to all elements
          SELECT CASE(nve)
          CASE(TRIA_NVETRI2D)
            WRITE(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                TRIM(sys_siL(iel,10))//''','''//&
                TRIM(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                TRIM(sys_siL(p_Imarker(iel),5))//''','''//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(0,10))//''','''//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(0,10))//''','''//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                TRIM(sys_siL(0,10))//''')"'
            WRITE(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'

          CASE(TRIA_NVEQUAD2D)
            WRITE(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                TRIM(sys_siL(iel,10))//''','''//&
                TRIM(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                TRIM(sys_siL(p_Imarker(iel),5))//''','''//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
            WRITE(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          END SELECT

        ELSE
          
          ! For all new elements there si no marker available
          WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
              '" fill="white" stroke="black" stroke-width="1"'
          
          ! Write data which is common to all elements
          SELECT CASE(nve)
          CASE(TRIA_NVETRI2D)
            WRITE(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                TRIM(sys_siL(iel,10))//''','''//&
                TRIM(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                TRIM(sys_siL(0,5))//''','''//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(0,10))//''','''//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(0,10))//''','''//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                TRIM(sys_siL(0,10))//''')"'
            WRITE(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'

          CASE(TRIA_NVEQUAD2D)
            WRITE(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                TRIM(sys_siL(iel,10))//''','''//&
                TRIM(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                TRIM(sys_siL(0,5))//''','''//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                TRIM(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
            WRITE(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          END SELECT
          
        END IF
               
        ! Each element is a polygon made up from 3/4 points
        WRITE(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        DO ive=1,nve
          xdim=qtree_getX(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(ive,iel))-x0
          ydim=qtree_getY(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(ive,iel))-y0
          WRITE(iunit,FMT='(A)',ADVANCE='NO') &
              TRIM(sys_siL(INT(dscale*xdim),10))//','//&
              TRIM(sys_siL(ysize-INT(dscale*ydim),10))//' '
        END DO
        
        ! And of course, the polygone must be closed !
        xdim=qtree_getX(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(1,iel))-x0
        ydim=qtree_getY(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(1,iel))-y0
        WRITE(iunit,FMT='(A)') &
            TRIM(sys_siL(INT(dscale*xdim),10))//','//&
            TRIM(sys_siL(ysize-INT(dscale*ydim),10))//'"/>'       

      END DO
      
    ELSE
      
      ! Output only elements without individual coloring
      DO iel=1,rhadapt%NEL

        ! Get number of vertices per element
        nve=hadapt_getNVE(rhadapt,iel)
        
        ! For all new elements there si no marker available
        WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
            '" fill="white" stroke="black" stroke-width="1"'

        ! Write data which is common to all elements
        SELECT CASE(nve)
        CASE(TRIA_NVETRI2D)
          WRITE(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
              TRIM(sys_siL(iel,10))//''','''//&
              TRIM(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
              TRIM(sys_siL(0,5))//''','''//&
              TRIM(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
              TRIM(sys_siL(0,10))//''','''//&
              TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
              TRIM(sys_siL(0,10))//''','''//&
              TRIM(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
              TRIM(sys_siL(0,10))//''')"'
          WRITE(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          
        CASE(TRIA_NVEQUAD2D)
          WRITE(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
              TRIM(sys_siL(iel,10))//''','''//&
              TRIM(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
              TRIM(sys_siL(0,5))//''','''//&
              TRIM(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
              TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
              TRIM(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
              TRIM(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
          WRITE(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
        END SELECT
        
        ! Each element is a polygon made up from 3/4 points
        WRITE(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        DO ive=1,nve
          xdim=qtree_getX(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(ive,iel))-x0
          ydim=qtree_getY(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(ive,iel))-y0
          WRITE(iunit,FMT='(A)',ADVANCE='NO') &
              TRIM(sys_siL(INT(dscale*xdim),10))//','//&
              TRIM(sys_siL(ysize-INT(dscale*ydim),10))//' '
        END DO
        
        ! And of course, the polygone must be closed !
        xdim=qtree_getX(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(1,iel))-x0
        ydim=qtree_getY(rhadapt%rVertexCoordinates2D,rhadapt%p_IverticesAtElement(1,iel))-y0
        WRITE(iunit,FMT='(A)') &
            TRIM(sys_siL(INT(dscale*xdim),10))//','//&
            TRIM(sys_siL(ysize-INT(dscale*ydim),10))//'"/>'       
        
      END DO
    END IF
         
    ! Loop over all vertices
    DO ivt=1,qtree_getsize(rhadapt%rVertexCoordinates2D)     
      xdim=qtree_getX(rhadapt%rVertexCoordinates2D,ivt)-x0
      ydim=qtree_getY(rhadapt%rVertexCoordinates2D,ivt)-y0
      
      ! Write vertices as points?
      IF (rhadapt%p_IvertexAge(ivt) > 0) THEN
        WRITE(iunit,FMT='(A)') '<circle id="vt'//TRIM(sys_siL(ivt,10))//'" cx="'//&
            TRIM(sys_siL(INT(dscale*xdim),10))//'" cy="'//&
            TRIM(sys_siL(ysize-INT(dscale*ydim),10))//&
            '" r="0pt" fill="white" stroke="black" stroke-width="1pt"'
      ELSE
        WRITE(iunit,FMT='(A)') '<circle id="vt'//TRIM(sys_siL(ivt,10))//'" cx="'//&
            TRIM(sys_siL(INT(dscale*xdim),10))//'" cy="'//&
            TRIM(sys_siL(ysize-INT(dscale*ydim),10))//&
            '" r="0pt" fill="black" stroke="black" stroke-width="1pt"'
      END IF
      
      ! Write data which is common to all vertices
      WRITE(iunit,FMT='(A)',ADVANCE='NO') ' onmousemove="ShowVertexInfo(evt,'''//&
          TRIM(sys_siL(ivt,10))//''','''//&
          TRIM(sys_siL(rhadapt%p_IvertexAge(ivt),5))//''','''
      
      ! Generate list of elements meeting at vertex
      ipos=arrlst_getNextInArrayList(rhadapt%rElementsAtVertex,ivt,.TRUE.)
      DO WHILE(ipos .GT. ARRLST_NULL)
        ! Get element number IEL
        iel=rhadapt%rElementsAtVertex%p_IData(ipos)
        
        ! Proceed to next entry in array list
        ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivt,.FALSE.)

        IF (ipos .GT. ARRLST_NULL) THEN
          WRITE(iunit,FMT='(A)',ADVANCE='NO') TRIM(sys_siL(iel,10))//','
        ELSE
          WRITE(iunit,FMT='(A)',ADVANCE='NO') TRIM(sys_siL(iel,10))
        END IF
      END DO
      WRITE(iunit,FMT='(A)') ''')"'
      WRITE(iunit,FMT='(A)') ' onmouseout="HideVertexInfo()" style="cursor: crosshair"/>'
    END DO


    ! Output info boxes
    WRITE(iunit,FMT='(A)') '<g id="vinfo" style="visibility: hidden">'
    WRITE(iunit,FMT='(A)') '  <rect id="vinfo_box"  x="0" y="0" width="200" height="110"'
    WRITE(iunit,FMT='(A)') '        style="fill: #FFFFCC; stroke: #000000; stroke-width: 0.5px;"/>'
    WRITE(iunit,FMT='(A)') '  <text id="vinfo_text1" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black">Vertex number:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="vinfo_text2" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black">Vertex age:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="vinfo_text3" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black">Elements meeting at vertex:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="vinfo_text4" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black">null</text>'
    WRITE(iunit,FMT='(A)') '  <text id="vinfo_text5" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black">null</text>'
    WRITE(iunit,FMT='(A)') '  <text id="vinfo_text6" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black">null</text>'
    WRITE(iunit,FMT='(A)') '</g>'

    WRITE(iunit,FMT='(A)') '<g id="einfo" style="visibility: hidden">'
    WRITE(iunit,FMT='(A)') '  <rect id="einfo_box"  x="0" y="0" width="200" height="210"'
    WRITE(iunit,FMT='(A)') '        style="fill: #FFFFCC; stroke: #000000; stroke-width: 0.5px;"/>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text1" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black">Element number:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text2" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">Element state:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text3" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">Element marker:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text4" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">Element neighbours:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text5" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">Element mid-neighbours:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text6" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">Element vertices:</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text7" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">null</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text8" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">null</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text9" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">null</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text10" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">null</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text11" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">null</text>'
    WRITE(iunit,FMT='(A)') '  <text id="einfo_text12" x="0" y="0" style="font-size:'//&
        TRIM(sys_siL(font_size,3))//'pt;font-family:'//TRIM(ADJUSTL(font_family))//&
        ';stroke:black;">null</text>'
    WRITE(iunit,FMT='(A)') '</g>'

    ! Close XML-file
    WRITE(iunit,FMT='(A)') '</svg>'
    CLOSE(iunit)
  END SUBROUTINE hadapt_writeGridSVG

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_writeGridGMV(rhadapt,coutputFile)

!<description>
    ! This subroutine outputs the current state of the adapted grid stored
    ! in the dynamic data structure to a given file in GMV format.
!</description>

!<input>
    ! Output file name w/o suffix .gmv
    CHARACTER(LEN=*), INTENT(IN)     :: coutputFile
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local parameters
    INTEGER(PREC_VERTEXIDX) :: ivt
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: iunit,nve
    INTEGER, SAVE :: iout=0

    ! Check if dynamic data structures generated
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA) THEN
      CALL output_line('Dynamic data structures are not generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridSVG')
      CALL sys_halt()
    END IF
    
    ! Increment the sample number
    iout=iout+1

    ! Open output file for writing
    CALL io_openFileForWriting(TRIM(ADJUSTL(coutputFile))//'.'//&
        TRIM(sys_siL(iout,5))//'.gmv',iunit,SYS_REPLACE,bformatted=.TRUE.)
    WRITE(UNIT=iunit,FMT='(A)') 'gmvinput ascii'

    ! Write vertices to output file
    WRITE(UNIT=iunit,FMT=*) 'nodes ', rhadapt%NVT
    DO ivt=1,qtree_getsize(rhadapt%rVertexCoordinates2D)
      WRITE(UNIT=iunit,FMT=10) qtree_getX(rhadapt%rVertexCoordinates2D,ivt)
    END DO
    DO ivt=1,qtree_getsize(rhadapt%rVertexCoordinates2D)
      WRITE(UNIT=iunit,FMT=10) qtree_getY(rhadapt%rVertexCoordinates2D,ivt)
    END DO
    DO ivt=1,qtree_getsize(rhadapt%rVertexCoordinates2D)
      WRITE(UNIT=iunit,FMT=10) 0._DP
    END DO

    ! Write cells to output file
    WRITE(UNIT=iunit,FMT=*) 'cells ', rhadapt%NEL
    DO iel=1,rhadapt%NEL
      nve=hadapt_getNVE(rhadapt,iel)
      SELECT CASE(nve)
      CASE(TRIA_NVETRI2D)
        WRITE(UNIT=iunit,FMT=*) 'tri 3'
        WRITE(UNIT=iunit,FMT=20) rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,iel)

      CASE(TRIA_NVEQUAD2D)
        WRITE(UNIT=iunit,FMT=*) 'quad 4'
        WRITE(UNIT=iunit,FMT=30) rhadapt%p_IverticesAtElement(1:TRIA_NVEQUAD2D,iel)
        
      CASE DEFAULT
        CALL output_line('Invalid element type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
        CALL sys_halt()
      END SELECT
    END DO

    ! Write velocity to output file
    WRITE(UNIT=iunit,FMT=*) 'velocity 1'
    DO ivt=1,rhadapt%NVT
      WRITE(UNIT=iunit,FMT=10) 0._DP
      WRITE(UNIT=iunit,FMT=10) 0._DP
      WRITE(UNIT=iunit,FMT=10) 0._DP
    END DO

    ! Write variable
    WRITE(UNIT=iunit,FMT=*) 'variable'
    WRITE(UNIT=iunit,FMT=*) 'vert_age 1'

    DO ivt=1,SIZE(rhadapt%p_IvertexAge,1)
      WRITE(UNIT=iunit,FMT=40) rhadapt%p_IvertexAge(ivt)
    END DO
    DO ivt=SIZE(rhadapt%p_IvertexAge,1)+1,rhadapt%NVT
      WRITE(UNIT=iunit,FMT=40) -99999
    END DO

    WRITE(UNIT=iunit,FMT=*) 'endvars'
    WRITE(UNIT=iunit,FMT=*) 'probtime ', 0._DP
    WRITE(UNIT=iunit,FMT=*) 'endgmv'

    ! Close output file
    CLOSE(iunit)

10  FORMAT(E15.6E3)
20  FORMAT(3(1X,I8))
30  FORMAT(4(1X,I8))
40  FORMAT(I8)
  END SUBROUTINE hadapt_writeGridGMV

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_checkConsistency(rhadapt)

!<description>
    ! This subroutine checks the internal consistency of the dynamic data structures.
    ! Note that this routine performs brute-force search, and hence, should not
    ! be called in a productive environment. In is meant for debugging purposes
    ! only. If an error occurs, it stop without possibility to resume.
!</description>

!<inputoutput>
    ! adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_VERTEXIDX)    :: ivt,idx
    INTEGER(PREC_ELEMENTIDX)   :: iel,jel,jelmid
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER :: h_IelementsAtVertexIdx,h_IelementsAtVertex
    INTEGER :: ive,jve,nve,mve
    LOGICAL :: btest,bfound

    ! Test #1: Consistency of element numbers
    btest=(rhadapt%NEL .EQ. SUM(rhadapt%InelOfType))
    CALL output_line('Test #1: Checking consistency of element numbers '//&
        MERGE('PASSED','FAILED',btest))

    ! Test #2: Vertex age must not exceed maximum refinement level
    btest=.TRUE.
    DO ivt=1,rhadapt%NVT
      btest=btest .OR. (rhadapt%p_IvertexAge(ivt) .GT. rhadapt%NSUBDIVIDEMAX)
    END DO
    CALL output_line('Test #2: Checking maximum vertex age '//&
        MERGE('PASSED','FAILED',btest))

    ! Test #3: Check consistency of element neighbours
    btest=.TRUE.
    DO iel=1,rhadapt%NEL

      ! Get number of vertices per element
      nve=hadapt_getNVE(rhadapt,iel)
      
      ! Loop over all adjacent elements
      DO ive=1,nve
        jel   =rhadapt%p_IneighboursAtElement(ive,iel)
        jelmid=rhadapt%p_ImidneighboursAtElement(ive,iel)

        ! Check that adjacent element number is not larger than the
        ! total number of elements present in the triangulation
        IF (jel > rhadapt%NEL .OR. jelmid > rhadapt%NEL) THEN
          btest=.FALSE.; CYCLE
        END IF

        ! Do nothing if we are adjacent to the boundary
        IF (jel .EQ. 0 .OR. jelmid .EQ. 0) CYCLE
        
        ! Is neighboring element subdivided?
        IF (jel .EQ. jelmid) THEN
          
          ! Get number of vertices per element
          mve=hadapt_getNVE(rhadapt,jel)
          
          ! Find element IEL in adjacency list of JEL
          bfound=.FALSE.
          DO jve=1,mve
            IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. &
                rhadapt%p_ImidneighboursAtElement(jve,jel)) THEN
              IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. iel) THEN
                bfound=.TRUE.; EXIT
              END IF
            ELSE
              IF (rhadapt%p_IneighboursAtElement(jve,jel)    .EQ. iel .OR.&
                  rhadapt%p_ImidneighboursAtElement(jve,jel) .EQ. iel) THEN
                bfound=.TRUE.; EXIT
              END IF
            END IF
          END DO
          btest=btest.AND.bfound
          
        ELSE

          ! Get number of vertices per element
          mve=hadapt_getNVE(rhadapt,jel)
          
          ! Find element IEL in adjacency list of JEL
          bfound=.FALSE.
          DO jve=1,mve
            IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. &
                rhadapt%p_ImidneighboursAtElement(jve,jel)) THEN
              IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. iel) THEN
                bfound=.TRUE.; EXIT
              END IF
            ELSE
              IF (rhadapt%p_IneighboursAtElement(jve,jel)    .EQ. iel .OR.&
                  rhadapt%p_ImidneighboursAtElement(jve,jel) .EQ. iel) THEN
                bfound=.TRUE.; EXIT
              END IF
            END IF
          END DO
          btest=btest.AND.bfound

          ! Get number of vertices per element
          mve=hadapt_getNVE(rhadapt,jelmid)
          
          ! Find element IEL in adjacency list of JELMID
          bfound=.FALSE.
          DO jve=1,mve
            IF (rhadapt%p_IneighboursAtElement(jve,jelmid) .EQ. &
                rhadapt%p_ImidneighboursAtElement(jve,jelmid)) THEN
              IF (rhadapt%p_IneighboursAtElement(jve,jelmid) .EQ. iel) THEN
                bfound=.TRUE.; EXIT
              END IF
            ELSE
              IF (rhadapt%p_IneighboursAtElement(jve,jelmid)    .EQ. iel .OR.&
                  rhadapt%p_ImidneighboursAtElement(jve,jelmid) .EQ. iel) THEN
                bfound=.TRUE.; EXIT
              END IF
            END IF
          END DO
          btest=btest.AND.bfound

        END IF
      END DO
    END DO
    CALL output_line('Test #3: Checking consistency of element neighbours '//&
        MERGE('PASSED','FAILED',btest))

    ! Test #4: Check consistency of common vertices between two edges
    btest=.TRUE.
    DO iel=1,rhadapt%NEL

      ! Get number of vertices per element
      nve=hadapt_getNVE(rhadapt,iel)
      
      ! Loop over all adjacent elements
      DO ive=1,nve
        jel   =rhadapt%p_IneighboursAtElement(ive,iel)
        jelmid=rhadapt%p_ImidneighboursAtElement(ive,iel)
        
        ! Check that adjacent element number is not larger than the
        ! total number of elements present in the triangulation
        IF (jel > rhadapt%NEL .OR. jelmid > rhadapt%NEL) THEN
          btest=.FALSE.; CYCLE
        END IF

        ! Do nothing if we are adjacent to the boundary
        IF (jel .EQ. 0 .OR. jelmid .EQ. 0) CYCLE
        
        ! Do nothing if there exists a temporal hanging node
        IF (jel .NE. jelmid) CYCLE

        ! Get number of vertices per element
        mve=hadapt_getNVE(rhadapt,jel)

        ! Find element IEL in adjacency list of JEL
        bfound=.FALSE.
        DO jve=1,mve
          IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. &
              rhadapt%p_ImidneighboursAtElement(jve,jel)) THEN
            IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. iel) THEN
              bfound=.TRUE.; EXIT
            END IF
          ELSE
            ! Do nothing if there exists a temporal hanging node
            CYCLE
          END IF
        END DO
        
        ! If the common edge has been found, check the two endpoints
        IF (bfound) THEN
          bfound=((rhadapt%p_IverticesAtElement(ive,iel) .EQ. &
                     rhadapt%p_IverticesAtElement(MODULO(jve,mve)+1,jel)) .AND. &
                     rhadapt%p_IverticesAtElement(MODULO(ive,nve)+1,iel) .EQ. &
                     rhadapt%p_IverticesAtElement(jve,jel))
        END IF
        btest=btest.AND.bfound
      END DO
    END DO
    CALL output_line('Test #4: Checking consistency of common vertices along edges '//&
        MERGE('PASSED','FAILED',btest))

    ! Test #5: Check consistency of element-meeting-at-vertex lists
    btest=(rhadapt%rElementsAtVertex%NTABLE .EQ. rhadapt%NVT)
    IF (btest) THEN
      
      ! Create index array
      CALL storage_new('hadapt_checkConsistency','IelementAtVertexIdx',rhadapt%NVT+1,&
          ST_INT,h_IelementsAtVertexIdx,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_IelementsAtVertexIdx,p_IelementsAtVertexIdx)

      ! Count number of elements meeting at vertex
      DO iel=1,rhadapt%NEL

        ! Get number of vertices per element
        nve=hadapt_getNVE(rhadapt,iel)

        ! Loop over corner vertices
        DO ive=1,nve
          ivt=rhadapt%p_IverticesAtElement(ive,iel)
          p_IelementsAtVertexIdx(ivt+1)=p_IelementsAtVertexIdx(ivt+1)+1
        END DO
      END DO

      ! Convert element couter into absolute position
      p_IelementsAtVertexIdx(1)=1
      DO ivt=2,rhadapt%NVT+1
        p_IelementsAtVertexIdx(ivt)=p_IelementsAtVertexIdx(ivt)+p_IelementsAtVertexIdx(ivt-1)
      END DO

      ! Create working array
      CALL storage_new('hadapt_checkConsistency','IelementAtVertex',&
          p_IelementsAtVertexIdx(rhadapt%NVT+1)-1,&
          ST_INT,h_IelementsAtVertex,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int(h_IelementsAtVertex,p_IelementsAtVertex)

      ! Retrieve the element numbers
      DO iel=1,rhadapt%NEL

        ! Get number of vertices per element
        nve=hadapt_getNVE(rhadapt,iel)

        ! Loop over corner vertices
        DO ive=1,nve
          ivt=rhadapt%p_IverticesAtElement(ive,iel)
          idx=p_IelementsAtVertexIdx(ivt)
          p_IelementsAtVertexIdx(ivt)=idx+1
          p_IelementsAtVertex(idx)   =iel
        END DO
      END DO
      
      ! Restore index array
      DO ivt=rhadapt%NVT+1,2,-1
        p_IelementsAtVertexIdx(ivt)=p_IelementsAtVertexIdx(ivt-1)
      END DO
      p_IelementsAtVertexIdx(1)=1

      ! Start to compare the temporal elements-meeting-at-vertex list
      ! and the dynamic data structure from the adaptivity structure
      DO ivt=1,rhadapt%NVT
        
        ! Get first entry in array list
        ipos=arrlst_getNextInArrayList(rhadapt%rElementsAtVertex,ivt,.TRUE.)
        
        ! Repeat until there is no entry left in the array list
        DO WHILE(ipos .GT. ARRLST_NULL)
          
          ! Get element number IEL
          iel=rhadapt%rElementsAtVertex%p_IData(ipos)

          ! Proceed to next entry in array list
          ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivt,.FALSE.)

          ! Look for element IEL in temporal elements-meeting-at-vertex list
          ! If it does exist, multiply its value by minus one so that it cannot
          ! be found twice. At the end, all positive entries in the temporal
          ! list are not present in the dynamic structure
          bfound=.FALSE.
          DO idx=p_IelementsAtVertexIdx(ivt),p_IelementsAtVertexIdx(ivt+1)-1
            IF (p_IelementsAtVertex(idx) .EQ. iel) THEN
              p_IelementsAtVertex(idx)=-p_IelementsAtVertex(idx)
              bfound=.TRUE.
              EXIT
            END IF
          END DO
          btest=btest.AND.bfound
        END DO
      END DO
      btest=btest.AND.ALL(p_IelementsAtVertex < 0)
      CALL output_line('Test #5: Checking consistency of elements meeting at vertices '//&
        MERGE('PASSED','FAILED',btest))

      ! Release auxiliary storage
      CALL storage_free(h_IelementsAtVertexIdx)
      CALL storage_free(h_IelementsAtVertex)
    ELSE

      CALL output_line('Test #5: Checking consistency of element-meeting-at-vertex list '//&
          MERGE('PASSED','FAILED',btest))
    END IF
  END SUBROUTINE hadapt_checkConsistency

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE hadapt_refreshAdaptation(rhadapt,rtriangulation)

!<description>
    ! This subroutine refreshes the pointers and internal structure of the 
    ! adaptation structure.
    ! NOTE: The triangulation structure must be compatible with the adaptivity 
    ! structure. This CANNOT be checked by the routine and must be guaranteed
    ! by the used. The philosophy behind this adaptivity structure is to not
    ! modify the triangulation externally. If external modifications are performed
    ! then all adaptivity structures must be released and rebuild from scratch.
    ! This is quite time-consuming, and moreover, the deletion of vertices is
    ! not possible since all vertices will be considered as initial vertices.
!</description>

!<input>
    ! Triangulation structure
    TYPE(t_triangulation), INTENT(IN) :: rtriangulation
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>


    ! Get pointer to InodalProperty
    IF (rhadapt%h_InodalProperty .EQ. &
        rtriangulation%h_InodalProperty) THEN
      CALL storage_getbase_int(rhadapt%h_InodalProperty,&
          rhadapt%p_InodalProperty)
    ELSE
      CALL output_line('Inconsistent handle h_InodalProperty',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      CALL sys_halt()
    END IF

    ! Get pointer to IverticesAtElement
    IF (rhadapt%h_IverticesAtElement .EQ. &
        rtriangulation%h_IverticesAtElement) THEN
      CALL storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
          rhadapt%p_IverticesAtElement)
    ELSE
      CALL output_line('Inconsistent handle h_IverticesAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      CALL sys_halt()
    END IF

    ! Get pointer to IneighboursAtElement
    IF (rhadapt%h_IneighboursAtElement .EQ. &
        rtriangulation%h_IneighboursAtElement) THEN
      CALL storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
          rhadapt%p_IneighboursAtElement)
    ELSE
      CALL output_line('Inconsistent handle h_IneighboursAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      CALL sys_halt()
    END IF

    ! Get pointer to ImidneighboursAtElement
    IF (rhadapt%h_ImidneighboursAtElement .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
          rhadapt%p_ImidneighboursAtElement)
    ELSE
      CALL output_line('Inconsistent handle h_ImidneighboursAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      CALL sys_halt()
    END IF

    ! Get pointer to IvertexAge
     IF (rhadapt%h_IvertexAge .NE. ST_NOHANDLE) THEN
      CALL storage_getbase_int(rhadapt%h_IvertexAge,&
          rhadapt%p_IvertexAge)
    ELSE
      CALL output_line('Inconsistent handle h_ImidneighboursAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      CALL sys_halt()
    END IF

  END SUBROUTINE hadapt_refreshAdaptation

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  ! Internal auxiliary routines and functions. No user should ever look at them

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE redgreen_refine(rhadapt,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine performs red-green refinement as proposed by R. Bank
!</description>

!<input>
    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</subroutine>
    
    ! local variables
    INTEGER, DIMENSION(:), POINTER :: p_Imarker
    INTEGER(PREC_ELEMENTIDX) :: iel
    
    ! Check if dynamic data structures are o.k. and if 
    ! cells are marked for refinement
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA .OR.&
        IAND(rhadapt%iSpec,HADAPT_MARKEDREFINE).NE.HADAPT_MARKEDREFINE) THEN
      CALL output_line('Dynamic data structures are not generated &
          &or no marker for refinement is available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'redgreen_refine')
      CALL sys_halt()
    END IF
    
    ! Set pointers
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
        
    ! Perform red-green refinement
    DO iel=1,SIZE(p_Imarker,1)
      
      SELECT CASE(p_Imarker(iel))
      CASE(MARK_ASIS_TRIA,MARK_ASIS_QUAD)
        ! Do nothing for elements that should be kept 'as is'

      CASE(:MARK_CRS_GENERIC)
        ! Do nothing for element that have been marked for coarsening

      CASE(MARK_REF_TRIA4TRIA)
        ! Red refinement triangle
        CALL refine_Tria4Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)
        p_Imarker(iel)=MARK_ASIS        

      CASE(MARK_REF_QUAD4QUAD)
        ! Red refinement quadrilateral
        CALL refine_Quad4Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)
        p_Imarker(iel)=MARK_ASIS
        
      CASE(MARK_REF_TRIA2TRIA_1,MARK_REF_TRIA2TRIA_2,MARK_REF_TRIA2TRIA_3)   
        ! Green refinement triangle
        CALL refine_Tria2Tria(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements=rhadapt%nGreenElements+2
        p_Imarker(iel)=MARK_ASIS
        
      CASE(MARK_REF_QUAD3TRIA_1,MARK_REF_QUAD3TRIA_2,&
           MARK_REF_QUAD3TRIA_3,MARK_REF_QUAD3TRIA_4)
        ! Green refinement quadrilateral
        CALL refine_Quad3Tria(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements=rhadapt%nGreenElements+3
        p_Imarker(iel)=MARK_ASIS
        
      CASE(MARK_REF_QUAD2QUAD_13,MARK_REF_QUAD2QUAD_24)
        ! Green refinement quadrilateral
        CALL refine_Quad2Quad(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements=rhadapt%nGreenElements+2
        p_Imarker(iel)=MARK_ASIS

      CASE(MARK_REF_QUAD4TRIA_12,MARK_REF_QUAD4TRIA_23,&
           MARK_REF_QUAD4TRIA_34,MARK_REF_QUAD4TRIA_14)
        ! Green refinement quadrilateral
        CALL refine_Quad4Tria(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements=rhadapt%nGreenElements+4
        p_Imarker(iel)=MARK_ASIS

      CASE DEFAULT
        CALL output_line('Invalid element refinement marker!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_refine')
        CALL sys_halt()
      END SELECT
    END DO

    ! Increase the number of refinement steps by one
    rhadapt%nRefinementSteps=rhadapt%nRefinementSteps+1
    
    ! Refinement has been performed.
    rhadapt%iSpec=IOR(rhadapt%ispec,HADAPT_REFINED)
    
    ! Hence, the markers are no longer valid
    rhadapt%iSpec=IAND(rhadapt%iSpec,NOT(HADAPT_MARKEDREFINE))
  END SUBROUTINE redgreen_refine

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE redgreen_coarsen(rhadapt,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine performs red-green coarsening as proposed by R. Bank
!</description>

!<input>
    ! callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt
    
    ! OPTIONAL Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(1) :: Ielements
    INTEGER(PREC_VERTEXIDX), DIMENSION(2)  :: Ivertices
    INTEGER,  DIMENSION(:), POINTER :: p_Imarker
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: iel,jel
    INTEGER(PREC_VERTEXIDX)  :: ivt,ivtReplace
    INTEGER :: ive

    ! Check if dynamic data structures are o.k. and 
    ! if  cells are marked for coarsening
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA .OR.&
        IAND(rhadapt%iSpec,HADAPT_MARKEDCOARSEN).NE.HADAPT_MARKEDCOARSEN) THEN
      CALL output_line('Dynamic data structures are not generated &
          &or no marker for coarsening is available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'redgreen_coarsen')
      CALL sys_halt()
    END IF
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
    
    ! Perform hierarchical red-green recoarsening
    DO iel=SIZE(p_Imarker,1),1,-1

      SELECT CASE(p_Imarker(iel))
      CASE(MARK_CRS_GENERIC:)
        ! Do nothing for elements ...
        ! - that should be kept 'as is'
        ! - that are only marked for generic recoarsening
        ! - that are marked for refinement.
        
      CASE(MARK_CRS_2TRIA1TRIA)
        CALL coarsen_2Tria1Tria(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)
        
      CASE(MARK_CRS_4TRIA1TRIA)
        CALL coarsen_4Tria1Tria(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_4TRIA2TRIA_1,&
           MARK_CRS_4TRIA2TRIA_2,&
           MARK_CRS_4TRIA2TRIA_3)
        CALL coarsen_4Tria2Tria(rhadapt,iel,p_Imarker(iel),&
            rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_3TRIA1QUAD)
        CALL coarsen_3Tria1Quad(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_4TRIA1QUAD)
        CALL coarsen_4Tria1Quad(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)
        
      CASE(MARK_CRS_4TRIA3TRIA_LEFT,&
           MARK_CRS_4TRIA3TRIA_RIGHT)
        CALL coarsen_4Tria3Tria(rhadapt,iel,p_Imarker(iel),&
            rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_2QUAD1QUAD)
        CALL coarsen_2Quad1Quad(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_2QUAD3TRIA)
        CALL coarsen_2Quad3Tria(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)
        
      CASE(MARK_CRS_4QUAD1QUAD)
        CALL coarsen_4Quad1Quad(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_4QUAD2QUAD)
        CALL coarsen_4Quad2Quad(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_4QUAD3TRIA)
        CALL coarsen_4Quad3Tria(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)
        
      CASE(MARK_CRS_4QUAD4TRIA)
        CALL coarsen_4Quad4Tria(rhadapt,iel,&
            rcollection,fcb_hadaptCallback)

      CASE DEFAULT
        CALL output_line('Invalid recoarsening marker!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_coarsen')
        CALL sys_halt()
      END SELECT
    END DO

    ! Loop over all vertices 1...NVT0 present in the triangulation before
    ! refinement and check if they are free for vertex removal.
    DO ivt=rhadapt%NVT0,1,-1
      
      ! If the vertex is locked, then skip this vertex
      IF (rhadapt%p_IvertexAge(ivt) .LE. 0) CYCLE
      
      ! Remove vertex physically. Note that this vertex is no longer associated
      ! to any element. All associations have been removed in the above element
      ! coarsening/conversion step. In order to prevent "holes" in the vertex list,
      ! vertex IVT is replaced by the last vertex if it is not the last one itself.
      CALL remove_vertex2D(rhadapt,ivt,ivtReplace)
      
      ! If vertex IVT was not the last one, update the "elements-meeting-at-vertex" list
      IF (ivtReplace .NE. 0) THEN
        
        ! Start with first element in "elements-meeting-at-vertex" list of the replaced vertex
        ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivtReplace,.TRUE.)
        update: DO WHILE(ipos .GT. ARRLST_NULL)
          
          ! Get element number JEL
          jel=rhadapt%rElementsAtVertex%p_IData(ipos)
          
          ! Proceed to next element
          ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivtReplace,.FALSE.)
          
          ! Look for vertex ivtReplace in element JEL and replace it by IVT
          DO ive=1,hadapt_getNVE(rhadapt,jel)
            IF (rhadapt%p_IverticesAtElement(ive,jel) .EQ. ivtReplace) THEN
              rhadapt%p_IverticesAtElement(ive,jel)=ivt
              CYCLE update
            END IF
          END DO
          
          ! If the replaced vertex ivtReplace could not be found in element JEL
          ! something is wrong and we stop the simulation
          CALL output_line('Unable to find replacement vertex in element',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_coarsen')
          CALL sys_halt()            
        END DO update
        
        ! Swap tables IVT and ivtReplace in arraylist and release table ivtReplace
        CALL arrlst_swapArrayList(rhadapt%rElementsAtVertex,ivt,ivtReplace)
        CALL arrlst_releaseArrayList(rhadapt%rElementsAtVertex,ivtReplace)
        
      ELSE
        
        ! Release table IVT
        CALL arrlst_releaseArrayList(rhadapt%rElementsAtVertex,ivt)
      END IF
            
      ! Optionally, invoke callback function
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/ivt,ivtReplace/); Ielements=(/0/)
        CALL fcb_hadaptCallback(rcollection,&
            HADAPT_OPR_REMOVEVERTEX,Ivertices,Ielements)
      END IF
    END DO
        
    ! Increase the number of recoarsening steps by one
    rhadapt%nCoarseningSteps=rhadapt%nCoarseningSteps+1

    ! Coarsening has been performed.
    rhadapt%iSpec=IOR(rhadapt%ispec,HADAPT_COARSENED)
    
    ! Hence, the markers are no longer valid
    rhadapt%iSpec=IAND(rhadapt%iSpec,NOT(HADAPT_MARKEDCOARSEN))
  END SUBROUTINE redgreen_coarsen
END MODULE hadaptivity

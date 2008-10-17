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
!#  7.) hadapt_setVertexCoords1D
!#      -> Set the coordinates of vertices to the adaptivity structure in 1D
!#
!#  8.) hadapt_getVertexCoords1D
!#      -> Get the coordinates of vertices from the adaptivity structure in 1D
!#
!#  9.) hadapt_setVertexCoords2D
!#      -> Set the coordinates of vertices to the adaptivity structure in 2D
!#
!# 10.) hadapt_getVertexCoords2D
!#      -> Get the coordinates of vertices from the adaptivity structure in 2D
!#
!# 11.) hadapt_setVertexCoords3D
!#      -> Set the coordinates of vertices to the adaptivity structure in 3D
!#
!# 12.) hadapt_getVertexCoords3D
!#      -> Get the coordinates of vertices from the adaptivity structure in 3D
!#
!# 13.) hadapt_setVerticesAtElement
!#      -> Set the "vertices-at-element" structure to the adaptivity structure
!#
!# 14.) hadapt_getVerticesAtElement
!#      -> Get the "vertices-at-element" structure from the adaptivity structure
!#
!# 15.) hadapt_setNeighboursAtElement
!#      -> Set the "neighbours-at-element" structure to the adaptivity structure
!#
!# 16.) hadapt_getNeighboursAtElement
!#      -> Get the "neighbours-at-element" structure from the adaptivity structure
!#
!# 17.) hadapt_setNelOfType
!#      -> Set the "InelOfType" array to the adaptivity structure
!#
!# 18.) hadapt_getNelOfType
!#      -> Get the "InelOfType" array from the adaptivity structure
!#
!# 19.) hadapt_setBoundary
!#      -> Set the boundary structure to the adaptivity structure
!#
!# 20.) hadapt_getBoundary
!#      -> Get the boundary structure form the adaptivity structure
!#
!# 21.) hadapt_setNodalProperty
!#      -> Set the "nodal property" list to the adaptivity structure
!#
!# 22.) hadapt_getNodalProperty
!#      -> Get the "nodal property" list from the adaptivity structure
!#
!# 23.) hadapt_genElementsAtVertex
!#      -> Generate the "elements-at-vertex" data structure
!#
!# 24.) hadapt_performAdaptation
!#      -> perform one step of grid adaptation
!#
!# 25.) hadapt_infoStatistics
!#      -> output information about the adaptivity structure
!#
!# 26.) hadapt_writeGridSVG
!#      -> write the adapted grid to file in SVG format
!#
!# 27.) hadapt_writeGridGMV
!#      -> write the adapted grid to file in GMV format
!#
!# 28.) hadapt_checkConsistency
!#      -> check the internal consistency of dynamic data structures
!#
!# 29.) hadapt_refreshAdaptation
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

module hadaptivity

  use arraylist
  use binarytree
  use boundary
  use collection
  use hadaptaux3d
  use io
  use linearsystemscalar
  use list
  use octree
  use paramlist
  use quadtree
  use sort
  use storage
  use triangulation
  use fsystem

  implicit none

  private
  public :: t_hadapt
  public :: hadapt_initFromParameterlist
  public :: hadapt_initFromTriangulation
  public :: hadapt_generateRawMesh
  public :: hadapt_releaseAdaptation
  public :: hadapt_duplicateAdaptation
  public :: hadapt_restoreAdaptation
  public :: hadapt_setVertexCoords2D
  public :: hadapt_getVertexCoords2D
  public :: hadapt_setVerticesAtElement
  public :: hadapt_getVerticesAtElement
  public :: hadapt_setNeighboursAtElement
  public :: hadapt_getNeighboursAtElement
  public :: hadapt_setNelOfType
  public :: hadapt_getNelOfType
  public :: hadapt_setBoundary
  public :: hadapt_getBoundary
  public :: hadapt_setNodalProperty
  public :: hadapt_getNodalProperty
  public :: hadapt_genElementsAtVertex
  public :: hadapt_performAdaptation
  public :: hadapt_infoStatistics
  public :: hadapt_writeGridSVG
  public :: hadapt_writeGridGMV
  public :: hadapt_checkConsistency
  public :: hadapt_refreshAdaptation

contains
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_initFromParameterlist(rhadapt, rparlist, ssection)

!<description>
    ! This subroutine initializes the adaptivity structure
    ! with the values supplied by the parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN)  :: rparlist

    ! name of the section
    character(LEN=*), intent(IN) :: ssection
!</input>

!<output>
    ! adaptivity structure
    type(t_hadapt), intent(OUT)  :: rhadapt
!</output>
!</subroutine>

    ! Get mandatory parameters from list
    call parlst_getvalue_int   (rparlist, ssection, "nsubdividemax",&
                                rhadapt%nsubdividemax)
    call parlst_getvalue_int   (rparlist, ssection, "iadaptationStrategy",&
                                rhadapt%iadaptationStrategy)
    call parlst_getvalue_double(rparlist, ssection, "drefinementTolerance",&
                                rhadapt%drefinementTolerance)
    call parlst_getvalue_double(rparlist, ssection, "dcoarseningTolerance",&
                                rhadapt%dcoarseningTolerance)

    ! Initialize data
    rhadapt%iSpec = HADAPT_HAS_PARAMETERS
  end subroutine hadapt_initFromParameterlist

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_initFromTriangulation(rhadapt, rtriangulation)

!<description>
    ! This subroutine initializes all required components of the adaptativit 
    ! structure from the triangulation structure rtriangulation.
!</description>

!<input>
    ! Triangulation structure
    type(t_triangulation), intent(IN) :: rtriangulation
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT)     :: rhadapt
!</inputoutput>
!</subroutine>

    ! Initialize duplication flag
    rhadapt%iduplicationFlag = 0

    ! Set dimension
    rhadapt%ndim = rtriangulation%ndim

    ! Set coordinates
    select case(rhadapt%ndim)
    case (NDIM1D)
      call hadapt_setVertexCoords1D(rhadapt,&
          rtriangulation%h_DvertexCoords, rtriangulation%NVT)

    case(NDIM2D)
      call hadapt_setVertexCoords2D(rhadapt,&
          rtriangulation%h_DvertexCoords, rtriangulation%NVT)

    case(NDIM3D)
      call hadapt_setVertexCoords3D(rhadapt,&
          rtriangulation%h_DvertexCoords, rtriangulation%NVT)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_initFromTriangulation')
      call sys_halt()
    end select

    ! Set nodal properties
    call hadapt_setNodalProperty(rhadapt,&
        rtriangulation%h_InodalProperty)

    ! Set element numbers
    call hadapt_setNelOfType(rhadapt,&
        rtriangulation%InelOfType)
    
    ! Set vertices at element
    call hadapt_setVerticesAtElement(rhadapt,&
        rtriangulation%h_IverticesAtElement, rtriangulation%NEL)

    ! Set elements adjacent to element
    call hadapt_setNeighboursAtElement(rhadapt,&
        rtriangulation%h_IneighboursAtElement)
    
    ! Set boundary for 2D
    if (rtriangulation%ndim .eq. NDIM2D) then
      call hadapt_setBoundary(rhadapt, rtriangulation%h_IboundaryCpIdx,&
                              rtriangulation%h_IverticesAtBoundary,&
                              rtriangulation%h_DvertexParameterValue,&
                              rtriangulation%NBCT, rtriangulation%NVBD)
    end if

print *, rtriangulation%ndim

print *, rhadapt%ndim
    ! Generate "elements-meeting-at-vertex" structure
    call hadapt_genElementsAtVertex(rhadapt)
    
    ! Create generation array and initialize all nodes with "age" 0
    call storage_new('hadapt_initFromTriangulation','p_IvertexAge',&
        rhadapt%NVT, ST_INT, rhadapt%h_IvertexAge, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rhadapt%h_IvertexAge, rhadapt%p_IvertexAge)
  end subroutine hadapt_initFromTriangulation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_generateRawMesh(rhadapt, rtriangulation)

!<description>
    ! This subroutine generates an (extended) raw mesh from the 
    ! adaptivity structure.
!</description>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT)        :: rhadapt
    
    ! Triangulation structure
    type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
!</subroutine>

    ! Get dimensions
    rtriangulation%ndim=rhadapt%ndim

    ! Get coordinates
    select case(rhadapt%ndim)
    case(NDIM2D)
      call hadapt_getVertexCoords2D(rhadapt,&
          rtriangulation%h_DvertexCoords, rtriangulation%NVT)

    case(NDIM3D)
      call hadapt_getVertexCoords3D(rhadapt,&
          rtriangulation%h_DvertexCoords, rtriangulation%NVT)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_generateRawMesh')
      call sys_halt()
    end select

    ! Get number of elements
    call hadapt_getNelOfType(rhadapt,&
        rtriangulation%InelOfType)

    ! Get vertices at element list
    call hadapt_getVerticesAtElement(rhadapt,&
        rtriangulation%h_IverticesAtElement, rtriangulation%NEL)

    ! Get element neighbours
    call hadapt_getNeighboursAtElement(rhadapt,&
        rtriangulation%h_IneighboursAtElement)

    ! Get nodal property list
    call hadapt_getNodalProperty(rhadapt,&
        rtriangulation%h_InodalProperty)

    ! Get boundary
    call hadapt_getBoundary(rhadapt, rtriangulation%h_IboundaryCpIdx,&
                            rtriangulation%h_IverticesAtBoundary,&
                            rtriangulation%h_DvertexParameterValue,&
                            rtriangulation%NVBD, rtriangulation%NBCT)
                            
    ! Generate extended raw mesh information for edge (and face)
    ! numbering.
    call tria_initExtendedRawMesh (rtriangulation)
    
  end subroutine hadapt_generateRawMesh
  
  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_releaseAdaptation(rhadapt)

!<description>
    ! This subroutine releases all internal structures of the
    ! adaptivity data structure rhadapt.
!</description>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    integer(I32) :: idupflag
    integer :: ibct

    idupflag = rhadapt%iduplicationFlag

    ! Check if quadtree exists
    if (iand(rhadapt%iSpec,HADAPT_HAS_COORDS) .eq. HADAPT_HAS_COORDS) then
      select case(rhadapt%ndim)
      case(NDIM2D)
        call qtree_releaseQuadtree(rhadapt%rVertexCoordinates2D)
        
      case(NDIM3D)
        call otree_releaseOctree(rhadapt%rVertexCoordinates3D)

      case DEFAULT
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_releaseAdaptation')
        call sys_halt()
      end select
    end if   

    ! Check if boundary structure exists
    if (associated(rhadapt%rBoundary)) then
      do ibct = 1, size(rhadapt%rBoundary, 1)
        call btree_releaseTree(rhadapt%rBoundary(ibct))
      end do
      deallocate(rhadapt%rBoundary)
      nullify(rhadapt%rBoundary)
    end if

    ! Release elements-meeting-at-vertex arraylist
    call arrlst_releaseArraylist(rhadapt%relementsAtVertex)

    ! Release storage which is no longer in use
    call checkAndRelease(idupflag, HADAPT_SHARE_IMARKER,&
                         rhadapt%h_Imarker)
    call checkAndRelease(idupflag, HADAPT_SHARE_IVERTEXAGE,&
                         rhadapt%h_IvertexAge)
    call checkAndRelease(idupflag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
                         rhadapt%h_ImidneighboursAtElement)
    
    ! Nullify "performance-pointers"
    nullify(rhadapt%p_IvertexAge)
    nullify(rhadapt%p_InodalProperty)
    nullify(rhadapt%p_IverticesAtElement)
    nullify(rhadapt%p_IneighboursAtElement)
    nullify(rhadapt%p_ImidneighboursAtElement)
    
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

  contains
    
    subroutine checkAndRelease (idupFlag, ibitfield, ihandle)
      integer(I32), intent(IN) :: ibitfield
      integer(I32), intent(IN) :: idupFlag
      integer, intent(INOUT) :: ihandle
      
      if (iand(idupFlag, ibitfield) .ne. ibitfield) then
        if (ihandle .ne. ST_NOHANDLE) call storage_free(ihandle)
      else
        ihandle = ST_NOHANDLE
      end if
      
    end subroutine checkAndRelease
  end subroutine hadapt_releaseAdaptation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_duplicateAdaptation(rhadapt, rhadaptBackup,&
                                        iduplicationFlag, bupdate)

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
    type(t_hadapt), intent(IN) :: rhadapt

    ! Bitfield that decides which handles are a copy of another structure, thus 
    ! which arrays are shared between the new and the old structure. 
    ! The bitfield contains a combination of HADAPT_SHARE_xxxx constants. Every 
    ! information whose flag is set in iduplicationFlag is shared between rhadapt
    ! and rhadaptBackup.
    ! Therefore e.g., iduplicationFlag=0 copies all arrays from rhadapt in memory,
    ! while HADAPT_SHARE_ALL will copy nothing, but will share everything between 
    ! rhadapt and rhadaptBackup.
    integer(I32), intent(IN)   :: iduplicationFlag
  
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
    logical, intent(IN), optional     :: bupdate
!</input>

!<inputoutput>
    ! The destination structure which receives the information
    type(t_hadapt), intent(INOUT) :: rhadaptBackup
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: idupFlag
    integer :: ibct
    logical :: bupd

    bupd = .false.
    if (present(bupdate)) bupd = bupdate
    
    if (.not. bupd) then
      ! Release any old data.
      call hadapt_releaseAdaptation(rhadaptBackup)

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

    else

      ! Create a bitfiled what to copy by ORing iduplicationFlag with what
      ! we have in rhadaptBackup. That way, only arrays that exist as real
      ! duplicates are copied from rhadapt to rhadaptBackup.

      idupFlag = ior(iduplicationFlag, rhadaptBackup%iduplicationFlag)

    end if

    ! Call checkAndCopy for all components. This will either copy the handle
    ! or allocate new memory and copy the content of the component.

    ! Bit   0: Imarker
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMARKER,&
                      rhadapt%h_Imarker, rhadaptBackup%h_Imarker)

    ! Bit   1: IvertexAge
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTEXAGE,&
                      rhadapt%h_IvertexAge, rhadaptBackup%h_IvertexAge)
    call storage_getbase_int(rhadaptBackup%h_IvertexAge,&
                             rhadaptBackup%p_IvertexAge)

    ! Bit   2: InodalProperty
    call checkAndCopy(idupFlag, HADAPT_SHARE_INODALPROPERTY,&
                      rhadapt%h_InodalProperty, rhadaptBackup%h_InodalProperty)
    call storage_getbase_int(rhadaptBackup%h_InodalProperty,&
                             rhadaptBackup%p_InodalProperty)

    ! Bit   3: IverticesAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTICESATELEMENT,&
                      rhadapt%h_IverticesAtElement,&
                      rhadaptBackup%h_IverticesAtElement)
    call storage_getbase_int2D(rhadaptBackup%h_IverticesAtElement,&
                               rhadaptBackup%p_IverticesAtElement)

    ! Bit   4: IneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_INEIGHATELEMENT,&
                      rhadapt%h_IneighboursAtElement,&
                      rhadaptBackup%h_IneighboursAtElement)
    call storage_getbase_int2D(rhadaptBackup%h_IneighboursAtElement,&
                               rhadaptBackup%p_IneighboursAtElement)

    ! Bit   5: ImidneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
                      rhadapt%h_ImidneighboursAtElement,&
                      rhadaptBackup%h_ImidneighboursAtElement)
    call storage_getbase_int2D(rhadaptBackup%h_ImidneighboursAtElement,&
                               rhadaptBackup%p_ImidneighboursAtElement)

    ! Bit   6: rVertexCoordinates
    if (iand(idupFlag, HADAPT_SHARE_RVERTEXCOORDINATES) .ne.&
                       HADAPT_SHARE_RVERTEXCOORDINATES) then
      
      select case(rhadaptBackup%ndim)
      case(NDIM2D)
        call qtree_duplicateQuadtree(rhadapt%rVertexCoordinates2D,&
                                     rhadaptBackup%rVertexCoordinates2D)
      case(NDIM3D)
        call otree_duplicateOctree(rhadapt%rVertexCoordinates3D,&
                                   rhadaptBackup%rVertexCoordinates3D)
      case DEFAULT
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_duplicateAdaptation')
        call sys_halt()
      end select
    end if
        
    ! Bit   7: rBoundary
    if (iand(idupFlag, HADAPT_SHARE_RBOUNDARY) .ne.&
                       HADAPT_SHARE_RBOUNDARY) then
      
      if (associated(rhadaptBackup%rBoundary)) then
        do ibct = 1, size(rhadaptBackup%rBoundary,1)
          call btree_releaseTree(rhadaptBackup%rBoundary(ibct))
        end do
        deallocate(rhadaptBackup%rBoundary)
      end if

      allocate(rhadaptBackup%rBoundary(size(rhadapt%rBoundary,1)))
      do ibct = 1, size(rhadapt%rBoundary,1)
        call btree_duplicateTree(rhadapt%rBoundary(ibct),&
                                 rhadaptBackup%rBoundary(ibct))
      end do
    end if
    
    ! Bit   8: rElementsAtVertex
    if (iand(idupFlag, HADAPT_SHARE_RELEMENTSATVERTEX) .ne.&
                       HADAPT_SHARE_RELEMENTSATVERTEX) then
      call arrlst_duplicateArrayList(rhadapt%rElementsAtVertex,&
                                     rhadaptBackup%rElementsAtVertex)
    end if

  contains
    
    subroutine checkAndCopy (idupFlag, ibitfield, isourcehandle, idesthandle)
      
      ! Checks if idupFlag has all bits ibitfield set.
      ! If yes, idesthandle is set to isourcehandle.
      ! Otherwise, the memory behind isourcehandle is duplicated in memory
      ! and idesthandle receives the handle to the new memory block.
      
      integer(I32), intent(IN) :: ibitfield
      integer(I32), intent(IN) :: idupFlag
      integer, intent(IN) :: isourcehandle
      integer, intent(INOUT) :: idesthandle
      
      if (iand(idupFlag, ibitfield) .ne. ibitfield) then
        if (isourcehandle .ne. ST_NOHANDLE) then
          call storage_copy(isourcehandle, idesthandle)
        end if
      else
        idesthandle = isourcehandle
      end if
      
    end subroutine checkAndCopy
  end subroutine hadapt_duplicateAdaptation

  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_restoreAdaptation(rhadaptBackup, rhadapt)

!<description>
    ! This subroutine restares data of an adaptivity structure. All information
    ! arrays which are not shared between rhadaptBackup and another adaptivity
    ! structure is copied (back) to rhadapt.
!</description>

!<input>
    ! Backup of an adaptivity structure
    type(t_hadapt), intent(IN) :: rhadaptBackup
!</input>

!<inputoutput>
    ! Destination adaptivity structure
    ! All components where a duplicates exist in rhadaptBackup are copied
    ! to rhadapt, overwriting the old information arrays.
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: idupFlag
    integer :: ibct
  
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
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMARKER,&
                      rhadapt%h_Imarker, rhadaptBackup%h_Imarker)
   
    ! Bit   1: IvertexAge
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTEXAGE,&
                      rhadapt%h_IvertexAge, rhadaptBackup%h_IvertexAge)
    call storage_getbase_int(rhadapt%h_IvertexAge,&
                             rhadapt%p_IvertexAge)

    ! Bit   2: InodalProperty
    call checkAndCopy(idupFlag, HADAPT_SHARE_INODALPROPERTY,&
                      rhadapt%h_InodalProperty, rhadaptBackup%h_InodalProperty)
    call storage_getbase_int(rhadapt%h_InodalProperty,&
                             rhadapt%p_InodalProperty)

    ! Bit   3: IverticesAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IVERTICESATELEMENT,&
                      rhadapt%h_IverticesAtElement,&
                      rhadaptBackup%h_IverticesAtElement)
    call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                               rhadapt%p_IverticesAtElement)

    ! Bit   4: IneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_INEIGHATELEMENT,&
                      rhadapt%h_IneighboursAtElement,&
                      rhadaptBackup%h_IneighboursAtElement)
    call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                               rhadapt%p_IneighboursAtElement)

    ! Bit   5: ImidneighboursAtElement
    call checkAndCopy(idupFlag, HADAPT_SHARE_IMIDNEIGHATELEMENT,&
                      rhadapt%h_ImidneighboursAtElement,&
                      rhadaptBackup%h_ImidneighboursAtElement)
    call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                               rhadapt%p_ImidneighboursAtElement)

    ! Bit   6: rVertexCoordinates
    if (iand(idupFlag, HADAPT_SHARE_RVERTEXCOORDINATES) .ne.&
                       HADAPT_SHARE_RVERTEXCOORDINATES) then
      
      select case(rhadaptBackup%ndim)
      case(NDIM2D)
        call qtree_restoreQuadtree(rhadaptBackup%rVertexCoordinates2D,&
                                   rhadapt%rVertexCoordinates2D)
      case(NDIM3D)
        call otree_restoreOctree(rhadaptBackup%rVertexCoordinates3D,&
                                 rhadapt%rVertexCoordinates3D)
      case DEFAULT
        call output_line('Invalid spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_restoreAdaptation')
        call sys_halt()
      end select
    end if

    ! Bit   7: rBoundary
    if (iand(idupFlag, HADAPT_SHARE_RBOUNDARY) .ne.&
                       HADAPT_SHARE_RBOUNDARY) then
      if (size(rhadaptBackup%rBoundary,1) .ne. &
          size(rhadapt%rBoundary,1)) then
        call output_line('Invalid number of boundary components!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_restoreAdaptation')
        call sys_halt()
      end if
      do ibct = 1, size(rhadapt%rBoundary,1)
        call btree_restoreTree(rhadaptBackup%rBoundary(ibct),&
                               rhadapt%rBoundary(ibct))
      end do
    end if
    

    ! Bit   8: rElementsAtVertex
    if (iand(idupFlag, HADAPT_SHARE_RELEMENTSATVERTEX) .ne.&
                       HADAPT_SHARE_RELEMENTSATVERTEX) then
      call arrlst_restoreArrayList(rhadaptBackup%rElementsAtVertex,&
                                   rhadapt%rElementsAtVertex)
    end if
  
  contains
    
    subroutine checkAndCopy (idupFlag, ibitfield, idesthandle, isourcehandle)
      
      ! Checks if idupFlag has all bits ibitfield set.
      ! If not, the memory behind isourcehandle is copied to idesthandle
      ! overwriting all previous information.
      
      integer(I32), intent(IN) :: ibitfield
      integer(I32), intent(IN) :: idupFlag
      integer, intent(IN) :: isourcehandle
      integer, intent(INOUT) :: idesthandle
      
      if (iand(idupFlag, ibitfield) .ne. ibitfield) then
        if (isourcehandle .ne. ST_NOHANDLE) then
          call storage_copy(isourcehandle, idesthandle)
        end if
      end if
      
    end subroutine checkAndCopy
  end subroutine hadapt_restoreAdaptation

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setVertexCoords1D(rhadapt, h_DvertexCoords, nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the two-dimensional array 
    ! p_DvertexCoords and stores them in the binary-search-tree.
!</description>

!<input>
    ! Handle to the vertex coordinates
    integer, intent(IN) :: h_DvertexCoords

    ! Number of vertices
    integer, intent(IN) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Check if handle is not empty
    if (h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords1D')
      call sys_halt()
    end if

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .eq. HADAPT_HAS_COORDS) then
      call qtree_releaseQuadtree(rhadapt%rVertexCoordinates2D)
    end if
    
    ! Set pointer
    call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords, nvt)
    
    ! Create quadtree for vertices
    call btree_createTree(rhadapt%rVertexCoordinates1D, 2*nvt,&
                          ST_DOUBLE, 0, 0, 0, 1.5_DP)
    
    ! Copy vertex coordinates to quadtree
    call btree_copyToTree(p_DvertexCoords(1,:), rhadapt%rVertexCoordinates1D)

    ! Set specifier for quadtree
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM2D
    rhadapt%NVT  = nvt
  end subroutine hadapt_setVertexCoords1D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVertexCoords1D(rhadapt, h_DvertexCoords, nvt, ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the binary tree
    ! and stores them in the two-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will 
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    integer, intent(INOUT) :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    integer, intent(OUT), optional :: nvt

    ! OPTIONAL: number of spatial dimensions
    integer, intent(OUT), optional :: ndim
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: h_Ddata

    ! Check if coordinates exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .ne. HADAPT_HAS_COORDS) then
      call output_line('Binary tree does not exist!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords1D')
      call sys_halt()
    end if

    ! Check if coordinates are given in 1D
    if (rhadapt%ndim .ne. NDIM1D) then
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords1D')
      call sys_halt()
    end if

    ! Copy binary tree to temporal handle 
    h_Ddata = ST_NOHANDLE
    call btree_copyFromTreeKey(rhadapt%rVertexCoordinates1D, h_Ddata)
    call storage_getbase_double(h_Ddata, p_Ddata, rhadapt%NVT)
    call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords)

    ! Check if vectors are compatible and reallocate
    ! the coordinate vector if required
    if (size(p_Ddata) .gt. size(p_DvertexCoords,2)) then
      call storage_realloc('hadapt_getVertexCoords1D', size(p_Ddata),&
                           h_DvertexCoords, ST_NEWBLOCK_NOINIT, .false.)
      call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords)
    end if
    
    ! Copy data to coordinate vector
    call lalg_copyVectorDble(p_Ddata, p_DvertexCoords(1,:), rhadapt%NVT)

    ! Release temporal memory
    call storage_free(h_Ddata)

    ! Set dimension
    if (present(ndim)) ndim = rhadapt%ndim
    if (present(nvt))  nvt  = rhadapt%NVT
  end subroutine hadapt_getVertexCoords1D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setVertexCoords2D(rhadapt, h_DvertexCoords, nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the two-dimensional array 
    ! p_DvertexCoords and stores them in the quadtree.
!</description>

!<input>
    ! Handle to the vertex coordinates
    integer, intent(IN) :: h_DvertexCoords

    ! Number of vertices
    integer, intent(IN) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP) :: xmin,xmax,ymin,ymax
    integer  :: nnode

    ! Check if handle is not empty
    if (h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords2D')
      call sys_halt()
    end if

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .eq. HADAPT_HAS_COORDS) then
      call qtree_releaseQuadtree(rhadapt%rVertexCoordinates2D)
    end if
    
    ! Set pointer
    call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords, nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin = minval(p_DvertexCoords(1,:))
    xmax = maxval(p_DvertexCoords(1,:))
    ymin = minval(p_DvertexCoords(2,:))
    ymax = maxval(p_DvertexCoords(2,:))
    
    ! Estimate number of initial quadrilaterals
    nnode = int(0.5_DP*nvt)

    ! Create quadtree for vertices
    call qtree_createQuadtree(rhadapt%rVertexCoordinates2D, nvt,&
                              nnode, xmin, ymin, xmax, ymax)
    
    ! Copy vertex coordinates to quadtree
    call qtree_copyToQuadtree(p_DvertexCoords, rhadapt%rVertexCoordinates2D)

    ! Set specifier for quadtree
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM2D
    rhadapt%NVT  = nvt
  end subroutine hadapt_setVertexCoords2D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVertexCoords2D(rhadapt, h_DvertexCoords, nvt, ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the quadtree
    ! and stores them in the two-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will 
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    integer, intent(INOUT) :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    integer, intent(OUT), optional :: nvt

    ! OPTIONAL: number of spatial dimensions
    integer, intent(OUT), optional :: ndim
!</output>
!</subroutine>

    ! Check if coordinates exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .ne. HADAPT_HAS_COORDS) then
      call output_line('Quadtree does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords2D')
      call sys_halt()
    end if

    ! Check if coordinates are given in 2D
    if (rhadapt%ndim .ne. NDIM2D) then
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords2D')
      call sys_halt()
    end if

    ! Copy quadtree to handle h_DvertexCoords
    call qtree_copyFromQuadtree(rhadapt%rVertexCoordinates2D, h_DvertexCoords)

    ! Set dimension
    if (present(ndim)) ndim = rhadapt%ndim
    if (present(nvt))  nvt  = rhadapt%NVT
  end subroutine hadapt_getVertexCoords2D
  
  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setVertexCoords3D(rhadapt, h_DvertexCoords, nvt)

!<description>
    ! This subroutine sets the vertex coordinates given by the handle
    ! h_DvertexCoords that points to the three-dimensional array 
    ! p_DvertexCoords and stores them in the octree.
!</description>

!<input>
    ! Handle to the vertex coordinates
    integer, intent(IN) :: h_DvertexCoords

    ! Number of vertices
    integer, intent(IN) :: nvt
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT)    :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP) :: xmin,xmax,ymin,ymax,zmin,zmax
    integer  :: nnode

    ! Check if handle is not empty
    if (h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVertexCoords3D')
      call sys_halt()
    end if

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .eq. HADAPT_HAS_COORDS) then
      call otree_releaseOctree(rhadapt%rVertexCoordinates3D)
    end if
    
    ! Set pointer
    call storage_getbase_double2D(h_DvertexCoords, p_DvertexCoords, nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin = minval(p_DvertexCoords(1,:))
    xmax = maxval(p_DvertexCoords(1,:))
    ymin = minval(p_DvertexCoords(2,:))
    ymax = maxval(p_DvertexCoords(2,:))
    zmin = minval(p_DvertexCoords(3,:))
    zmax = maxval(p_DvertexCoords(3,:))
    
    ! Estimate number of initial quadrilaterals
    nnode = int(0.5_DP*nvt)

    ! Create octree for vertices
    call otree_createOctree(rhadapt%rVertexCoordinates3D, nvt,&
                            nnode, xmin, ymin, zmin, xmax, ymax, zmax)
    
    ! Copy vertex coordinates to octree
    call otree_copyToOctree(p_DvertexCoords, rhadapt%rVertexCoordinates3D)

    ! Set specifier for quadtree
    rhadapt%iSpec=ior(rhadapt%iSpec, HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM3D
    rhadapt%NVT  = nvt
  end subroutine hadapt_setVertexCoords3D

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVertexCoords3D(rhadapt, h_DvertexCoords, nvt, ndim)

!<description>
    ! This subroutine gets the vertex coordinates from the octree
    ! and stores them in the three-dimensional array associated to the
    ! handle h_DvertexCoords. Note that the allocated memory will 
    ! be reallocated if is does not provide the correct dimensions.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Handlt to the vertex coordinate vector
    integer, intent(INOUT) :: h_DvertexCoords
!</inputoutput>

!<output>
    ! OPTIONAL: number of vertices
    integer, intent(OUT), optional :: nvt

    ! OPTIONAL: number of spatial dimensions
    integer, intent(OUT), optional :: ndim
!</output>
!</subroutine>

    ! Check if coordinates exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_COORDS) .ne. HADAPT_HAS_COORDS) then
      call output_line('Octree does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords3D')
      call sys_halt()
    end if

    ! Check if coordinates are given in 2D
    if (rhadapt%ndim .ne. NDIM3D) then
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVertexCoords3D')
      call sys_halt()
    end if

    ! Copy octree to handle h_DvertexCoords
    call otree_copyFromOctree(rhadapt%rVertexCoordinates3D, h_DvertexCoords)

    ! Set dimension
    if (present(ndim)) ndim = rhadapt%ndim
    if (present(nvt))  nvt  = rhadapt%NVT
  end subroutine hadapt_getVertexCoords3D

  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_setVerticesAtElement(rhadapt, h_IverticesAtElement, nel)

!<description>
    ! This routine assigns the handle to the "vertices-at-element" array
    ! to the adaptivity structure and sets up internal links.
!</description>

!<input>
    ! Handle to p_IverticesAtElement
    integer, intent(IN) :: h_IverticesAtElement

    ! Total number of elements
    integer, intent(IN) :: nel
!</input>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    if (h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setVerticesAtElement')
      call sys_halt()
    end if

    rhadapt%h_IverticesAtElement = h_IverticesAtElement
    call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                               rhadapt%p_IverticesAtElement)
    
    ! Set specifier for IverticesAtElement
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_VERTATELEM)

    ! Set dimensions
    rhadapt%NEL = nel
  end subroutine hadapt_setVerticesAtElement

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getVerticesAtElement(rhadapt, h_IverticesAtElement, nel)

!<description>
    ! This routine assigns the "vertices-at-element" array from the
    ! adaptivity structure to the given handle.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_IverticesAtElement
    integer, intent(INOUT) :: h_IverticesAtElement
!</inputoutput>

!<output>
    ! OPTIONAL: number of elements
    integer, intent(OUT), optional :: nel
!</output>
!</subroutine>

    ! Check if "vertices-at-element" array exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_VERTATELEM) .ne. HADAPT_HAS_VERTATELEM) then
      call output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getVerticesAtElement')
      call sys_halt()
    end if

    ! Check if handle needs to be freed first
    if (h_IverticesAtElement .ne. ST_NOHANDLE .and.&
        h_IverticesAtElement .ne. rhadapt%h_IverticesAtElement)&
        call storage_free(h_IverticesAtElement)

    ! Assign handle
    h_IverticesAtElement = rhadapt%h_IverticesAtElement

    ! Set dimensions
    if(present(nel))   nel = rhadapt%NEL
  end subroutine hadapt_getVerticesAtElement
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_setNeighboursAtElement(rhadapt, h_IneighboursAtElement)

!<description>
    ! This routine assigns the handle to the "neighbours-at-element" array
    ! to the adaptivity structure and sets up internal links.
!</description>

!<input>
    ! Handle to p_IneighboursAtElement
    integer, intent(IN) :: h_IneighboursAtElement
!</input>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    if (h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setNeighboursAtElement')
      call sys_halt()
    end if

    rhadapt%h_IneighboursAtElement = h_IneighboursAtElement
    call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                               rhadapt%p_IneighboursAtElement)
    
    ! Create structure of "mid-adjacent" neighbours
    call storage_copy(rhadapt%h_IneighboursAtElement,&
                      rhadapt%h_ImidneighboursAtElement)
    call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                               rhadapt%p_ImidneighboursAtElement)
    
    ! Set specifier for IverticesAtElement
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_NEIGHATELEM)
  end subroutine hadapt_setNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getNeighboursAtElement(rhadapt,&
                                           h_IneighboursAtElement, nel)

!<description>
    ! This routine assigns the "neighbours-at-element" array from the
    ! adaptivity structure to the given handle
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_IneighboursAtElement
    integer, intent(INOUT) :: h_IneighboursAtElement
!</inputoutput>

!<output>
    ! OPTIONAL: number of elements
    integer, intent(OUT), optional :: nel
!</output>
!</subroutine>

    ! Check if "neighbours-at-element" array exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_NEIGHATELEM) .ne. HADAPT_HAS_NEIGHATELEM) then
      call output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNeighboursAtElement')
      call sys_halt()
    end if

    ! Check if handle needs to be freed first
    if (h_IneighboursAtElement .ne. ST_NOHANDLE .and.&
        h_IneighboursAtElement .ne. rhadapt%h_IneighboursAtElement)&
        call storage_free(h_IneighboursAtElement)

    ! Assign handle
    h_IneighboursAtElement = rhadapt%h_IneighboursAtElement

    ! Set dimensions
    if(present(nel))   nel = rhadapt%NEL
  end subroutine hadapt_getNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setNelOfType(rhadapt, InelOfType)

!<description>
    ! This subroutine sets the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</description>

!<input>
    ! Number of elements with a defined number of vertices per element.
    integer, dimension(TRIA_MAXNVE), intent(IN) :: InelOfType
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    rhadapt%InelOfType = InelOfType

    ! Set specifier for InelOfType
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_NELOFTYPE)
  end subroutine hadapt_setNelOfType

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getNelOfType(rhadapt, InelOfType)

!<description>
    ! This subroutine returns the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<output>
    ! Number of elements with a defined number of vertices per element.
    integer, dimension(TRIA_MAXNVE), intent(OUT) :: InelOfType
!</output>
!</subroutine>

    ! Check if "InelOfType" array exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_NELOFTYPE) .ne. HADAPT_HAS_NELOFTYPE) then
      call output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNelOfType')
      call sys_halt()
    end if
    
    InelOfType = rhadapt%InelOfType
  end subroutine hadapt_getNelOfType

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setBoundary(rhadapt, h_IboundaryCpIdx, h_IverticesAtBoundary,&
                                h_DvertexParameterValue, nbct, nvbd)

!<description>
    ! This subroutine sets the boundary structure to the
    ! adaptivity structure and initializes internal data
!</description>

!<input>
    ! Handle to p_IboundaryCpIdx
    integer, intent(IN) :: h_IboundaryCpIdx

    ! Handle to p_IverticesAtBoundary
    integer, intent(IN) :: h_IverticesAtBoundary

    ! Handle to p_DvertexParameterValue
    integer, intent(IN) :: h_DvertexParameterValue

    ! Number of boundary components
    integer, intent(IN) :: nbct

    ! Number of vertices at the boundary 
    integer, intent(IN) :: nvbd
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexParameterValue2D
    real(DP), dimension(:), pointer   :: p_DvertexParameterValue
    integer, dimension(:), pointer    :: p_IboundaryCpIdx
    integer, dimension(:), pointer    :: p_IverticesAtBoundary
    integer, dimension(:,:), pointer  :: p_IneighboursAtBoundary
    integer(I32), dimension(2) :: Isize
    integer :: ioff,lvbd,ivbdStart,ivbdEnd,ibct,h_IneighboursAtBoundary
    
    ! Check if handle are not empty
    if ((h_IboundaryCpIdx .eq. ST_NOHANDLE) .or.&
        (h_IverticesAtBoundary .eq. ST_NOHANDLE) .or.&
        (h_DvertexParameterValue .eq. ST_NOHANDLE)) then
      call output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setBoundary')
      call sys_halt()
    end if

    ! Check if boundary structure exists and remove it
    if (associated(rhadapt%rBoundary)) then
      do ibct = 1, size(rhadapt%rBoundary,1)
        call btree_releaseTree(rhadapt%rBoundary(ibct))
      end do
      deallocate(rhadapt%rBoundary)
    end if

    ! Set pointers
    call storage_getbase_double(h_DvertexParameterValue, p_DvertexParameterValue, nvbd)
    call convert_pointer(nvbd,p_DvertexParameterValue, p_DvertexParameterValue2D)
    call storage_getbase_int(h_IboundaryCpIdx, p_IboundaryCpIdx, nbct+1)
    call storage_getbase_int(h_IverticesAtBoundary, p_IverticesAtBoundary, nvbd)
    
    ! Allocate array of search trees
    allocate(rhadapt%rBoundary(nbct))
    
    ! Create auxiliary array
    Isize = (/2,nvbd/)
    call storage_new('hadapt_setBoundary', 'IData', Isize, ST_INT,&
                     h_IneighboursAtBoundary, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int2D(h_IneighboursAtBoundary, p_IneighboursAtBoundary)

    ! Initialization
    ioff = 0

    do ibct = 1, nbct

      ! Create a separate search tree for each boundary component
      lvbd = p_IboundaryCpIdx(ibct+1)-p_IboundaryCpIdx(ibct)
      call btree_createTree(rhadapt%rBoundary(ibct), lvbd, ST_INT, 1, 0, 2)

      ! Set subdimensions
      ivbdStart = ioff+1
      ivbdEnd   = ioff+lvbd
      ioff      = ioff+lvbd

      ! Fill auxiliary working array p_IneighboursAtBoundary
      ! For each item ivbd we have:
      !   p_IneighboursAtBoundary(BdrPrev,ivbd) -> previous neighbour of ivbd
      !   p_IneighboursAtBoundary(BdrNext,ivbd) -> following neighbour of ivbd
      p_IneighboursAtBoundary(BdrPrev,ivbdStart)           = p_IverticesAtBoundary(ivbdEnd)
      p_IneighboursAtBoundary(BdrPrev,ivbdStart+1:ivbdEnd) = p_IverticesAtBoundary(ivbdStart:ivbdEnd-1)
      p_IneighboursAtBoundary(BdrNext,ivbdStart:ivbdEnd-1) = p_IverticesAtBoundary(ivbdStart+1:ivbdEnd)
      p_IneighboursAtBoundary(BdrPrev,ivbdEnd)             = p_IverticesAtBoundary(ivbdStart)

      ! Fill search tree
      call btree_copyToTree(p_IverticesAtBoundary(ivbdStart:ivbdEnd),&
                            rhadapt%rBoundary(ibct),&
                            p_DData=p_DvertexParameterValue2D(:,ivbdStart:ivbdEnd),&
                            p_IData=p_IneighboursAtBoundary(:,ivbdStart:ivbdEnd))
    end do
    
    ! Set dimensions
    rhadapt%NBCT  = nbct
    rhadapt%NVBD  = nvbd
    rhadapt%NVBD0 = nvbd

    ! Free auxiliary storage
    call storage_free(h_IneighboursAtBoundary)

    ! Set specifier for boundary
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_BOUNDARY)

  contains
    
    ! ****************************************
    ! The convertion routine 1D -> 2D
    
    subroutine convert_pointer(isize, ptr_1d, ptr_2d)
      integer(I32), intent(IN)                 :: isize
      real(DP), dimension(1,isize), intent(IN), target :: ptr_1d
      real(DP), dimension(:,:), pointer                :: ptr_2d
      
      ptr_2d => ptr_1d
    end subroutine convert_pointer
  end subroutine hadapt_setBoundary
  
  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getBoundary(rhadapt, h_IboundaryCpIdx, h_IverticesAtBoundary,&
                                h_DvertexParameterValue, nvbd, nbct)

!<description>
    ! This subroutine extracts the boundary data from the adaptivity structure
    ! and generates the the arrays p_IboundaryCpIdx, p_IverticesAtBoundary and 
    ! p_DvertexParameterValue. If the optional parameter nvbd is given
    ! then the number of boundary vertices is assigned.
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Handle to p_IboundaryCpIdx
    integer, intent(INOUT) :: h_IboundaryCpIdx
    
    ! Handle to p_IverticesAtBoundary
    integer, intent(INOUT) :: h_IverticesAtBoundary

    ! Handle to p_DvertexParameterValue
    integer, intent(INOUT) :: h_DvertexParameterValue

    ! OPTIONAL: number of vertices at the boundary
    integer, intent(OUT), optional :: nvbd

    ! OPTIONAL: number of boundary components
    integer, intent(OUT), optional :: nbct
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer  :: p_IboundaryCpIdx
    integer, dimension(:), pointer  :: p_IverticesAtBoundary
    integer(I32) :: isize
    integer :: ioff,lvbd,ivbdStart,ivbdEnd,ibct

    ! Check if boundary data exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_BOUNDARY) .ne. HADAPT_HAS_BOUNDARY) then
      call output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getBoundary')
      call sys_halt()
    end if

    ! Due to the fact, that the boundary data may be stored in multiple 
    ! binary search trees (one for each component) the static arrays 
    ! are pre-allocated and filled step-by-step from the dynamic data
    ! in the boundary search tree(s).

    ! Check if the arrays p_IboundaryCpIdx, p_IverticesAtBoundary and 
    ! p_DvertexParameterValue have correct dimension
    if (h_IboundaryCpIdx .eq. ST_NOHANDLE) then
      call storage_new('hadapt_getBoundary', 'p_IboundaryCpIdx',&
                       rhadapt%NBCT+1, ST_INT, h_IboundaryCpIdx,&
                       ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_IboundaryCpIdx, isize)
      if (isize < rhadapt%NBCT+1) then
        call storage_realloc('hadapt_getBoundary', rhadapt%NBCT+1,&
                             h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    if (h_IverticesAtBoundary .eq. ST_NOHANDLE) then
      call storage_new('hadapt_getBoundary', 'p_IverticesAtBoundary',&
                       rhadapt%NVBD, ST_INT, h_IverticesAtBoundary,&
                       ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_IverticesAtBoundary, isize)
      if (isize < rhadapt%NVBD) then
        call storage_realloc('hadapt_getBoundary', rhadapt%NVBD,&
                             h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    if (h_DvertexParameterValue .eq. ST_NOHANDLE) then
      call storage_new('hadapt_getBoundary', 'p_DvertexParameterValue',&
                       rhadapt%NVBD, ST_DOUBLE, h_DvertexParameterValue,&
                       ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_DvertexParameterValue, isize)
      if (isize < rhadapt%NVBD) then
        call storage_realloc('hadapt_getBoundary', rhadapt%NVBD,&
                             h_DvertexParameterValue, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! Set pointers
    call storage_getbase_int(h_IboundaryCpIdx,&
                             p_IboundaryCpIdx, rhadapt%NBCT+1)
    call storage_getbase_int(h_IverticesAtBoundary,&
                             p_IverticesAtBoundary, rhadapt%NVBD)
    call storage_getbase_double(h_DvertexParameterValue,&
                                p_DvertexParameterValue, rhadapt%NVBD)

    ! Initialization
    p_IboundaryCpIdx(1) = 1
    ioff = 0

    do ibct = 1, rhadapt%NBCT
      
      ! Set subdimensions
      lvbd      = rhadapt%rBoundary(ibct)%NA
      ivbdStart = ioff+1
      ivbdEnd   = ioff+lvbd
      ioff      = ioff+lvbd

      ! Set index for next boundary component
      p_IboundaryCpIdx(ibct+1) = ioff+1

      ! Get static boundary data from tree
      call btree_copyFromTreeKey(rhadapt%rBoundary(ibct),&
                                 p_IverticesAtBoundary(ivbdStart:ivbdEnd))
      call btree_copyFromTreeDble(rhadapt%rBoundary(ibct),&
                                  p_DvertexParameterValue(ivbdStart:ivbdEnd), 1)

      ! Sort array w.r.t. increasing parameter values
      call sort_dp(p_DvertexParameterValue(ivbdStart:ivbdEnd),&
                   SORT_QUICK, p_IverticesAtBoundary(ivbdStart:ivbdEnd))
    end do

    ! Set dimension
    if(present(nvbd)) nvbd = rhadapt%NVBD
    if(present(nbct)) nbct = rhadapt%NBCT
  end subroutine hadapt_getBoundary

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_setNodalProperty(rhadapt, h_InodalProperty)

!<description>
    ! This subroutine sets the nodal property list to the adaptivity structure
!</description>

!<input>
    ! Handle to p_InodalProperty
    integer, intent(IN) :: h_InodalProperty
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! Check if handle is not empty
    if (h_InodalProperty .eq. ST_NOHANDLE) then
      call output_line('Invalid handle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_setNodalProperty')
      call sys_halt()
    end if

    rhadapt%h_InodalProperty = h_InodalProperty
    call storage_getbase_int(rhadapt%h_InodalProperty,&
                             rhadapt%p_InodalProperty)
    
    ! Set specifier for InodalProperty
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_NODALPROP)
  end subroutine hadapt_setNodalProperty

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_getNodalProperty(rhadapt, h_InodalProperty)

!<description>
    ! This subroutine assignes the "nodal property" array from the
    ! adaptivity structure to the given handle
!</description>

!<input>
    ! Adaptivity structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>

!<inputoutput>
    ! Hande to p_InodalProperty
    integer, intent(INOUT) :: h_InodalProperty
!</inputoutput>
!</subroutine>

    ! Check if "nodal property" list exits
    if (iand(rhadapt%iSpec, HADAPT_HAS_NODALPROP) .ne. HADAPT_HAS_NODALPROP) then
      call output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_getNodalProperty')
      call sys_halt()
    end if

    ! Check if handle needs to be freed first
    if (h_InodalProperty .ne. ST_NOHANDLE .and.&
        h_InodalProperty .ne. rhadapt%h_InodalProperty)&
        call storage_free(h_InodalProperty)

    ! Assign handle
    h_InodalProperty = rhadapt%h_InodalProperty
  end subroutine hadapt_getNodalProperty

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_genElementsAtVertex(rhadapt)

!<description>
    ! This subroutine generates a list of arrays for the "elements-at-vertex"
    ! structure. This is an internal structure of the adaptivity structure
    ! which is used, e.g., in the removel of vertices.
!</description>

!<inputoutput>
    ! Adaptivity structur
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    integer :: ipos,iel,ive

    ! Check if "vertices-at-element" list exists
    if (iand(rhadapt%iSpec, HADAPT_HAS_VERTATELEM) .ne. HADAPT_HAS_VERTATELEM) then
      call output_line('Structure does not exist!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_genElementsAtVertex')
      call sys_halt()
    end if

    ! Create arraylist
    call arrlst_createArraylist(rhadapt%relementsAtVertex,&
                                2*rhadapt%NVT, 8*rhadapt%NEL,&
                                ST_INT, ARRAYLIST_UNORDERED)

    ! Fill arraylist
    do iel = 1, rhadapt%NEL
      do ive = 1, hadapt_getNVE(rhadapt,iel)
        call arrlst_appendToArraylist(rhadapt%relementsAtVertex,&
                                      rhadapt%p_IverticesAtElement(ive,iel),&
                                      iel, ipos)
      end do
    end do
    
    ! Set specifier for relementsAtVertex
    rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_ELEMATVERTEX)
  end subroutine hadapt_genElementsAtVertex

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_performAdaptation(rhadapt, rindicator,&
                                      rcollection, fcb_hadaptCallback)

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
    type(t_vectorScalar), intent(IN)            :: rindicator

    ! callback routines
    include 'intf_hadaptcallback.inc'
    optional                                    :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(1) :: Ielements, Ivertices
    integer(I32), dimension(2) :: Isize
    integer :: nvt,nel

    ! Check if dynamic data structures are available
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA) .ne. HADAPT_HAS_DYNAMICDATA) then
      call output_line('Dynamic data structures are not generated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
      call sys_halt()
    end if

    ! Initialize initial dimensions
    call storage_getsize(rhadapt%h_IverticesAtElement, Isize)
    rhadapt%NELMAX = Isize(2)
    call storage_getsize(rhadapt%h_IneighboursAtElement, Isize)
    rhadapt%NELMAX = min(rhadapt%NELMAX,Isize(2))

    rhadapt%InelOfType0 = rhadapt%InelOfType
    rhadapt%NVT0        = rhadapt%NVT
    rhadapt%NEL0        = rhadapt%NEL
    rhadapt%NVBD0       = rhadapt%NVBD
    rhadapt%increaseNVT = 0
    
    ! What kind of grid refinement should be performed
    select case(rhadapt%iadaptationStrategy)
    case (HADAPT_NOADAPTATION)   ! No grid refinement

      
    case (HADAPT_REDGREEN)       ! Red-green grid refinement

      ! Which spatial dimensions are we?
      select case(rhadapt%ndim)
        
      case (NDIM2D)
        ! Mark elements for refinement based on indicator function
        call mark_refinement2D(rhadapt, rindicator)
        
        ! Mark additional elements to restore conformity
        call redgreen_mark_refinement2D(rhadapt, rcollection, fcb_hadaptCallback)
        
        ! Mark element for recoarsening based on indicator function
        call redgreen_mark_coarsening2D(rhadapt, rindicator)
        
      case DEFAULT
        call output_line('Unsupported spatial dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
        call sys_halt()
      end select

      ! Compute new dimensions
      nvt = rhadapt%NVT+rhadapt%increaseNVT
      nel = NumberOfElements(rhadapt)

      ! Adjust vertex age array and nodal property array
      call storage_realloc('hadapt_performAdaptation', nvt,&
                           rhadapt%h_IvertexAge, ST_NEWBLOCK_NOINIT, .true.)
      call storage_realloc('hadapt_performAdaptation', nvt,&
                           rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)
     
      ! Adjust elemental arrays
      call storage_realloc('hadapt_performAdaptation', nel,&
                           rhadapt%h_IverticesAtElement, ST_NEWBLOCK_NOINIT, .true.)
      call storage_realloc('hadapt_performAdaptation', nel,&
                           rhadapt%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)
      call storage_realloc('hadapt_performAdaptation', nel,&
                           rhadapt%h_ImidneighboursAtElement, ST_NEWBLOCK_NOINIT, .true.)

      ! Reset pointers
      call storage_getbase_int(rhadapt%h_IvertexAge,&
                               rhadapt%p_IvertexAge)
      call storage_getbase_int(rhadapt%h_InodalProperty,&
                               rhadapt%p_InodalProperty)
      call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                                 rhadapt%p_IverticesAtElement)
      call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                                 rhadapt%p_IneighboursAtElement)
      call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                                 rhadapt%p_ImidneighboursAtElement)

      ! Adjust dimension of solution vector
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        Ivertices = (/nvt/)
        Ielements = (/0/)
        call fcb_hadaptCallback(rcollection, HADAPT_OPR_ADJUSTVERTEXDIM,&
                                Ivertices, Ielements)
      end if

      ! Perform refinement
      call redgreen_refine(rhadapt, rcollection, fcb_hadaptCallback)

      ! Perform coarsening
      call redgreen_coarsen(rhadapt, rcollection, fcb_hadaptCallback)

      ! Adjust nodal property array
      call storage_realloc('hadapt_performAdaptation', rhadapt%NVT,&
                           rhadapt%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)
            
    case DEFAULT
      call output_line('Unsupported refinement strategy!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_performAdaptation')
      call sys_halt()
    end select

  contains

    ! Here, some auxiliary routines follow
    
    !*************************************************************
    ! This function computes the number of elements after refinement

    function NumberOfElements(rhadapt) result(nel)
      
      type(t_hadapt), intent(IN) :: rhadapt
      integer :: nel
      
      ! local variables
      integer, dimension(:), pointer :: p_Imarker
      integer :: iel
      
      call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)
      
      ! Initialize number of elements by current number
      nel = rhadapt%NEL0
      
      ! Loop over all elements and check marker
      do iel = 1, rhadapt%NEL0
        
        select case(p_Imarker(iel))
        case(MARK_REF_TRIA2TRIA_1,&
             MARK_REF_TRIA2TRIA_2,&
             MARK_REF_TRIA2TRIA_3,&
             MARK_REF_QUAD2QUAD_13,&
             MARK_REF_QUAD2QUAD_24,&
             MARK_CRS_2QUAD3TRIA)
          ! Interestingly enought, the coarsening of 2 quadrilaterals into three
          ! triangles which reduces the number of vertices by one also increases
          ! the number of elements by one which has to be taken into account here.
          nel = nel+1

        case(MARK_REF_QUAD3TRIA_1,&
             MARK_REF_QUAD3TRIA_2,&
             MARK_REF_QUAD3TRIA_3,&
             MARK_REF_QUAD3TRIA_4,&
             MARK_REF_TRIA3TRIA_12,&
             MARK_REF_TRIA3TRIA_23,&
             MARK_REF_TRIA3TRIA_13)
          nel = nel+2
          
        case(MARK_REF_TRIA4TRIA,&
            MARK_REF_QUAD4QUAD,&
            MARK_REF_QUAD4TRIA_12,&
             MARK_REF_QUAD4TRIA_23,&
             MARK_REF_QUAD4TRIA_34,&
             MARK_REF_QUAD4TRIA_14)
          nel = nel+3
        end select
      end do
    end function NumberOfElements
  end subroutine hadapt_performAdaptation

  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_infoStatistics(rhadapt)

!<description>
    ! This subroutine outputs statistical info about the adaptivity data structure
!</description>

!<input>
    ! adaptivity data structure
    type(t_hadapt), intent(IN) :: rhadapt
!</input>
!</subroutine>

    ! local variables
    integer :: ibct
    
    call output_line('Adaptivity statistics:')
    call output_line('----------------------')
    call output_line('Total number of grid refinement steps:           '//&
        trim(sys_siL(rhadapt%nRefinementSteps,3)))
    call output_line('Total number of grid coarsening steps:           '//&
        trim(sys_siL(rhadapt%nCoarseningSteps,3)))
    call output_line('Total number of grid smoothing  steps:           '//&
        trim(sys_siL(rhadapt%nSmoothingSteps,3)))
    call output_line('Strategy for grid refinement and coarsening:     '//&
        trim(sys_siL(rhadapt%iadaptationStrategy,3)))
    call output_line('Total number of elements (initially):            '//&
        trim(sys_siL(rhadapt%NEL,15))//"("//trim(sys_siL(rhadapt%NEL0,15))//")")
    call output_line('Total number of vertices (initially):            '//&
        trim(sys_siL(rhadapt%NVT,15))//"("//trim(sys_siL(rhadapt%NVT0,15))//")")
    call output_line('Total number of vertices at boundary(initially): '//&
        trim(sys_siL(rhadapt%NVBD,15))//"("//trim(sys_siL(rhadapt%NVBD0,15))//")")
    call output_lbrk

    call output_line('Handles:')
    call output_line('--------')
    call output_line('h_Imarker:                 '//&
        trim(sys_siL(rhadapt%h_Imarker,15)))
    call output_line('h_IvertexAge:              '//&
        trim(sys_siL(rhadapt%h_IvertexAge,15)))
    call output_line('h_InodalProperty:          '//&
        trim(sys_siL(rhadapt%h_InodalProperty,15)))
    call output_line('h_IverticesAtElement:      '//&
        trim(sys_siL(rhadapt%h_IverticesAtElement,15)))
    call output_line('h_IneighboursAtElement:    '//&
        trim(sys_siL(rhadapt%h_IneighboursAtElement,15)))
    call output_line('h_ImidneighboursAtElement: '//&
        trim(sys_siL(rhadapt%h_ImidneighboursAtElement,15)))
    call output_lbrk

    call output_line('Coordinates:')
    call output_line('------------')
    select case(rhadapt%ndim)
    case(NDIM2D)
      call qtree_infoQuadtree(rhadapt%rVertexCoordinates2D)
    case(NDIM3D)
      call otree_infoOctree(rhadapt%rVertexCoordinates3D)
    end select
    call output_lbrk
    
    call output_line('Boundary:')
    call output_line('---------')
    do ibct=1,rhadapt%NBCT
      call output_line('Boundary component '//trim(sys_siL(ibct,3))//":")
      call output_line('------------------------')
      call btree_infoTree(rhadapt%rBoundary(ibct))
    end do
    call output_lbrk

    call output_line('Element at vertex structure:')
    call output_line('----------------------------')
    call arrlst_infoArraylist(rhadapt%relementsAtVertex)
  end subroutine hadapt_infoStatistics

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_writeGridSVG(rhadapt, coutputFile, width, height)

!<description>
    ! This subroutine outputs the current state of the adapted grid stored
    ! in the dynamic data structure to a given file in SVG format.
!</description>

!<input>
    ! Output file name w/o suffix .svg
    character(LEN=*), intent(IN)     :: coutputFile

    ! OPTIONAL: Width of the generated file
    integer, intent(IN), optional :: width

    ! OPTIONAL: Heigh of the generated file
    integer, intent(IN), optional :: height
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local parameters
    integer, parameter :: defaultwidth  = 1280
    integer, parameter :: defaultheight = 1024
    integer, parameter :: xoffset       = 100
    integer, parameter :: yoffset       = 100
    integer, parameter :: font_size     = 18
    integer, parameter :: colsep        = 10
    integer, parameter :: linesep       = 32
    character(LEN=*), parameter :: font_family = "Arial"
    
    ! local variables
    real(DP), dimension(2*NDIM2D) :: bbox
    real(DP) :: x0,y0,xdim,ydim,dscale
    
    integer, dimension(:), pointer :: p_Imarker
    integer :: ivt,iel,ipos,xsize,ysize,iunit,ive,nve
    integer, save :: iout=0
    
    ! Check if dynamic data structures generated
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA) .ne. HADAPT_HAS_DYNAMICDATA) then
      call output_line('Dynamic data structures are not generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridSVG')
      call sys_halt()
    end if
    
    ! Increment the sample number
    iout=iout+1
    
    ! Open output file for writing
    call io_openFileForWriting(trim(adjustl(coutputFile))//'.'//&
        trim(sys_siL(iout,5))//'.svg', iunit, SYS_REPLACE, bformatted=.true.)
    
    ! Write prolog, XML-declaration, document type declaration)
    write(iunit,FMT='(A)') '<?xml version="1.0" encoding="utf-8" standalone="yes"?>'
    write(iunit,FMT='(A)') '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN"'
    write(iunit,FMT='(A)') ' "http://www.w3.org/Graphics/SVG/1.0/DTD/svg11.dtd">'
    write(iunit,FMT='(A)')
    write(iunit,FMT='(A)') '<!-- Created with Featflow2 (http://www.featflow.de/) -->'
    write(iunit,FMT='(A)')
    write(iunit,FMT='(A)') '<svg version="1.0"'
    write(iunit,FMT='(A)') ' xmlns="http://www.w3.org/2000/svg"'
    write(iunit,FMT='(A)') ' xmlns:xlink="http://www.w3.org/1999/xlink"'
    
    if (present(width)) then
      xsize=width
    else
      xsize=defaultwidth
    end if

    if (present(height)) then
      ysize=height
    else
      ysize=defaultheight
    end if

    ! Determine bounding box
    bbox   = qtree_getBoundingBox(rhadapt%rVertexCoordinates2D)
    xdim   = bbox(3)-bbox(1); x0=bbox(1)
    ydim   = bbox(4)-bbox(2); y0=bbox(2)
    dscale = min(xsize/xdim,ysize/ydim)

    ! Set height and width of image
    write(iunit,FMT='(A)') ' width="100%" height="100%" xml:space="preserve"'
    write(iunit,FMT='(A)') ' viewBox="'//&
        trim(sys_siL(-xoffset,10))//' '//&
        trim(sys_siL(-yoffset,10))//' '//&
        trim(sys_siL(2*xoffset+xsize,10))//' '//&
        trim(sys_siL(2*yoffset+ysize,10))//'">'

    ! Write embedded java-scripts to file
    write(iunit,FMT='(A)') '<defs>'
    write(iunit,FMT='(A)') '<script type="text/javascript">'
    write(iunit,FMT='(A)') '<![CDATA['

    write(iunit,FMT='(A)') '  function ShowVertexInfo(evt,ivt,age,elems) {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var x=evt.clientX; var y=evt.clientY;'
    write(iunit,FMT='(A)') '    var m=svgdoc.getScreenCTM();'
    write(iunit,FMT='(A)') '    var p=svgdoc.createSVGPoint();'
    write(iunit,FMT='(A)') '    p.x=evt.clientX; p.y=evt.clientY;'
    write(iunit,FMT='(A)') '    p=p.matrixTransform(m.inverse());'
    write(iunit,FMT='(A)') '    p.x=p.x+15; p.y=p.y+15;'

    write(iunit,FMT='(A)') '    var vinfo=svgdoc.getElementById("vinfo");'
    write(iunit,FMT='(A)') '    var vinfo_box=svgdoc.getElementById("vinfo_box");'
    write(iunit,FMT='(A)') '    var vinfo_text1=svgdoc.getElementById("vinfo_text1");'
    write(iunit,FMT='(A)') '    var vinfo_text2=svgdoc.getElementById("vinfo_text2");'
    write(iunit,FMT='(A)') '    var vinfo_text3=svgdoc.getElementById("vinfo_text3");'
    write(iunit,FMT='(A)') '    var vinfo_text4=svgdoc.getElementById("vinfo_text4");'
    write(iunit,FMT='(A)') '    var vinfo_text5=svgdoc.getElementById("vinfo_text5");'
    write(iunit,FMT='(A)') '    var vinfo_text6=svgdoc.getElementById("vinfo_text6");'
    
    write(iunit,FMT='(A)') '    vinfo_box.setAttribute("x",p.x);'
    write(iunit,FMT='(A)') '    vinfo_box.setAttribute("y",p.y);'
    write(iunit,FMT='(A)') '    vinfo_text1.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text1.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text2.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text2.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text3.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text3.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text4.setAttribute("x",p.x+vinfo_text1.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    vinfo_text4.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text4.firstChild.nodeValue = ivt;'
    write(iunit,FMT='(A)') '    vinfo_text5.setAttribute("x",p.x+vinfo_text2.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    vinfo_text5.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text5.firstChild.nodeValue = age;'
    write(iunit,FMT='(A)') '    vinfo_text6.setAttribute("x",p.x+vinfo_text3.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    vinfo_text6.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    vinfo_text6.firstChild.nodeValue = elems;'

    write(iunit,FMT='(A)') '    var textlen = vinfo_text1.getComputedTextLength()+&
        &vinfo_text4.getComputedTextLength();'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,vinfo_text2.getComputedTextLength()&
        &+vinfo_text5.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,vinfo_text3.getComputedTextLength()&
        &+vinfo_text6.getComputedTextLength());'
    write(iunit,FMT='(A)') '    vinfo_box.setAttribute("width",textlen+30);'

    write(iunit,FMT='(A)') '    vinfo.setAttribute("style","visibility: visible");'
    write(iunit,FMT='(A)') '  }'
    
    write(iunit,FMT='(A)') '  function HideVertexInfo() {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var vinfo=svgdoc.getElementById("vinfo");'
    write(iunit,FMT='(A)') '    vinfo.setAttribute("style","visibility: hidden");'
    write(iunit,FMT='(A)') '  }'

    write(iunit,FMT='(A)') '  function ShowElementInfo(evt,iel,state,marker,kadj,kmidadj,kvert) {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var x=evt.clientX; var y=evt.clientY;'
    write(iunit,FMT='(A)') '    var m=svgdoc.getScreenCTM();'
    write(iunit,FMT='(A)') '    var p=svgdoc.createSVGPoint();'
    write(iunit,FMT='(A)') '    p.x=evt.clientX; p.y=evt.clientY;'
    write(iunit,FMT='(A)') '    p=p.matrixTransform(m.inverse());'
    write(iunit,FMT='(A)') '    p.x=p.x+15; p.y=p.y+15;'

    write(iunit,FMT='(A)') '    var einfo=svgdoc.getElementById("einfo");'
    write(iunit,FMT='(A)') '    var einfo_box=svgdoc.getElementById("einfo_box");'
    write(iunit,FMT='(A)') '    var einfo_text1=svgdoc.getElementById("einfo_text1");'
    write(iunit,FMT='(A)') '    var einfo_text2=svgdoc.getElementById("einfo_text2");'
    write(iunit,FMT='(A)') '    var einfo_text3=svgdoc.getElementById("einfo_text3");'
    write(iunit,FMT='(A)') '    var einfo_text4=svgdoc.getElementById("einfo_text4");'
    write(iunit,FMT='(A)') '    var einfo_text5=svgdoc.getElementById("einfo_text5");'
    write(iunit,FMT='(A)') '    var einfo_text6=svgdoc.getElementById("einfo_text6");'
    write(iunit,FMT='(A)') '    var einfo_text7=svgdoc.getElementById("einfo_text7");'
    write(iunit,FMT='(A)') '    var einfo_text8=svgdoc.getElementById("einfo_text8");'
    write(iunit,FMT='(A)') '    var einfo_text9=svgdoc.getElementById("einfo_text9");'
    write(iunit,FMT='(A)') '    var einfo_text10=svgdoc.getElementById("einfo_text10");'
    write(iunit,FMT='(A)') '    var einfo_text11=svgdoc.getElementById("einfo_text11");'
    write(iunit,FMT='(A)') '    var einfo_text12=svgdoc.getElementById("einfo_text12");'
    
    write(iunit,FMT='(A)') '    einfo_box.setAttribute("x",p.x);'
    write(iunit,FMT='(A)') '    einfo_box.setAttribute("y",p.y);'
    write(iunit,FMT='(A)') '    einfo_text1.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text1.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text2.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text2.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text3.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text3.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text4.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text4.setAttribute("y",p.y+4*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text5.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text5.setAttribute("y",p.y+5*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text6.setAttribute("x",p.x+'//trim(sys_siL(colsep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text6.setAttribute("y",p.y+6*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text7.setAttribute("x",p.x+einfo_text1.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text7.setAttribute("y",p.y+'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text7.firstChild.nodeValue = iel;'
    write(iunit,FMT='(A)') '    einfo_text8.setAttribute("x",p.x+einfo_text2.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text8.setAttribute("y",p.y+2*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text8.firstChild.nodeValue = state;'
    write(iunit,FMT='(A)') '    einfo_text9.setAttribute("x",p.x+einfo_text3.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text9.setAttribute("y",p.y+3*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text9.firstChild.nodeValue = marker;'
    write(iunit,FMT='(A)') '    einfo_text10.setAttribute("x",p.x+einfo_text4.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text10.setAttribute("y",p.y+4*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text10.firstChild.nodeValue = kadj;'
    write(iunit,FMT='(A)') '    einfo_text11.setAttribute("x",p.x+einfo_text5.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text11.setAttribute("y",p.y+5*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text11.firstChild.nodeValue = kmidadj;'
    write(iunit,FMT='(A)') '    einfo_text12.setAttribute("x",p.x+einfo_text6.getComputedTextLength()+20);'
    write(iunit,FMT='(A)') '    einfo_text12.setAttribute("y",p.y+6*'//trim(sys_siL(linesep,10))//');'
    write(iunit,FMT='(A)') '    einfo_text12.firstChild.nodeValue = kvert;'

    write(iunit,FMT='(A)') '    var textlen = einfo_text1.getComputedTextLength()+&
        &einfo_text7.getComputedTextLength();'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text2.getComputedTextLength()+&
        &einfo_text8.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text3.getComputedTextLength()+&
        &einfo_text9.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text4.getComputedTextLength()+&
        &einfo_text10.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text5.getComputedTextLength()+&
        &einfo_text11.getComputedTextLength());'
    write(iunit,FMT='(A)') '    textlen = Math.max(textlen,einfo_text6.getComputedTextLength()+&
        &einfo_text12.getComputedTextLength());'
    write(iunit,FMT='(A)') '    einfo_box.setAttribute("width",textlen+30);'

    write(iunit,FMT='(A)') '    einfo.setAttribute("style","visibility: visible");'
    write(iunit,FMT='(A)') '  }'

    write(iunit,FMT='(A)') '  function HideElementInfo() {'
    write(iunit,FMT='(A)') '    var svgdoc=document.documentElement;'
    write(iunit,FMT='(A)') '    var einfo=svgdoc.getElementById("einfo");'
    write(iunit,FMT='(A)') '    einfo.setAttribute("style","visibility: hidden");'
    write(iunit,FMT='(A)') '  }'
    write(iunit,FMT='(A)') ']]>'
    write(iunit,FMT='(A)') '</script>'
    write(iunit,FMT='(A)') '</defs>'
       
    !---------------------------------------------------------------------------
    ! Output all elements
    !---------------------------------------------------------------------------
    if (iand(rhadapt%iSpec, HADAPT_MARKEDREFINE)  .eq. HADAPT_MARKEDREFINE .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDCOARSEN) .eq. HADAPT_MARKEDCOARSEN) then

      ! Set pointer to marker
      call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

      ! Output elements and color those which are marked for refinement
      do iel = 1, rhadapt%NEL
        
        ! Get number of vertices per elements
        nve = hadapt_getNVE(rhadapt, iel)
        
        if (iel .le. rhadapt%NEL0) then
          
          ! What kind of element are we?
          select case(p_Imarker(iel))
            
          ! Element is neither marked for refinement nor coarsening
          case(MARK_ASIS_TRIA,&
               MARK_ASIS_QUAD)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="white" stroke="black" stroke-width="1"'
            
          ! Element is marked for green refinement
          case(MARK_REF_TRIA2TRIA_1,&
               MARK_REF_TRIA2TRIA_2,&
               MARK_REF_TRIA2TRIA_3,&
               MARK_REF_QUAD3TRIA_1,&
               MARK_REF_QUAD3TRIA_2,&
               MARK_REF_QUAD3TRIA_3,&
               MARK_REF_QUAD3TRIA_4,&
               MARK_REF_QUAD4TRIA_12,&
               MARK_REF_QUAD4TRIA_23,&
               MARK_REF_QUAD4TRIA_34,&
               MARK_REF_QUAD4TRIA_14,&
               MARK_REF_QUAD2QUAD_13,&
               MARK_REF_QUAD2QUAD_24)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="green" stroke="black" stroke-width="1"'

          ! Element is marked for blue refinement
          case(MARK_REF_TRIA3TRIA_12,&
               MARK_REF_TRIA3TRIA_23,&
               MARK_REF_TRIA3TRIA_13,&
               MARK_REF_QUADBLUE_412,&
               MARK_REF_QUADBLUE_234,&
               MARK_REF_QUADBLUE_123,&
               MARK_REF_QUADBLUE_341)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="blue" stroke="black" stroke-width="1"'

          ! Element is marked for red refinement
          case(MARK_REF_TRIA4TRIA,&
               MARK_REF_QUAD4QUAD)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="red" stroke="black" stroke-width="1"'

          ! Element is marked for generic coarsening
          case(MARK_CRS_GENERIC)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="lightgray" stroke="black" stroke-width="1"'
            
          ! Element is marked for specific coarsening
          case(MARK_CRS_4QUAD4TRIA:MARK_CRS_4TRIA1TRIA)
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="yellow" stroke="black" stroke-width="1"'

          ! Unknown element marker
          case DEFAULT
            write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
                '" fill="hotpink" stroke="black" stroke-width="1"'
          end select
          
          ! Write data which is common to all elements
          select case(nve)
          case(TRIA_NVETRI2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(p_Imarker(iel),5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'

          case(TRIA_NVEQUAD2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(p_Imarker(iel),5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          end select

        else
          
          ! For all new elements there si no marker available
          write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
              '" fill="white" stroke="black" stroke-width="1"'
          
          ! Write data which is common to all elements
          select case(nve)
          case(TRIA_NVETRI2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(0,5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(0,10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'

          case(TRIA_NVEQUAD2D)
            write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
                trim(sys_siL(iel,10))//''','''//&
                trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
                trim(sys_siL(0,5))//''','''//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
                trim(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
            write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          end select
          
        end if
               
        ! Each element is a polygon made up from 3/4 points
        write(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        do ive =1 , nve
          xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-x0
          ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-y0
          write(iunit,FMT='(A)',ADVANCE='NO') &
              trim(sys_siL(int(dscale*xdim),10))//','//&
              trim(sys_siL(ysize-int(dscale*ydim),10))//' '
        end do
        
        ! And of course, the polygone must be closed !
        xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-x0
        ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-y0
        write(iunit,FMT='(A)') &
            trim(sys_siL(int(dscale*xdim),10))//','//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//'"/>'       

      end do
      
    else
      
      ! Output only elements without individual coloring
      do iel = 1, rhadapt%NEL

        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)
        
        ! For all new elements there si no marker available
        write(iunit,FMT='(A)') '<polygon id="el'//trim(sys_siL(iel,9))//&
            '" fill="white" stroke="black" stroke-width="1"'

        ! Write data which is common to all elements
        select case(nve)
        case(TRIA_NVETRI2D)
          write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
              trim(sys_siL(iel,10))//''','''//&
              trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
              trim(sys_siL(0,5))//''','''//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(0,10))//''','''//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(0,10))//''','''//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
              trim(sys_siL(0,10))//''')"'
          write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
          
        case(TRIA_NVEQUAD2D)
          write(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
              trim(sys_siL(iel,10))//''','''//&
              trim(sys_siL(redgreen_getState2D(rhadapt,iel),4))//''','''//&
              trim(sys_siL(0,5))//''','''//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IneighboursAtElement(4,iel),10))//''','''//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(3,iel),10))//','//&
              trim(sys_siL(rhadapt%p_ImidneighboursAtElement(4,iel),10))//''','''//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(1,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(2,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(3,iel),10))//','//&
              trim(sys_siL(rhadapt%p_IverticesAtElement(4,iel),10))//''')"'
          write(iunit,FMT='(A)') ' onmouseout="HideElementInfo()" style="cursor: help"'
        end select
        
        ! Each element is a polygon made up from 3/4 points
        write(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        do ive = 1, nve
          xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-x0
          ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                            rhadapt%p_IverticesAtElement(ive,iel))-y0
          write(iunit,FMT='(A)',ADVANCE='NO') &
              trim(sys_siL(int(dscale*xdim),10))//','//&
              trim(sys_siL(ysize-int(dscale*ydim),10))//' '
        end do
        
        ! And of course, the polygone must be closed !
        xdim = qtree_getX(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-x0
        ydim = qtree_getY(rhadapt%rVertexCoordinates2D,&
                          rhadapt%p_IverticesAtElement(1,iel))-y0
        write(iunit,FMT='(A)') &
            trim(sys_siL(int(dscale*xdim),10))//','//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//'"/>'       
        
      end do
    end if
         
    ! Loop over all vertices
    do ivt = 1, qtree_getsize(rhadapt%rVertexCoordinates2D)     
      xdim = qtree_getX(rhadapt%rVertexCoordinates2D,ivt)-x0
      ydim = qtree_getY(rhadapt%rVertexCoordinates2D,ivt)-y0
      
      ! Write vertices as points?
      if (rhadapt%p_IvertexAge(ivt) > 0) then
        write(iunit,FMT='(A)') '<circle id="vt'//trim(sys_siL(ivt,10))//'" cx="'//&
            trim(sys_siL(int(dscale*xdim),10))//'" cy="'//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//&
            '" r="0pt" fill="white" stroke="black" stroke-width="1pt"'
      else
        write(iunit,FMT='(A)') '<circle id="vt'//trim(sys_siL(ivt,10))//'" cx="'//&
            trim(sys_siL(int(dscale*xdim),10))//'" cy="'//&
            trim(sys_siL(ysize-int(dscale*ydim),10))//&
            '" r="0pt" fill="black" stroke="black" stroke-width="1pt"'
      end if
      
      ! Write data which is common to all vertices
      write(iunit,FMT='(A)',ADVANCE='NO') ' onmousemove="ShowVertexInfo(evt,'''//&
          trim(sys_siL(ivt,10))//''','''//&
          trim(sys_siL(rhadapt%p_IvertexAge(ivt),5))//''','''
      
      ! Generate list of elements meeting at vertex
      ipos = arrlst_getNextInArrayList(rhadapt%rElementsAtVertex, ivt, .true.)
      do while(ipos .gt. ARRLST_NULL)
        ! Get element number IEL
        iel = rhadapt%rElementsAtVertex%p_IData(ipos)
        
        ! Proceed to next entry in array list
        ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex, ivt, .false.)

        if (ipos .gt. ARRLST_NULL) then
          write(iunit,FMT='(A)',ADVANCE='NO') trim(sys_siL(iel,10))//','
        else
          write(iunit,FMT='(A)',ADVANCE='NO') trim(sys_siL(iel,10))
        end if
      end do
      write(iunit,FMT='(A)') ''')"'
      write(iunit,FMT='(A)') ' onmouseout="HideVertexInfo()" style="cursor: crosshair"/>'
    end do


    ! Output info boxes
    write(iunit,FMT='(A)') '<g id="vinfo" style="visibility: hidden">'
    write(iunit,FMT='(A)') '  <rect id="vinfo_box"  x="0" y="0" width="200" height="110"'
    write(iunit,FMT='(A)') '        style="fill: #FFFFCC; stroke: #000000; stroke-width: 0.5px;"/>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text1" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Vertex number:</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text2" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Vertex age:</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text3" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Elements meeting at vertex:</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text4" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">null</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text5" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">null</text>'
    write(iunit,FMT='(A)') '  <text id="vinfo_text6" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">null</text>'
    write(iunit,FMT='(A)') '</g>'

    write(iunit,FMT='(A)') '<g id="einfo" style="visibility: hidden">'
    write(iunit,FMT='(A)') '  <rect id="einfo_box"  x="0" y="0" width="200" height="210"'
    write(iunit,FMT='(A)') '        style="fill: #FFFFCC; stroke: #000000; stroke-width: 0.5px;"/>'
    write(iunit,FMT='(A)') '  <text id="einfo_text1" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black">Element number:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text2" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element state:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text3" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element marker:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text4" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element neighbours:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text5" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element mid-neighbours:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text6" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">Element vertices:</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text7" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text8" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text9" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text10" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text11" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '  <text id="einfo_text12" x="0" y="0" style="font-size:'//&
        trim(sys_siL(font_size,3))//'pt;font-family:'//trim(adjustl(font_family))//&
        ';stroke:black;">null</text>'
    write(iunit,FMT='(A)') '</g>'

    ! Close XML-file
    write(iunit,FMT='(A)') '</svg>'
    close(iunit)
  end subroutine hadapt_writeGridSVG

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_writeGridGMV(rhadapt, coutputFile)

!<description>
    ! This subroutine outputs the current state of the adapted grid stored
    ! in the dynamic data structure to a given file in GMV format.
!</description>

!<input>
    ! Output file name w/o suffix .gmv
    character(LEN=*), intent(IN)     :: coutputFile
!</input>

!<inputoutput>
    ! Adaptive data structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local parameters
    integer,  dimension(:), pointer :: p_Imarker
    integer :: ivt,iel,iunit,nve
    integer, save :: iout=0

    ! Check if dynamic data structures generated
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA) .ne. HADAPT_HAS_DYNAMICDATA) then
      call output_line('Dynamic data structures are not generated',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
      call sys_halt()
    end if
    
    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)

    ! Increment the sample number
    iout=iout+1

    ! Open output file for writing
    call io_openFileForWriting(trim(adjustl(coutputFile))//'.'//&
        trim(sys_siL(iout,5))//'.gmv', iunit, SYS_REPLACE, bformatted=.true.)
    write(UNIT=iunit,FMT='(A)') 'gmvinput ascii'

    ! Write vertices to output file
    select case(rhadapt%ndim)
    case (NDIM2D)
      write(UNIT=iunit,FMT=*) 'nodes ', rhadapt%NVT
      do ivt = 1, qtree_getsize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) qtree_getX(rhadapt%rVertexCoordinates2D, ivt)
      end do
      do ivt = 1, qtree_getsize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) qtree_getY(rhadapt%rVertexCoordinates2D, ivt)
      end do
      do ivt = 1, qtree_getsize(rhadapt%rVertexCoordinates2D)
        write(UNIT=iunit,FMT=10) 0._DP
      end do

    case (NDIM3D)
      write(UNIT=iunit,FMT=*) 'nodes ', rhadapt%NVT
      do ivt = 1, otree_getsize(rhadapt%rVertexCoordinates3D)
        write(UNIT=iunit,FMT=10) otree_getX(rhadapt%rVertexCoordinates3D, ivt)
      end do
      do ivt = 1, otree_getsize(rhadapt%rVertexCoordinates3D)
        write(UNIT=iunit,FMT=10) otree_getY(rhadapt%rVertexCoordinates3D, ivt)
      end do
      do ivt = 1, otree_getsize(rhadapt%rVertexCoordinates3D)
        write(UNIT=iunit,FMT=10) otree_getZ(rhadapt%rVertexCoordinates3D, ivt)
      end do

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
      call sys_halt()
    end select
      

    ! Write cells to output file
    write(UNIT=iunit,FMT=*) 'cells ', rhadapt%NEL
    do iel = 1, rhadapt%NEL
      nve = hadapt_getNVE(rhadapt,iel)
      select case(nve)
      case(TRIA_NVETRI2D)
        write(UNIT=iunit,FMT=*) 'tri 3'
        write(UNIT=iunit,FMT=20) rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D, iel)

      case(TRIA_NVEQUAD2D)
        write(UNIT=iunit,FMT=*) 'quad 4'
        write(UNIT=iunit,FMT=30) rhadapt%p_IverticesAtElement(1:TRIA_NVEQUAD2D, iel)
        
      case DEFAULT
        call output_line('Invalid element type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hadapt_writeGridGMV')
        call sys_halt()
      end select
    end do

    ! Write velocity to output file
    write(UNIT=iunit,FMT=*) 'velocity 1'
    do ivt = 1, rhadapt%NVT
      write(UNIT=iunit,FMT=10) 0._DP
      write(UNIT=iunit,FMT=10) 0._DP
      write(UNIT=iunit,FMT=10) 0._DP
    end do

    ! Write variable for vertex age
    write(UNIT=iunit,FMT=*) 'variable'
    write(UNIT=iunit,FMT=*) 'vert_age 1'

    do ivt = 1, min(rhadapt%NVT, size(rhadapt%p_IvertexAge, 1))
      write(UNIT=iunit,FMT=40) rhadapt%p_IvertexAge(ivt)
    end do
    do ivt = min(rhadapt%NVT, size(rhadapt%p_IvertexAge, 1))+1, rhadapt%NVT
      write(UNIT=iunit,FMT=40) -99
    end do

    ! Write variable for element marker
    write(UNIT=iunit,FMT=*) 'elem_mark 0'
    do iel = 1, min(rhadapt%NEL, size(p_Imarker, 1))
      write(UNIT=iunit,FMT=40) p_Imarker(iel)
    end do
    do iel = min(rhadapt%NEL, size(p_Imarker, 1))+1, rhadapt%NEL
      write(UNIT=iunit,FMT=40) -99
    end do

    write(UNIT=iunit,FMT=*) 'endvars'
    write(UNIT=iunit,FMT=*) 'probtime ', 0._DP
    write(UNIT=iunit,FMT=*) 'endgmv'

    ! Close output file
    close(iunit)

10  format(E15.6E3)
20  format(3(1X,I8))
30  format(4(1X,I8))
40  format(I8)
  end subroutine hadapt_writeGridGMV

  ! ***************************************************************************

!<subroutine>

  subroutine hadapt_checkConsistency(rhadapt)

!<description>
    ! This subroutine checks the internal consistency of the dynamic data structures.
    ! Note that this routine performs brute-force search, and hence, should not
    ! be called in a productive environment. In is meant for debugging purposes
    ! only. If an error occurs, it stop without possibility to resume.
!</description>

!<inputoutput>
    ! adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx,p_IelementsAtVertex
    integer :: ipos,ivt,idx,iel,jel,jelmid,ive,jve,nve,mve
    integer :: h_IelementsAtVertexIdx,h_IelementsAtVertex
    logical :: btest,bfound

    ! Test #1: Consistency of element numbers
    btest = (rhadapt%NEL .eq. sum(rhadapt%InelOfType))
    call output_line('Test #1: Checking consistency of element numbers '//&
        merge('PASSED','FAILED',btest))

    ! Test #2: Vertex age must not exceed maximum refinement level
    btest = .true.
    do ivt = 1, rhadapt%NVT
      btest = btest .or. (rhadapt%p_IvertexAge(ivt) .gt. rhadapt%NSUBDIVIDEMAX)
    end do
    call output_line('Test #2: Checking maximum vertex age '//&
        merge('PASSED','FAILED',btest))

    ! Test #3: Check consistency of element neighbours
    btest = .true.
    do iel = 1, rhadapt%NEL

      ! Get number of vertices per element
      nve = hadapt_getNVE(rhadapt, iel)
      
      ! Loop over all adjacent elements
      do ive = 1, nve
        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)

        ! Check that adjacent element number is not larger than the
        ! total number of elements present in the triangulation
        if (jel > rhadapt%NEL .or. jelmid > rhadapt%NEL) then
          btest = .false.; cycle
        end if

        ! Do nothing if we are adjacent to the boundary
        if (jel .eq. 0 .or. jelmid .eq. 0) cycle
        
        ! Is neighboring element subdivided?
        if (jel .eq. jelmid) then
          
          ! Get number of vertices per element
          mve = hadapt_getNVE(rhadapt,jel)
          
          ! Find element IEL in adjacency list of JEL
          bfound = .false.
          do jve = 1, mve
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. &
                rhadapt%p_ImidneighboursAtElement(jve, jel)) then
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            else
              if (rhadapt%p_IneighboursAtElement(jve, jel)    .eq. iel .or.&
                  rhadapt%p_ImidneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            end if
          end do
          btest = btest .and. bfound
          
        else

          ! Get number of vertices per element
          mve = hadapt_getNVE(rhadapt, jel)
          
          ! Find element IEL in adjacency list of JEL
          bfound = .false.
          do jve = 1, mve
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. &
                rhadapt%p_ImidneighboursAtElement(jve, jel)) then
              if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            else
              if (rhadapt%p_IneighboursAtElement(jve, jel)    .eq. iel .or.&
                  rhadapt%p_ImidneighboursAtElement(jve, jel) .eq. iel) then
                bfound = .true.; exit
              end if
            end if
          end do
          btest = btest .and. bfound

          ! Get number of vertices per element
          mve = hadapt_getNVE(rhadapt, jelmid)
          
          ! Find element IEL in adjacency list of JELMID
          bfound = .false.
          do jve = 1, mve
            if (rhadapt%p_IneighboursAtElement(jve, jelmid) .eq. &
                rhadapt%p_ImidneighboursAtElement(jve, jelmid)) then
              if (rhadapt%p_IneighboursAtElement(jve, jelmid) .eq. iel) then
                bfound = .true.; exit
              end if
            else
              if (rhadapt%p_IneighboursAtElement(jve, jelmid)    .eq. iel .or.&
                  rhadapt%p_ImidneighboursAtElement(jve, jelmid) .eq. iel) then
                bfound = .true.; exit
              end if
            end if
          end do
          btest = btest .and. bfound

        end if
      end do
    end do
    call output_line('Test #3: Checking consistency of element neighbours '//&
        merge('PASSED','FAILED',btest))

    ! Test #4: Check consistency of common vertices between two edges
    btest = .true.
    do iel = 1, rhadapt%NEL

      ! Get number of vertices per element
      nve = hadapt_getNVE(rhadapt, iel)
      
      ! Loop over all adjacent elements
      do ive = 1, nve
        jel    = rhadapt%p_IneighboursAtElement(ive, iel)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive, iel)
        
        ! Check that adjacent element number is not larger than the
        ! total number of elements present in the triangulation
        if (jel > rhadapt%NEL .or. jelmid > rhadapt%NEL) then
          btest = .false.; cycle
        end if

        ! Do nothing if we are adjacent to the boundary
        if (jel .eq. 0 .or. jelmid .eq. 0) cycle
        
        ! Do nothing if there exists a temporal hanging node
        if (jel .ne. jelmid) cycle

        ! Get number of vertices per element
        mve = hadapt_getNVE(rhadapt, jel)

        ! Find element IEL in adjacency list of JEL
        bfound = .false.
        do jve = 1, mve
          if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. &
              rhadapt%p_ImidneighboursAtElement(jve, jel)) then
            if (rhadapt%p_IneighboursAtElement(jve, jel) .eq. iel) then
              bfound = .true.; exit
            end if
          else
            ! Do nothing if there exists a temporal hanging node
            cycle
          end if
        end do
        
        ! If the common edge has been found, check the two endpoints
        if (bfound) then
          bfound = ((rhadapt%p_IverticesAtElement(ive, iel) .eq. &
                     rhadapt%p_IverticesAtElement(modulo(jve,mve)+1, jel)) .and. &
                     rhadapt%p_IverticesAtElement(modulo(ive,nve)+1, iel) .eq. &
                     rhadapt%p_IverticesAtElement(jve, jel))
        end if
        btest = btest .and. bfound
      end do
    end do
    call output_line('Test #4: Checking consistency of common vertices along edges '//&
        merge('PASSED','FAILED',btest))

    ! Test #5: Check consistency of element-meeting-at-vertex lists
    btest = (rhadapt%rElementsAtVertex%NTABLE .eq. rhadapt%NVT)
    if (btest) then
      
      ! Create index array
      call storage_new('hadapt_checkConsistency', 'IelementAtVertexIdx', rhadapt%NVT+1,&
                       ST_INT, h_IelementsAtVertexIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)

      ! Count number of elements meeting at vertex
      do iel = 1, rhadapt%NEL

        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)

        ! Loop over corner vertices
        do ive = 1, nve
          ivt = rhadapt%p_IverticesAtElement(ive, iel)
          p_IelementsAtVertexIdx(ivt+1) = p_IelementsAtVertexIdx(ivt+1)+1
        end do
      end do

      ! Convert element couter into absolute position
      p_IelementsAtVertexIdx(1) = 1
      do ivt = 2, rhadapt%NVT+1
        p_IelementsAtVertexIdx(ivt) = p_IelementsAtVertexIdx(ivt)+&
                                      p_IelementsAtVertexIdx(ivt-1)
      end do

      ! Create working array
      call storage_new('hadapt_checkConsistency', 'IelementAtVertex',&
                       p_IelementsAtVertexIdx(rhadapt%NVT+1)-1,&
                       ST_INT, h_IelementsAtVertex, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(h_IelementsAtVertex, p_IelementsAtVertex)

      ! Retrieve the element numbers
      do iel = 1, rhadapt%NEL

        ! Get number of vertices per element
        nve = hadapt_getNVE(rhadapt, iel)

        ! Loop over corner vertices
        do ive = 1, nve
          ivt = rhadapt%p_IverticesAtElement(ive, iel)
          idx = p_IelementsAtVertexIdx(ivt)
          p_IelementsAtVertexIdx(ivt) = idx+1
          p_IelementsAtVertex(idx)    = iel
        end do
      end do
      
      ! Restore index array
      do ivt = rhadapt%NVT+1, 2, -1
        p_IelementsAtVertexIdx(ivt) = p_IelementsAtVertexIdx(ivt-1)
      end do
      p_IelementsAtVertexIdx(1) = 1

      ! Start to compare the temporal elements-meeting-at-vertex list
      ! and the dynamic data structure from the adaptivity structure
      do ivt = 1, rhadapt%NVT
        
        ! Get first entry in array list
        ipos = arrlst_getNextInArrayList(rhadapt%rElementsAtVertex, ivt, .true.)
        
        ! Repeat until there is no entry left in the array list
        do while(ipos .gt. ARRLST_NULL)
          
          ! Get element number IEL
          iel = rhadapt%rElementsAtVertex%p_IData(ipos)

          ! Proceed to next entry in array list
          ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex, ivt, .false.)

          ! Look for element IEL in temporal elements-meeting-at-vertex list
          ! If it does exist, multiply its value by minus one so that it cannot
          ! be found twice. At the end, all positive entries in the temporal
          ! list are not present in the dynamic structure
          bfound = .false.
          do idx = p_IelementsAtVertexIdx(ivt), p_IelementsAtVertexIdx(ivt+1)-1
            if (p_IelementsAtVertex(idx) .eq. iel) then
              p_IelementsAtVertex(idx) = -p_IelementsAtVertex(idx)
              bfound = .true.
              exit
            end if
          end do
          btest = btest .and. bfound
        end do
      end do
      btest = btest .and. all(p_IelementsAtVertex < 0)
      call output_line('Test #5: Checking consistency of elements meeting at vertices '//&
        merge('PASSED','FAILED',btest))

      ! Release auxiliary storage
      call storage_free(h_IelementsAtVertexIdx)
      call storage_free(h_IelementsAtVertex)
    else

      call output_line('Test #5: Checking consistency of element-meeting-at-vertex list '//&
          merge('PASSED','FAILED',btest))
    end if
  end subroutine hadapt_checkConsistency

  ! ***************************************************************************

!<subroutine>
  
  subroutine hadapt_refreshAdaptation(rhadapt, rtriangulation)

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
    type(t_triangulation), intent(IN) :: rtriangulation
!</input>

!<inputoutput>
    ! Adaptivity structure
    type(t_hadapt), intent(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>


    ! Get pointer to InodalProperty
    if (rhadapt%h_InodalProperty .eq. &
        rtriangulation%h_InodalProperty) then
      call storage_getbase_int(rhadapt%h_InodalProperty,&
                               rhadapt%p_InodalProperty)
    else
      call output_line('Inconsistent handle h_InodalProperty',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      call sys_halt()
    end if

    ! Get pointer to IverticesAtElement
    if (rhadapt%h_IverticesAtElement .eq. &
        rtriangulation%h_IverticesAtElement) then
      call storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
                                 rhadapt%p_IverticesAtElement)
    else
      call output_line('Inconsistent handle h_IverticesAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      call sys_halt()
    end if

    ! Get pointer to IneighboursAtElement
    if (rhadapt%h_IneighboursAtElement .eq. &
        rtriangulation%h_IneighboursAtElement) then
      call storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
                                 rhadapt%p_IneighboursAtElement)
    else
      call output_line('Inconsistent handle h_IneighboursAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      call sys_halt()
    end if

    ! Get pointer to ImidneighboursAtElement
    if (rhadapt%h_ImidneighboursAtElement .ne. ST_NOHANDLE) then
      call storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
                                 rhadapt%p_ImidneighboursAtElement)
    else
      call output_line('Inconsistent handle h_ImidneighboursAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      call sys_halt()
    end if

    ! Get pointer to IvertexAge
     if (rhadapt%h_IvertexAge .ne. ST_NOHANDLE) then
      call storage_getbase_int(rhadapt%h_IvertexAge,&
                               rhadapt%p_IvertexAge)
    else
      call output_line('Inconsistent handle h_ImidneighboursAtElement',&
          OU_CLASS_ERROR,OU_MODE_STD,'hadapt_refreshAdaptation')
      call sys_halt()
    end if

  end subroutine hadapt_refreshAdaptation

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  ! Internal auxiliary routines and functions. No user should ever look at them

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine redgreen_refine(rhadapt, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine performs red-green refinement as proposed by R. Bank
!</description>

!<input>
    ! Callback function
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adaptive data structure
    type(t_hadapt), intent(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    type(t_collection), intent(INOUT), optional :: rcollection
!</subroutine>
    
    ! local variables
    integer, dimension(:), pointer :: p_Imarker
    integer :: iel
    
    ! Check if dynamic data structures are o.k. and if 
    ! cells are marked for refinement
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA) .ne. HADAPT_HAS_DYNAMICDATA .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDREFINE)    .ne. HADAPT_MARKEDREFINE) then
      call output_line('Dynamic data structures are not generated &
          &or no marker for refinement is available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'redgreen_refine')
      call sys_halt()
    end if
    
    ! Set pointers
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)
        
    ! Perform red-green refinement
    do iel = 1, size(p_Imarker, 1)
      
      select case(p_Imarker(iel))
      case(MARK_ASIS_TRIA,&
           MARK_ASIS_QUAD)

        ! Do nothing for elements that should be kept 'as is'

      case(:MARK_CRS_GENERIC)

        ! Do nothing for element that have been marked for coarsening

      case(MARK_REF_TRIA4TRIA)

        ! Red refinement triangle
        call refine_Tria4Tria(rhadapt, iel,&
                              rcollection, fcb_hadaptCallback)
        p_Imarker(iel) = MARK_ASIS        

      case(MARK_REF_QUAD4QUAD)

        ! Red refinement quadrilateral
        call refine_Quad4Quad(rhadapt, iel,&
                              rcollection, fcb_hadaptCallback)
        p_Imarker(iel) = MARK_ASIS
        
      case(MARK_REF_TRIA2TRIA_1,&
           MARK_REF_TRIA2TRIA_2,&
           MARK_REF_TRIA2TRIA_3)   

        ! Green refinement triangle
        call refine_Tria2Tria(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_REF_QUAD3TRIA_1,&
           MARK_REF_QUAD3TRIA_2,&
           MARK_REF_QUAD3TRIA_3,&
           MARK_REF_QUAD3TRIA_4)

        ! Green refinement quadrilateral
        call refine_Quad3Tria(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+3
        p_Imarker(iel)         = MARK_ASIS
        
      case(MARK_REF_QUAD2QUAD_13,&
           MARK_REF_QUAD2QUAD_24)

        ! Green refinement quadrilateral
        call refine_Quad2Quad(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2
        p_Imarker(iel)         = MARK_ASIS

      case(MARK_REF_QUAD4TRIA_12,&
           MARK_REF_QUAD4TRIA_23,&
           MARK_REF_QUAD4TRIA_34,&
           MARK_REF_QUAD4TRIA_14)

        ! Green refinement quadrilateral
        call refine_Quad4Tria(rhadapt, iel, p_Imarker(iel),&
                              rcollection, fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+4
        p_Imarker(iel)         = MARK_ASIS

      case DEFAULT
        call output_line('Invalid element refinement marker!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_refine')
        call sys_halt()
      end select
    end do

    ! Increase the number of refinement steps by one
    rhadapt%nRefinementSteps = rhadapt%nRefinementSteps+1
    
    ! Refinement has been performed.
    rhadapt%iSpec = ior(rhadapt%ispec, HADAPT_REFINED)
    
    ! Hence, the markers are no longer valid
    rhadapt%iSpec = iand(rhadapt%iSpec, not(HADAPT_MARKEDREFINE))
  end subroutine redgreen_refine

  ! ***************************************************************************

!<subroutine>

  subroutine redgreen_coarsen(rhadapt, rcollection, fcb_hadaptCallback)

!<description>
    ! This subroutine performs red-green coarsening as proposed by R. Bank
!</description>

!<input>
    ! callback routines
    include 'intf_hadaptcallback.inc'
    optional :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    type(t_hadapt), intent(INOUT)               :: rhadapt
    
    ! OPTIONAL Collection
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer,  dimension(:), pointer :: p_Imarker
    integer, dimension(1) :: Ielements
    integer, dimension(2)  :: Ivertices
    integer :: ipos,iel,jel,ivt,ivtReplace,ive

    ! Check if dynamic data structures are o.k. and 
    ! if  cells are marked for coarsening
    if (iand(rhadapt%iSpec, HADAPT_HAS_DYNAMICDATA) .ne. HADAPT_HAS_DYNAMICDATA .or.&
        iand(rhadapt%iSpec, HADAPT_MARKEDCOARSEN)   .ne. HADAPT_MARKEDCOARSEN) then
      call output_line('Dynamic data structures are not generated &
          &or no marker for coarsening is available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'redgreen_coarsen')
      call sys_halt()
    end if
    call storage_getbase_int(rhadapt%h_Imarker, p_Imarker)
    

    ! Perform hierarchical red-green recoarsening
    element: do iel = size(p_Imarker,1), 1, -1

      select case(p_Imarker(iel))
      case(MARK_CRS_GENERIC:)

        ! Do nothing for elements that should be kept 'as is'
        ! and those which are marked for refinement.
        
      case(MARK_CRS_2TRIA1TRIA)
        call coarsen_2Tria1Tria(rhadapt,iel,&
                                rcollection, fcb_hadaptCallback)
        
      case(MARK_CRS_4TRIA1TRIA)
        call coarsen_4Tria1Tria(rhadapt,iel,&
                                rcollection, fcb_hadaptCallback)

      case(MARK_CRS_4TRIA2TRIA_1,&
           MARK_CRS_4TRIA2TRIA_2,&
           MARK_CRS_4TRIA2TRIA_3)
        call coarsen_4Tria2Tria(rhadapt, iel, p_Imarker(iel),&
                                rcollection, fcb_hadaptCallback)

      case(MARK_CRS_3TRIA1QUAD)
        call coarsen_3Tria1Quad(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)

      case(MARK_CRS_4TRIA1QUAD)
        call coarsen_4Tria1Quad(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)
        
      case(MARK_CRS_4TRIA3TRIA_LEFT,&
           MARK_CRS_4TRIA3TRIA_RIGHT)
        call coarsen_4Tria3Tria(rhadapt, iel, p_Imarker(iel),&
                                rcollection, fcb_hadaptCallback)

      case(MARK_CRS_2QUAD1QUAD)
        call coarsen_2Quad1Quad(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)

      case(MARK_CRS_2QUAD3TRIA)
        call coarsen_2Quad3Tria(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)
        
      case(MARK_CRS_4QUAD1QUAD)
        call coarsen_4Quad1Quad(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)

      case(MARK_CRS_4QUAD2QUAD)
        call coarsen_4Quad2Quad(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)

      case(MARK_CRS_4QUAD3TRIA)
        call coarsen_4Quad3Tria(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)
        
      case(MARK_CRS_4QUAD4TRIA)
        call coarsen_4Quad4Tria(rhadapt, iel,&
                                rcollection, fcb_hadaptCallback)

      case DEFAULT
        call output_line('Invalid recoarsening marker!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_coarsen')
        call sys_halt()
      end select

    end do element


    ! Loop over all vertices 1...NVT0 present in the triangulation before
    ! refinement and check if they are free for vertex removal.
    vertex: do ivt = rhadapt%NVT0, 1, -1
      
      ! If the vertex is locked, then skip this vertex
      if (rhadapt%p_IvertexAge(ivt) .le. 0) cycle vertex
      
      ! Remove vertex physically. Note that this vertex is no longer associated
      ! to any element. All associations have been removed in the above element
      ! coarsening/conversion step. In order to prevent "holes" in the vertex list,
      ! vertex IVT is replaced by the last vertex if it is not the last one itself.
      call remove_vertex2D(rhadapt, ivt, ivtReplace)
      
      ! If vertex IVT was not the last one, update the "elements-meeting-at-vertex" list
      if (ivtReplace .ne. 0) then
        
        ! Start with first element in "elements-meeting-at-vertex" list of the replaced vertex
        ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,&
                                         ivtReplace, .true.)

        update: do while(ipos .gt. ARRLST_NULL)
          
          ! Get element number JEL
          jel = rhadapt%rElementsAtVertex%p_IData(ipos)
          
          ! Proceed to next element
          ipos = arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,&
                                           ivtReplace, .false.)
          
          ! Look for vertex ivtReplace in element JEL and replace it by IVT
          do ive = 1, hadapt_getNVE(rhadapt, jel)
            if (rhadapt%p_IverticesAtElement(ive, jel) .eq. ivtReplace) then
              rhadapt%p_IverticesAtElement(ive, jel) = ivt
              cycle update
            end if
          end do
          
          ! If the replaced vertex ivtReplace could not be found in element JEL
          ! something is wrong and we stop the simulation
          call output_line('Unable to find replacement vertex in element',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_coarsen')
          call sys_halt()            
        end do update
        
        ! Swap tables IVT and ivtReplace in arraylist and release table ivtReplace
        call arrlst_swapArrayList(rhadapt%rElementsAtVertex, ivt, ivtReplace)
        call arrlst_releaseArrayList(rhadapt%rElementsAtVertex, ivtReplace)
        
      else
        
        ! Release table IVT
        call arrlst_releaseArrayList(rhadapt%rElementsAtVertex, ivt)
      end if
            
      ! Optionally, invoke callback function
      if (present(fcb_hadaptCallback) .and. present(rcollection)) then
        Ivertices = (/ivt,ivtReplace/); Ielements = (/0/)
        call fcb_hadaptCallback(rcollection, HADAPT_OPR_REMOVEVERTEX,&
                                Ivertices, Ielements)
      end if
    end do vertex
        
    ! Increase the number of recoarsening steps by one
    rhadapt%nCoarseningSteps = rhadapt%nCoarseningSteps+1

    ! Coarsening has been performed.
    rhadapt%iSpec = ior(rhadapt%ispec, HADAPT_COARSENED)
    
    ! Hence, the markers are no longer valid
    rhadapt%iSpec = iand(rhadapt%iSpec, not(HADAPT_MARKEDCOARSEN))
  end subroutine redgreen_coarsen
end module hadaptivity

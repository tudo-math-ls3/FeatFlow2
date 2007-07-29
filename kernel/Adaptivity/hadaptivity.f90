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
!#  5.) hadapt_setVertexCoords2D
!#      -> Set the coordinates of vertices to the adaptivity structure in 2D
!#
!#  6.) hadapt_getVertexCoords2D
!#      -> Get the coordinates of vertices from the adaptivity structure in 2D
!#
!#  7.) hadapt_setVertexCoords3D
!#      -> Set the coordinates of vertices to the adaptivity structure in 3D
!#
!#  8.) hadapt_getVertexCoords3D
!#      -> Get the coordinates of vertices from the adaptivity structure in 3D
!#
!#  9.) hadapt_setVerticesAtElement
!#      -> Set the "vertices-at-element" structure to the adaptivity structure
!#
!# 10.) hadapt_getVerticesAtElement
!#      -> Get the "vertices-at-element" structure from the adaptivity structure
!#
!# 11.) hadapt_setNeighboursAtElement
!#      -> Set the "neighbours-at-element" structure to the adaptivity structure
!#
!# 12.) hadapt_getNeighboursAtElement
!#      -> Get the "neighbours-at-element" structure from the adaptivity structure
!#
!# 13.) hadapt_setNelOfType
!#      -> Set the "InelOfType" array to the adaptivity structure
!#
!# 14.) hadapt_getNelOfType
!#      -> Get the "InelOfType" array from the adaptivity structure
!#
!# 15.) hadapt_setBoundary
!#      -> Set the boundary structure to the adaptivity structure
!#
!# 16.) hadapt_getBoundary
!#      -> Get the boundary structure form the adaptivity structure
!#
!# 17.) hadapt_setNodalProperty
!#      -> Set the "nodal property" list to the adaptivity structure
!#
!# 18.) hadapt_getNodalProperty
!#      -> Get the "nodal property" list from the adaptivity structure
!#
!# 19.) hadapt_genElementsAtVertex
!#      -> Generate the "elements-at-vertex" data structure
!#
!# 20.) hadapt_performAdaptation
!#      -> perform one step of grid adaptation
!#
!# 21.) hadapt_info
!#      -> output information about the adaptivity structure
!#
!# 22.) hadapt_writeGridSVG
!#      -> write the adapted grid to file in SVG format
!#
!# The following internal routines are available:
!#
!#  1.) add_vertex2D = add_vertex_atEdgeMidpoint2D /
!#                     add_vertex_atElementCenter2D
!#      -> add a new vertex to the adaptation data structure in 2D
!#
!#  2.) remove_vertex2D
!#      -> remove an existing vertex from the adaptation data structure in 2D
!#
!#  3.) replace_element2D = replace_elementTria /
!#                          replace_elementQuad
!#      -> replace an existing element by another element of he same type in 2D
!#
!#  4.) add_element2D = add_elementTria /
!#                      add_elementQuad
!#      -> add a new element to the adaptation data structure in 2D
!#
!#  5.) remove_element2D
!#      -> remove an existing element from the adaptation data structure in 2D
!#
!#  6.) update_ElementNeighbors2D = update_ElemEdgeNeighbors2D /
!#                                  update_ElemNeighbors2D
!#      -> update the list of neighboring elements along a given edge in 2D
!#
!#  7.) refine_Tria2Tria
!#      -> refine a triangle by subdivision into two triangles
!#
!#  8.) refine_Tria3Tria
!#      -> refine a triangle by subdivision into three triangles
!#
!#  9.) refine_Tria4Tria
!#      -> refine a triangle by subdivision into four triangles
!#
!# 10.) refine_Quad2Quad
!#      -> refine a quadrilateral by subdivision into two quadrilaterals
!#
!# 11.) refine_Quad3Tria
!#      -> refine a quadrilateral by subdivision into three triangles
!#
!# 12.) refine_Quad4Tria
!#      -> refine a quadrilateral by subdivision into four triangles
!#
!# 13.) refine_Quad4Quad
!#      -> refine a quadrilateral by subdivision into four quadrilaterals
!#
!# 14.) convert_Tria2Tria
!#      -> convert two neighboring triangles into four similar triangle
!#
!# 15.) convert_Quad2Quad
!#      -> convert two neighboring quadrilaterals into four similar quadrilaterals
!#
!# 16.) convert_Quad3Tria
!#      -> convert three neighboring triangles into four similar quadrilaterals
!#
!# 17.) convert_Quad4Tria
!#      -> convert four neighboring triangles into four similar quadrilaterals
!#
!# 19.) coarsen_2Tria1Tria
!#      -> coarsen two green triangles by combining them to the macro element
!#
!# 20.) coarsen_4Tria1Tria
!#      -> coarsen four red triangles by combining them to the macro element
!#
!# 21.) coarsen_4Tria2Tria
!#      -> coarsen four red triangles by combining them to two green elements
!#
!# 22.) mark_refinement2D
!#      -> mark elements for refinement in 2D
!#
!# 23.) redgreen_mark_coarsening2D
!#      -> mark elements for coarsening in 2D
!#
!# 24.) redgreen_mark_refinement2D
!#      -> mark elements for refinement due to Red-Green strategy in 2D
!#
!# 25.) redgreen_refine
!#      -> perform red-green refinement for marked elements
!#
!# 26.) redgreen_coarsen
!#      -> perform red-green coarsening for marked elements
!#
!# 27.) redgreen_getState
!#      -> return the state of an element
!#
!# 28.) redgreen_getStateTria
!#      -> return the state of a triangle
!#
!# 29.) redgreen_getStateQuad
!#      -> return the state of a quadrilateral
!#
!# 30.) redgreen_rotateState
!#      -> compute the state of an element after rotation
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

  USE fsystem
  USE storage
  USE paramlist
  USE collection
  USE triangulation
  USE linearsystemscalar
  USE quadtree
  USE octree
  USE arraylist
  USE binarytree
  USE list
  USE sort
  USE io
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: hadapt_initFromParameterlist
  PUBLIC :: hadapt_initFromTriangulation
  PUBLIC :: hadapt_generateRawMesh
  PUBLIC :: hadapt_releaseAdaptation
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
  PUBLIC :: hadapt_info
  PUBLIC :: hadapt_writeGridSVG


!<constants>

!<constantblock description="Global constants for grid modification operations">

  ! Operation identifier for initialization of callback function
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_INITCALLBACK     = -1

  ! Operation identifier for finalization of callback function
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_DONECALLBACK     = -2

  ! Operation identifier for adjustment of vertex dimension
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_ADJUSTVERTEXDIM  = 1

  ! Operation identifier for vertex insertion at edge midpoint
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_INSERTVERTEXEDGE = 2

  ! Operation identifier for vertex insertion at element centroid
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_INSERTVERTEXCENTR= 3

  ! Operation identifier for vertex removal
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REMOVEVERTEX     = 4

  ! Operation identifier for refinment: 1-tria : 2-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_TRIA2TRIA    = 5

  ! Operation identifier for refinment: 1-tria : 3-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_TRIA3TRIA12  = 6
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_TRIA3TRIA23  = 7
 
  ! Operation identifier for refinment: 1-tria : 4-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_TRIA4TRIA    = 8

  ! Operation identifier for refinment: 1-quad : 2-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_QUAD2QUAD    = 9
  
  ! Operation identifier for refinment: 1-quad : 3-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_QUAD3TRIA    = 10

  ! Operation identifier for refinment: 1-quad : 4-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_QUAD4TRIA    = 11

  ! Operation identifier for refinment: 1-quad : 4-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_REF_QUAD4QUAD    = 12

  ! Operation identifier for conversion: 2-tria : 4-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CVT_TRIA2TRIA    = 13

  ! Operation identifier for conversion: 2-quad : 4-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CVT_QUAD2QUAD    = 14

  ! Operation identifier for conversion: 3-tria : 4-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CVT_QUAD3TRIA    = 15

  ! Operation identifier for conversion: 4-tria : 4-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CVT_QUAD4TRIA    = 16

  ! Operation identifier fo coarsening: 2-tria : 1-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_2TRIA1TRIA   = 17

  ! Operation identifier fo coarsening: 4-tria : 1-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA1TRIA   = 18

  ! Operation identifier fo coarsening: 4-tria : 2-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA2TRIA1  = 19
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA2TRIA2  = 20
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA2TRIA3  = 21

!</constantblock>

!<constantblock description="Global flags for grid refinement/coarsening">

  ! No refinement
  INTEGER, PARAMETER :: HADAPT_NOREFINEMENT          = 0

  ! No coarsening
  INTEGER, PARAMETER :: HADAPT_NOCOARSENING          = 0

  ! Red-Green refinement strategy (R. Banks)
  INTEGER, PARAMETER :: HADAPT_REDGREEN              = 1

  ! Longest edge bisection strategy (M. Rivara)
  INTEGER, PARAMETER :: HADAPT_LONGESTEDGE           = 2

!</constantblock>

!<constantblock description="Bitfield identifiers for state of adaptation">

  ! Adaptation is undefined
  INTEGER, PARAMETER, PUBLIC :: HADAPT_UNDEFINED     = 2**0

  ! Parameters of adaptivity structure are initialized
  INTEGER, PARAMETER :: HADAPT_HAS_PARAMETERS        = 2**1

  ! Quadtree/octree for vertex coordinates is generated
  INTEGER, PARAMETER :: HADAPT_HAS_COORDS            = 2**2

  ! Array for IverticesAtElement is generated
  INTEGER, PARAMETER :: HADAPT_HAS_VERTATELEM        = 2**3

  ! Array for IneighboursAtElement is generated
  INTEGER, PARAMETER :: HADAPT_HAS_NEIGHATELEM       = 2**4

  ! Boundary data is generated
  INTEGER, PARAMETER :: HADAPT_HAS_BOUNDARY          = 2**5

  ! Nodal property is generated
  INTEGER, PARAMETER :: HADAPT_HAS_NODALPROP         = 2**6

  ! Number of elements for predefined type
  INTEGER, PARAMETER :: HADAPT_HAS_NELOFTYPE         = 2**7

  ! Array for IelementsAtVertex is generated
  INTEGER, PARAMETER :: HADAPT_HAS_ELEMATVERTEX      = 2**8

  ! Dynamic data structures are all generated
  INTEGER, PARAMETER :: HADAPT_HAS_DYNAMICDATA       = HADAPT_HAS_PARAMETERS+&
                                                       HADAPT_HAS_COORDS+&
                                                       HADAPT_HAS_VERTATELEM+&
                                                       HADAPT_HAS_NEIGHATELEM+&
                                                       HADAPT_HAS_BOUNDARY+&
                                                       HADAPT_HAS_NODALPROP+&
                                                       HADAPT_HAS_NELOFTYPE+&
                                                       HADAPT_HAS_ELEMATVERTEX

  ! Cells are marked for refinement
  INTEGER, PARAMETER :: HADAPT_MARKEDREFINE          = 2**9

  ! Cells are marked for coarsening
  INTEGER, PARAMETER :: HADAPT_MARKEDCOARSEN         = 2**10

  ! Cells are marked
  INTEGER, PARAMETER :: HADAPT_MARKED                = HADAPT_MARKEDREFINE+&
                                                       HADAPT_MARKEDCOARSEN
  
  ! Grid has been refined
  INTEGER, PARAMETER :: HADAPT_REFINED               = 2**11
  
  ! Grid has been coarsened
  INTEGER, PARAMETER :: HADAPT_COARSENED             = 2**12

!</constantblock>

!<constantblock description="Constants for element marker">

  ! Mark sub-element for generic recoarsening. Note that the detailed
  ! recoarsening strategy must be determined form the element state.
  INTEGER, PARAMETER :: MARK_CRS_GENERIC            = -1

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into the macro element
  INTEGER, PARAMETER :: MARK_CRS_4TRIA1TRIA         = -2

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into two green triangles, whereby the first vertex
  ! if the inner triangle is connected to the opposite vertex
  INTEGER, PARAMETER :: MARK_CRS_4TRIA2TRIA_1       = -3

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into two green triangles, whereby the second vertex
  ! if the inner triangle is connected to the opposite vertex
  INTEGER, PARAMETER :: MARK_CRS_4TRIA2TRIA_2       = -4

  ! Mark inner triangle of a 1-tria : 4-tria refinement for
  ! recoarsening into two green triangles, whereby the third vertex
  ! if the inner triangle is connected to the opposite vertex
  INTEGER, PARAMETER :: MARK_CRS_4TRIA2TRIA_3       = -5

  ! Mark inner green triangle of a 1-quad : 3-tria refinement for
  ! recoarsening into the macro quadrilateral together with is neighbors
  INTEGER, PARAMETER :: MARK_CRS_3TRIA1QUAD         = -6

  ! Mark left green triangle of a 1-tria : 2-tia refinement for 
  ! recoarsening into the macro triangle together with its right neighbor
  INTEGER, PARAMETER :: MARK_CRS_2TRIA1TRIA_LEFT    = -7

  ! Mark right green triangle of a 1-tria : 2-tia refinement for 
  ! recoarsening into the macro triangle together with its left neighbor
  INTEGER, PARAMETER :: MARK_CRS_2TRIA1TRIA_RIGHT   = -8

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! recoarsening into the macro quadrilateral together with its three neighbors
  INTEGER, PARAMETER :: MARK_CRS_4TRIA1QUAD         = -9

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! conversion into three triangles keeping the left neighboring triangle
  INTEGER, PARAMETER :: MARK_CRS_4TRIA3TRIA_LEFT    = -10

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! conversion into three triangles keeping the right neighboring triangle
  INTEGER, PARAMETER :: MARK_CRS_4TRIA3TRIA_RIGHT   = -11

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for recoarsening
  INTEGER, PARAMETER :: MARK_CRS_4QUAD1QUAD         = -12

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! three triangles keeping the midpoint of the first edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD3TRIA_1       = -13

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! three triangles keeping the midpoint of the second edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD3TRIA_2       = -14

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! three triangles keeping the midpoint of the third edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD3TRIA_3       = -15

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! three triangles keeping the midpoint of the fourth edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD3TRIA_4       = -16

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! four triangles keeping the midpoints of the first and second edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD4TRIA_12      = -17

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! four triangles keeping the midpoints of the second and third edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD4TRIA_23      = -18

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! four triangles keeping the midpoints of the third and fourth edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD4TRIA_34      = -19

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! four triangles keeping the midpoints of the first and fourth edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD4TRIA_14      = -20

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! two quadrilaterals keeping the midpoints of the first and third edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD2QUAD_13      = -21

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for conversion into
  ! two quadrilaterals keeping the midpoints of the second and fourth edge
  INTEGER, PARAMETER :: MARK_CRS_4QUAD2QUAD_24      = -22

  ! Mark for keeping element 'as is'
  INTEGER, PARAMETER :: MARK_ASIS                   = 0
  INTEGER, PARAMETER :: MARK_ASIS_TRIA              = 0
  INTEGER, PARAMETER :: MARK_ASIS_QUAD              = 1

  ! Mark element for 1-tria : 2-tria refinement along first edge
  INTEGER, PARAMETER :: MARK_REF_TRIA2TRIA_1        = 2

  ! Mark element for 1-tria : 2-tria refinement along second edge
  INTEGER, PARAMETER :: MARK_REF_TRIA2TRIA_2        = 4

  ! Mark element for 1-tria : 2-tria refinement along third edge
  INTEGER, PARAMETER :: MARK_REF_TRIA2TRIA_3        = 8

  ! Mark element for 1-tria : 3-tria refinement along first and second edge
  INTEGER, PARAMETER :: MARK_REF_TRIA3TRIA_12       = 6

  ! Mark element for 1-tria : 3-tria refinement along second and third edge
  INTEGER, PARAMETER :: MARK_REF_TRIA3TRIA_23       = 12

  ! Mark element for 1-tria : 3-tria refinement along first and third edge
  INTEGER, PARAMETER :: MARK_REF_TRIA3TRIA_13       = 10

  ! Mark element for 1-tria : 4-tria red refinement
  INTEGER, PARAMETER :: MARK_REF_TRIA4TRIA          = 14

  ! Mark element for 1-quad : 3-tria refinement along first edge
  INTEGER, PARAMETER :: MARK_REF_QUAD3TRIA_1        = 3
  
  ! Mark element for 1-quad : 3-tria refinement along second edge
  INTEGER, PARAMETER :: MARK_REF_QUAD3TRIA_2        = 5
  
  ! Mark element for 1-quad : 3-tria refinement along third edge
  INTEGER, PARAMETER :: MARK_REF_QUAD3TRIA_3        = 9
  
  ! Mark element for 1-quad : 3-tria refinement along fourth edge
  INTEGER, PARAMETER :: MARK_REF_QUAD3TRIA_4        = 17

  ! Mark element for 1-quad : 4-tria refinement along first and second edge
  INTEGER, PARAMETER :: MARK_REF_QUAD4TRIA_12       = 7

  ! Mark element for 1-quad : 4-tria refinement along second and third edge
  INTEGER, PARAMETER :: MARK_REF_QUAD4TRIA_23       = 13

  ! Mark element for 1-quad : 4-tria refinement along third and fourth edge
  INTEGER, PARAMETER :: MARK_REF_QUAD4TRIA_34       = 25

  ! Mark element for 1-quad : 4-tria refinement along first and fourth edge
  INTEGER, PARAMETER :: MARK_REF_QUAD4TRIA_14       = 19

  ! Mark element for 1-quad : 4-quad red refinement
  INTEGER, PARAMETER :: MARK_REF_QUAD4QUAD          = 31

  ! Mark element for 1-quad : 2-quad refinement along first and third edge
  INTEGER, PARAMETER :: MARK_REF_QUAD2QUAD_13       = 11

  ! Mark element for 1-quad : 2-quad refinement along second and fourth edge
  INTEGER, PARAMETER :: MARK_REF_QUAD2QUAD_24       = 21

  ! Mark element for 1-quad blue refinement
  INTEGER, PARAMETER :: MARK_REF_QUADBLUE_412       = 23
  INTEGER, PARAMETER :: MARK_REF_QUADBLUE_234       = 29
  INTEGER, PARAMETER :: MARK_REF_QUADBLUE_123       = 15
  INTEGER, PARAMETER :: MARK_REF_QUADBLUE_341       = 27
  
!</constant>

!<constantblock description="Constants for element states">

  ! Triangle from the root triangulation
  INTEGER, PARAMETER :: STATE_TRIA_ROOT             = 0

  ! Inner triangle of a 1-tria : 4-tria red refinement
  INTEGER, PARAMETER :: STATE_TRIA_REDINNER         = 14

  ! Outer triangle of a 1-tria : 4-tria red refinement or
  ! outer/inner triangle of 1-quad : 4-tria refinement
  INTEGER, PARAMETER :: STATE_TRIA_OUTERINNER       = 4
  INTEGER, PARAMETER :: STATE_TRIA_OUTERINNER1      = 2   ! only theoretically
  INTEGER, PARAMETER :: STATE_TRIA_OUTERINNER2      = 8   ! only theoretically

  ! Inner triangle of a 1-quad : 3-tria refinement
  INTEGER, PARAMETER :: STATE_TRIA_GREENINNER       = -2

  ! Outer trignale of a green refinement
  INTEGER, PARAMETER :: STATE_TRIA_GREENOUTER_LEFT  = -4
  INTEGER, PARAMETER :: STATE_TRIA_GREENOUTER_RIGHT = -8

  ! Quadrilateral form the root triangulation
  INTEGER, PARAMETER :: STATE_QUAD_ROOT             = 1

  ! Quadrilateral of a 1-quad : 4-quad red refinement
  INTEGER, PARAMETER :: STATE_QUAD_RED1             = 25
  INTEGER, PARAMETER :: STATE_QUAD_RED2             = 19
  INTEGER, PARAMETER :: STATE_QUAD_RED3             = 7
  INTEGER, PARAMETER :: STATE_QUAD_RED4             = 13

  ! Quadrilateral of a 1-quad : 2-quad refinement
  INTEGER, PARAMETER :: STATE_QUAD_HALF1            = 21
  INTEGER, PARAMETER :: STATE_QUAD_HALF2            = 11
  

!</constantblock>

!<constantblock description="Constants for grid adaptation">
  
  ! Array position of the boundary
  INTEGER, PARAMETER :: BdrValue = 1
  
  ! Array position of the previous boundary vertex
  INTEGER, PARAMETER :: BdrPrev  = 1

  ! Array position of the next boundary vertex
  INTEGER, PARAMETER :: BdrNext  = 2

!</constantblock>

!</constants>  

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<types>

  !<typeblock>
  
  ! This type contains all data structures to handle 
  ! adaptive grid refinement and grid coarsening.
  TYPE, PUBLIC :: t_hadapt
    ! Format Tag: Specifies the state of adaptation
    INTEGER :: iSpec                                 = HADAPT_UNDEFINED

    ! Tag: Specified the strategy for grid refinement
    INTEGER :: irefinementStrategy                   = HADAPT_NOREFINEMENT

    ! Tag: Specifies the coarsening strategy
    INTEGER :: icoarseningStrategy                   = HADAPT_NOCOARSENING

    ! Maximum number of subdivisions from the original mesh
    INTEGER :: NSUBDIVIDEMAX                         = 0

    ! Total number of grid refinement steps
    INTEGER :: nRefinementSteps                      = 0

    ! Total number of grid coarsening steps
    INTEGER :: nCoarseningSteps                      = 0

    ! Total number of grid smoothing steps
    INTEGER :: nSmoothingSteps                       = 0

    ! Tolerance for refinement
    REAL(DP) :: drefinementTolerance                 = 0

    ! Tolerance for coarsening
    REAL(DP) :: dcoarseningTolerance                 = 0

    ! Dimension of the triangulation
    INTEGER :: ndim                                  = 0

    ! Total number of vertices (initially)
    INTEGER(PREC_VERTEXIDX) :: NVT0                  = 0
    
    ! Total number of vertices
    INTEGER(PREC_VERTEXIDX) :: NVT                   = 0

    ! Increment of vertices
    INTEGER(PREC_VERTEXIDX) :: increaseNVT           = 0

    ! Total number of boundary vertives (initially)
    INTEGER :: NVBD0                                 = 0

    ! Total number of boundary vertives
    INTEGER :: NVBD                                  = 0
    
    ! Total number of boundary components (should not change)
    INTEGER :: NBCT                                  = 0
    
    ! Total number of elements (initially)
    INTEGER(PREC_ELEMENTIDX) :: NEL0                 = 0
    
    ! Total number of elements
    INTEGER(PREC_ELEMENTIDX) :: NEL                  = 0

    ! Maximum number of elements (before reallocation)
    INTEGER(PREC_ELEMENTIDX) :: NELMAX               = 0
    
    ! Total number of green elements (required internally)
    INTEGER(PREC_ELEMENTIDX) :: nGreenElements       = 0
    
    ! Nuber of elements with a defined number of vertices per element.
    ! InelOfType(TRIA_NVETRI2D)  = number of triangles in the mesh (2D).
    ! InelOfType(TRIA_NVEQUAD2D) = number of quadrilaterals in the mesh (2D).
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE) :: InelOfType = 0

    ! Same as InelOfType but this array stores the number of elements
    ! which are initially present in the mesh
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE) :: InelOfType0 = 0
    
    ! Element marker array.
    ! Handle to
    !       p_Imarker = array [1..NEL] of integer.
    ! For each element there is one identifier [0,..,31] in the marker.
    ! BIT0: Is 0 for triangle and 1 for quadrilateral
    ! BIT1: Is 0 if first edge is not subdivided, 1 otherwise.
    ! BIT2: Is 0 if second edge is not subdivided, 1 otherwise.
    ! BIT3: Is 0 if third edge is not subdivided, 1 otherwise.
    ! BIT4: Is 0 if fourth edge is not subdivided, 1 otherwise.
    INTEGER :: h_Imarker = ST_NOHANDLE
    
    ! Vertex age array.
    ! Handle to
    !       p_IvertexAge = array [1..NVT] of integer
    ! Each node has an individual age which is used to determine
    ! the type of the element (green/red triangle/quadrilateral).
    ! The nodes of the initial grid are assigned age 0 by default.
    ! Each node which is introduced at the midpoint of one edge (I,J)
    ! is given the age MAX(AGE(I),AGE(J))+1. For vertices which are
    ! inserted at the center of a quadrilateral, the maximum age of
    ! all four surrounding vertices is adopted and increased by one.
    ! A negative age means, the the vertex is lockes, that is, it
    ! cannot be removed. By definition, vertices of the initial
    ! triangulation cannot be removed and are always locked.
    INTEGER :: h_IvertexAge = ST_NOHANDLE

    ! Pointer to h_IvertexAge.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER, DIMENSION(:), POINTER :: p_IvertexAge => NULL()

    ! Nodal property array.
    ! Handle to
    !       p_InodalProperty = array [1..NVT+NMT] of integer
    ! p_InodalProperty(i) defines for each vertex i=(1..NVT) 
    ! and each edge i=(NVT+1..NVT+NMT) its function inside of the
    ! geometry. Generally said, the range of the p_InodalProperty-array 
    ! characterizes the type of the node (=vertex/edge):
    ! = 0    : The vertex/edge is an inner vertex/edge
    ! > 0    : The vertex/edge is a boundary vertex/edge on the real
    !           boundary. KNPR(.) defines the number of the boundary
    !           component.
    INTEGER :: h_InodalProperty = ST_NOHANDLE

    ! Pointer to h_InodalProperty.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER, DIMENSION(:), POINTER :: p_InodalProperty => NULL()
    
    ! Vertices adjacent to an element.
    ! Handle to 
    !       p_IverticesAtElement = array [1..TRIA_MAXNVE2D,1..NEL] of integer.
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    INTEGER :: h_IverticesAtElement = ST_NOHANDLE

    ! Pointer to h_IverticesAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement => NULL ()

    ! Neighbour elements adjacent to an element.
    ! Handle to
    !       p_IneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element, the numbers of adjacent elements in mathematically
    ! positive sense, metting the element in an edge.
    INTEGER :: h_IneighboursAtElement = ST_NOHANDLE

    ! Pointer to h_IneighboursAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement => NULL ()

    ! Midneighbour elements adjacent to an element.
    ! Handle to
    !       p_ImidneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! H-adaptivity is performed element-by-element. If one element is
    ! refined, then its neighbors need to be informed that there
    ! are possibly two elements adjacent along one edge. This is a
    ! nonconforming state. When time comes to process the other
    ! element (with hanging node), the element knows which elements
    ! are adjacent along the first and the second half of the edge.
    INTEGER :: h_ImidneighboursAtElement = ST_NOHANDLE

    ! Pointer to h_ImidneighboursAtElement.
    ! This array is introduced to increase performance and must
    ! not be modified by the user
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_ImidneighboursAtElement => NULL ()
    
    ! Quadtree storing the nodal coordinates in 2D
    TYPE(t_quadtree) :: rVertexCoordinates2D

    ! Octree storing the nodal coordinates in 2D
    TYPE(t_octree) :: rVertexCoordinates3D
    
    ! Array of binary search trees storing the boundary data
    ! p_IboundaryCpIdx and p_IverticesAtBoundary
    TYPE(t_btree), DIMENSION(:), POINTER :: rBoundary => NULL()

    ! Arraylist for elements-meeting-at-vertex structure
    TYPE(t_arraylist) :: relementsAtVertex
  END TYPE t_hadapt

  !</typeblock> 

!</Types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  INTERFACE add_vertex2D
    MODULE PROCEDURE add_vertex_atEdgeMidpoint2D
    MODULE PROCEDURE add_vertex_atElementCenter2D
  END INTERFACE

  INTERFACE replace_element2D
    MODULE PROCEDURE replace_elementTria
    MODULE PROCEDURE replace_elementQuad
  END INTERFACE

  INTERFACE add_element2D
    MODULE PROCEDURE add_elementTria
    MODULE PROCEDURE add_elementQuad
  END INTERFACE

  INTERFACE update_ElementNeighbors2D
    MODULE PROCEDURE update_ElemEdgeNeighbors2D
    MODULE PROCEDURE update_ElemNeighbors2D
  END INTERFACE

CONTAINS
  
  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_initFromParameterlist(rhadapt,rparlist,ssection)

!<description>
    ! This subroutine initializes the adaptivity structure
    ! with the values supplied by the parameter list
!</descrpition>

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
    CALL parlst_getvalue_int   (rparlist,ssection,"nsubdividemax",rhadapt%nsubdividemax)
    CALL parlst_getvalue_int   (rparlist,ssection,"irefinementStrategy",rhadapt%irefinementStrategy)
    CALL parlst_getvalue_int   (rparlist,ssection,"icoarseningStrategy",rhadapt%icoarseningStrategy)
    CALL parlst_getvalue_double(rparlist,ssection,"drefinementTolerance",rhadapt%drefinementTolerance)
    CALL parlst_getvalue_double(rparlist,ssection,"dcoarseningTolerance",rhadapt%dcoarseningTolerance)

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

    ! Set dimension
    rhadapt%ndim = rtriangulation%ndim

    ! Set coordinates
    SELECT CASE(rhadapt%ndim)
    CASE(NDIM2D)
      CALL hadapt_setVertexCoords2D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE(NDIM3D)
      CALL hadapt_setVertexCoords3D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE DEFAULT
      PRINT *, "hadapt_initFromTriangulation: Invalid spatial dimension!"
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
        rtriangulation%h_IverticesAtBoundary,rtriangulation%h_DvertexParameterValue,&
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
    rtriangulation%ndim = rhadapt%ndim

    ! Get coordinates
    SELECT CASE(rhadapt%ndim)
    CASE(NDIM2D)
      CALL hadapt_getVertexCoords2D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE(NDIM3D)
      CALL hadapt_getVertexCoords3D(rhadapt,&
          rtriangulation%h_DvertexCoords,rtriangulation%NVT)

    CASE DEFAULT
      PRINT *, "hadapt_generateRawMesh: Invalid spatial dimension!"
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
    INTEGER :: ibct

    ! Check if quadtree exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).EQ.HADAPT_HAS_COORDS) THEN
      SELECT CASE(rhadapt%ndim)
      CASE(NDIM2D)
        CALL qtree_releaseQuadtree(rhadapt%rVertexCoordinates2D)
        
      CASE(NDIM3D)
        CALL otree_releaseOctree(rhadapt%rVertexCoordinates3D)

      CASE DEFAULT
        PRINT *, "hadapt_releaseAdaptation: Invalid spatial dimension!"
        CALL sys_halt()
      END SELECT
    END IF   

    ! Check if boundary structure exists
    IF (ASSOCIATED(rhadapt%rBoundary)) THEN
      DO ibct=1,SIZE(rhadapt%rBoundary)
        CALL btree_releaseTree(rhadapt%rBoundary(ibct))
      END DO
      DEALLOCATE(rhadapt%rBoundary)
      NULLIFY(rhadapt%rBoundary)
    END IF

    ! Release elements-meeting-at-vertex arraylist
    CALL arrlst_releaseArraylist(rhadapt%relementsAtVertex)

    ! Release storage which is no longer in use
    IF (rhadapt%h_Imarker /= ST_NOHANDLE)    CALL storage_free(rhadapt%h_Imarker)
    IF (rhadapt%h_IvertexAge /= ST_NOHANDLE) CALL storage_free(rhadapt%h_IvertexAge)
    IF (rhadapt%h_ImidneighboursAtElement /= ST_NOHANDLE)&
        CALL storage_free(rhadapt%h_ImidneighboursAtElement)
    
    ! Nullify "performance-pointers"
    NULLIFY(rhadapt%p_IvertexAge)
    NULLIFY(rhadapt%p_InodalProperty)
    NULLIFY(rhadapt%p_IverticesAtElement)
    NULLIFY(rhadapt%p_IneighboursAtElement)
    NULLIFY(rhadapt%p_ImidneighboursAtElement)
    
    ! Clear parameters
    rhadapt%irefinementStrategy  = HADAPT_NOREFINEMENT
    rhadapt%icoarseningStrategy  = HADAPT_NOCOARSENING
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
  END SUBROUTINE hadapt_releaseAdaptation

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
!</inputoutput
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_QTREEIDX)            :: nnode
    REAL(DP)                          :: xmin,xmax,ymin,ymax

    ! Check if handle is not empty
    IF (h_DvertexCoords == ST_NOHANDLE) THEN
      PRINT *, "hadapt_setVertexCoords2D: Invalid handle!"
      CALL sys_halt()
    END IF

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).EQ.HADAPT_HAS_COORDS) THEN
      CALL qtree_releaseQuadtree(rhadapt%rVertexCoordinates2D)
    END IF
    
    ! Set pointer
    CALL storage_getbase_double2D(h_DvertexCoords,p_DvertexCoords,nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin = MINVAL(p_DvertexCoords(1,:))
    xmax = MAXVAL(p_DvertexCoords(1,:))
    ymin = MINVAL(p_DvertexCoords(2,:))
    ymax = MAXVAL(p_DvertexCoords(2,:))
    
    ! Estimate number of initial quadrilaterals
    nnode = INT(0.5_DP*nvt)

    ! Create quadtree for vertices
    CALL qtree_createQuadtree(rhadapt%rVertexCoordinates2D,nvt,&
        nnode,xmin,ymin,xmax,ymax)
    
    ! Copy vertex coordinates to quadtree
    CALL qtree_copyToQuadtree(p_DvertexCoords,rhadapt%rVertexCoordinates2D)

    ! Set specifier for quadtree
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM2D
    rhadapt%NVT  = nvt
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
!</subroutine>

    ! Check if coordinates exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).NE.HADAPT_HAS_COORDS) THEN
      PRINT *, "hadapt_getVertexCoords2D: quadtree does not exist!"
      CALL sys_halt()
    END IF

    ! Check if coordinates are given in 2D
    IF (rhadapt%ndim /= NDIM2D) THEN
      PRINT *, "hadapt_getVertexCoords2D: invalid spatial dimension!"
      CALL sys_halt()
    END IF

    ! Copy quadtree to handle h_DvertexCoords
    CALL qtree_copyFromQuadtree(rhadapt%rVertexCoordinates2D,h_DvertexCoords)

    ! Set dimension
    IF (PRESENT(ndim)) ndim=rhadapt%ndim
    IF (PRESENT(nvt)) nvt=rhadapt%NVT
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
!</inputoutput
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    INTEGER(PREC_QTREEIDX)            :: nnode
    REAL(DP)                          :: xmin,xmax,ymin,ymax,zmin,zmax

    ! Check if handle is not empty
    IF (h_DvertexCoords == ST_NOHANDLE) THEN
      PRINT *, "hadapt_setVertexCoords3D: Invalid handle!"
      CALL sys_halt()
    END IF

    ! Check if quadtree is already generated, then remove old quadtree/octree first
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).EQ.HADAPT_HAS_COORDS) THEN
      CALL otree_releaseOctree(rhadapt%rVertexCoordinates3D)
    END IF
    
    ! Set pointer
    CALL storage_getbase_double2D(h_DvertexCoords,p_DvertexCoords,nvt)

    ! Get outer bounding-box of vertex coordinates
    xmin = MINVAL(p_DvertexCoords(1,:))
    xmax = MAXVAL(p_DvertexCoords(1,:))
    ymin = MINVAL(p_DvertexCoords(2,:))
    ymax = MAXVAL(p_DvertexCoords(2,:))
    zmin = MINVAL(p_DvertexCoords(3,:))
    zmax = MAXVAL(p_DvertexCoords(3,:))
    
    ! Estimate number of initial quadrilaterals
    nnode = INT(0.5_DP*nvt)

    ! Create octree for vertices
    CALL otree_createOctree(rhadapt%rVertexCoordinates3D,nvt,&
        nnode,xmin,ymin,zmin,xmax,ymax,zmax)
    
    ! Copy vertex coordinates to octree
    CALL otree_copyToOctree(p_DvertexCoords,rhadapt%rVertexCoordinates3D)

    ! Set specifier for quadtree
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_COORDS)

    ! Set dimensions
    rhadapt%ndim = NDIM3D
    rhadapt%NVT  = nvt
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
!</subroutine>

    ! Check if coordinates exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_COORDS).NE.HADAPT_HAS_COORDS) THEN
      PRINT *, "hadapt_getVertexCoords3D: quadtree does not exist!"
      CALL sys_halt()
    END IF

    ! Check if coordinates are given in 2D
    IF (rhadapt%ndim /= NDIM3D) THEN
      PRINT *, "hadapt_getVertexCoords3D: invalid spatial dimension!"
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
    IF (h_IverticesAtElement == ST_NOHANDLE) THEN
      PRINT *, "hadapt_setVerticesAtElement: Invalid handle!"
      CALL sys_halt()
    END IF

    rhadapt%h_IverticesAtElement = h_IverticesAtElement
    CALL storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
        rhadapt%p_IverticesAtElement)
    
    ! Set specifier for IverticesAtElement
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_VERTATELEM)

    ! Set dimensions
    rhadapt%NEL  = nel
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
!</subroutine>

    ! Check if "vertices-at-element" array exists
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_VERTATELEM).NE.HADAPT_HAS_VERTATELEM) THEN
      PRINT *, "hadapt_getVerticesAtElement: structure does not exist!"
      CALL sys_halt()
    END IF

    ! Check if handle needs to be freed first
    IF (h_IverticesAtElement /= ST_NOHANDLE .AND.&
        h_IverticesAtElement /= rhadapt%h_IverticesAtElement)&
        CALL storage_free(h_IverticesAtElement)

    ! Assign handle
    h_IverticesAtElement = rhadapt%h_IverticesAtElement

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
    IF (h_IneighboursAtElement == ST_NOHANDLE) THEN
      PRINT *, "hadapt_setNeighboursAtElement: Invalid handle!"
      CALL sys_halt()
    END IF

    rhadapt%h_IneighboursAtElement = h_IneighboursAtElement
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
      PRINT *, "hadapt_getNeighboursAtElement: structure does not exist!"
      CALL sys_halt()
    END IF

    ! Check if handle needs to be freed first
    IF (h_IneighboursAtElement /= ST_NOHANDLE .AND.&
        h_IneighboursAtElement /= rhadapt%h_IneighboursAtElement)&
        CALL storage_free(h_IneighboursAtElement)

    ! Assign handle
    h_IneighboursAtElement = rhadapt%h_IneighboursAtElement

    ! Set dimensions
    IF(PRESENT(nel))   nel=rhadapt%NEL
  END SUBROUTINE hadapt_getNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_setNelOfType(rhadapt,InelOfType)

!<description>
    ! This subroutine sets the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</desciption>

!<input>
    ! Number of elements with a defined number of vertices per element.
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE), INTENT(IN) :: InelOfType
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    rhadapt%InelOfType = InelOfType

    ! Set specifier for InelOfType
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_NELOFTYPE)
  END SUBROUTINE hadapt_setNelOfType

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_getNelOfType(rhadapt,InelOfType)

!<description>
    ! This subroutine returns the number of elements with a defined number
    ! of vertices per elements, that is, InelOfType(TRIA_MAXNVE)
!</desciption>

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
      PRINT *, "hadapt_getNelOfType: structure does not exist!"
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
    IF ((h_IboundaryCpIdx == ST_NOHANDLE) .OR.&
        (h_IverticesAtBoundary == ST_NOHANDLE) .OR.&
        (h_DvertexParameterValue == ST_NOHANDLE)) THEN
      PRINT *, "hadapt_setNodalProperty: Invalid handle!"
      CALL sys_halt()
    END IF

    ! Check if boundary structure exists and remove it
    IF (ASSOCIATED(rhadapt%rBoundary)) THEN
      DO ibct=1,SIZE(rhadapt%rBoundary)
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
      lvbd = p_IboundaryCpIdx(ibct+1)-p_IboundaryCpIdx(ibct)
      CALL btree_createTree(rhadapt%rBoundary(ibct),lvbd,ST_INT,1,0,2)

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
      CALL btree_copyToTree(p_IverticesAtBoundary(ivbdStart:ivbdEnd),rhadapt%rBoundary(ibct),&
          p_DData = p_DvertexParameterValue2D(:,ivbdStart:ivbdEnd),&
          p_IData = p_IneighboursAtBoundary(:,ivbdStart:ivbdEnd))
    END DO
    
    ! Set dimensions
    rhadapt%NBCT  = nbct
    rhadapt%NVBD  = nvbd
    rhadapt%NVBD0 = nvbd

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
      PRINT *, "hadapt_getBoundary: structure does not exist!"
      CALL sys_halt()
    END IF

    ! Due to the fact, that the boundary data may be stored in multiple 
    ! binary search trees (one for each component) the static arrays 
    ! are pre-allocated and filled step-by-step from the dynamic data
    ! in the boundary search tree(s).

    ! Check if the arrays p_IboundaryCpIdx, p_IverticesAtBoundary and 
    ! p_DvertexParameterValue have correct dimension
    IF (h_IboundaryCpIdx == ST_NOHANDLE) THEN
      CALL storage_new('hadapt_getBoundary','p_IboundaryCpIdx',&
          rhadapt%NBCT+1,ST_INT,h_IboundaryCpIdx,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_IboundaryCpIdx,isize)
      IF (isize < rhadapt%NBCT+1) THEN
        CALL storage_realloc('hadapt_getBoundary',rhadapt%NBCT+1,&
            h_IboundaryCpIdx,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF

    IF (h_IverticesAtBoundary == ST_NOHANDLE) THEN
      CALL storage_new('hadapt_getBoundary','p_IverticesAtBoundary',&
          rhadapt%NVBD,ST_INT,h_IverticesAtBoundary,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_IverticesAtBoundary,isize)
      IF (isize < rhadapt%NVBD) THEN
        CALL storage_realloc('hadapt_getBoundary',rhadapt%NVBD,&
            h_IverticesAtBoundary,ST_NEWBLOCK_NOINIT,.FALSE.)
      END IF
    END IF

    IF (h_DvertexParameterValue == ST_NOHANDLE) THEN
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
    p_IboundaryCpIdx(1) = 1
    ioff = 0

    DO ibct=1,rhadapt%NBCT
      
      ! Set subdimensions
      lvbd      = rhadapt%rBoundary(ibct)%NA
      ivbdStart = ioff+1
      ivbdEnd   = ioff+lvbd
      ioff      = ioff+lvbd

      ! Set index for next boundary component
      p_IboundaryCpIdx(ibct+1) = ioff+1

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
    IF (h_InodalProperty == ST_NOHANDLE) THEN
      PRINT *, "hadapt_setNodalProperty: Invalid handle!"
      CALL sys_halt()
    END IF

    rhadapt%h_InodalProperty = h_InodalProperty
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
      PRINT *, "hadapt_getNodalProperty: structure does not exist!"
      CALL sys_halt()
    END IF

    ! Check if handle needs to be freed first
    IF (h_InodalProperty /= ST_NOHANDLE .AND.&
        h_InodalProperty /= rhadapt%h_InodalProperty)&
        CALL storage_free(h_InodalProperty)

    ! Assign handle
    h_InodalProperty = rhadapt%h_InodalProperty
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
      PRINT *, "hadapt_genElementsAtVertex: structure does not exist!"
      CALL sys_halt()
    END IF

    ! Create arraylist
    CALL arrlst_createArraylist(rhadapt%relementsAtVertex,&
        2*rhadapt%NVT,8*rhadapt%NEL,ST_INT,ARRAYLIST_UNORDERED)

    ! Fill arraylist
    DO iel=1,rhadapt%NEL
      DO ive=1,get_NVE(rhadapt,iel)
        CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,&
            rhadapt%p_IverticesAtElement(ive,iel),iel,ipos)
      END DO
    END DO
    
    ! Set specifier for relementsAtVertex
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_HAS_ELEMATVERTEX)
  END SUBROUTINE hadapt_genElementsAtVertex

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_performAdaptation(rhadapt,rindicator,rcollection,fcb_hadaptCallback)

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
    TYPE(t_vectorScalar), INTENT(IN)         :: rindicator

    ! callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)            :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_VERTEXIDX)    :: nvt
    INTEGER(PREC_ELEMENTIDX)   :: nel
    
    ! Check if dynamic data structures are available
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA) THEN
      PRINT *, "hadapt_performAdaptation: dynamic data structures are not generated!"
      CALL sys_halt()
    END IF

    ! Initialize initial dimensions
    CALL storage_getsize(rhadapt%h_IverticesAtElement,Isize)
    rhadapt%NELMAX      = Isize(2)
    CALL storage_getsize(rhadapt%h_IneighboursAtElement,Isize)
    rhadapt%NELMAX      = MIN(rhadapt%NELMAX,Isize(2))

    rhadapt%InelOfType0 = rhadapt%InelOfType
    rhadapt%NVT0        = rhadapt%NVT
    rhadapt%NEL0        = rhadapt%NEL
    rhadapt%NVBD0       = rhadapt%NVBD
    rhadapt%increaseNVT = 0
    
    ! What kind of grid refinement should be performed
    SELECT CASE(rhadapt%irefinementStrategy)
    CASE (HADAPT_NOREFINEMENT)   ! No grid refinement

      
    CASE (HADAPT_REDGREEN)   ! Red-green grid refinement

      ! Mark elements for refinement based on indicator function
      CALL mark_refinement2D(rhadapt,rindicator)

      ! Mark additional elements to restore conformity
      CALL redgreen_mark_refinement2D(rhadapt,rcollection,fcb_hadaptCallback)

      ! Mark element for recoarsening based on indicator function
      CALL redgreen_mark_coarsening2D(rhadapt,rindicator)

      ! Compute new dimensions
      nvt = rhadapt%NVT+rhadapt%increaseNVT
      nel = NumberOfElements(rhadapt)

      ! Adjust vertex age array and nodal property array
      CALL storage_realloc('hadapt_performAdaptation',nvt,&
          rhadapt%h_IvertexAge,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_realloc('hadapt_performAdaptation',nvt,&
          rhadapt%h_InodalProperty,ST_NEWBLOCK_NOINIT,.TRUE.)
     
      ! Adjust elemental arrays
      CALL storage_realloc('hadapt_performAdaptation',nel,&
          rhadapt%h_IverticesAtElement,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_realloc('hadapt_performAdaptation',nel,&
          rhadapt%h_IneighboursAtElement,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_realloc('hadapt_performAdaptation',nel,&
          rhadapt%h_ImidneighboursAtElement,ST_NEWBLOCK_NOINIT,.TRUE.)

      ! Reset pointers
      CALL storage_getbase_int(rhadapt%h_IvertexAge,rhadapt%p_IvertexAge)
      CALL storage_getbase_int(rhadapt%h_InodalProperty,&
          rhadapt%p_InodalProperty)
      CALL storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
          rhadapt%p_IverticesAtElement)
      CALL storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
          rhadapt%p_IneighboursAtElement)
      CALL storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
          rhadapt%p_ImidneighboursAtElement)

      ! Adjust dimension of solution vector
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
          CALL fcb_hadaptCallback(rcollection,&
          HADAPT_OPR_ADJUSTVERTEXDIM,(/nvt/),(/0/))
      
      CALL hadapt_writegridsvg(rhadapt,"mygrid",ibset(0,0))
      PAUSE

      ! Perform refinement
      CALL redgreen_refine(rhadapt,rcollection,fcb_hadaptCallback)

      ! Perform coarsening
      CALL redgreen_coarsen(rhadapt,rcollection,fcb_hadaptCallback)

      CALL hadapt_writegridsvg(rhadapt,"mygrid",ibset(0,0))
      PAUSE

      ! Adjust nodal property array
      CALL storage_realloc('hadapt_performAdaptation',rhadapt%NVT,&
          rhadapt%h_InodalProperty,ST_NEWBLOCK_NOINIT,.TRUE.)
      

    CASE (HADAPT_LONGESTEDGE)   ! Bisection of Longest edge

      ! Mark elements for refinement based on indicator function
      CALL mark_refinement2D(rhadapt,rindicator)
      
      ! Perform refinement
      CALL bisection_refine(rhadapt)

      
    CASE DEFAULT
      PRINT *, "hadapt_performAdaptation: Unsupported refinement strategy!"
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
      
      CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
      
      ! Initialize number of elements by current number
      nel = rhadapt%NEL0
      
      ! Loop over all elements and check marker
      DO iel=1,rhadapt%NEL0
        
        SELECT CASE(p_Imarker(iel))
        CASE(MARK_REF_TRIA2TRIA_1,MARK_REF_TRIA2TRIA_2,MARK_REF_TRIA2TRIA_3,&
             MARK_REF_QUAD2QUAD_13,MARK_REF_QUAD2QUAD_24)
          nel=nel+1

        CASE(MARK_REF_QUAD3TRIA_1,MARK_REF_QUAD3TRIA_2,MARK_REF_QUAD3TRIA_3,&
             MARK_REF_QUAD3TRIA_4,MARK_REF_TRIA3TRIA_12,MARK_REF_TRIA3TRIA_23,&
             MARK_REF_TRIA3TRIA_13)
          nel=nel+2
          
        CASE(MARK_REF_TRIA4TRIA,MARK_REF_QUAD4QUAD,MARK_REF_QUAD4TRIA_12,&
             MARK_REF_QUAD4TRIA_23,MARK_REF_QUAD4TRIA_34,MARK_REF_QUAD4TRIA_14)
          nel=nel+3
        END SELECT
      END DO
    END FUNCTION NumberOfElements
  END SUBROUTINE hadapt_performAdaptation

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE hadapt_info(rhadapt)

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
    
    WRITE(*,FMT='(2X,72("+"))')
    WRITE(*,FMT='(2X,A)') 'Grid adaptivity:'
    WRITE(*,FMT='(2X,72("-"))')
    WRITE(*,FMT='(2X,A,T52,I15)') 'Total number of grid refinement steps:', rhadapt%nRefinementSteps
    WRITE(*,FMT='(2X,A,T52,I15)') 'Total number of grid coarsening steps:', rhadapt%nCoarseningSteps
    WRITE(*,FMT='(2X,A,T52,I15)') 'Total number of grid smoothing steps:',  rhadapt%nSmoothingSteps
    WRITE(*,*)
    WRITE(*,FMT='(2X,A,T52,I3)')  'Strategy for grid refinement:', rhadapt%irefinementStrategy
    WRITE(*,FMT='(2X,A,T52,I3)')  'Strategy for grid coarsening:', rhadapt%icoarseningStrategy
    WRITE(*,*)
    WRITE(*,FMT='(2X,A,T40,I8,1X,A,1X,I8)') 'Number of elements:',          rhadapt%NEL,'initially',rhadapt%NEL0
    WRITE(*,FMT='(2X,A,T40,I8,1X,A,1X,I8)') 'Number of vertices:',          rhadapt%NVT,'initially',rhadapt%NVT0
    WRITE(*,FMT='(2X,A,T40,I8,1X,A,1X,I8)') 'Number of boundary vertices:', rhadapt%NVBD,'initially',rhadapt%NVBD0
    WRITE(*,*)
    
    WRITE(*,FMT='(2X,A)') 'Handles:'
    WRITE(*,FMT='(2X,72("-"))')
    WRITE(*,FMT='(2X,A,T52,I15)') 'h_Imarker',                rhadapt%h_Imarker
    WRITE(*,FMT='(2X,A,T52,I15)') 'h_IvertexAge',             rhadapt%h_IvertexAge
    WRITE(*,FMT='(2X,A,T52,I15)') 'h_InodalProperty',         rhadapt%h_InodalProperty
    WRITE(*,FMT='(2X,A,T52,I15)') 'h_IverticesAtElement',     rhadapt%h_IverticesAtElement
    WRITE(*,FMT='(2X,A,T52,I15)') 'h_IneighboursAtElement',   rhadapt%h_IneighboursAtElement
    WRITE(*,FMT='(2X,A,T52,I15)') 'h_ImidneighboursAtElement',rhadapt%h_ImidneighboursAtElement
    WRITE(*,*)

    WRITE(*,FMT='(2X,A)') 'Coordinates:'
    WRITE(*,FMT='(2X,72("-"))')
    SELECT CASE(rhadapt%ndim)
    CASE(NDIM2D)
      CALL qtree_infoQuadtree(rhadapt%rVertexCoordinates2D)
      WRITE(*,*)
    CASE(NDIM3D)
      CALL otree_infoOctree(rhadapt%rVertexCoordinates3D)
      WRITE(*,*)
    END SELECT

    WRITE(*,FMT='(2X,A)') 'Boundary:'
    WRITE(*,FMT='(2X,72("-"))')
    DO ibct=1,rhadapt%NBCT
      WRITE(*,FMT='(2X,A,I2,A)') 'Boundary component',ibct,':'
      WRITE(*,FMT='(2X,72("-"))')
      CALL btree_infoTree(rhadapt%rBoundary(ibct))
    END DO
    WRITE(*,*)

    WRITE(*,FMT='(2X,A)') 'Element-at-Vertex:'
    WRITE(*,FMT='(2X,72("-"))')
    CALL arrlst_infoArraylist(rhadapt%relementsAtVertex)
    
    WRITE(*,FMT='(2X,72("+"))')
    WRITE(*,*)
  END SUBROUTINE hadapt_info

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hadapt_writeGridSVG(rhadapt,coutputFile,cinfo,width,height)

!<description>
    ! This subroutine outputs the current state of the adapted grid stored
    ! in the dynamic data structure to a given file in SVG format.
!</description>

!<input>
    ! Output file name w/o suffix .svg
    CHARACTER(LEN=*), INTENT(IN)     :: coutputFile

    ! Format bitvector: determined additional output
    ! BIT0 : output vertices
    ! BIT1 : output vertex numbers
    ! BIT2 : output vertex age
    ! BIT3 : output element numbers
    ! BIT4 : output element state
    ! BIT5 : output element marker
    ! BIT6 : output adjacency numbers
    INTEGER, INTENT(IN) :: cinfo

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
    INTEGER, PARAMETER :: defaultwidth  = 3000
    INTEGER, PARAMETER :: defaultheight = 3000
    INTEGER, PARAMETER :: xoff = 100
    INTEGER, PARAMETER :: yoff = 100

    CHARACTER(LEN=*), PARAMETER :: font_family     = "Verdana"
    CHARACTER(LEN=*), PARAMETER :: font_size       = "64"
    CHARACTER(LEN=*), PARAMETER :: font_size_small = "14"
    
    ! local variables
    REAL(DP), DIMENSION(2*NDIM2D) :: bbox
    REAL(DP) :: x0,y0,xdim,ydim,xscale,yscale,xmid,ymid
    INTEGER :: xsize,ysize
    INTEGER, DIMENSION(:), POINTER                     :: p_Imarker
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE2D) :: Kadj,Kmidadj
    INTEGER(PREC_VERTEXIDX),  DIMENSION(TRIA_MAXNVE2D) :: Kvert
    INTEGER(PREC_VERTEXIDX)  :: ivt
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: iunit,ive,nve,istate
    INTEGER, SAVE :: iout=0
    
    ! Check if dynamic data structures generated
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA) THEN
      PRINT *, "adapt_output_svg: dynamic data structures are not generated"
      CALL sys_halt()
    END IF
    
    ! Increment the sample number
    iout=iout+1
    
    ! Open output file for writing
    CALL io_openFileForWriting(TRIM(ADJUSTL(coutputFile))//'.'//&
        TRIM(sys_siL(iout,5))//'.svg',iunit,SYS_REPLACE,bformatted=.TRUE.)
    
    ! Write XML-header to file
    WRITE(iunit,FMT='(A)') '<?xml version="1.0" encoding="utf-8" standalone="no"?>'
    WRITE(iunit,FMT='(A)') '<!-- Created with Featflow2 (http://www.featflow.de/) -->'
    WRITE(iunit,FMT='(A)') '<svg'
    WRITE(iunit,FMT='(A)') ' xmlns:svg="http://www.w3.org/2000/svg"'
    WRITE(iunit,FMT='(A)') ' xmlns="http://www.w3.org/2000/svg"'
    WRITE(iunit,FMT='(A)') ' version="1.0"'
    
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
    bbox   = qtree_getBoundingBox(rhadapt%rVertexCoordinates2D)
    xdim   = bbox(3)-bbox(1); x0 = bbox(1)
    ydim   = bbox(4)-bbox(2); y0 = bbox(2)
    xscale = xsize*MIN(1._DP,xdim/ydim)/xdim
    yscale = ysize*MIN(1._DP,xdim/ydim)/ydim

    ! Set height and width of image
    WRITE(iunit,FMT='(A)') ' width="100%" height="100%" viewBox="0 0 '//&
        TRIM(sys_siL(INT(xscale*xdim)+2*xoff,10))//' '//&
        TRIM(sys_siL(INT(yscale*ydim)+2*yoff,10))//'" xml:space="preserve">'
       
    !---------------------------------------------------------------------------
    ! Output all elements
    !---------------------------------------------------------------------------
    IF (IAND(rhadapt%iSpec,HADAPT_MARKEDREFINE).EQ.HADAPT_MARKEDREFINE .OR.&
        IAND(rhadapt%iSpec,HADAPT_MARKEDCOARSEN).EQ.HADAPT_MARKEDCOARSEN) THEN

      ! Set pointer to marker
      CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)

      ! Output elements and color those which are marked for refinement
      DO iel=1,rhadapt%NEL
        
        ! Get number of vertices per elements
        nve = get_NVE(rhadapt,iel)

        ! Get local data
        Kvert(1:nve)   = rhadapt%p_IverticesAtElement(1:nve,iel)
        Kadj(1:nve)    = rhadapt%p_IneighboursAtElement(1:nve,iel)
        Kmidadj(1:nve) = rhadapt%p_ImidneighboursAtElement(1:nve,iel)
        
        IF (iel <= rhadapt%NEL0) THEN

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
          CASE(MARK_CRS_4QUAD2QUAD_24:MARK_CRS_4TRIA1TRIA)
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="yellow" stroke="black" stroke-width="1"'

          ! Unknown element marker
          CASE DEFAULT
            WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
                '" fill="hotpink" stroke="black" stroke-width="1"'
          END SELECT

        ELSE

          ! For all new elements no marker is available
          WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
              '" fill="white" stroke="black" stroke-width="1"'
        END IF
               
        ! Each element is a polygon made up from 3/4 points
        WRITE(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        DO ive=1,nve
          xdim = qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(ive))-x0
          ydim = qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(ive))-y0
          WRITE(iunit,FMT='(A)',ADVANCE='NO') &
              TRIM(sys_siL(xoff+INT(xscale*xdim),9))//','//&
              TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//' '
        END DO
        
        ! And of course, the polygone must be closed !
        xdim = qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1))-x0
        ydim = qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1))-y0
        WRITE(iunit,FMT='(A)') &
            TRIM(sys_siL(xoff+INT(xscale*xdim),9))//','//&
            TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//'"/>'
        
        ! Output elemental data
        IF (btest(cinfo,3) .OR.&
            btest(cinfo,4) .OR.&
            btest(cinfo,5)) THEN
          
          ! Determine state of element
          SELECT CASE(nve)
          CASE(TRIA_NVETRI2D)
            xdim   = SUM(qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1:3)))/3._DP-x0
            ydim   = SUM(qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1:3)))/3._DP-y0
            istate = redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))

          CASE(TRIA_NVEQUAD2D)
            xdim   = SUM(qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1:4)))/4._DP-x0
            ydim   = SUM(qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1:4)))/4._DP-y0
            istate = redgreen_getStateQuad(rhadapt%p_IvertexAge(Kvert(1:4)))
          END SELECT
          
          ! What data should be written?
          IF (btest(cinfo,3)) THEN
            
            ! Write element number
            WRITE(iunit,FMT='(A)') '<text id="nel'//TRIM(sys_siL(iel,9))//&
                '" x="'//TRIM(sys_siL( xoff+INT(xscale*xdim),9))//&
                '" y="'//TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size//'" fill="black">'//&
                TRIM(sys_siL(iel,9))//'</text>'
            
          ELSEIF (btest(cinfo,4)) THEN
            
            ! Write element state
            WRITE(iunit,FMT='(A)') '<text id="nel'//TRIM(sys_siL(iel,9))//&
                '" x="'//TRIM(sys_siL( xoff+INT(xscale*xdim),9))//&
                '" y="'//TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size//'" fill="black">'//&
                TRIM(sys_siL(istate,3))//'</text>'
            
          ELSEIF (btest(cinfo,5)) THEN
            
            ! Write element marker
            WRITE(iunit,FMT='(A)') '<text id="nel'//TRIM(sys_siL(iel,9))//&
                '" x="'//TRIM(sys_siL( xoff+INT(xscale*xdim),9))//&
                '" y="'//TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size//'" fill="black">'//&
                TRIM(sys_siL(p_Imarker(iel),3))//'</text>'
          END IF
        END IF
        
        ! Output adjacency information?
        IF (btest(cinfo,6)) THEN
          xmid = SUM(qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1:nve)))/REAL(nve,DP)-x0
          ymid = SUM(qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1:nve)))/REAL(nve,DP)-y0
          
          DO ive=1,nve
            xdim = (2*qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(ive))+xmid)/3._DP
            ydim = (2*qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(ive))+ymid)/3._DP
            WRITE(iunit,FMT='(A)') ' <text x="'//TRIM(sys_siL(xoff+INT(xscale*xdim),9))//&
                '" y="'//TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size_small//&
                '" text-anchor="middle" fill="sienne">'//&
                TRIM(sys_siL(Kadj(ive),9))//","//TRIM(sys_siL(Kmidadj(ive),9))//'</text>'
          END DO
        END IF
      END DO
      
    ELSE
      
      ! Output only elements without individual coloring
      DO iel=1,rhadapt%NEL

        ! Get number of vertices per element
        nve = get_NVE(rhadapt,iel)
        
        ! Get local data
        Kvert(1:nve)   = rhadapt%p_IverticesAtElement(1:nve,iel)
        Kadj(1:nve)    = rhadapt%p_IneighboursAtElement(1:nve,iel)
        Kmidadj(1:nve) = rhadapt%p_ImidneighboursAtElement(1:nve,iel)
        
        ! Start new element ...
        WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
            '" fill="white" stroke="black" stroke-width="1" '
        WRITE(iunit,FMT='(A)',ADVANCE='NO') ' points="'
        
        ! Output vertices of element if required
        IF (btest(cinfo,0)) THEN
          DO ive=1,nve
            xdim = qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(ive))-x0
            ydim = qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(ive))-y0
            WRITE(iunit,FMT='(A)',ADVANCE='NO') &
                TRIM(sys_siL(xoff+INT(xscale*xdim),9))//','//&
                TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//' '
          END DO
          
          ! The first vertex is also the last vertex
          xdim = qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1))-x0
          ydim = qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1))-y0
          WRITE(iunit,FMT='(A)') &
              TRIM(sys_siL(xoff+INT(xscale*xdim),9))//','//&
              TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//'"/>'
        END IF

        ! Output elemental data
        IF (btest(cinfo,3) .OR.&
            btest(cinfo,4)) THEN
          
          ! Determine state of element
          SELECT CASE(nve)
          CASE(TRIA_NVETRI2D)
            xdim   = SUM(qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1:3)))/3._DP-x0
            ydim   = SUM(qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1:3)))/3._DP-y0
            istate = redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))
          
          CASE(TRIA_NVEQUAD2D)
            xdim   = SUM(qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1:4)))/4._DP-x0
            ydim   = SUM(qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1:4)))/4._DP-y0
            istate = redgreen_getStateQuad(rhadapt%p_IvertexAge(Kvert(1:4)))
          END SELECT
          
          ! What data should be written?
          IF (btest(cinfo,3)) THEN
            
            ! Write element number
            WRITE(iunit,FMT='(A)') '<text id="nel'//TRIM(sys_siL(iel,9))//&
                '" x="'//TRIM(sys_siL( xoff+INT(xscale*xdim),9))//&
                '" y="'//TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size//'" fill="black">'//&
                TRIM(sys_siL(iel,9))//'</text>'
            
          ELSEIF(btest(cinfo,4)) THEN

            ! Write element state
            WRITE(iunit,FMT='(A)') '<text id="nel'//TRIM(sys_siL(iel,9))//&
                '" x="'//TRIM(sys_siL( xoff+INT(xscale*xdim),9))//&
                '" y="'//TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size//'" fill="black">'//&
                TRIM(sys_siL(istate,3))//'</text>'
          END IF
        END IF
        
        ! Output adjacency information?
        IF (btest(cinfo,6)) THEN
          xmid = SUM(qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(1:nve)))/REAL(nve,DP)-x0
          ymid = SUM(qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(1:nve)))/REAL(nve,DP)-y0
          
          DO ive=1,nve
            xdim = (2*qtree_getX(rhadapt%rVertexCoordinates2D,Kvert(ive))+xmid)/3._DP
            ydim = (2*qtree_getY(rhadapt%rVertexCoordinates2D,Kvert(ive))+ymid)/3._DP
            WRITE(iunit,FMT='(A)') ' <text x="'//TRIM(sys_siL(xoff+INT(xscale*xdim),9))//&
                '" y="'//TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size_small//&
                '" text-anchor="middle" fill="sienne">'//&
                TRIM(sys_siL(Kadj(ive),9))//","//TRIM(sys_siL(Kmidadj(ive),9))//'</text>'
          END DO
        END IF
      END DO
    END IF
    
    ! Output nodal data
    IF (btest(cinfo,0) .OR.&
        btest(cinfo,1) .OR.&
        btest(cinfo,2)) THEN
      
      ! Loop over all vertices
      DO ivt=1,qtree_getsize(rhadapt%rVertexCoordinates2D)
       
        xdim = qtree_getX(rhadapt%rVertexCoordinates2D,ivt)-x0
        ydim = qtree_getY(rhadapt%rVertexCoordinates2D,ivt)-y0
        
        ! Write vertices as points?
        IF (btest(cinfo,0)) THEN
          IF (rhadapt%p_IvertexAge(ivt) > 0) THEN
            WRITE(iunit,FMT='(A)') '<circle cx="'//&
                TRIM(sys_siL(xoff+INT(xscale*xdim),9))//'" cy="'//&
                TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" r="15"  fill="white" stroke="black" stroke-width="3"/>'
          ELSE
            WRITE(iunit,FMT='(A)') '<circle cx="'//&
                TRIM(sys_siL(xoff+INT(xscale*xdim),9))//'" cy="'//&
                TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" r="15"  fill="black" stroke="black" stroke-width="1"/>'
          END IF
        END IF
        
        IF (btest(cinfo,1)) THEN
          ! Write vertex number?
          WRITE(iunit,FMT='(A)') '<text x="'//&
              TRIM(sys_siL(xoff+INT(xscale*xdim),9))//'" y="'//&
              TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
              '" font-family="'//font_family//'" font-size="'//font_size//'" fill="black">'//&
              TRIM(sys_siL(ivt,9))//'</text>'
          
        ELSEIF (btest(cinfo,2)) THEN
          ! Write vertex age?
          IF (rhadapt%p_IvertexAge(ivt) > 0) THEN
            WRITE(iunit,FMT='(A)') '<text x="'//&
                TRIM(sys_siL(xoff+INT(xscale*xdim),9))//'" y="'//&
                TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size//'" fill="black">'//&
                TRIM(sys_siL(rhadapt%p_IvertexAge(ivt),3))//'</text>'
          ELSE
            WRITE(iunit,FMT='(A)') '<text x="'//&
                TRIM(sys_siL(xoff+INT(xscale*xdim),9))//'" y="'//&
                TRIM(sys_siL(ysize+yoff-INT(yscale*ydim),9))//&
                '" font-family="'//font_family//'" font-size="'//font_size//'" fill="fuchsia">'//&
                TRIM(sys_siL(-rhadapt%p_IvertexAge(ivt),3))//'</text>'
          END IF
        END IF
      END DO
    END IF
    
    ! Close XML-file
    WRITE(iunit,FMT='(A)') '</svg>'
    CLOSE(iunit)
  END SUBROUTINE hadapt_writeGridSVG

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

  ! Internal auxiliary routines and functions. No user should ever look at them

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<function>

  PURE FUNCTION get_NVE(rhadapt,iel) RESULT(nve)

!<description>
    ! This function returns the number of vertices present in the given element
!</description>

!<input>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(IN)           :: rhadapt

    ! Number of the element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
!</input>

!<result>
    ! number of vertices per element
    INTEGER :: nve
!</result>
!</function

    ! Are we in 2D or 3D?
    IF (rhadapt%ndim == NDIM2D) THEN
      
      ! Do we have quadrilaterals in the triangulation?
      IF (rhadapt%InelOfType(TRIA_NVEQUAD2D).EQ.0) THEN

        ! There are no quadrilaterals in the current triangulation.
        ! Hence, return TRIA_NVETRI2D by default.
        nve = TRIA_NVETRI2D
        
      ELSE
        
        ! There are quadrilaterals and possible also triangles in
        ! the current triangulatin. If the last entry of the
        ! vertices-at-element list is nonzero then TRIA_NVEQUAD2D vertices
        ! are present in the current element. Otherwise return TRIA_NVETRI2D.
        IF (rhadapt%p_IverticesAtElement(TRIA_NVEQUAD2D,iel).EQ.0) THEN
          nve = TRIA_NVETRI2D
        ELSE
          nve = TRIA_NVEQUAD2D
        END IF
        
      END IF

    ELSE
      
      
    END IF
  END FUNCTION get_NVE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE add_vertex_atEdgeMidpoint2D(rhadapt,i1,i2,e1,i12,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine adds a new vertex at the midpoint of a given egde. 
    ! First, the coordinates of the new vertex are computed and added to
    ! the quadtree. If the new vertex is located at the boundary then its
    ! parametrization is computed and its normal vector is stored in
    ! the boundary data structure.
!</description>

!<input>
    ! First point of edge on which new vertex will be added
    INTEGER(PREC_VERTEXIDX), INTENT(IN)  :: i1

    ! Second point of edge on which new vertex will be added
    INTEGER(PREC_VERTEXIDX), INTENT(IN)  :: i2

    ! Number of the right-adjacent element w.r.t. to the oriented edge (I1,I2)
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e1

    ! Callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>

!<output>
    ! Number of the new vertex located between i1 and i2
    INTEGER(PREC_VERTEXIDX), INTENT(OUT)        :: i12
!</output>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(NDIM2D) :: Dcoord
    REAL(DP)                    :: x1,y1,x2,y2,dvbdp1,dvbdp2
    INTEGER(PREC_QTREEIDX)      :: inode
    INTEGER(PREC_TREEIDX)       :: ipred,ipos
    INTEGER                     :: ibct
    
    ! Get coordinates of edge vertices
    x1 = qtree_getX(rhadapt%rVertexCoordinates2D,i1)
    y1 = qtree_getY(rhadapt%rVertexCoordinates2D,i1)
    x2 = qtree_getX(rhadapt%rVertexCoordinates2D,i2)
    y2 = qtree_getY(rhadapt%rVertexCoordinates2D,i2)

    ! Compute coordinates of new vertex 
    Dcoord = 0.5_DP*(/x1+x2,y1+y2/)

    ! Search for vertex coordinates in quadtree: 
    ! If the vertex already exists, e.g., it was added when the adjacent element
    ! was refined, then nothing needs to be done for this vertex
    IF (qtree_searchInQuadtree(rhadapt%rVertexCoordinates2D,&
        Dcoord,inode,ipos,i12) == QTREE_FOUND) RETURN
    
    ! Otherwise, update number of vertices
    rhadapt%NVT = rhadapt%NVT+1
    i12         = rhadapt%NVT
    
    ! Set age of vertex
    rhadapt%p_IvertexAge(i12) = &
        1+MAX(ABS(rhadapt%p_IvertexAge(i1)),ABS(rhadapt%p_IvertexAge(i2)))
    
    ! Set nodal property
    IF (e1 == 0) THEN
      rhadapt%p_InodalProperty(i12) = rhadapt%p_InodalProperty(i1)
    ELSE
      rhadapt%p_InodalProperty(i12) = 0
    END IF
    
    ! Add new entry to vertex coordinates
    CALL qtree_insertIntoQuadtree(rhadapt%rVertexCoordinates2D,i12,Dcoord,inode)
    
    ! Are we at the boundary?
    IF (e1 == 0) THEN
      ! Increment number of boundary nodes
      rhadapt%NVBD = rhadapt%NVBD+1
      
      ! Get number of boundary component
      ibct = rhadapt%p_InodalProperty(i1)
      
      ! Get parameter values of the boundary nodes
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i1,ipred) == BTREE_NOT_FOUND) THEN
        PRINT *, "add_vertex: Unable to find vertex in boudary data structure!"
        CALL sys_halt()
      END IF
      ipos   = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      dvbdp1 = rhadapt%rBoundary(ibct)%DData(BdrValue,ipos)
      
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i2,ipred) == BTREE_NOT_FOUND) THEN
        PRINT *, "add_vertex: Unable to find vertex in boudary data structure!"
        CALL sys_halt()
      END IF
      ipos   = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      dvbdp2 = rhadapt%rBoundary(ibct)%DData(BdrValue,ipos)
      
      ! If I2 is last(=first) node on boundary component IBCT round DVBDP2 to next integer
      IF (dvbdp2 <= dvbdp1) dvbdp2=CEILING(dvbdp1)
      
      ! Add new entry to boundary structure
      CALL btree_insertIntoTree(rhadapt%rBoundary(ibct),i12,&
          Idata=(/i1,i2/),Ddata=(/0.5_DP*(dvbdp1+dvbdp2)/))
    END IF
      
    ! Optionally, invoke callback function
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,&
        HADAPT_OPR_INSERTVERTEXEDGE,(/i12,i1,i2/),(/0/))
  END SUBROUTINE add_vertex_atEdgeMidpoint2D

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE add_vertex_atElementCenter2D(rhadapt,i1,i2,i3,i4,i5,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine adds a new vertex at the center of a given quadtrilateral.
    ! First, the coordinates of the new vertex computed and added to the
    ! quadtree. The new vertex cannot be located at the boundary.
!</description>

!<input>
    ! Four corners of the quadrilateral
    INTEGER(PREC_VERTEXIDX), INTENT(IN)  :: i1,i2,i3,i4

    ! Callback routine
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>

!<output>
    ! Number of the new vertex
    INTEGER(PREC_VERTEXIDX), INTENT(OUT)        :: i5
!</output>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(NDIM2D) :: Dcoord
    REAL(DP) :: x1,y1,x2,y2,x3,y3,x4,y4,x21,y21,x31,y31,x24,y24,alpha
    INTEGER(PREC_QTREEIDX) :: inode,ipos
    
    ! Compute coordinates of new vertex
    x1 = qtree_getX(rhadapt%rVertexCoordinates2D,i1)
    y1 = qtree_getY(rhadapt%rVertexCoordinates2D,i1)
    x2 = qtree_getX(rhadapt%rVertexCoordinates2D,i2)
    y2 = qtree_getY(rhadapt%rVertexCoordinates2D,i2)
    x3 = qtree_getX(rhadapt%rVertexCoordinates2D,i3)
    y3 = qtree_getY(rhadapt%rVertexCoordinates2D,i3)
    x4 = qtree_getX(rhadapt%rVertexCoordinates2D,i4)
    y4 = qtree_getY(rhadapt%rVertexCoordinates2D,i4)

    x21 = x2-x1; x31 = x3-x1; x24 = x2-x4
    y21 = y2-y1; y31 = y3-y1; y24 = y2-y4
    alpha = (x21*y24-y21*x24)/(x31*y24-y31*x24)
    
    Dcoord = (/x1+alpha*x31,y1+alpha*y31/)

    ! Search for vertex coordinates in quadtree
    IF (qtree_searchInQuadtree(rhadapt%rVertexCoordinates2D,&
        Dcoord,inode,ipos,i5) == QTREE_NOT_FOUND) THEN
      
      ! Update number of vertices
      rhadapt%NVT = rhadapt%NVT+1
      i5             = rhadapt%NVT
      
      ! Set age of vertex
      rhadapt%p_IvertexAge(I5) = &
          1+MAX(ABS(rhadapt%p_IvertexAge(i1)),&
                ABS(rhadapt%p_IvertexAge(i2)),&
                ABS(rhadapt%p_IvertexAge(i3)),&
                ABS(rhadapt%p_IvertexAge(i4)))

      ! Set nodal property
      rhadapt%p_InodalProperty(i5) = 0
      
      ! Add new entry to DCORVG
      CALL qtree_insertIntoQuadtree(rhadapt%rVertexCoordinates2D,i5,Dcoord,inode)
    END IF
    
    ! Optionally, invoke callback function
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,&
        HADAPT_OPR_INSERTVERTEXCENTR,(/i5,i1,i2,i3,i4/),(/0/))
  END SUBROUTINE add_vertex_atElementCenter2D

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE remove_vertex2D(rhadapt,ivt,ivtReplace)
  
!<description>
    ! This subroutine removes an existing vertex from the adaptivity structure
    ! and moves the last vertex at its position. The number of the replacement
    ! vertex is returned as ivtReplace. If the vertex to be replace is the last
    ! vertex then ivtReplace=0 is returned on output.
!</description>

!<input>
    ! Number of the vertex to be deleted
    INTEGER(PREC_VERTEXIDX), INTENT(IN) :: ivt

    
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>

!<output>
    ! Number of the vertex to replace the deleted one
    INTEGER(PREC_VERTEXIDX), INTENT(OUT) :: ivtReplace
!</output>
!</subroutine>

    ! local variables
    INTEGER(PREC_VERTEXIDX) :: i1,i2
    INTEGER(PREC_TREEIDX)   :: ipred,ipos
    INTEGER :: ibct

    ! Get last vertex and decrease number of vertices by one
    ivtReplace  = rhadapt%NVT
    rhadapt%NVT = rhadapt%NVT-1
    
    ! If IVT is a boundary node remove it from the boundary and
    ! connect its boundary neighbors with each other
    ibct = rhadapt%p_InodalProperty(ivt)
    
    ! Are we at the boundary?
    IF (ibct .NE. 0) THEN
      
      ! Find position of vertex IVT in boundary array
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),ivt,ipred) .EQ.&
          BTREE_NOT_FOUND) THEN
        PRINT *, "remove_vertex: Unable to find vertex in boundary data structure"
        CALL sys_halt()
      END IF
      ipos = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      
      ! Get the two boundary neighbors: I1 <- IVT -> I2
      i1 = rhadapt%rBoundary(ibct)%IData(BdrPrev,ipos)
      i2 = rhadapt%rBoundary(ibct)%IData(BdrNext,ipos)
      
      ! Connect the boundary neighbors with each other: I1 <=> I2
      ! First, set I2 as next neighboring of I1
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i1,ipred) .EQ.&
          BTREE_NOT_FOUND) THEN
        PRINT *, "remove_vertex: Unable to find vertex in boundary data structure"
        CALL sys_halt()
      END IF
      ipos = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      rhadapt%rBoundary(ibct)%IData(BdrNext,ipos) = i2
      
      ! Second, set I1 as previous neighbor of I2
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i2,ipred) .EQ.&
          BTREE_NOT_FOUND) THEN
        PRINT *, "remove_vertex: Unable to find vertex in boundary data structure"
        CALL sys_halt()
      END IF
      ipos = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      rhadapt%rBoundary(ibct)%IData(BdrPrev,ipos) = i1
      
      ! And finally, delete IVT from the boundary
      IF (btree_deleteFromTree(rhadapt%rBoundary(ibct),ivt) .EQ.&
          BTREE_NOT_FOUND) THEN
        PRINT *, "remove_vertex: Unable to delete vertex from boundary data structure"
        CALL sys_halt()
      END IF
    END IF
    
    ! If IVT is not the last node then copy the data for vertex IVTREPLACE
    ! to IVT and prepare IVTREPLACE for elimination
    IF (ivt < ivtReplace) THEN
      
      ! If IVTREPLACE is a boundary node then remove IVTREPLACE from the boundary
      ! vector, insert IVT into the boundary vector instead, and connect the
      ! boundary neighbors of IVTREPLACE with IVT
      ibct = rhadapt%p_InodalProperty(ivtReplace)

      ! Are we at the boundary?
      IF (ibct .NE. 0) THEN

        IF (btree_searchInTree(rhadapt%rBoundary(ibct),ivtReplace,ipred) .EQ.&
            BTREE_NOT_FOUND) THEN
          PRINT *, "remove_vertex: Unable to find vertex in boundary data structure"
          CALL sys_halt()
        END IF
        ipos = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        
        ! Insert IVT into the boundary vector
        CALL btree_insertIntoTree(rhadapt%rBoundary(ibct),ivt,&
            Idata=rhadapt%rBoundary(ibct)%IData(:,ipos),&
            Ddata=rhadapt%rBoundary(ibct)%DData(:,ipos))
        
        ! Get the two boundary neighbors: I1 <- IVTREPLACE -> I2
        i1 = rhadapt%rBoundary(ibct)%IData(BdrPrev,ipos)
        i2 = rhadapt%rBoundary(ibct)%IData(BdrNext,ipos)
        
        ! Connect the boundary neighbors with IVT: I1 <- IVT -> I2
        ! First, set IVT as next neighbor of I1
        IF (btree_searchInTree(rhadapt%rBoundary(ibct),i1,ipred) .EQ.&
            BTREE_NOT_FOUND) THEN
          PRINT *, "remove_vertex: Unable to find vertex in boundary data structure"
          CALL sys_halt()
        END IF
        ipos = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        rhadapt%rBoundary(ibct)%IData(BdrNext,ipos) = ivt
        
        ! Second, set IVT as previous neighbor of I2
        IF (btree_searchInTree(rhadapt%rBoundary(ibct),i2,ipred) .EQ.&
            BTREE_NOT_FOUND) THEN
          PRINT *, "remove_vertex: Unable to find vertex in boundary data structure"
          CALL sys_halt()
        END IF
        ipos = rhadapt%rBoundary(ibct)%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        rhadapt%rBoundary(ibct)%IData(BdrPrev,ipos) = ivt
        
        ! Finally, delete IVTREPLACE from the boundary
        IF (btree_deleteFromTree(rhadapt%rBoundary(ibct),ivtReplace) .EQ.&
            BTREE_NOT_FOUND) THEN
          PRINT *, "remove_vertex: Unable to delete vertex from the boundary data structure"
          CALL sys_halt()
        END IF
      END IF
      
      ! Copy data from node IVTREPLACE to node IVT
      rhadapt%p_InodalProperty(ivt) = rhadapt%p_InodalProperty(ivtReplace)
      rhadapt%p_IvertexAge(ivt)     = rhadapt%p_IvertexAge(ivtReplace)
    ELSE

      ! IVT is the last vertex of the adaptivity structure
      ivtReplace = 0
    END IF
  END SUBROUTINE remove_vertex2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE replace_elementTria(rhadapt,ipos,i1,i2,i3,e1,e2,e3,e4,e5,e6)
  
!<description>
    ! This subroutine replaces the vertices and adjacent elements for
    ! a given element
!</description>

!<input>
    ! position number of the element in dynamic data structure
    INTEGER(PREC_TREEIDX), INTENT(IN)    :: ipos

    ! numbers of the element nodes
    INTEGER(PREC_VERTEXIDX), INTENT(IN)  :: i1,i2,i3

    ! numbers of the surrounding elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e1,e2,e3

    ! numbers of the surrounding mid-elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e4,e5,e6
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)  :: rhadapt
!</inputoutput>
!</subroutine

    ! Replace triangular element
    rhadapt%p_IverticesAtElement(:,ipos)      = (/i1,i2,i3,0/)
    rhadapt%p_IneighboursAtElement(:,ipos)    = (/e1,e2,e3,0/)
    rhadapt%p_ImidneighboursAtElement(:,ipos) = (/e4,e5,e6,0/)    
  END SUBROUTINE replace_elementTria

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE replace_elementQuad(rhadapt,ipos,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)
  
!<description>
    ! This subroutine replaces the vertices and adjacent elements for
    ! a given element
!</description>

!<input>
    ! position number of the element in dynamic data structure
    INTEGER(PREC_TREEIDX), INTENT(IN)    :: ipos

    ! numbers of the element nodes
    INTEGER(PREC_VERTEXIDX), INTENT(IN)  :: i1,i2,i3,i4

    ! numbers of the surrounding elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e1,e2,e3,e4

    ! numbers of the surrounding mid-elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e5,e6,e7,e8
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)     :: rhadapt
!</inputoutput>
!</subroutine

    ! Replace quadrilateral element
    rhadapt%p_IverticesAtElement(:,ipos)      = (/i1,i2,i3,i4/)
    rhadapt%p_IneighboursAtElement(:,ipos)    = (/e1,e2,e3,e4/)
    rhadapt%p_ImidneighboursAtElement(:,ipos) = (/e5,e6,e7,e8/)
  END SUBROUTINE replace_elementQuad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE add_elementTria(rhadapt,i1,i2,i3,e1,e2,e3,e4,e5,e6)

!<description>
    ! This subroutine adds a new element connected to three vertices 
    ! and surrounded by three adjacent elements
!</description>

!<input>
    ! numbers of the element nodes
    INTEGER(PREC_VERTEXIDX), INTENT(IN)  :: i1,i2,i3

    ! numbers of the surrounding elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e1,e2,e3

    ! numbers of the surrounding mid-elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e4,e5,e6
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)     :: rhadapt
!</inputoutput>
!</subroutine
    
    ! Increase number of elements and number of triangles
    rhadapt%NEL = rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)+1

    rhadapt%p_IverticesAtElement(:,rhadapt%NEL)      = (/i1,i2,i3,0/)
    rhadapt%p_IneighboursAtElement(:,rhadapt%NEL)    = (/e1,e2,e3,0/)
    rhadapt%p_ImidneighboursAtElement(:,rhadapt%NEL) = (/e4,e5,e6,0/)
  END SUBROUTINE add_elementTria

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE add_elementQuad(rhadapt,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)

!<description>
    ! This subroutine adds a new element connected to four vertices 
    ! and surrounded by four adjacent elements
!</description>

!<input>
    ! numbers of the element nodes
    INTEGER(PREC_VERTEXIDX), INTENT(IN)  :: i1,i2,i3,i4

    ! numbers of the surrounding elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e1,e2,e3,e4

    ! numbers of the surrounding mid-elements
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e5,e6,e7,e8
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)     :: rhadapt
!</inputoutput>
!</subroutine
    
    ! Increase number of elements and number of quadrilaterals
    rhadapt%NEL   = rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+1

    rhadapt%p_IverticesAtElement(:,rhadapt%NEL)      = (/i1,i2,i3,i4/)
    rhadapt%p_IneighboursAtElement(:,rhadapt%NEL)    = (/e1,e2,e3,e4/)
    rhadapt%p_ImidneighboursAtElement(:,rhadapt%NEL) = (/e5,e6,e7,e8/)
  END SUBROUTINE add_elementQuad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE remove_element2D(rhadapt,iel,ielReplace)
  
!<description>
    ! This subroutine removes an existing element and moves the last
    ! element of the adaptation data structure to its position.
    ! The routine returns the former number ielReplace of the last 
    ! element. If iel is the last element, then ielReplace=0 is returned.
!</description>

!<input>
    ! Number of the element that should be removed
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)  :: iel
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)      :: rhadapt
!</inputoutput>

!<output>
    ! Former number of the replacement element
    INTEGER(PREC_ELEMENTIDX), INTENT(OUT) :: ielReplace
!</output>
!</subroutine>

    ! local variables
    INTEGER(PREC_ELEMENTIDX) :: jel
    INTEGER :: ive,jve
    
    ! Replace element by the last element and delete last element
    ielReplace = rhadapt%NEL

    ! Which kind of element are we?
    SELECT CASE(get_NVE(rhadapt,iel))
    CASE(TRIA_NVETRI2D)
      rhadapt%InelOfType(TRIA_NVETRI2D) = rhadapt%InelOfType(TRIA_NVETRI2D)-1

    CASE(TRIA_NVEQUAD2D)
      rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)-1

    CASE DEFAULT
      PRINT *, "remove_element2D: Invalid element type!"
      CALL sys_halt()
    END SELECT

    ! Check if we are the last element
    IF (iel .NE. ielReplace) THEN
      
      ! Element is not the last one.
      IF (iel > ielReplace) THEN
        PRINT *, "remove_element2D: Invalid element!"
        CALL sys_halt()
      END IF

      ! Copy data from element ielReplaced to element iel
      rhadapt%p_IverticesAtElement(:,iel)      = rhadapt%p_IverticesAtElement(:,ielReplace)
      rhadapt%p_IneighboursAtElement(:,iel)    = rhadapt%p_IneighboursAtElement(:,ielReplace)
      rhadapt%p_ImidneighboursAtElement(:,iel) = rhadapt%p_ImidneighboursAtElement(:,ielReplace)

      ! Update the adjacency list of surrounding elements
      DO ive=1,get_NVE(rhadapt,iel)
        
        ! Get element number of adjacent element jel
        jel = rhadapt%p_IneighboursAtElement(ive,iel)
        
        ! Are we at the boundary?
        IF (jel .EQ. 0) CYCLE

        ! Find replacement element in adjacency list of element 
        ! iel and update entry to new element 
        DO jve=1,get_NVE(rhadapt,jel)
          IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. ielReplace) THEN
            rhadapt%p_IneighboursAtElement(jve,jel)    = iel
            rhadapt%p_ImidneighboursAtElement(jve,jel) = iel
            EXIT
          END IF
        END DO
      END DO
      
      ! Delete replacement element
      rhadapt%p_IverticesAtElement(:,ielReplace)      = 0
      rhadapt%p_IneighboursAtElement(:,ielReplace)    = 0
      rhadapt%p_ImidneighboursAtElement(:,ielReplace) = 0
      
    ELSE

      ! Delete element
      rhadapt%p_IverticesAtElement(:,iel)      = 0
      rhadapt%p_IneighboursAtElement(:,iel)    = 0
      rhadapt%p_ImidneighboursAtElement(:,iel) = 0

      ! Element iel is the last element
      ielReplace = 0
    END IF

    ! Decrease number of elements
    rhadapt%NEL = rhadapt%NEL-1
  END SUBROUTINE remove_element2D

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE update_ElemEdgeNeighbors2D(rhadapt,e0adj,e0midadj,e0,eadj,emidadj)

!<description>
    ! This subroutine updates the list of elements adjecent to another 
    ! element and the list of elements mid-adjacent to another element.
    !
    ! The situation is as follows:
    !
    !  +---------------------+            +---------------------+
    !  |                     |            |          .          |
    !  |          E0         |            | EMIDADJ  .   EADJ   |
    !  |                     |            |          .          |
    !  |                     |            |          |          |
    !  |                     |            |          |          |
    !  +---------------------+    --->    +----------+----------+
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  |  E0ADJ   . E0MIDADJ |            | E0ADJ    . E0MIDADJ |
    !  |          .          |            |          .          |
    !  +---------------------+            +---------------------+
    !
    ! "Sitting" on element E0 we want to update the elements lists:
    !
    ! 1) If E0 is located at the boundary then nothing needs to be done.
    ! 2) If the adjacent element has not been subdivided, that is,
    !    E0ADJ == E0MIDADJ then it suffices to update its adjacency
    !    list for the entry E0 adopting the values EMIDADJ and EADJ.
    ! 3) If the adjacent element has been subdivided, that is,
    !    E0ADJ /= E0MIDADJ then the adjacency list if each of these two 
    !    elements is updated by the value EMIDADJ and EADJ, respectively.
!</description>

!<input>
    ! Number of the neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e0adj

    ! Number of the mid-neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e0midadj

    ! Number of the updated macro-element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: e0

    ! Number of the new neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: eadj

    ! Number of the new mid-neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: emidadj
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER :: ive,nve
    
    ! Do nothing for elements adjacent to the boundary
    IF (e0adj .EQ. 0 .OR. e0midadj .EQ. 0) RETURN

    ! Check if adjacent and mid-adjacent elements are the same.
    IF (e0adj .EQ. e0midadj) THEN

      ! Case 1: Adjacent element has not been subdivided, that is, the 
      !         current edge contains a hanging node for the adjacent element.
      
      ! What kind of element is neighboring element?
      nve = get_NVE(rhadapt,e0adj)
      
      ! Perform operations for elements/vertices?
      IF (e0 < 0) THEN
        
        ! Check for vertices. Hence, loop over all entries in the
        ! list of vertices for element E0ADJ and check if the value E0
        ! is present. If this is the case, then update the corrseponding
        ! entries in the lists of (mid-)adjacent element neighbors.
        DO ive=1,nve
          IF (rhadapt%p_IverticesAtElement(ive,e0adj) .EQ. e0) THEN
            rhadapt%p_IneighboursAtElement(ive,e0adj)    = eadj
            rhadapt%p_ImidneighboursAtElement(ive,e0adj) = emidadj
            EXIT
          END IF
        END DO

      ELSE
        
        ! Check for elements. Hence, loop over all entries in the
        ! list of adjacent elements for element E0ADJ and check if the value E0
        ! is present. If this is the case, then update the corrseponding
        ! entries in the lists of (mid-)adjacent element neighbors.
        DO ive=1,nve
          IF (rhadapt%p_IneighboursAtElement(ive,e0adj).EQ.e0) THEN
            rhadapt%p_IneighboursAtElement(ive,e0adj)    = eadj
            rhadapt%p_ImidneighboursAtElement(ive,e0adj) = emidadj
            EXIT
          END IF
        END DO
        
      END IF
                  
    ELSE
      
      ! Case 2: Adjacent element has already been subdivided.
      
      ! What kind of element is neighboring element 
      nve = get_NVE(rhadapt,e0adj)

      ! Perform operations for elements/vertices?
      IF (e0 < 0) THEN
        
        ! Check for vertices. Hence, loop over all entries in the
        ! list of vertices for element E0ADJ and check if the value E0
        ! is present. If this is the case, then update the corrseponding
        ! entries in the lists of (mid-)adjacent element neighbors.
        DO ive=1,nve
          IF (rhadapt%p_IverticesAtElement(ive,e0adj) .EQ. e0) THEN
            rhadapt%p_IneighboursAtElement(ive,e0adj)    = emidadj
            rhadapt%p_ImidneighboursAtElement(ive,e0adj) = emidadj
            EXIT
          END IF
        END DO

      ELSE
        
        ! Check for elements. Hence, loop over all entries in the
        ! list of adjacent elements for element E0ADJ and check if the value E0
        ! is present. If this is the case, then update the corrseponding
        ! entries in the lists of (mid-)adjacent element neighbors.
        DO ive=1,nve
          IF (rhadapt%p_IneighboursAtElement(ive,e0adj).EQ.e0) THEN
            rhadapt%p_IneighboursAtElement(ive,e0adj)    = emidadj
            rhadapt%p_ImidneighboursAtElement(ive,e0adj) = emidadj
            EXIT
          END IF
        END DO
        
      END IF
      
      ! What kind of element is neighboring element 
      nve = get_NVE(rhadapt,e0midadj)
      
      ! Perform operations for elements/vertices?
      IF (e0 < 0) THEN
        
        ! Check for vertices. Hence, loop over all entries in the
        ! list of vertices for element E0MIDADJ and check if the value E0
        ! is present. If this is the case, then update the corrseponding
        ! entries in the lists of (mid-)adjacent element neighbors.
        DO ive=1,nve
          IF (rhadapt%p_IverticesAtElement(ive,e0midadj) .EQ. e0) THEN
            rhadapt%p_IneighboursAtElement(ive,e0midadj)    = eadj
            rhadapt%p_ImidneighboursAtElement(ive,e0midadj) = eadj
            EXIT
          END IF
        END DO

      ELSE
        
        ! Check for elements. Hence, loop over all entries in the
        ! list of adjacent elements for element E0MIDADJ and check if the value E0
        ! is present. If this is the case, then update the corrseponding
        ! entries in the lists of (mid-)adjacent element neighbors.
        DO ive=1,nve
          IF (rhadapt%p_IneighboursAtElement(ive,e0midadj).EQ.e0) THEN
            rhadapt%p_IneighboursAtElement(ive,e0midadj)    = eadj
            rhadapt%p_ImidneighboursAtElement(ive,e0midadj) = eadj
            EXIT
          END IF
        END DO
        
      END IF

    END IF
  END SUBROUTINE update_ElemEdgeNeighbors2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE update_ElemNeighbors2D(rhadapt,iel0,iel)

!<description>
    ! This subroutine updates the list of elements adjacent to another elements.
    ! For all elements jel which are adjacent to the old element iel0 the new 
    ! value iel is stored in the neighbours-at-element structure.
!</description>

!<input>
    ! Number of the element to be updated
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel0

    ! New value of the element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ELEMENTIDX) :: jel
    INTEGER :: ive,jve

    ! Loop over adjacent elements
    adjacent: DO ive=1,get_NVE(rhadapt,iel0)
      ! Get number of adjacent element
      jel = rhadapt%p_IneighboursAtElement(ive,iel0)

      ! Are we at the boundary?
      IF (jel .EQ. 0) CYCLE adjacent
      
      ! Find position of element IEL0 in adjacent element JEL
      DO jve=1,get_NVE(rhadapt,jel)
        IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. iel0) THEN
          rhadapt%p_IneighboursAtElement(jve,jel) = iel
          CYCLE adjacent
        END IF
      END DO
      
      ! If the old element number was not found in adjacent element JEL
      ! then something went wrong and we should not proceed.
      PRINT *, "update_ElemNeighbors2D: Inconsistent adjacency lists!"
      CALL sys_halt()
    END DO adjacent
  END SUBROUTINE update_ElemNeighbors2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Tria2Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into two triangular 
    ! elements by subdividing one edge.
    !
    ! In the illustration, i1,i2,i3 and i4 denote the vertices whereby i4 
    ! is the new vertex. Moreover, (e1)-(e6) stand for the element numbers
    ! which are adjecent and mid-adjacent to current element.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                            i3
    !            +                             +
    !           / \                           /|\
    !     (e3) /   \ (e5)               (e3) / | \ (e5)
    !         /     \                       /  |  \
    !        /       \          ->         /   |   \
    !       /   iel   \                   /    |    \
    ! (e6) /           \ (e2)       (e6) / iel |nel+1\ (e2)
    !     /*            \               /*     |     *\
    !    +---------------+             +-------+-------+
    !   i1 (e1)     (e4) i2            i1 (e1) i4 (e4) i2
    !
!</description>

!<input>
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
    
    ! Identifiert for element marker
    INTEGER, INTENT(IN)                  :: imarker

    ! Callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4
    INTEGER :: loc1,loc2,loc3
    
    ! Find local position of edge to be subdivided
    SELECT CASE(imarker)
    CASE(2)
      loc1=1; loc2=2; loc3=3

    CASE(4)
      loc1=2; loc2=3; loc3=1

    CASE(8)
      loc1=3; loc2=1; loc3=2

    CASE DEFAULT
      PRINT *, "refine_Tria2Tria2: Invalid marker",imarker
      CALL sys_halt()
    END SELECT
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1,iel)
    i2 = rhadapt%p_IverticesAtElement(loc2,iel)
    i3 = rhadapt%p_IverticesAtElement(loc3,iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1,iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2,iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3,iel)
    e4 = rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e5 = rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc3,iel)
    
    ! Store values before refinement
    nel0 = rhadapt%NEL
    
    ! Add one new vertex I4 at the midpoint of edge (I1,I2).
    CALL add_vertex2D(rhadapt,i1,i2,e1,i4,rcollection,fcb_hadaptCallback)
    
    ! Replace element IEL and add one new element E4
    CALL replace_element2D(rhadapt,iel,i1,i4,i3,e1,nel0+1,e3,e1,nel0+1,e6)
    CALL add_element2D(rhadapt,i2,i3,i4,e2,iel,e4,e5,iel,e4)
    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,nel0+1,nel0+1)
    
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Tria2Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel,ipos)

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA2TRIA,&
        (/i1,i2,i3,i4/),(/e1,e2,e3,e4,e5,e6/))
  END SUBROUTINE refine_Tria2Tria

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE refine_Tria3Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into three triangular 
    ! elements by subdividing the longest edge and connecting the new vertex 
    ! with the opposite midpoint.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           / \                       /|\
    !     (e3) /   \ (e5)           (e3) / |*\ (e5)
    !         /     \                   /  |ne\
    !        /       \          ->     /   |l+2+
    !       /   iel   \               /    |  / \
    ! (e6) /           \ (e2)   (e6) / iel | /nel\ (e2)
    !     /*            \           /*     |/ +1 *\
    !    +---------------+         +-------+-------+
    !   i1 (e1)     (e4) i2       i1 (e1)     (e4) i2
    !
!</description>

!<input> 
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
    
    ! Identifier for element marker
    INTEGER, INTENT(IN)                  :: imarker

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP) :: dlen12,dlen23,x,y
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5
    INTEGER :: loc1,loc2,loc3
    
    ! Find corresponding edges according to convention to be subdivided
    SELECT CASE(imarker)
    CASE(6)
      loc1=1; loc2=2; loc3=3

    CASE(12)
      loc1=2; loc2=3; loc3=1

    CASE(10)
      loc1=3; loc2=1; loc3=2

    CASE DEFAULT
      PRINT *, "refine_Tria3Tria: Invalid marker",imarker
      CALL sys_halt()
    END SELECT

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1,iel)
    i2 = rhadapt%p_IverticesAtElement(loc2,iel)
    i3 = rhadapt%p_IverticesAtElement(loc3,iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1,iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2,iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3,iel)
    e4 = rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e5 = rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc3,iel)
    
    ! Store values before refinement
    nel0 = rhadapt%NEL
    
    ! Add two new vertices I4 and I5 at the midpoint of edges (I1,I2) and (I2,I3).
    CALL add_vertex2D(rhadapt,i1,i2,e1,i4,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i2,i3,e2,i5,rcollection,fcb_hadaptCallback)
    
    ! Compute the length of edges (I1,I2) and (I2,I3)
    x = qtree_getX(rhadapt%rVertexCoordinates2D,i1)-&
        qtree_getX(rhadapt%rVertexCoordinates2D,i2)
    y = qtree_getY(rhadapt%rVertexCoordinates2D,i1)-&
        qtree_getY(rhadapt%rVertexCoordinates2D,i2)
    dlen12=SQRT(x*x+y*y)

    x = qtree_getX(rhadapt%rVertexCoordinates2D,i2)-&
        qtree_getX(rhadapt%rVertexCoordinates2D,i3)
    y = qtree_getY(rhadapt%rVertexCoordinates2D,i2)-&
        qtree_getY(rhadapt%rVertexCoordinates2D,i3)
    dlen23=SQRT(x*x+y*y)
    
    IF (dlen12 > dlen23) THEN
      
      ! 1st CASE: longest edge is (I1,I2)
      
      ! Replace element IEL and add two new elements E4, and E5
      CALL replace_element2D(rhadapt,iel,i1,i4,i3,e1,nel0+2,e3,e1,nel0+2,e6)
      CALL add_element2D(rhadapt,i2,i5,i4,e2,nel0+2,e4,e2,nel0+2,e4)
      CALL add_element2D(rhadapt,i3,i4,i5,iel,nel0+1,e5,iel,nel0+1,e5)
      
      ! Update list of nieghboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+1,iel)
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,nel0+2,nel0+1)
            
      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "refine_Tria3Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel,ipos)

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)

      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
          CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA3TRIA12,&
          (/i1,i2,i3,i4,i5/),(/e1,e2,e3,e4,e5/))
      
    ELSE
      
      ! 2nd CASE: longest edge is (I2,I3)
      
      ! Replace element IEL and add two new elements E4 and E5
      CALL replace_element2D(rhadapt,iel,i1,i5,i3,nel0+2,e5,e3,nel0+2,e5,e6)
      CALL add_element2D(rhadapt,i2,i5,i4,e2,nel0+2,e4,e2,nel0+2,e4)
      CALL add_element2D(rhadapt,i1,i4,i5,e1,nel0+1,iel,e1,nel0+1,iel)
      
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+1,nel0+2)
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,iel,nel0+1)
      
      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "refine_Tria3Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
      
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,nel0+2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)

      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
          CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA3TRIA23,&
          (/i1,i2,i3,i4,i5/),(/e1,e2,e3,e4,e5/))
    END IF
  END SUBROUTINE refine_Tria3Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Tria4Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into four similar 
    ! triangular elements by connecting the three edge midpoints.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           / \                       /*\
    !     (e3) /   \ (e5)           (e3) /   \ (e5)
    !         /     \                   /nel+2\
    !        /       \          ->  i6 +-------+ i5
    !       /   iel   \               / \nel+3/ \
    ! (e6) /           \ (e2)   (e6) /   \   /nel\ (e2)
    !     /*            \           /*iel \*/ +1 *\
    !    +---------------+         +-------+-------+
    !   i1 (e1)     (e4) i2       i1 (e1)  i4 (e4) i2
    !
!</description>

!<input>
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i2 = rhadapt%p_IverticesAtElement(2,iel)
    i3 = rhadapt%p_IverticesAtElement(3,iel)

    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e2 = rhadapt%p_IneighboursAtElement(2,iel)
    e3 = rhadapt%p_IneighboursAtElement(3,iel)
    e4 = rhadapt%p_ImidneighboursAtElement(1,iel)
    e5 = rhadapt%p_ImidneighboursAtElement(2,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(3,iel)
        
    ! Store values before refinement
    nel0 = rhadapt%NEL
    
    ! Add three new vertices I4,I5 and I6 at the midpoint of edges (I1,I2),
    ! (I2,I3) and (I1,I3), respectively. 
    CALL add_vertex2D(rhadapt,i1,i2,e1,i4,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i2,i3,e2,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i1,e3,i6,rcollection,fcb_hadaptCallback)
    
    ! Replace element IEL and add three new elements E4, E5, and E6
    CALL replace_element2D(rhadapt,iel,i1,i4,i6,e1,nel0+3,e6,e1,nel0+3,e6)
    CALL add_element2D(rhadapt,i2,i5,i4,e2,nel0+3,e4,e2,nel0+3,e4)
    CALL add_element2D(rhadapt,i3,i6,i5,e3,nel0+3,e5,e3,nel0+3,e5)
    CALL add_element2D(rhadapt,i4,i5,i6,nel0+1,nel0+2,iel,nel0+1,nel0+2,iel)
    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,nel0+2,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel,iel,nel0+2)

    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Tria4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Tria4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+2,ipos)

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+3,ipos)

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA4TRIA,&
        (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6/))
  END SUBROUTINE refine_Tria4Tria

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE refine_Quad2Quad(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element 
    ! into two quadrilateral elements
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i3          i4       i6      i3
    !      +---------------+            +-------+-------+
    !      |               |            |       |      *|
    ! (e4) |               | (e6)  (e4) |       |       | (e6)
    !      |               |            |       |       |
    !      |     iel       |    ->      |  iel  | nel+1 |
    !      |               |            |       |       |
    ! (e8) |               | (e2)  (e8) |       |       | (e2)
    !      |*              |            |*      |       |
    !      +---------------+            +-------+-------+
    !     i1 (e1)     (e5) i2          i1       i5      i2
    ! 
!</description>

!<input>
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Identifier for element marker
    INTEGER, INTENT(IN)                  :: imarker

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>    

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6
    INTEGER :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    SELECT CASE(imarker)
    CASE(11)
      loc1=1; loc2=2; loc3=3; loc4=4

    CASE(21)
      loc1=2; loc2=3; loc3=4; loc4=1

    CASE DEFAULT
      PRINT *, "refine_Quad2Quad: Invalid marker",imarker
      CALL sys_halt()
    END SELECT
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1,iel)
    i2 = rhadapt%p_IverticesAtElement(loc2,iel)
    i3 = rhadapt%p_IverticesAtElement(loc3,iel)
    i4 = rhadapt%p_IverticesAtElement(loc4,iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1,iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2,iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3,iel)
    e4 = rhadapt%p_IneighboursAtElement(loc4,iel)
    e5 = rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e7 = rhadapt%p_ImidneighboursAtElement(loc3,iel)
    e8 = rhadapt%p_ImidneighboursAtElement(loc4,iel)
    
    ! Store values before refinement
    nel0 = rhadapt%NEL
    
    ! Add two new vertices I5 and I6 at the midpoint of edges (I1,I2) and
    ! (I3,I4), respectively.
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i4,e3,i6,rcollection,fcb_hadaptCallback)
    
    ! Replace element IEL and add one new element E5
    CALL replace_element2D(rhadapt,iel,i1,i5,i6,i4,e1,nel0+1,e7,e4,e1,nel0+1,e7,e8)
    CALL add_element2D(rhadapt,i3,i6,i5,i2,e3,iel,e5,e2,e3,iel,e5,e6)
    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel,nel0+1,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel,nel0+1)
    
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad2Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad2Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    
    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD2QUAD,&
        (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6,e7,e8/))
  END SUBROUTINE refine_Quad2Quad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Quad3Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element
    ! into three triangular elements
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7)     (e3) i3
    !      +----------------+            +---------------+
    !      |                |            |\             /|
    ! (e4) |                | (e6)  (e4) | \           / | (e6)
    !      |                |            |  \  nel+2  /  |
    !      |      iel       |     ->     |   \       /   |
    !      |                |            |    \     /    |
    ! (e8) |                | (e2)  (e8) |     \   / nel | (e2)
    !      |*               |            |* iel \*/  +1 *| 
    !      +----------------+            +-------+-------+
    !     i1 (e1)      (e5) i2          i1 (e1)  i5 (e5) i2
    !
!</description>

!<input>
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Identifier for element marker
    INTEGER, INTENT(IN)                  :: imarker

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)     :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5
    INTEGER :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    SELECT CASE(imarker)
    CASE(3)
      loc1=1; loc2=2; loc3=3; loc4=4

    CASE(5)
      loc1=2; loc2=3; loc3=4; loc4=1

    CASE(9)
      loc1=3; loc2=4; loc3=1; loc4=2

    CASE(17)
      loc1=4; loc2=1; loc3=2; loc4=3

    CASE DEFAULT
      PRINT *, "refine_Quad3Tria: Invalid marker",imarker
      CALL sys_halt()
    END SELECT
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1,iel)
    i2 = rhadapt%p_IverticesAtElement(loc2,iel)
    i3 = rhadapt%p_IverticesAtElement(loc3,iel)
    i4 = rhadapt%p_IverticesAtElement(loc4,iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1,iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2,iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3,iel)
    e4 = rhadapt%p_IneighboursAtElement(loc4,iel)
    e5 = rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e7 = rhadapt%p_ImidneighboursAtElement(loc3,iel)
    e8 = rhadapt%p_ImidneighboursAtElement(loc4,iel)

    ! Store values before refinement
    nel0 = rhadapt%NEL
    
    ! Add one new vertex I5 it the midpoint of edge (I1,I2)
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)
    
    ! Replace element IEL and add two new elements E5 and E6
    CALL replace_element2D(rhadapt,iel,i1,i5,i4,e1,nel0+2,e4,e1,nel0+2,e8)
    CALL add_element2D(rhadapt,i2,i3,i5,e2,nel0+2,e5,e6,nel0+2,e5)
    CALL add_element2D(rhadapt,i5,i3,i4,nel0+1,e3,iel,nel0+1,e7,iel)    
    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel,nel0+1,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,nel0+2,nel0+2)
    
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad3Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad3Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+2,ipos)

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD3TRIA,&
        (/i1,i2,i3,i4,i5/),(/e1,e2,e3,e4,e5,e6,e7,e8/))
  END SUBROUTINE refine_Quad3Tria
  
! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Quad4Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element
    ! into four triangular elements
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e6)      (e3) i3          i4 (e6)     (e3) i3
    !      +----------------+            +---------------+
    !      |                |            |\\\\          *|
    ! (e4) |                | (e6)  (e4) | \*  \\\ nel+2 | (e6)
    !      |                |            |  \      \\\   | 
    !      |      iel       |     ->     |   \         \\+i6
    !      |                |            |    \ nel+3  / |
    ! (e8) |                | (e2)  (e8) | iel \     /nel| (e2)
    !      |*               |            |*     \  / +1 *|
    !      +----------------+            +-------+-------+
    !     i1 (e1)      (e5) i2          i1 (e1)  i5 (e5) i2
    !
!</description>

!<input>
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Identifier for element marker
    INTEGER, INTENT(IN)                  :: imarker

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6
    INTEGER :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    SELECT CASE(imarker)
    CASE(7)
      loc1=1; loc2=2; loc3=3; loc4=4

    CASE(13)
      loc1=2; loc2=3; loc3=4; loc4=1

    CASE(25)
      loc1=3; loc2=4; loc3=1; loc4=2

    CASE(19)
      loc1=4; loc2=1; loc3=2; loc4=3

    CASE DEFAULT
      PRINT *, "refine_Quad4Tria: Invalid marker",imarker
      CALL sys_halt()
    END SELECT

    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(loc1,iel)
    i2 = rhadapt%p_IverticesAtElement(loc2,iel)
    i3 = rhadapt%p_IverticesAtElement(loc3,iel)
    i4 = rhadapt%p_IverticesAtElement(loc4,iel)

    e1 = rhadapt%p_IneighboursAtElement(loc1,iel)
    e2 = rhadapt%p_IneighboursAtElement(loc2,iel)
    e3 = rhadapt%p_IneighboursAtElement(loc3,iel)
    e4 = rhadapt%p_IneighboursAtElement(loc4,iel)
    e5 = rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e7 = rhadapt%p_ImidneighboursAtElement(loc3,iel)
    e8 = rhadapt%p_ImidneighboursAtElement(loc4,iel)
    
    ! Store values before refinement
    nel0 = rhadapt%NEL
    
    ! Add new vertices I5 and I6 at the midpoint of edges (I1,I2) and
    ! (I2,I3), respectively.
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i2,i3,e2,i6,rcollection,fcb_hadaptCallback)
    
    ! Replace element IEL and add three new elements E5,E6 and E7
    CALL replace_element2D(rhadapt,iel,i1,i5,i4,e1,nel0+3,e4,e1,nel0+3,e8)
    CALL add_element2D(rhadapt,i2,i6,i5,e2,nel0+3,e5,e2,nel0+3,e5)
    CALL add_element2D(rhadapt,i3,i4,i6,e3,nel0+3,e6,e7,nel0+3,e6)
    CALL add_element2D(rhadapt,i4,i5,i6,iel,nel0+1,nel0+2,iel,nel0+1,nel0+2)
    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel,nel0+2,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,nel0+2,nel0+2)
    
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+2,ipos)

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+3,ipos)

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD4TRIA,&
        (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6,e7,e8/))
  END SUBROUTINE refine_Quad4Tria
    
    ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE refine_Quad4Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element
    ! four similar quadrilateral elements
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7) i7  (e3) i3
    !      +----------------+            +-------+-------+
    !      |                |            |*      |      *|
    ! (e4) |                | (e6)  (e4) | nel+3 | nel+2 | (e6)
    !      |                |            |       |       |
    !      |       iel      |     ->   i8+-------+-------+i6
    !      |                |            |     i9|       |
    ! (e8) |                | (e2)  (e8) | iel   | nel+1 | (e2)
    !      |*               |            |*      |      *|
    !      +----------------+            +-------+-------+
    !     i1 (e1)      (e5) i2          i1 (e1) i5  (e5) i2
    !
!</description>

!<input>
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Store vertex- and element-values of the current element
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i2 = rhadapt%p_IverticesAtElement(2,iel)
    i3 = rhadapt%p_IverticesAtElement(3,iel)
    i4 = rhadapt%p_IverticesAtElement(4,iel)

    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e2 = rhadapt%p_IneighboursAtElement(2,iel)
    e3 = rhadapt%p_IneighboursAtElement(3,iel)
    e4 = rhadapt%p_IneighboursAtElement(4,iel)
    e5 = rhadapt%p_ImidneighboursAtElement(1,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(2,iel)
    e7 = rhadapt%p_ImidneighboursAtElement(3,iel)
    e8 = rhadapt%p_ImidneighboursAtElement(4,iel)

    ! Store values before refinement
    nel0 = rhadapt%NEL
    
    ! Add five new vertices I5,I6,I7,I8 and I9 at the midpoint of edges
    ! (I1,I2), (I2,I3), (I3,I4) and (I1,I4) and at the center of element IEL
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i2,i3,e2,i6,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i4,e3,i7,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i4,i1,e4,i8,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i1,i2,i3,i4,i9,rcollection,fcb_hadaptCallback)
    
    ! Replace element IEL and add three new elements E5, E6, and E7
    CALL replace_element2D(rhadapt,iel,i1,i5,i9,i8,e1,nel0+1,nel0+3,e8,e1,nel0+1,nel0+3,e8)
    CALL add_element2D(rhadapt,i2,i6,i9,i5,e2,nel0+2,iel,e5,e2,nel0+2,iel,e5)
    CALL add_element2D(rhadapt,i3,i7,i9,i6,e3,nel0+3,nel0+1,e6,e3,nel0+3,nel0+1,e6)
    CALL add_element2D(rhadapt,i4,i8,i9,i7,e4,iel,nel0+2,e7,e4,iel,nel0+2,e7)
    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel,nel0+2,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,nel0+3,nel0+2)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel,nel0+3)
        
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad4Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad4Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+2,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "refine_Quad4Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+3,ipos)

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+3,ipos)

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) &
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD4QUAD,&
        (/i1,i2,i3,i4,i5,i6,i7,i8,i9/),(/e1,e2,e3,e4,e5,e6,e7,e8/))
  END SUBROUTINE refine_Quad4Quad
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE convert_Tria2Tria(rhadapt,iel1,iel2,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines two neighboring triangles into one triangle 
    ! and performs regular refinement into four similar triangles afterwards.
    ! The local orientation of both elements can be uniquely determined
    ! from the elemental states so that the first node of each triangle
    ! is located at the midpoint of the bisected edge.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                         i3
    !            +                          +
    !           /|\                        /*\
    !     (e3) / | \ (e5)            (e3) /   \ (e5)
    !         /  |  \                    /nel+1\
    !        /   |   \        ->      i7+-------+i6
    !       /    |    \                / \nel+2/ \
    ! (e6) / iel | jel \ (e2)    (e6) /   \   /   \ (e2)
    !     /*     |     *\            /* iel\*/jel *\
    !    +-------+-------+          +-------+-------+
    !   i1(e1,e7)i5(e4,e8)i2       i1(e1,e7)i5(e4,e8)i2
    !
!</description>

!<input>
    ! Number of first element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel1

    ! Number of second element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel2

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,iel,jel
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6
    INTEGER(PREC_ELEMENTIDX) :: e1,e2,e3,e4,e5,e6,e7,e8

    ! Make sure that IEL < JEL
    IF (iel1 < iel2) THEN
      iel = iel1; jel = iel2
    ELSE
      iel = iel2; jel = iel1
    END IF

    ! Find local positions of element iel from is state
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i4 = rhadapt%p_IverticesAtElement(2,iel)
    i3 = rhadapt%p_IverticesAtElement(3,iel)
    
    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e3 = rhadapt%p_IneighboursAtElement(3,iel)
    
    e7 = rhadapt%p_ImidneighboursAtElement(1,iel)
    e6 = rhadapt%p_ImidneighboursAtElement(3,iel)

    ! Find local positions of element jel from is state
    i2 = rhadapt%p_IverticesAtElement(1,jel)
    
    e2 = rhadapt%p_IneighboursAtElement(1,jel)
    e4 = rhadapt%p_IneighboursAtElement(3,jel)
    
    e5 = rhadapt%p_ImidneighboursAtElement(1,jel)
    e8 = rhadapt%p_ImidneighboursAtElement(3,jel)
    
    ! Store values before conversion
    nel0 = rhadapt%NEL

    ! Add two new vertices I5 and I6 at the midpoint of edges (I2,I3)
    ! and (I1,I3), respectively.
    CALL add_vertex2D(rhadapt,i2,i3,e2,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i1,e3,i6,rcollection,fcb_hadaptCallback)

    ! Replace elements IEL and JEL and add two new elements
    CALL replace_element2D(rhadapt,iel,i1,i4,i6,e1,nel0+2,e6,e7,nel0+2,e6)
    CALL replace_element2D(rhadapt,jel,i2,i5,i4,e2,nel0+2,e4,e2,nel0+2,e8)
    CALL add_element2D(rhadapt,i3,i6,i5,e3,nel0+2,e5,e3,nel0+2,e5)
    CALL add_element2D(rhadapt,i4,i5,i6,jel,nel0+1,iel,jel,nel0+1,iel)

    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e2,e5,jel,nel0+1,jel)
    CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel,iel,nel0+1)

    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Tria2Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Tria2Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+2,ipos)

    ! "Lock" all vertices
    rhadapt%p_IvertexAge(i1) = -ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2) = -ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3) = -ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4) = -ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5) = -ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6) = -ABS(rhadapt%p_IvertexAge(i6))
    
    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_TRIA2TRIA,&
        (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6,e7,e8/))
  END SUBROUTINE convert_Tria2Tria

! ***************************************************************************

!<subroutine>

  SUBROUTINE convert_Quad2Quad(rhadapt,iel1,iel2,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines two neighboring quadrilaterals into one
    ! and performs regular refinement into four similar quadrilaterals.
    ! The local orientation of both elements can be uniquely determined
    ! from the elemental states so that the first node of each
    ! Depending on the given state of the element, the corresponding
    ! neighboring element is given explicitly.
    !
    !     initial quadrilateral      subdivided quadrilateral
    !
    !     i4(e7,e3)i7(f5,f5)i3         i4(e7,e3)i7(f5,f1)i3
    !      +-------+-------+            +-------+-------+
    !      |       |      *|            |*      |      *|
    ! (e4) |       |       | (f8)  (e4) | nel+2 | nel+1 | (f8)
    !      |       |       |            |       |       |
    !      |  iel  | jel   |     ->   i8+-------+-------+i6
    !      |       |       |            |       |i9     |
    ! (e8) |       |       | (f4)  (e8) |  iel  |  jel  | (f4)
    !      |*      |       |            |*      |      *|
    !      +-------+-------+            +-------+-------+
    !     i1(e1,e5)i5(f3,f7)i2         i1(e1,e5)i5(f3,f7)i2
    ! 
!</description>

!<input>
    ! Number of first element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel1

    ! Number of second element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel2
    
    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)     :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,iel,jel,e1,e3,e4,e5,e7,e8,f1,f3,f4,f5,f7,f8
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Make sure that IEL < JEL
    IF (iel1 < iel2) THEN
      iel = iel1; jel = iel2
    ELSE
      iel = iel2; jel = iel1
    END IF

    ! Find local positions of element iel
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i7 = rhadapt%p_IverticesAtElement(3,iel)
    i4 = rhadapt%p_IverticesAtElement(4,iel)

    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e3 = rhadapt%p_IneighboursAtElement(3,iel)
    e4 = rhadapt%p_IneighboursAtElement(4,iel)
    
    e5 = rhadapt%p_ImidneighboursAtElement(1,iel)
    e7 = rhadapt%p_ImidneighboursAtElement(3,iel)
    e8 = rhadapt%p_ImidneighboursAtElement(4,iel)

    ! Find local positions of element jel
    i2 = rhadapt%p_IverticesAtElement(4,jel)
    i3 = rhadapt%p_IverticesAtElement(1,jel)

    f1 = rhadapt%p_IneighboursAtElement(1,jel)   
    f3 = rhadapt%p_IneighboursAtElement(3,jel)
    f4 = rhadapt%p_IneighboursAtElement(4,jel)
    
    f5 = rhadapt%p_ImidneighboursAtElement(1,jel)
    f7 = rhadapt%p_ImidneighboursAtElement(3,jel)
    f8 = rhadapt%p_ImidneighboursAtElement(4,jel)

    ! Store values before conversion
    nel0 = rhadapt%NEL

    ! Add two new vertices I6, I8, and I9 at the midpoint of edges (I2,I3),
    ! (I1,I4) and (I1,I2), respectively.
    CALL add_vertex2D(rhadapt,i2,i3,f4,i6,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i4,i1,e4,i8,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i1,i2,i3,i4,i9,rcollection,fcb_hadaptCallback)

    ! Replace element IEL and JEL and add two new elements KEL and LEL
    CALL replace_element2D(rhadapt,iel,i1,i5,i9,i8,e1,jel,nel0+2,e8,e5,jel,nel0+2,e8)
    CALL replace_element2D(rhadapt,jel,i2,i6,i9,i5,f4,nel0+1,iel,f3,f4,nel0+1,iel,f7)
    CALL add_element2D(rhadapt,i3,i7,i9,i6,f1,nel0+2,jel,f8,f5,nel0+2,jel,f8)
    CALL add_element2D(rhadapt,i4,i8,i9,i7,e4,iel,nel0+1,e3,e4,iel,nel0+1,e7)

    ! Update list of neighboring elements
    CALL update_elementNeighbors2D(rhadapt,f4,f8,jel,nel0+1,jel)
    CALL update_elementNeighbors2D(rhadapt,e4,e4,iel,iel,nel0+2)
    CALL update_elementNeighbors2D(rhadapt,f1,f5,jel,nel0+1,nel0+1)
    CALL update_elementNeighbors2D(rhadapt,e3,e7,iel,nel0+2,nel0+2)

    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad2Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad2Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad2Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad2Quad: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+2,ipos)
    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,nel0+2,ipos)

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+2,ipos)

    ! "Lock" all vertices
    rhadapt%p_IvertexAge(i1) = -ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2) = -ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3) = -ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4) = -ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5) = -ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6) = -ABS(rhadapt%p_IvertexAge(i6))
    rhadapt%p_IvertexAge(i7) = -ABS(rhadapt%p_IvertexAge(i7))
    rhadapt%p_IvertexAge(i8) = -ABS(rhadapt%p_IvertexAge(i8))
    rhadapt%p_IvertexAge(i9) = -ABS(rhadapt%p_IvertexAge(i9))

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_QUAD2QUAD,&
        (/i1,i2,i3,i4,i5,i6,i7,i8,i9/),(/e1,f4,f1,e4,f3,f8,e3,e8/))
  END SUBROUTINE convert_Quad2Quad

   ! ***************************************************************************

!<subroutine>

  SUBROUTINE convert_Quad3Tria(rhadapt,iel1,iel2,iel3,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines three neighboring triangles which result
    ! from a Quad4Tria refinement into one quadrilateral and performs
    ! regular refinement into four quadrilaterals afterwards.
    ! This subroutine is based on the convention that IEL1 denotes the
    ! left element, IEL2 denots the right element and IEL3 stands for 
    ! the triangle which connects IEL1 and IEL2.
    !
    ! initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i2          i4 (e7)  i7 (e3)  i3
    !      +---------------+            +-------+-------+
    !      |\             /|            |*      |      *|
    ! (e4) | \           / | (e6)  (e4) | nel+1 | iel3  | (e6)
    !      |  \   iel3  /  |            |       |       |
    !      |   \       /   |     ->   i8+-------+-------+i6
    !      |    \     /    |            |       |i9     |
    ! (e8) | iel1\   /iel2 | (e2)  (e8) | iel1  | iel2  | (e2)
    !      |*     \*/     *|            |*      |      *|
    !      +-------+-------+            +-------+-------+
    !     i1(e1,e9)i5(e5,e10)i2        i1(e1,e9)i5(e5,e10)i2
    !
!</description>

!<input>
    ! Number of first triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel1

    ! Number of second triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel2
    
    ! Number of third triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel3

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)     :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: nel0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Get local data from element iel1
    i1 = rhadapt%p_IverticesAtElement(1,iel1)
    i5 = rhadapt%p_IverticesAtElement(2,iel1)
    i4 = rhadapt%p_IverticesAtElement(3,iel1)

    e1 = rhadapt%p_IneighboursAtElement(1,iel1)
    e4 = rhadapt%p_IneighboursAtElement(3,iel1)

    e9 = rhadapt%p_ImidneighboursAtElement(1,iel1)
    e8 = rhadapt%p_ImidneighboursAtElement(3,iel1)
    
    ! Get local data from element iel2
    i2 = rhadapt%p_IverticesAtElement(1,iel2)
    i3 = rhadapt%p_IverticesAtElement(2,iel2)
    
    e2 = rhadapt%p_IneighboursAtElement(1,iel2)
    e5 = rhadapt%p_IneighboursAtElement(3,iel2)

    e6 = rhadapt%p_ImidneighboursAtElement(1,iel2)
    e10= rhadapt%p_ImidneighboursAtElement(3,iel2)

    ! Get local data from element iel3
    e3 = rhadapt%p_IneighboursAtElement(2,iel3)
    e7 = rhadapt%p_ImidneighboursAtElement(2,iel3)
    
    ! Store values before conversion
    nel0 = rhadapt%NEL

    ! Add four new vertices I6,I7,I8 and I9 at the midpoint of edges 
    ! (I2,I3), (I3,I4) and (I1,I4) and at the center of element IEL
    CALL add_vertex2D(rhadapt,i2,i3,e2,i6,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i4,e3,i7,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i4,i1,e4,i8,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i1,i2,i3,i4,i9,rcollection,fcb_hadaptCallback)

    ! Replace elements IEL1, IEL2 and IEL3 and add one new element
    CALL replace_element2D(rhadapt,iel1,i1,i5,i9,i8,e1,iel2,nel0+1,e8,e9,iel2,nel0+1,e8)
    CALL replace_element2D(rhadapt,iel2,i2,i6,i9,i5,e2,iel3,iel1,e5,e2,iel3,iel1,e10)
    CALL replace_element2D(rhadapt,iel3,i3,i7,i9,i6,e3,nel0+1,iel2,e6,e3,nel0+1,iel2,e6)
    CALL add_element2D(rhadapt,i4,i8,i9,i7,e4,iel1,iel3,e7,e4,iel1,iel3,e7)

    ! Update element neighbors
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel3,iel2)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,nel0+1,iel3)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel1,iel1,nel0+1)

    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad3Tria: Unable to delete element from vertex list i3!",i3
      CALL sys_halt()
    END IF

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad3Tria: Unable to delete element from vertex list i4!",i4
      CALL sys_halt()
    END IF

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad3Tria: Unable to delete element from vertex list i4!",i4
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad3Tria: Unable to delete element from vertex list i5!",i5
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,iel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel1,ipos)    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+1,ipos)

    ! Finally, adjust numbers of triangles/quadrilaterals
    rhadapt%InelOfType(TRIA_NVETRI2D)  = rhadapt%InelOfType(TRIA_NVETRI2D)-3
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+3

    ! "Lock" all vertices
    rhadapt%p_IvertexAge(i1) = -ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2) = -ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3) = -ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4) = -ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5) = -ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6) = -ABS(rhadapt%p_IvertexAge(i6))
    rhadapt%p_IvertexAge(i7) = -ABS(rhadapt%p_IvertexAge(i7))
    rhadapt%p_IvertexAge(i8) = -ABS(rhadapt%p_IvertexAge(i8))
    rhadapt%p_IvertexAge(i9) = -ABS(rhadapt%p_IvertexAge(i9))
    
    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_QUAD3TRIA,&
        (/i1,i2,i3,i4,i5,i6,i7,i8,i9/),(/e1,e2,e3,e4,e5,e6,e7,e8/))
  END SUBROUTINE convert_Quad3Tria

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE convert_Quad4Tria(rhadapt,iel1,iel2,iel3,iel4,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four neighboring triangles which result
    ! from a Quad4Tria refinement into one quadrilateral and performs
    ! regular refinement into four quadrilaterals afterwards.
    ! This subroutine is based on the convention, that all four elements
    ! are given in couterclockwise order. More precisely, IEL2 and IEL3
    ! make up the inner "diamond" of the refinement, whereas IEL1 and IEL4
    ! are the right and left outer triangles, respectively.
    !
    !    initial quadrilateral         subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i3           i4 (e7)  i7 (e3) i3
    !      +---------------+             +-------+-------+
    !      |\\\\    iel3  *| (e12)       |*      |      *| (e12)
    ! (e4) | \*  \\\       | (e6)   (e4) | iel4  | iel3  | (e6)
    !      |  \      \\\   |             |       |       | 
    !      |   \  iel4   \\+i6    ->   i8+-------+-------+i6
    !      |    \        / |             |       |i9     |
    ! (e8) | iel1\     /   | (e11)  (e8) | iel1  | iel2  | (e11)
    !      |*     \  /iel2*| (e2)        |*      |      *| (e2)
    !      +-------+-------+             +-------+-------+
    !     i1(e1,e9)i5(e5,e10)i2         i1(e1,e9)i5(e5,e10)i2
    !
!</description>

!<input>
    ! Number of first triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel1

    ! Number of second triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel2

    ! Number of third triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel3

    ! Number of fourth triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel4

    ! Callback function
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Get local data from element iel1
    i1 = rhadapt%p_IverticesAtElement(1,iel1)
    i5 = rhadapt%p_IverticesAtElement(2,iel1)
    i4 = rhadapt%p_IverticesAtElement(3,iel1)

    e1 = rhadapt%p_IneighboursAtElement(1,iel1)
    e4 = rhadapt%p_IneighboursAtElement(3,iel1)

    e9 = rhadapt%p_ImidneighboursAtElement(1,iel1)
    e8 = rhadapt%p_ImidneighboursAtElement(3,iel1)
    
    ! Get local data from element iel2
    i2 = rhadapt%p_IverticesAtElement(1,iel2)
    i6 = rhadapt%p_IverticesAtElement(2,iel2)

    e2 = rhadapt%p_IneighboursAtElement(1,iel2)
    e5 = rhadapt%p_IneighboursAtElement(3,iel2)

    e11= rhadapt%p_ImidneighboursAtElement(1,iel2)
    e10= rhadapt%p_ImidneighboursAtElement(3,iel2)

    ! Get local data from element iel3
    i3 = rhadapt%p_IverticesAtElement(1,iel3)

    e3 = rhadapt%p_IneighboursAtElement(1,iel3)
    e6 = rhadapt%p_IneighboursAtElement(3,iel3)

    e7 = rhadapt%p_ImidneighboursAtElement(1,iel3)
    e12= rhadapt%p_ImidneighboursAtElement(3,iel3)
    
    ! Add three new vertices I7, I8 and I9 at the midpoint of edges
    ! (I3,I4) and (I1,I4) and at the center of element IEL
    CALL add_vertex2D(rhadapt,i3,i4,e3,i7,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i4,i1,e4,i8,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i1,i2,i3,i4,i9,rcollection,fcb_hadaptCallback)

    ! Replace all four elements
    CALL replace_element2D(rhadapt,iel1,i1,i5,i9,i8,e1,iel2,iel4,e8,e9,iel2,iel4,e8)
    CALL replace_element2D(rhadapt,iel2,i2,i6,i9,i5,e2,iel3,iel1,e5,e11,iel3,iel1,e10)
    CALL replace_element2D(rhadapt,iel3,i3,i7,i9,i6,e3,iel4,iel2,e6,e3,iel4,iel2,e12)
    CALL replace_element2D(rhadapt,iel4,i4,i8,i9,i7,e4,iel1,iel3,e7,e4,iel1,iel3,e7)

    ! Update element neighbors
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel4,iel3)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel1,iel1,iel4)

    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel4).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel4).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF

    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "convert_Quad4Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,iel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,iel4,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel4,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel4,ipos)

    ! Finally, adjust numbers of triangles/quadrilaterals
    rhadapt%InelOfType(TRIA_NVETRI2D)  = rhadapt%InelOfType(TRIA_NVETRI2D)-4
    rhadapt%InelOfType(TRIA_NVEQUAD2D) = rhadapt%InelOfType(TRIA_NVEQUAD2D)+4

    ! "Lock" all vertices
    rhadapt%p_IvertexAge(i1) = -ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2) = -ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3) = -ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4) = -ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5) = -ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6) = -ABS(rhadapt%p_IvertexAge(i6))
    rhadapt%p_IvertexAge(i7) = -ABS(rhadapt%p_IvertexAge(i7))
    rhadapt%p_IvertexAge(i8) = -ABS(rhadapt%p_IvertexAge(i8))
    rhadapt%p_IvertexAge(i9) = -ABS(rhadapt%p_IvertexAge(i9))

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_QUAD4TRIA,&
        (/i1,i2,i3,i4,i5,i6,i7,i8,i9/),(/e1,e2,e3,e4,e5,e6,e7,e8/))
  END SUBROUTINE convert_Quad4Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_2Tria1Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines two triangles resulting from a 1-tria : 2-tria
    ! green refinement into the original macro triangle. By definition, iel is
    ! the number of the triangle with smaller element number.
    ! Be warned, this numbering convention is not checked since the
    ! routine is only called internally and cannot be used from outside.
    !
    ! Due to the numbering convention used in the refinement procedure
    ! the "second" element neighbor is the other green triangle that is removed.
    !
    !     initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /|\                       / \
    !          / | \                     /   \
    !         /  |  \                   /     \
    !  (e3)  /   |   \  (e2)  ->  (e3) /       \ (e2)
    !       /    |    \               /   iel   \
    !      / iel | jel \             /           \
    !     /*     |      \           /*            \
    !    +-------+-------+         +---------------+
    !    i1 (e1) i4 (e4) i2        i1     (e1)     i2
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection    
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: jel,e1,e2,e3,e4,ielReplace
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4

    ! Get adjacent green elements
    jel = rhadapt%p_IneighboursAtElement(2,iel)

    ! Store vertex- and element values of the two elements
    i1 = rhadapt%p_IverticesAtElement(1,iel)
    i4 = rhadapt%p_IverticesAtElement(2,iel)
    i3 = rhadapt%p_IverticesAtElement(3,iel)
    i2 = rhadapt%p_IverticesAtElement(1,jel)
    
    e1 = rhadapt%p_IneighboursAtElement(1,iel)
    e3 = rhadapt%p_IneighboursAtElement(3,iel)
    e2 = rhadapt%p_IneighboursAtElement(1,jel)
    e4 = rhadapt%p_IneighboursAtElement(3,jel)

    ! Delete element JEL
    CALL remove_element2D(rhadapt,jel,ielReplace)
    IF (ielReplace.NE.0) CALL update_ElementNeighbors2D(rhadapt,ielReplace,jel)
    
    ! Update element IEL
    CALL replace_element2D(rhadapt,iel,i1,i2,i3,e1,e2,e3,e4,e2,e3)

    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e4,e4,jel,iel,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e2,jel,iel,iel)

    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "coarsen_2Tria1Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "coarsen_2Tria1Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel,ipos)

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_2TRIA1TRIA,&
        (/i1,i2,i3,i4/),(/e1,e2,e3,e4/))
  END SUBROUTINE coarsen_2Tria1Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Tria1Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-tria : 4-tria refinement into the original macro triangle.
    ! By definition, iel is the number of the inner red triangle.
    ! Be warned, this numbering convention is not checked since the
    ! routine is only called internally and cannot be used from outside.
    !
    ! Due to the numbering convention used in the refinement procedure
    ! the "third" element neighbor of the inner triangle has the number
    ! of the original macro element which will be used for resulting triangle.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /*\                       / \
    !     (e3) /   \ (e5)           (e3) /   \ (e5)
    !         / iel2\                   /     \
    !      i6+-------+i5      ->       /       \
    !       / \ iel / \               /  iel0   \
    ! (e6) /   \   /   \ (e2)   (e6) /           \ (e2)
    !     /*iel0\*/iel1*\           /             \
    !    +-------+-------+         +---------------+
    !   i1 (e1)  i4 (e4) i2       i1  (e1)   (e4)  i2
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: iel0,iel1,iel2,e1,e2,e3,e4,e5,e6,ielReplace
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6

    ! Store vertex- and element-values of the three neighboring elements
    i4 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i6 = rhadapt%p_IverticesAtElement(3,iel)
    
    iel1 = rhadapt%p_IneighboursAtElement(1,iel)
    iel2 = rhadapt%p_IneighboursAtElement(2,iel)
    iel0 = rhadapt%p_IneighboursAtElement(3,iel)
    
    i1 = rhadapt%p_IverticesAtElement(1,iel0)
    i2 = rhadapt%p_IverticesAtElement(1,iel1)
    i3 = rhadapt%p_IverticesAtElement(1,iel2)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1,iel0)
    e6 = rhadapt%p_IneighboursAtElement(3,iel0)

    e2 = rhadapt%p_IneighboursAtElement(1,iel1)
    e4 = rhadapt%p_IneighboursAtElement(3,iel1)

    e3 = rhadapt%p_IneighboursAtElement(1,iel2)
    e5 = rhadapt%p_IneighboursAtElement(3,iel2)

    ! Delete elements IEL, IEL1 and IEL2
    CALL remove_element2D(rhadapt,iel,ielReplace)
    IF (ielReplace.NE.0) CALL update_ElementNeighbors2D(rhadapt,ielReplace,iel)

    CALL remove_element2D(rhadapt,iel2,ielReplace)
    IF (ielReplace.NE.0) CALL update_ElementNeighbors2D(rhadapt,ielReplace,iel2)

    CALL remove_element2D(rhadapt,iel1,ielReplace)
    IF (ielReplace.NE.0) CALL update_ElementNeighbors2D(rhadapt,ielReplace,iel1)

    ! Update element IEL0
    CALL replace_element2D(rhadapt,iel0,i1,i2,i3,e1,e2,e3,e4,e5,e6)

    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e2,e2,iel1,iel0,iel0)
    CALL update_ElementNeighbors2D(rhadapt,e3,e3,iel2,iel0,iel0)
    CALL update_ElementNeighbors2D(rhadapt,e4,e4,iel1,iel0,iel0)
    CALL update_ElementNeighbors2D(rhadapt,e5,e5,iel2,iel0,iel0)

    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "coarsen_4Tria1Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      PRINT *, "coarsen_4Tria1Tria: Unable to delete element from vertex list!"
      CALL sys_halt()
    END IF
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel0,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel0,ipos)

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA1TRIA,&
        (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6/))
  END SUBROUTINE coarsen_4Tria1Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Tria2Tria(rhadapt,iel,iloc,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-tria : 4-tria refinement into two green triangle.
    ! The local position (node 1,2,3 of the interior element) of the midpoint 
    ! vertex that should be  kept is denoted by iloc. 
    ! By definition, iel is the number of the inner red triangle.
    ! Be warned, this numbering convention is not checked since the
    ! routine is only called internally and cannot be used from outside.
    !
    ! Due to the numbering convention used in the refinement procedure
    ! the "third" element neighbor of the inner triangle has the number
    ! of the original macro element which will be used for the first resulting 
    ! triangle. The number of the second triangle will be the number of the 
    ! "first" element neighbor of the inner triangle.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /*\                       /|\
    !    (e3)  /   \ (e5)           (e3) / | \ (e5)
    !         /iel2 \                   /  |  \
    !     i6 +-------+ i5     ->       /   |   \
    !       / \ iel / \               /    |    \
    ! (e6) /   \   /   \ (e2)   (e6) /     |     \ (e2)
    !     /*iel0\*/iel1*\           /*iel0 | iel1*\
    !    +-------+-------+         +-------+-------+
    !    i1 (e1) i4 (e4) i2        i1 (e1) i4 (e4)  i2
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Relative position of the midpoint vertex that is kept
    INTEGER, INTENT(IN) :: iloc

    ! callback routines
    include 'intf_hadaptcallback.inc'
    OPTIONAL :: fcb_hadaptCallback
!</input>

!<inputoutput>
    ! Adaptivity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: iel0,iel1,iel2,e1,e2,e3,e4,e5,e6,ielReplace
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6
    
    ! Store vertex- and element-values of the three neighboring elements
    i4 = rhadapt%p_IverticesAtElement(1,iel)
    i5 = rhadapt%p_IverticesAtElement(2,iel)
    i6 = rhadapt%p_IverticesAtElement(3,iel)
    
    iel1 = rhadapt%p_IneighboursAtElement(1,iel)
    iel2 = rhadapt%p_IneighboursAtElement(2,iel)
    iel0 = rhadapt%p_IneighboursAtElement(3,iel)
    
    i1 = rhadapt%p_IverticesAtElement(1,iel0)
    i2 = rhadapt%p_IverticesAtElement(1,iel1)
    i3 = rhadapt%p_IverticesAtElement(1,iel2)

    ! Store values of the elements adjacent to the resulting macro element
    e1 = rhadapt%p_IneighboursAtElement(1,iel0)
    e6 = rhadapt%p_IneighboursAtElement(3,iel0)

    e2 = rhadapt%p_IneighboursAtElement(1,iel1)
    e4 = rhadapt%p_IneighboursAtElement(3,iel1)

    e3 = rhadapt%p_IneighboursAtElement(1,iel2)
    e5 = rhadapt%p_IneighboursAtElement(3,iel2)

    ! Delete elements IEL, and IEL2
    CALL remove_element2D(rhadapt,iel,ielReplace)
    IF (ielReplace.NE.0) CALL update_ElementNeighbors2D(rhadapt,ielReplace,iel)

    CALL remove_element2D(rhadapt,iel2,ielReplace)
    IF (ielReplace.NE.0) CALL update_ElementNeighbors2D(rhadapt,ielReplace,iel2)
    
    ! Which midpoint vertex should be kept?
    SELECT CASE(iloc)
    CASE(1)
      ! Update elements IEL and IEL1
      CALL replace_element2D(rhadapt,iel0,i1,i4,i3,e1,iel1,e3,e1,iel1,e6)
      CALL replace_element2D(rhadapt,iel1,i2,i3,i4,e2,iel0,e4,e5,iel0,e4)

      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e3,e3,iel2,iel0,iel0)
      CALL update_ElementNeighbors2D(rhadapt,e5,e5,iel2,iel1,iel1)

      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel0,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel1,ipos)

      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
          CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA2TRIA1,&
          (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6/))

    CASE(2)
      ! Update elements IEL and IEL1
      CALL replace_element2D(rhadapt,iel0,i2,i5,i1,e2,iel1,e1,e2,iel1,e4)
      CALL replace_element2D(rhadapt,iel1,i3,i1,i5,e3,iel0,e5,e6,iel0,e5)

      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e4,e4,iel1,iel0,iel0)
      CALL update_ElementNeighbors2D(rhadapt,e2,e2,iel1,iel0,iel0)
      CALL update_ElementNeighbors2D(rhadapt,e5,e5,iel2,iel1,iel1)
      CALL update_ElementNeighbors2D(rhadapt,e3,e3,iel2,iel1,iel1)
      CALL update_ElementNeighbors2D(rhadapt,e6,e6,iel0,iel1,iel1)

      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel0,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel0,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,iel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel1,ipos)

      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
          CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA2TRIA2,&
          (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6/))

    CASE(3)
      ! Update elements IEL and IEL1
      CALL replace_element2D(rhadapt,iel0,i1,i2,i6,e1,iel1,e6,e4,iel1,e6)
      CALL replace_element2D(rhadapt,iel1,i3,i6,i2,e3,iel0,e2,e3,iel0,e5)

      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e4,e4,iel1,iel0,iel0)
      CALL update_ElementNeighbors2D(rhadapt,e5,e5,iel2,iel1,iel1)
      CALL update_ElementNeighbors2D(rhadapt,e6,e6,iel2,iel1,iel1)

      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        PRINT *, "coarsen_4Tria2Tria: Unable to delete element from vertex list!"
        CALL sys_halt()
      END IF
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel0,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel1,ipos)

      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
          CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA2TRIA3,&
          (/i1,i2,i3,i4,i5,i6/),(/e1,e2,e3,e4,e5,e6/))

    CASE DEFAULT
      PRINT *, "coarsen_4Tria2Tria: Invalid position of midpoint vertex!"
      CALL sys_halt()
    END SELECT
  END SUBROUTINE coarsen_4Tria2Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mark_refinement2D(rhadapt,rindicator)

!<description>
    ! This subroutine marks all elements that should be refined due
    ! to accuracy reasons. The decision is based on some indicator
    ! vector which must be given elementwise. Moreover, this subroutine 
    ! filters all triangular elements and sets their fourth marker 
    ! to -1. Finally, the age of vertices connected to elements 
    ! which are marked for refinement are negated.
!</description>

!<input>
    ! Indicator vector for refinement
    TYPE(t_vectorScalar), INTENT(IN) :: rindicator
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER    :: p_Dindicator
    INTEGER, DIMENSION(:), POINTER     :: p_Imarker
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D)  :: Kvert
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE2D) :: Kadj
    INTEGER(PREC_VERTEXIDX)  :: ivt
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: ive

    ! Check if dynamic data structures are generated and contain data
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA) THEN
      PRINT *, "mark_refinement2D: dynamic data structures are not generated"
      CALL sys_halt()
    END IF

    ! Initialize marker structure for NEL0 elements
    IF (rhadapt%h_Imarker /= ST_NOHANDLE) &
        CALL storage_free(rhadapt%h_Imarker)
    CALL storage_new('mark_refinement2D','Imarker',rhadapt%NEL0,&
        ST_INT,rhadapt%h_Imarker,ST_NEWBLOCK_ZERO)
       
    ! Set pointers
    CALL lsyssc_getbase_double(rindicator,p_Dindicator)
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)

    ! Set state of all vertices to "free". Note that vertices of the
    ! initial triangulation are always "locked", i.e. have no positive age.
    DO ivt=1,SIZE(rhadapt%p_IvertexAge)
      rhadapt%p_IvertexAge(ivt) = ABS(rhadapt%p_IvertexAge(ivt))
    END DO

    ! Loop over all elements and mark those for which
    ! the indicator is greater than the prescribed treshold
    mark: DO iel=1,SIZE(p_Dindicator)

      IF (p_Dindicator(iel) .GT. rhadapt%drefinementTolerance) THEN
        
        ! Mark element for refinement
        Kvert = rhadapt%p_IverticesAtElement(:,iel)
        Kadj  = rhadapt%p_IneighboursAtElement(:,iel)

        ! An element can only be refined, if all of its vertices do
        ! not exceed the number of admissible subdivision steps.
        ! So, check the element type and loop over the 3 or 4 vertices.
        IF (Kvert(4) == 0) THEN
          
          ! If triangle has reached maximum number of refinement levels,
          ! then enforce no further refinement of this element
          IF (ANY(ABS(rhadapt%p_IvertexAge(Kvert(1:3))).EQ.&
              rhadapt%NSUBDIVIDEMAX)) THEN
            p_Imarker(iel) = MARK_ASIS
            
            ! According to the indicator, this element should be refined. Since the 
            ! maximum admissible refinement level has been reached no refinement 
            ! was performed. At the same time, all vertices of the element should
            ! be "locked" to prevent this element from coarsening
            DO ive=1,3
              rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive))) 
            END DO
            
            CYCLE mark
          END IF
          
          ! Otherwise, we can mark the triangle for refinement
          p_Imarker(iel) = MARK_REF_TRIA4TRIA
          
          ! Moreover, we have to "lock" its vertices from recoarsening
          DO ive=1,3
            rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
          END DO
          
          ! Update number of new vertices. In principle, we can increase the number of 
          ! new vertices by 3, i.e., one for each each. However, some care must be taken
          ! if the edge belongs to some adjacent element which has been marked for
          ! refinement previously. Hence, a new vertex is only added, if the edge
          ! is connected to the boundary or if the adjacent element has not been marked.
          DO ive=1,3
            IF (Kadj(ive)==0) THEN
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
            ELSEIF((p_Imarker(Kadj(ive)) .NE. MARK_REF_TRIA4TRIA) .AND.&
                   (p_Imarker(Kadj(ive)) .NE. MARK_REF_QUAD4QUAD)) THEN
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
            END IF
          END DO
        ELSE
          
          ! If quadrilateral has reached maximum number of refinement levels,
          ! then enforce no further refinement of this element
          IF (ANY(ABS(rhadapt%p_IvertexAge(Kvert(1:4))).EQ.&
              rhadapt%NSUBDIVIDEMAX)) THEN
            p_Imarker(iel) = MARK_ASIS

            ! According to the indicator, this element should be refined. Since the 
            ! maximum admissible refinement level has been reached no refinement 
            ! was performed. At the same time, all vertices of the element should
            ! be "locked" to prevent this element from coarsening
            DO ive=1,4
              rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive))) 
            END DO

            CYCLE mark
          END IF
          
          ! Otherwise, we can mark the quadrilateral for refinement
          p_Imarker(iel) = MARK_REF_QUAD4QUAD
          
          ! Moreover, we have to "lock" its vertices from recoarsening
          DO ive=1,4
            rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
          END DO
          
          ! Update number of new vertices
          DO ive=1,4
            IF (Kadj(ive)==0) THEN
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
            ELSEIF((p_Imarker(Kadj(ive)) .NE. MARK_REF_TRIA4TRIA) .AND.&
                   (p_Imarker(Kadj(ive)) .NE. MARK_REF_QUAD4QUAD)) THEN
              rhadapt%increaseNVT = rhadapt%increaseNVT+1
            END IF
          END DO

          ! And don't forget the new vertex in the interior of the quadrilateral
          rhadapt%increaseNVT = rhadapt%increaseNVT+1
        END IF
        
      ELSE
        ! Unmark element for refinement
        p_Imarker(iel) = MARK_ASIS
      END IF
    END DO mark

    ! Set specifier to "marked for refinement"
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_MARKEDREFINE)
  END SUBROUTINE mark_refinement2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE redgreen_mark_coarsening2D(rhadapt,rindicator)

!<description>
    ! This routine marks all elements that sould be recoarsened due to accuracy
    ! reasons. The decision is based on some indicator vector which must be 
    ! given elementwise. The recoarsening strategy is as follows: A subset of
    ! elements can only be coarsened if this subset results from a previous
    ! refinement step. In other words, the grid cannot become coarser than the
    ! initial triangulation. Moreover, one recoarsening step can only "undo"
    ! what one refinement step can "do". This is in contrast to other "node-
    ! removal" techniques, which remove all needless vertices and retriangulates
    ! the generated "holes". However, such algorithms cannot guarantee a grid
    ! hierarchy between refined and recoarsened grids.
    !
    ! In the context of hierarchical red-green mesh adaptivity, each recoarsening
    ! step is the inverse operation of a previous refinement step. In other words,
    ! the result of a complete recoarsening step could also be generated by locally
    ! refining the initial triangulation. The algorithm is as follows:
    !
    ! 1) All vertices present in the initial triangulation cannot be removed. A 
    ! vertex belongs to the initial triangulation if and only if its age is zero.
    !
    ! 2) All other nodes which do not belong to elements that are marked for 
    ! refinement are potential candidates for removal. In order to prevent the
    ! removel of such vertices, the nodes must be "locked" by some mechanism.
    ! 
    ! 3) All three/four nodes of a red triangle/quadrilateral are "locked" if 
    ! the element should not be coarsened due to accuracy reasons. Otherwise, if
    ! a red element should be coarsened then proceed as follows:
    !  a) For a red quadrilateral "lock" the single vertex of the macro element
    !  b) For a red outer triangle "lock" the single vertex of the macro element
    !  c) For a red inner tringle do nothing
    !
    ! 4) For a green element "lock" the vertices of the macro element.
    !
    ! As a consequence, all vertices of macro elements are locked and in addition
    ! all vertices belonging to red elements which should not be coarsened due to
    ! accuracy reasons are removed from the list of removable nodes.
    ! It remains to prevent the generaiton of blue elements to to nodal removal.
    ! This must be done iteratively since "locking" of on vertex may alter the
    ! situation for all adjacent elements. The generaiton of blue elements can
    ! be prevented by applying the following rules:
    !
    ! 1) If an inner red triangle has two locked vertices, then the third vertex
    ! is locked, too.
    !
    ! 2) If the youngest vertex of an outer green triangle (i.e. the vertex which
    ! is responsible for the green refinement) is locked, then lock all nodes
    ! connected to the green triangle. There is one exception to this rule. For
    ! the "inner" green triangles resulting from a 1-quad : 4-tria refinement, 
    ! this rule should not be applied since the four triangles can be converted
    ! into three triangles if only one of the two "green vertices" are locked.
    ! (apply 4-tria : 1-quad coarsening and eliminate the handing node by suitable
    ! 1-quad : 3-tria refinement = 4-tria : 3-tria conversion).
!</description>

!<input>
    ! Indicator vector for refinement
    TYPE(t_vectorScalar), INTENT(IN) :: rindicator
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER                    :: p_Dindicator
    INTEGER,  DIMENSION(:), POINTER                    :: p_Imarker
    INTEGER(PREC_VERTEXIDX),  DIMENSION(TRIA_MAXNVE2D) :: Kvert,KvertJ
    INTEGER(PREC_VERTEXIDX)  :: i,j
    INTEGER(PREC_ELEMENTIDX) :: iel,jel,kel
    INTEGER :: ive,nve,istate,jstate,icount,ishift,iVertexLocked
    LOGICAL :: isModified

    ! Check if dynamic data structures are generated and contain data
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA .OR.&
        IAND(rhadapt%iSpec,HADAPT_MARKEDREFINE).NE.HADAPT_MARKEDREFINE) THEN
      PRINT *, "redgreen_mark_coarsening2D: dynamic data structures are not &
          & generated or no marker for grid refinement exists!"
      CALL sys_halt()
    END IF
    
    ! Set pointers
    IF (rhadapt%h_Imarker == ST_NOHANDLE) THEN
      PRINT *, "redgreen_mark_coarsening2D: marker array is not available"
      CALL sys_halt()
    END IF
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
    CALL lsyssc_getbase_double(rindicator,p_Dindicator)

    ! All nodes of the initial triangulation have age-0 and, hence, will
    ! never be deleted. Loop over all elements and "lock" nodes which 
    ! should not be removed from the triangulation due to the indicator.
    DO iel=1,SIZE(p_Dindicator)
      
      ! Check if the current element is marked for refinement. Then
      ! all vertices are already locked by the refinement procedure
      IF (p_Imarker(iel) .NE. MARK_ASIS_TRIA .AND.&
          p_Imarker(iel) .NE. MARK_ASIS_QUAD) CYCLE

      ! Get number of vertices per element
      nve = get_NVE(rhadapt,iel)

      ! Get local data for element iel
      Kvert(1:nve) = rhadapt%p_IverticesAtElement(1:nve,iel)

      ! Get state of current element
      SELECT CASE(nve)
      CASE(TRIA_NVETRI2D)
        istate=redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))

      CASE(TRIA_NVEQUAD2D)
        istate=redgreen_getStateQuad(rhadapt%p_IvertexAge(Kvert(1:4)))

      CASE DEFAULT
        PRINT *, "redgreen_mark_coarsening2D: Invalid number of vertices per element!"
        CALL sys_halt()
      END SELECT
      
      ! Phase 1: Depending on the state if the element IEL and the indicator
      ! "lock" some or all of its nodes so that they will not be removed.
      SELECT CASE(istate)

        !-----------------------------------------------------------------------
        ! The starting point of the re-coarsening algorithm are the red elements.
        ! Strictly speaking, if all four subelements of a red-refinement should
        ! be coarsened, then the aim is to convert them back into their original
        ! macro element. On the other hand, if one or more red triangle cannot
        ! be combined, then all of its nodes must be "locked".

      CASE(STATE_TRIA_ROOT,STATE_QUAD_ROOT)
        ! The vertices of a root triangle/quadrilateral have zero age by
        ! definition. Hence, they are never deleted from the triangulation.

        !-----------------------------------------------------------------------
        ! Non-green elements: If the element is not marked for coarsening,
        ! then "lock" all of its vertices. For quadrilateral elements, this
        ! is sufficient to determine those elements which can be combined.
        ! This algorithm works also for inner red triangles, but for outer
        ! red triangles some more care must be taken. An outer red triangle
        ! has the same state as an inner triangle of a 1-quad : 4-tria 
        ! refinement. Hence, we must check, if the current triangle is 
        ! adjacent to an inner red triangle in order to apply the algorithm.

      CASE(STATE_QUAD_RED1,STATE_QUAD_RED2,STATE_QUAD_RED3)
        PRINT *, "redgreen_mark_coarsening2D: Theoretically, this state should &
            & not appear. Please inform the author of this routine."
        CALL sys_halt()

      CASE(STATE_QUAD_RED4)
        ! Element IEL is red sub-element of a Quad4Quad refinement.

        ! Should the element be coarsened?
        IF (p_Dindicator(iel) .GE. rhadapt%dcoarseningTolerance) THEN

          ! If this is not the case, then "lock" all of its four nodes
          DO ive=1,TRIA_NVEQUAD2D
            i = rhadapt%p_IverticesAtElement(ive,iel)
            rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
          END DO
        ELSE

          ! Otherwise, "lock" only the node from the macro element.
          ! Due to the refinement strategy, this is the first vertex.
          i = rhadapt%p_IverticesAtElement(1,iel)
          rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

          ! Provisionally, mark element for generic coarsening
          p_Imarker(iel) = MARK_CRS_GENERIC
        END IF

      CASE(STATE_TRIA_REDINNER)
        ! Element IEL is inner red triangle of a Tria4Tria refinement.
        
        ! Should the element be coarsened?
        IF (p_Dindicator(iel) .GE. rhadapt%dcoarseningTolerance) THEN
          
          ! If this is not the case, then "lock" all of its three nodes
          DO ive=1,TRIA_NVETRI2D
            i = rhadapt%p_IverticesAtElement(ive,iel)
            rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
          END DO

        ELSE
          ! Otherwise, mark element for generic coarsening
          p_Imarker(iel) = MARK_CRS_GENERIC
        END IF

      CASE(STATE_TRIA_OUTERINNER)
        ! Element IEL is either an outer red triangle of a 1-tria : 4-tria re-
        ! finement or one of the two inner triangles of a 1-quad : 4-tria refinement
        
        ! Are we outer red triangle that was not created during the 
        ! red-green marking routine as a result of element conversion?
        jel    = rhadapt%p_IneighboursAtElement(2,iel)
        jstate = redgreen_getStateTria(rhadapt%p_IvertexAge(&
                 rhadapt%p_IverticesAtElement(1:3,jel)))
        
        IF (jstate .EQ. STATE_TRIA_REDINNER) THEN
          
          ! Should the element be coarsened?
          IF (p_Dindicator(iel) .GE. rhadapt%dcoarseningTolerance) THEN

            ! If this is not the case, then "lock" all of its three nodes
            DO ive=1,TRIA_NVETRI2D
              i = rhadapt%p_IverticesAtElement(ive,iel)
              rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
            END DO
          ELSE
            
            ! Otherwise, "lock" only the node from the macro element.
            ! Due to the refinement strategy, this is the first vertex.
            i = rhadapt%p_IverticesAtElement(1,iel)
            rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

            ! Provisionally, mark element for generic coarsening
            p_Imarker(iel) = MARK_CRS_GENERIC
          END IF
        ELSE
          
          ! We are an inner green triangle resulting from a 1-quad : 4-tria
          ! refinement. By construction, the first vertex belongs to the macro
          ! element and, consequently, has to be "locked" from removal
          i = rhadapt%p_IverticesAtElement(1,iel)
          rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

          ! Provisionally, mark element for generic coarsening
          p_Imarker(iel) = MARK_CRS_GENERIC
        END IF

        !-----------------------------------------------------------------------
        ! Green elements: "Lock" the vertices of the macro element.
        ! This sounds rather complicated but it is an easy task.
        ! All non-green elements have been processed before, so we can be sure
        ! to deal only with green elements. For a green element that results 
        ! from a 1-tria : 2-tria refinement, we know that the vertices of the
        ! macro element are located at the first position. The other "old" node
        ! is the second/third vertex depending on the position of the green
        ! triangle, i.e. right/left. By construction, this is also true for
        ! the outer green triangles resulting from a 1-quad : 3/4-tria
        ! refinement. It suffices to check for inner green triangles of a
        ! 1-quad : 3-tria refinement, since the inner triangles of a 1-quad :
        ! 4-tria refinement are processed in (STATE_TRIA_OUTERINNER)
        !-----------------------------------------------------------------------
        
      CASE(STATE_TRIA_GREENOUTER_RIGHT)
        ! Element IEL is right green triangle resulting from a 1-tria : 2-tria
        ! refinement. Here, the third vertex is the one which was last inserted.
        ! Hence, the first vertex is older then the third one by construction.
        i = rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

        i = rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

        ! Provisionally, mark element for generic coarsening
        p_Imarker(iel) = MARK_CRS_GENERIC
        
      CASE(STATE_TRIA_GREENOUTER_LEFT)
        ! Element IEL is left green triangle resulting from a 1-tria : 2-tria
        ! refinement. Here, the second vertex is the one which was last insrted.
        ! Hence, the first vertex is older than the second one by construction.
        i = rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

        i = rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

        ! Provisionally, mark element for generic coarsening
        p_Imarker(iel) = MARK_CRS_GENERIC

      CASE(STATE_TRIA_GREENINNER)
        ! Element IEL is the inner green triangle resulting from a 1-quad :
        ! 3-tria refinement. By construction, the first vertex is the newly
        ! introduced one which is the youngest vertex of the element. Hence,
        ! "lock" both the second and the third vertex which belong to the
        ! original macro element which was a quadrilateral
        i = rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
        
        j = rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(j) = -ABS(rhadapt%p_IvertexAge(j))

        ! Provisionally, mark element for generic coarsening
        p_Imarker(iel) = MARK_CRS_GENERIC

      CASE(STATE_QUAD_HALF1,STATE_QUAD_HALF2)
        PRINT *, "redgree_mark_coarsening2D: Not implemented yet."
        CALL sys_halt()

      CASE DEFAULT
        PRINT *, "redgreen_mark_coarsening2D: Invalid state",istate
        CALL sys_halt()
      END SELECT
    END DO

    ! Phase 2: Prevent the recoarsening algorithm from generating blue elements.
    ! The above loop "locked" all vertices, which should not be removed either
    ! based on the indicator vector or due to the fact, that the vertex 
    ! correcponds to the outer macro element. However, we have to keep in mind
    ! that blue elements could be generated during recoarsening. In short,
    ! inner red triangles with two "locked" vertices and red quadrilaterals
    ! for which the macro element has only one or two "unlocked" nodes must
    ! be locked completely. Moreover, all nodes of a green element must be
    ! "locked" if one of its youngest vertices (which have been created during
    ! the last green-refinement) must not be removed, i.e., it is locked.
    ! Keep in mind, that the "locking" of a vertex requires the recursive check
    ! of all surrounding elements. In order to prevent this recursion, we make 
    ! use of a do-while loop which is terminated of no more vertices are touched.

    isModified = .TRUE.
    DO WHILE(isModified)
      isModified = .FALSE.

      ! Loop over all elements (from the initial triangulation + those which
      ! have been created during the conversion of green into red elements)
      DO iel=1,rhadapt%NEL

        ! Only process those elements which are candidates for removal
        IF (p_Imarker(iel) .NE. MARK_CRS_GENERIC) CYCLE

        ! Get number of vertices per element
        nve = get_NVE(rhadapt,iel)
        
        ! Get local data for element iel
        Kvert(1:nve) = rhadapt%p_IverticesAtElement(1:nve,iel)

        ! Check if all vertices of the element are locked then delete element
        ! from the list of removable elements. This could also be done in the
        ! third phase when the locked vertices are translated into coarsening
        ! rules. However, this check is relatively cheap as compared to the
        ! computation of the exact element state for some elements. Heuristically,
        ! a large number of elements can be filtered out by this check in advance
        ! so that no detailed investigation of their vertices is required.
        IF (ALL(rhadapt%p_IvertexAge(Kvert(1:nve)).LE.0)) THEN
          p_Imarker(iel) = MERGE(MARK_ASIS_TRIA,MARK_ASIS_QUAD,nve == 3)
          CYCLE
        END IF
        
        ! Get state of current element
        SELECT CASE(nve)
        CASE(TRIA_NVETRI2D)
          istate=redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))

        CASE(TRIA_NVEQUAD2D)
          istate=redgreen_getStateQuad(rhadapt%p_IvertexAge(Kvert(1:4)))

        CASE DEFAULT
          PRINT *, "redgreen_mark_coarsening2D: Invalid number of vertices per element!"
          CALL sys_halt()
        END SELECT
        
        ! Do we have to "lock" some vertices?
        SELECT CASE(istate)

        CASE(STATE_QUAD_RED4)
          ! Element IEL is one of four red quadrilaterals of a 1-quad : 4-quad refinement.

          ! If the interior vertex is locked then lock all vertices of the element
          IF (rhadapt%p_IvertexAge(Kvert(3)) .LE. 0) THEN

            ! Lock all other vertices. Note that the first vertiex and the third
            ! one are already locked, so that we have to consider vertices 2 and 4
            DO ive=2,4,2
              rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
            END DO

            ! Delete element from list of removable elements
            p_Imarker(iel) = MARK_ASIS_QUAD
            
            ! We modified some vertex in this iteration
            isModified=.TRUE.

          ! If three midpoint vertices of adjacent quadrilaterals are locked
          ! then all vertices of the four red sub-quadrilaterals must be locked
          ! in order to prevent the generation of blue quadrilaterals
          ELSEIF(rhadapt%p_IvertexAge(Kvert(2)) .LE. 0 .AND.&
                 rhadapt%p_IvertexAge(Kvert(4)) .LE. 0) THEN
            
            ! Check the midpoint of the counterclockwise neighboring element
            jel = rhadapt%p_IneighboursAtElement(2,iel)
            i   = rhadapt%p_IverticesAtElement(2,jel)
            IF (rhadapt%p_IvertexAge(i) .LE. 0) THEN
              
              ! Lock all vertices of element IEL. Due to the fact that the interior
              ! vertex is locked, the nodes of all other quadrilaterals will be
              ! locked in the next loop of phase 2.
              DO ive=1,TRIA_NVEQUAD2D
                rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
              END DO

              ! Delete element from list of removable elements
              p_Imarker(iel) = MARK_ASIS_QUAD
              
              ! We modified some vertex in this iteration
              isModified=.TRUE.
              CYCLE
            END IF

            ! Check the midpoint of the clockwise neighboring element
            jel = rhadapt%p_IneighboursAtElement(3,iel)
            i   = rhadapt%p_IverticesAtElement(4,jel)
            IF (rhadapt%p_IvertexAge(i) .LE. 0) THEN
              
              ! Lock all vertices of element IEL. Due to the fact that the interior
              ! vertex is locked, the nodes of all other quadrilaterals will be
              ! locked in the next loop of phase 2.
              DO ive=1,TRIA_NVEQUAD2D
                rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
              END DO

              ! Delete element from list of removable elements
              p_Imarker(iel) = MARK_ASIS_QUAD
              
              ! We modified some vertex in this iteration
              isModified=.TRUE.
              CYCLE
            END IF
          END IF

        CASE(STATE_TRIA_REDINNER)
          ! Element IEL is inner red triangle of a 1-tria : 4-tria refinement.
          
          ! Determine number of locked vertices of the inner triangle.
          icount = 0
          DO ive=1,TRIA_NVETRI2D
            IF (rhadapt%p_IvertexAge(Kvert(ive)) .LE. 0) icount=icount+1
          END DO
          
          ! If exactly two vertices are locked, then "lock" the third one, too.
          IF (icount .EQ. 2) THEN
            DO ive=1,TRIA_NVETRI2D
              rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
            END DO
            
            ! Delete element from list of removable elements
            p_Imarker(iel) = MARK_ASIS
            
            ! We modified some vertex in this iteration
            isModified=.TRUE.
            
          ELSEIF(icount .EQ. 3) THEN
            ! Delete element from list of removable elements
            p_Imarker(iel) = MARK_ASIS
          END IF

        CASE(STATE_TRIA_GREENOUTER_RIGHT)
          ! Element IEL is right green triangle of a 1-tria : 2-tria refinement.
          
          ! If the third vertex is locked, then "lock" all other vertices, too.
          IF (rhadapt%p_IvertexAge(Kvert(3)) .LE. 0) THEN
            
            ! Lock all other vertices
            DO ive=1,2
              rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
            END DO
            
            ! Delete element from list of removable elements              
            p_Imarker(iel) = MARK_ASIS
            
            ! We modified some vertex in this iteration
            isModified=.TRUE.
          END IF
          
        CASE(STATE_TRIA_GREENOUTER_LEFT)
          ! Element IEL is left green triangle of a 1-tria : 2-tria refinement.
          
          ! If the second vertex is locked, then "lock" all other vertices, too.
          IF (rhadapt%p_IvertexAge(Kvert(2)) .LE. 0) THEN
            
            ! Lock all other vertices
            DO ive=1,3,2
              rhadapt%p_IvertexAge(Kvert(ive)) = -ABS(rhadapt%p_IvertexAge(Kvert(ive)))
            END DO
            
            ! Delete element from list of removable elements              
            p_Imarker(iel) = MARK_ASIS
            
            ! We modified some vertex in this iteration
            isModified=.TRUE.
          END IF

        END SELECT
      END DO
    END DO

    ! Phase 3: All vertices that are still "free" can be removed from the
    ! triangulation. At the moment, elements are only marked for generic
    ! removal and some of them cannot be removed, e.g., if all of its vertices
    ! are locked. In a loop over all elements which are marked for generic
    ! removal we determine how to remove the element or delete it from the
    ! list of removable elements if all of its vertices are locked.
    DO iel=1,rhadapt%NEL
      
      ! Only process those elements which are candidates for removal
      IF (p_Imarker(iel) .NE. MARK_CRS_GENERIC) CYCLE

      ! Get number of vertices per element
      nve = get_NVE(rhadapt,iel)

      ! Get local data for element iel
      Kvert(1:nve) = rhadapt%p_IverticesAtElement(1:nve,iel)
      
      ! Get state of current element
      SELECT CASE(nve)
      CASE(TRIA_NVETRI2D)
        istate=redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))

      CASE(TRIA_NVEQUAD2D)
        istate=redgreen_getStateQuad(rhadapt%p_IvertexAge(Kvert(1:4)))

      CASE DEFAULT
        PRINT *, "redgreen_mark_coarsening2D: Invalid number of vertices per element!"
        CALL sys_halt()
      END SELECT
      
      ! What state are we?
      SELECT CASE(istate)
        
      CASE(STATE_TRIA_REDINNER)
        ! If all vertices of the element are free, then the inner red triangle
        ! can be combined with its three nieghboring red triangles so as to
        ! recover the original macro triangle. If one vertex of the inner
        ! triangle is locked, then the inner red tiangle together with its three
        ! neighboring elements can be convertec into two green triangles.
        IF (rhadapt%p_IvertexAge(Kvert(1)) .LE. 0) THEN
          p_Imarker(iel) = MARK_CRS_4TRIA2TRIA_1
        ELSEIF(rhadapt%p_IvertexAge(Kvert(2)) .LE. 0) THEN
          p_Imarker(iel) = MARK_CRS_4TRIA2TRIA_2
        ELSEIF(rhadapt%p_IvertexAge(Kvert(3)) .LE. 0) THEN
          p_Imarker(iel) = MARK_CRS_4TRIA2TRIA_3
        ELSE
          p_Imarker(iel) = MARK_CRS_4TRIA1TRIA
        END IF
        
      CASE(STATE_TRIA_GREENINNER)
        ! If all vertices of the element are locked, then we can delete the
        ! element from the list of removable elements. Recall that both 
        ! vertices of the macro element are locked by definition so that it
        ! suffices to check the remaining (first) vertex.
        IF (rhadapt%p_IvertexAge(Kvert(1)) .LE. 0) THEN
          
          ! Delete element from list of removable elements              
          p_Imarker(iel) = MARK_ASIS
        ELSE

          ! Mark element for recoarsening together with its two neighbors
          p_Imarker(iel) = MARK_CRS_3TRIA1QUAD
        END IF
        
      CASE(STATE_TRIA_GREENOUTER_LEFT)
        ! If this is the left green triangle of a 1-tria : 2-tria refinement then
        ! get the element number of the corresponding right triangle and mark the
        ! element with the smaller number for refinement
        jel = rhadapt%p_IneighboursAtElement(2,iel)
        Kvert(1:3)  = rhadapt%p_IverticesAtElement(1:3,jel)
        jstate = redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))
        IF (jstate .EQ. STATE_TRIA_GREENOUTER_RIGHT) THEN
          
          ! Which is the element with smaller number?
          IF (iel < jel) THEN
            p_Imarker(iel) = MARK_CRS_2TRIA1TRIA_LEFT
          ELSE
            p_Imarker(iel) = MARK_CRS_2TRIA1TRIA_RIGHT
          END IF
        END IF
        
!!$      CASE(STATE_TRIA_GREENOUTER_RIGHT)
!!$        ! If this is the right green triangle of a 1-tria : 2-tria refinement then
!!$        ! get the element number of the corresponding left triangle and mark the
!!$        ! element with the smaller number for refinement
!!$        jel = rhadapt%p_IneighboursAtElement(2,iel)
!!$        Kvert(1:3)  = rhadapt%p_IverticesAtElement(1:3,jel)
!!$        jstate = redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))
!!$        IF (jstate .EQ. STATE_TRIA_GREENOUTER_LEFT) THEN
!!$          
!!$          ! Which is the element with smaller number?
!!$          IF (iel < jel) THEN
!!$            p_Imarker(iel) = MARK_CRS_2TRIA1TRIA_RIGHT
!!$          ELSE
!!$            p_Imarker(iel) = MARK_CRS_2TRIA1TRIA_LEFT
!!$          END IF
!!$        END IF
        
      CASE(STATE_TRIA_OUTERINNER)
        ! If this is the inner green tiangle of a 1-quad : 4-tria refinement then
        ! mark the element for 4-tria : 1-quad or 4-tria : 3-tria coarsening.
        jel = rhadapt%p_IneighboursAtElement(2,iel)
        KvertJ(1:3) = rhadapt%p_IverticesAtElement(1:3,jel)
        jstate = redgreen_getStateTria(rhadapt%p_IvertexAge(KvertJ(1:3)))
        
        ! Are we green triangle of a 1-quad : 4-tria refinement?
        IF (jstate .EQ. STATE_TRIA_OUTERINNER) THEN

          ! Which of elements IEL and JEL is the "most inner" one?
          kel = rhadapt%p_IneighboursAtElement(1,iel)

          IF (redgreen_getStateTria(rhadapt%p_IvertexAge(&
              rhadapt%p_IverticesAtElement(1:3,kel))).EQ.&
              STATE_TRIA_GREENOUTER_LEFT) THEN
            
            ! Element IEL is the "most inner" element
            IF (rhadapt%p_IvertexAge(Kvert(2)).LE.0) THEN
              p_Imarker(iel) = MARK_CRS_4TRIA3TRIA_RIGHT
            ELSEIF(rhadapt%p_IvertexAge(Kvert(3)).LE.0) THEN
              p_Imarker(iel) = MARK_CRS_4TRIA3TRIA_LEFT
            ELSE
              p_Imarker(iel) = MARK_CRS_4TRIA1QUAD
            END IF

          ELSE
            ! Element JEL is the "most inner" element
            IF (rhadapt%p_IvertexAge(KvertJ(2)).LE.0) THEN
              p_Imarker(jel) = MARK_CRS_4TRIA3TRIA_RIGHT
            ELSEIF(rhadapt%p_IvertexAge(KvertJ(3)).LE.0) THEN
              p_Imarker(jel) = MARK_CRS_4TRIA3TRIA_LEFT
            ELSE
              p_Imarker(jel) = MARK_CRS_4TRIA1QUAD
            END IF   
          END IF
        END IF

      CASE(STATE_QUAD_RED4)
        ! This is one of four quadrilaterals which result from a 1-quad : 4-quad
        ! refinement. Here, the situation is slightly more difficult because 
        ! multiple rotationally-symmetric situations have to be considered.
                
        ! As a first step, determine the element with the smallest number that
        ! originates from the resulting macro element. Starting from element
        ! IEL visit the "right" neighbor and compare its element number to
        ! the smallest element number currently available. Start with IEL.
        
        jel = iel
        kel = iel
        ishift = 0
        
        ! In the same loop, determine those vertices (exlcuding the interior
        ! vertex) which are connected to two adjacent quadilaterals that are locked
        IF (rhadapt%p_IvertexAge(Kvert(2)).LE.0) THEN
          iVertexLocked = ibset(0,0)
          icount        = 1
        ELSE
          iVertexLocked = 0
          icount        = 0
        END IF

        ! Ok, loop over all remaining elements (counterclockwise)
        DO ive=1,3
          
          ! Get counterclockwise neighboring element
          jel = rhadapt%p_IneighboursAtElement(2,jel)
          
          ! Is its second vertex locked?
          IF (rhadapt%p_IvertexAge(rhadapt%p_IverticesAtElement(2,jel)).LE.0) THEN
            iVertexLocked = ibset(iVertexLocked,ive)
            icount        = icount+1
          END IF
          
          ! Do we have a new smallest element number?
          IF (jel < kel) THEN
            kel    = jel
            ishift = ive
          END IF
        END DO

        ! How many vertices are locked?
        SELECT CASE(icount)
        CASE(0)
          ! Simply, mark element with smallest number for 4-quad : 1-quad coarsening
          p_Imarker(kel) = MARK_CRS_4QUAD1QUAD

        CASE(1)
          ! Ok, the four red quadrilaterals have to be converted into three triangles.
          ! Now, we have to find out the exact conversion routine. The bits of 
          ! iVertexLocked are set if and only if the vertex is locked. Moreover,
          ! the first bit corresponds to the second vertex of element IEL, the second
          ! bit corresponds to the second vertex of  the right neighbor of element 
          ! IEL and so on. The value ishift determines the relative position of the
          ! element with minimum number w.r.t. to IEL. That is, if IEL is the 
          ! element with the smallest number, then ishift=0. Let us perform a 
          ! cyclic bit-shift by the value ishift of iVertexLocked so that afterwards
          ! the first bit corresponds to the element with the minimum number.
          IF (ishift.NE.0) iVertexLocked = ISHFTC(iVertexLocked,ishift,TRIA_MAXNVE2D)

          ! Which vertex is locked?
          SELECT CASE(iVertexLocked)
          CASE(1)
            p_Imarker(kel) = MARK_CRS_4QUAD3TRIA_1

          CASE(2)
            p_Imarker(kel) = MARK_CRS_4QUAD3TRIA_2

          CASE(4)
            p_Imarker(kel) = MARK_CRS_4QUAD3TRIA_3

          CASE(8)
            p_Imarker(kel) = MARK_CRS_4QUAD3TRIA_4

          CASE DEFAULT
            PRINT *, "redgreen_mark_coarsening2D: Invalid locking."
            CALL sys_halt()
          END SELECT
          
        CASE(2)
          ! Ok, the four red quadrilaterals have to be converted into three triangles.
          ! Now, we have to find out the exact conversion routine. The bits of 
          ! iVertexLocked are set if and only if the vertex is locked. Moreover,
          ! the first bit corresponds to the second vertex of element IEL, the second
          ! bit corresponds to the second vertex of  the right neighbor of element 
          ! IEL and so on. The value ishift determines the relative position of the
          ! element with minimum number w.r.t. to IEL. That is, if IEL is the 
          ! element with the smallest number, then ishift=0. Let us perform a 
          ! cyclic bit-shift by the value ishift of iVertexLocked so that afterwards
          ! the first bit corresponds to the element with the minimum number.
          IF (ishift.NE.0) iVertexLocked = ISHFTC(iVertexLocked,ishift,TRIA_MAXNVE2D)
          
          ! Which vertices are locked?
          SELECT CASE(iVertexLocked)
          CASE(3)
            p_Imarker(kel) = MARK_CRS_4QUAD4TRIA_12
            
          CASE(6)
            p_Imarker(kel) = MARK_CRS_4QUAD4TRIA_23

          CASE(9)
            p_Imarker(kel) = MARK_CRS_4QUAD4TRIA_34

          CASE(12)
            p_Imarker(kel) = MARK_CRS_4QUAD4TRIA_14

          CASE(5)
            p_Imarker(kel) = MARK_CRS_4QUAD2QUAD_13

          CASE(10)
            p_Imarker(kel) = MARK_CRS_4QUAD2QUAD_24

          CASE DEFAULT
            PRINT *, "redgreen_mark_coarsening2D: Invalid locking."
            CALL sys_halt()
          END SELECT
          
        CASE DEFAULT
          PRINT *, "redgreen_mark_coarsening2D: Invalid number of locked vertices!"
          CALL sys_halt()
          
        END SELECT
        
      END SELECT
    END DO

10  rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_MARKEDCOARSEN)
  END SUBROUTINE redgreen_mark_coarsening2D

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE redgreen_mark_refinement2D(rhadapt,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine initializes tha adaptive data structure for red-green refinement.
    ! Starting from the marker array the neighbors of elements marked for refinement 
    ! are also marked for refinement until the resulting mesh satiesfies global
    ! conformity. This subroutine is implemented in an iterative fashion rather than
    ! making use of recursive subroutine calls. 
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
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER, DIMENSION(:), POINTER :: p_Imarker,p_Imodified
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D)  :: Kvert
    INTEGER(PREC_VERTEXIDX)  :: i,nvt
    INTEGER(PREC_ELEMENTIDX) :: nel,iel,jel,kel,lel,iel1,iel2
    INTEGER :: ive,jve,nve,mve,istate,jstate
    INTEGER :: h_Imodified,imodifier
    LOGICAL :: isConform
    
    !--------------------------------------------------------------------------
    ! At the moment, only those cells are marked for regular refinementfor which 
    ! the cell indicator does not satisfy some prescribed treshold. In order to 
    ! globally restore the conformity of the final grid, repeatedly loop over all
    ! elements and check if each of the three or four adjacent edges satisfies 
    ! conformity. Some care must be taken for green elements. In order to retain
    ! the shape quality of the elements, green elements need to be converted into
    ! red ones before further refinement is allowed.
    !--------------------------------------------------------------------------

    ! In the course of element marking, some green elements may have to be converted
    ! into regularly refined ones. In the worst case (Quad2Quad-refinement) each
    ! green element gives rise to 1.5 vertices and 2 new elements.
    ! Hence, the nodal arrays p_IvertexAge and p_InodalProperty as well as the
    ! element arrays p_IverticesAtElement, p_IneighboursAtElement,
    ! p_ImidneighboursAtElement and p_Imarker are enlarged precautionary. 
    IF (rhadapt%nGreenElements > 0) THEN

      ! Predict new dimensions for the worst case
      nvt = CEILING(rhadapt%NVT0+1.5*rhadapt%nGreenElements)
      nel = rhadapt%NEL0+2*rhadapt%nGreenElements

      ! Adjust nodal/elemental arrays
      CALL storage_realloc('redgreen_mark_refinement2D',nvt,&
          rhadapt%h_IvertexAge,ST_NEWBLOCK_ZERO,.TRUE.)
      CALL storage_realloc('redgreen_mark_refinement2D',nvt,&
          rhadapt%h_InodalProperty,ST_NEWBLOCK_ZERO,.TRUE.)
      CALL storage_realloc('redgreen_mark_refinement2D',nel,&
          rhadapt%h_Imarker,ST_NEWBLOCK_ZERO,.TRUE.)
      CALL storage_realloc('redgreen_mark_refinement2D',nel,&
          rhadapt%h_IverticesAtElement,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_realloc('redgreen_mark_refinement2D',nel,&
          rhadapt%h_IneighboursAtElement,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_realloc('redgreen_mark_refinement2D',nel,&
          rhadapt%h_ImidneighboursAtElement,ST_NEWBLOCK_NOINIT,.TRUE.)
      
      ! Reset pointers
      CALL storage_getbase_int(rhadapt%h_IvertexAge,rhadapt%p_IvertexAge)
      CALL storage_getbase_int(rhadapt%h_InodalProperty,&
          rhadapt%p_InodalProperty)
      CALL storage_getbase_int2D(rhadapt%h_IverticesAtElement,&
          rhadapt%p_IverticesAtElement)
      CALL storage_getbase_int2D(rhadapt%h_IneighboursAtElement,&
          rhadapt%p_IneighboursAtElement)
      CALL storage_getbase_int2D(rhadapt%h_ImidneighboursAtElement,&
          rhadapt%p_ImidneighboursAtElement)

      ! Adjust dimension of solution vector
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection))&
          CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_ADJUSTVERTEXDIM,(/nvt/),(/0/))

      ! Create new array for modifier
      CALL storage_new ('redgreen_mark_refinement2D','p_Imodified',nel,&
          ST_INT,h_Imodified,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_Imodified,p_Imodified)
    
    ELSE
      
      ! No green elements have to be considered, hence use NEL0      
      CALL storage_new('redgreen_mark_refinement2D','p_Imodified',rhadapt%NEL0,&
          ST_INT,h_Imodified,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_Imodified,p_Imodified)

    END IF

    ! The task of the two arrays p_Imarker and p_Imodified are as follows:
    ! p_Imarker stores for each element its individual state from which the task of
    ! refinement can be uniquely determined (see above; data type t_hadapt).
    ! If the state of one element is modified, then potentionally all neighboring
    ! elements have to be checked, if conformity is violated. In order to prevent this
    ! (expensive) check to be performed for all elements again and again, only those
    ! elements are investigated for which the modifier p_Imodified is active, i.e.
    ! p_Imodified(iel) == imodifier. In order to pre-select elements for investigation
    ! in the next loop, p_Imodified(iel) is set to -imodifier and the modifier is reversed
    ! after each iteration. This complicated step is necessary for the following reason:
    ! If one element with mediate number initially is not marked but gets marked for
    ! subdivision at one edge, then it will be considered in the same iteration and
    ! green elements will eventually converted. However, it is possible that the same
    ! element would be marked for subdivision at another edge due to an adjacent element
    ! with large number. Hence, it is advisable to first process all elements which are
    ! activated and pre-select their neighbors for the next iteration.

    ! Ok, so let's start. Initially, indicate all elements as modified which are marked
    ! for refinement due to accuracy reasons.
    imodifier=1
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
    DO iel=1,rhadapt%NEL0
      IF (p_Imarker(iel) .NE. MARK_ASIS) p_Imodified(iel) = imodifier
    END DO

    isConform = .FALSE.
    conformity: DO
      
      ! If conformity is guaranteed for all cells, then exit
      IF (isConform) EXIT conformity
      isConform = .TRUE.
      
      ! Otherwise, loop over all elements present in the initial grid (IEL <= NEL0)
      ! which are modified and mark their neighbors for further refinement 
      ! if conformity is violated for some edge     
      DO iel=1,rhadapt%NEL0
        
        ! Skip those element which have not been modified
        IF (p_Imodified(iel) .NE. imodifier) CYCLE

        ! Get number of vertices per element
        nve = get_NVE(rhadapt,iel)
        
        ! Get local data for element iel
        Kvert(1:nve)   = rhadapt%p_IverticesAtElement(1:nve,iel)
      
        ! Are we triangular or quadrilateral element?
        SELECT CASE(nve)
        CASE(TRIA_NVETRI2D)
          p_Imarker(iel) = ibclr(p_Imarker(iel),0)
          istate         = redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))

        CASE(TRIA_NVEQUAD2D)
          p_Imarker(iel) = ibset(p_Imarker(iel),0)
          istate         = redgreen_getStateQuad(rhadapt%p_IvertexAge(Kvert(1:4)))

        CASE DEFAULT
          PRINT *, "redgreen_mark_refinement2D: Invalid number of vertices per element!"
          CALL sys_halt()
        END SELECT

        !------------------------------------------------------------
        ! Check if the state of the current element IEL allows direct
        ! refinement or if the elements needs to be "converted" first
        !------------------------------------------------------------
        SELECT CASE(istate)
        CASE (STATE_TRIA_ROOT,STATE_QUAD_ROOT,STATE_TRIA_REDINNER,&
              STATE_QUAD_RED1,STATE_QUAD_RED2,STATE_QUAD_RED3,STATE_QUAD_RED4)
          ! States which can be directly accepted

        CASE(STATE_TRIA_OUTERINNER1,STATE_TRIA_OUTERINNER2)
          ! Theoretically, these states may occure and have to be treated 
          ! like CASE(4), see below. Due to the fact, that these states can
          ! only be generated by means of local red-green refinement and are
          ! not present in the initial grid we make use of additional know-
          ! ledge: Triangles which have too youngest vertices can either
          ! be outer red elements resulting from a Tria4Tria refinement 
          ! or one of the opposite triangles resulting from Quad4Tria refinement.
          ! In all cases, the edge which connects the two nodes is opposite
          ! to the first local vertex. Hence, work must only be done for CASE(4)
          PRINT *, "redgreen_mark_refinement2D: State 2 or 8 must not occur"
          CALL sys_halt()
          
        CASE(STATE_TRIA_OUTERINNER)
          ! The triangle can either be an outer red element resulting from a Tria4Tria 
          ! refinement or one of the two opposite triangles resulting from Quad4Tria 
          ! refinement. This can be easily checked. If the opposite element is not 
          ! an inner red triangle (14), then the element IEL and its opposite neighbor 
          ! make up the inner "diamond" resulting from a Quad4Tria refinement. 
          jel = rhadapt%p_IneighboursAtElement(2,iel)
       
          ! First, we need to check a special case. If the second edge of element IEL
          ! has two different neighbors, then the original neighbor along this edge
          ! was a green triangle that has been converted into an (outer) red one. In
          ! this case, the current element IEL does not need further modifications.
          IF (jel .NE. rhadapt%p_ImidneighboursAtElement(2,iel)) GOTO 100
          
          ! Otherwise, determine the state of the edge neighbor JEL
          jstate = redgreen_getState(rhadapt,jel)
          
          IF (jstate .EQ. STATE_TRIA_OUTERINNER) THEN
            ! We know that element IEL and JEL make up the inner "diamond" resulting
            ! from a Quad4Tria refinement. Hence, the element with larger number must
            ! be the inner one so that conversion into four quadrilaterals can be done.
            IF (iel > jel) THEN
              
              ! To begin with, we need to find the two missing triangles KEL and LEL
              ! which make up the original quadrilateral
              kel = rhadapt%p_IneighboursAtElement(1,iel)
              lel = rhadapt%p_IneighboursAtElement(3,iel)

              ! Mark the edge of the element adjacent to KEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(3,kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3,kel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)
              
              ! Mark the edge of the element adjacent to LEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(1,lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,lel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

              ! Now, we can physically convert the four triangles into four quadrilaterals
              CALL convert_Quad4Tria(rhadapt,kel,jel,lel,iel,rcollection,fcb_hadaptCallback)
              isConform=.FALSE.

              ! All four elements have to be converted from triangles to quadrilaterals.
              p_Imarker(iel) = ibset(0,0)

              ! The second and third edge of KEL must be unmarked. Moreover, the first edge is
              ! marked for refinement if and only if it is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(1,kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
              IF (ismarked_edge(iel1,iel2,kel)) THEN
                p_Imarker(kel) = ibset(ibset(0,1),0)
              ELSE
                p_Imarker(kel) = ibset(0,0)
              END IF

              ! The second and third edge of JEL must be unmarked. 
              p_Imarker(jel) = ibset(0,0)

              ! The fourth edge is only marked if it is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(4,jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),4)
              
              ! The first edge is only marked if it is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(1,jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)
              
              ! The first, second and third edge of LEL must be unmarked.
              p_Imarker(lel) = ibset(0,0)

              ! The fourth edge is only marked if it is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(4,lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4,lel)
              IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)
              
            ELSE

              ! To begin with, we need to find the two missing triangles KEL and LEL
              ! which make up the original quadrilateral
              kel  = rhadapt%p_IneighboursAtElement(1,jel)
              lel  = rhadapt%p_IneighboursAtElement(3,jel)

              ! Mark the edge of the element adjacent to KEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(3,kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3,kel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)

              ! Mark the edge of the element adjacent to LEL for subdivision
              iel1 = rhadapt%p_IneighboursAtElement(1,lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,lel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

              ! Now, we can physically convert the four triangles into four quadrilaterals
              CALL convert_Quad4Tria(rhadapt,kel,iel,lel,jel,rcollection,fcb_hadaptCallback)
              isConform=.FALSE.

              ! All four elements have to be converted from triangles to quadrilaterals.
              p_Imarker(jel) = ibset(0,0)

              ! The second and third edge of KEL must be unmarked. Moreover, the first edge is
              ! marked for refinement if and only if it is also marked from the adjacent element.
              iel1 = rhadapt%p_IneighboursAtElement(1,kel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
              IF (ismarked_edge(iel1,iel2,kel)) THEN
                p_Imarker(kel) = ibset(ibset(0,1),0)
              ELSE
                p_Imarker(kel) = ibset(0,0)
              END IF

              ! The second and third edge of IEL must be unmarked.
              p_Imarker(iel) = ibset(0,0)

              ! The fourth edge is only marked if it is also marked from the adjacent element
              iel1 = rhadapt%p_IneighboursAtElement(4,iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),4)
              
              ! The first edge is only marked if it is also marked from the adjacent element              
              iel1 = rhadapt%p_IneighboursAtElement(1,iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)
                      
              ! The first, second and third edge of LEL must be unmarked.
              p_Imarker(lel) = ibset(0,0)
              
              ! The fourth edge is only marked if it is also marked from the adjacent element              
              iel1 = rhadapt%p_IneighboursAtElement(4,lel)
              iel2 = rhadapt%p_ImidneighboursAtElement(4,lel)
              IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)

            END IF
          END IF
                    
        CASE(STATE_QUAD_HALF1)
          ! Element is quadrilateral that results from a Quad2Quad refinement.
          ! Due to our orientation convention the other "halved" element is 
          ! adjacent to the second edge. 
          jel = rhadapt%p_IneighboursAtElement(2,iel)

          ! Mark the fourth edge of elements IEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(4,iel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4,iel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

          ! Mark the fourth edge of elements JEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(4,jel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4,jel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
          
          ! Now, we can physically convert the two quadrilaterals into four quadrilaterals
          CALL convert_Quad2Quad(rhadapt,iel,jel,rcollection,fcb_hadaptCallback)
          isConform=.FALSE.

          ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
          ! The conversion routine checks element IEL and JEL and keeps the element
          ! with smaller number at its "corner". That is, the "first" vertex of the
          ! element with smaller number remains the "first" vertex whereas the 
          ! "first" vertex of the element with larger number will be converted to
          ! the second "vertex". This is due to the fact, that all corners vertices
          ! are the "first" vertex if the correcponding element. Hence, we have to
          ! check if IEL<JEL or vice versa and transfer some markers from elements
          ! IEL, JEL to the new elements NEL0+1, NEL0+2. In addition, we have to
          ! remove some markers from elements IEL and JEL afterwards(!!!)
          p_Imarker(iel) = ibset(0,0)
          p_Imarker(jel) = ibset(0,0)
          p_Imarker(kel) = ibset(0,0)
          p_Imarker(lel) = ibset(0,0)
          
          ! Which element is smaller?
          IF (iel < jel) THEN
            kel =  rhadapt%p_IneighboursAtElement(2,jel)
            lel =  rhadapt%p_IneighboursAtElement(3,iel)
            
            iel1 = rhadapt%p_IneighboursAtElement(1,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)

            iel1 = rhadapt%p_IneighboursAtElement(4,jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,jel)
            IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),4)

            iel1 = rhadapt%p_IneighboursAtElement(1,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            iel1 = rhadapt%p_IneighboursAtElement(4,lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,lel)
            IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)

          ELSE
            kel =  rhadapt%p_IneighboursAtElement(2,iel)
            lel =  rhadapt%p_IneighboursAtElement(3,jel)

            iel1 = rhadapt%p_IneighboursAtElement(1,jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,jel)
            IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)

            iel1 = rhadapt%p_IneighboursAtElement(4,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),4)

            iel1 = rhadapt%p_IneighboursAtElement(1,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            iel1 = rhadapt%p_IneighboursAtElement(4,lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,lel)
            IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)
          END IF
          
          ! Finally, mark as quadrilaterals
          p_Imarker(lel) = ibset(p_Imarker(lel),0)
          p_Imarker(kel) = ibset(p_Imarker(kel),0)
         
        CASE(STATE_QUAD_HALF2)
          ! Theoretically, this state is not possible. In general, it has to be treated
          ! like state 21 = STATE_QUAD_HALF1
          PRINT *, "redgreen_mark_refinement2D: State 11 must not occur"
          CALL sys_halt()
          
        CASE(STATE_TRIA_GREENINNER)
          ! We are processing a green triangle. Due to our refinement convention, element IEL
          ! can only be the inner triangle resulting from a Quad3Tria refinement.
          jel = rhadapt%p_IneighboursAtElement(3,iel)
          kel = rhadapt%p_IneighboursAtElement(1,iel)

          ! To begin with, we need to mark the edges of the adjacent elements for subdivision.
          iel1 = rhadapt%p_IneighboursAtElement(2,iel)
          iel2 = rhadapt%p_ImidneighboursAtElement(2,iel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
          
          ! Mark the edge of the element adjacent to JEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(3,jel)
          iel2 = rhadapt%p_ImidneighboursAtElement(3,jel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)

          ! Mark the edge of the element adjacent to KEL for subdivision
          iel1 = rhadapt%p_IneighboursAtElement(1,kel)
          iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)

          ! Now, we can physically convert the three elements IEL,JEL and KEL
          ! into four similar quadrilaterals
          CALL convert_Quad3Tria(rhadapt,jel,kel,iel,rcollection,fcb_hadaptCallback)
          isConform=.FALSE.
          
          ! The new element NEL0+1 has zero marker by construction but that of IEL must
          ! be nullified. The markers for the modified triangles also need to be adjusted.
          p_Imarker(iel) = ibset(0,0)
          p_Imarker(jel) = ibset(0,0)
          p_Imarker(kel) = ibset(0,0)

          ! The first edge is only marked if it is also marked from the adjacent element              
          iel1 = rhadapt%p_IneighboursAtElement(1,jel)
          iel2 = rhadapt%p_ImidneighboursAtElement(1,jel)
          IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)
          
          ! The fourth edge is only marked if it is also marked from the adjacent element              
          iel1 = rhadapt%p_IneighboursAtElement(4,kel)
          iel2 = rhadapt%p_ImidneighboursAtElement(4,kel)
          IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)
          
        CASE(STATE_TRIA_GREENOUTER_LEFT)
          ! We are processing a green triangle. Here, we have to consider several cases.
          ! First, let us find out the state of the neighboring element JEL.
          jel    = rhadapt%p_IneighboursAtElement(2,iel)
          jstate = redgreen_getState(rhadapt,jel)

          ! What state is element JEL
          SELECT CASE(jstate)
          CASE(STATE_TRIA_GREENOUTER_RIGHT)
            ! Element IEL and JEL are the result of a Tria2Tria refinement, whereby
            ! triangle IEL is located left to element JEL. We can safely convert both
            ! elements into one and perform regular refinement afterwards.

            ! To begin with, we need to mark the edges of the elements adjacent 
            ! to IEL and JEL for subdivision. 
            iel1 = rhadapt%p_IneighboursAtElement(3,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
            
            ! The same procedure must be applied to the neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(1,jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
            
            ! Now, we can physically convert the two elements IEL and JEL into four similar triangles
            CALL convert_Tria2Tria(rhadapt,iel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
            ! The markers for the modified elements IEL and JEL need to be adjusted.
            p_Imarker(iel) = 0
            p_Imarker(jel) = 0
            
            ! Which is the smaller elements?
            IF (iel < jel) THEN
              iel1 = rhadapt%p_IneighboursAtElement(1,iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)
              
              iel1 = rhadapt%p_IneighboursAtElement(3,jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),3)
            ELSE
              iel1 = rhadapt%p_IneighboursAtElement(1,jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)
              
              iel1 = rhadapt%p_IneighboursAtElement(3,iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),3)
            END IF
          
          CASE(STATE_TRIA_GREENINNER)
            ! Element IEL and JEL are the result of a Quad3Tria refinement, whereby
            ! element JEL is the inner triangle and IEL is its left neighbor.

            ! To begin with, we need to mark the edges of the elements adjacent to IEL,
            ! JEL and KEL for subdivision. Let us start with the third neighbor of IEL.
            iel1 = rhadapt%p_IneighboursAtElement(3,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

            ! The same procedure must be applied to the second neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(2,jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(2,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
            
            ! And again, the same procedure must be applied to the first neighbor of 
            ! element KEL which is the first neighbor of JEL, whereby KEL needs to 
            ! be found in the dynamic data structure first
            kel = rhadapt%p_IneighboursAtElement(1,jel)
            
            ! Ok, now we can proceed to its first neighbor
            iel1 = rhadapt%p_IneighboursAtElement(1,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)
            
            ! Now, we can physically convert the three elements IEL,JEL and KEL
            ! into four similar quadrilaterals
            CALL convert_Quad3Tria(rhadapt,iel,kel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new element NEL0+1 has zero marker by construction but that of JEL must 
            ! be nullified. The markers for the modified triangles also need to be adjusted.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)
           
          CASE(STATE_TRIA_OUTERINNER)
            ! Element IEL and JEL are the result of a Quad4Tria refinement, whereby
            ! element JEL is the inner triangles  and IEL is its right neighbor.
            kel = rhadapt%p_IneighboursAtElement(2,jel)
            lel = rhadapt%p_IneighboursAtElement(3,jel)

            ! Mark the edge of the element adjacent to IEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(3,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
            
            ! Mark the edge of the element adjacent to LEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(1,lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,lel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

            ! Now, we can physically convert the four triangles into four quadrilaterals
            CALL convert_Quad4Tria(rhadapt,iel,kel,lel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! All four elements have to be converted from triangles to quadrilaterals.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)
            p_Imarker(lel) = ibset(0,0)
            
            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4,lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,lel)
            IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel) = ibset(p_Imarker(lel),4)

          CASE DEFAULT
            PRINT *, "redgreen_mark_refinement2D: Invalid state",istate,iel
            CALL sys_halt()
          END SELECT

        CASE(STATE_TRIA_GREENOUTER_RIGHT)
          ! We are processing a green triangle. Here, we have to consider several cases.
          ! First, let us find out the state of the neighboring element JEL.
          jel    = rhadapt%p_IneighboursAtElement(2,iel)
          jstate = redgreen_getState(rhadapt,jel)
          
          ! What state is element JEL
          SELECT CASE(jstate)
          CASE(STATE_TRIA_GREENOUTER_LEFT)
            ! Element IEL and JEL are the result of a Tria2Tria refinement, whereby
            ! triangle IEL is located right to element JEL. We can safely convert both
            ! elements into one and perform regular refinement afterwards.
            
            ! To begin with, we need to mark the edges of the elements adjacent
            ! to IEL and JEL for subdivision.
            iel1 = rhadapt%p_IneighboursAtElement(1,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
            
            ! The same procedure must be applied to the neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(3,jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
            
            ! Now, we can physically convert the two elements IEL and JEL into four similar triangles
            CALL convert_Tria2Tria(rhadapt,jel,iel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
            ! The markers for the modified elements IEL and JEL need to be adjusted
            p_Imarker(iel) = 0
            p_Imarker(jel) = 0
            
            ! Whis is the smaller element?
            IF (iel < jel) THEN
              iel1 = rhadapt%p_IneighboursAtElement(1,iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),1)
              
              iel1 = rhadapt%p_IneighboursAtElement(3,jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),3)
            ELSE
              iel1 = rhadapt%p_IneighboursAtElement(1,jel)
              iel2 = rhadapt%p_ImidneighboursAtElement(1,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel) = ibset(p_Imarker(jel),1)
              
              iel1 = rhadapt%p_IneighboursAtElement(3,iel)
              iel2 = rhadapt%p_ImidneighboursAtElement(3,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),3)
            END IF


          CASE(STATE_TRIA_GREENINNER)
            ! Element IEL and JEL are the result of a Quad3Tria refinement, whereby
            ! element JEL is the inner triangle and IEL is its left neighbor.

            ! To begin with, we need to mark the edges of the elements adjacent to IEL,
            ! JEL and KEL for subdivision. Let us start with the firth neighbor of IEL.
            iel1 = rhadapt%p_IneighboursAtElement(1,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

            ! The same procedure must be applied to the second neighbor of element JEL
            iel1 = rhadapt%p_IneighboursAtElement(2,jel)
            iel2 = rhadapt%p_ImidneighboursAtElement(2,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)

            ! And again, the same procedure must be applied to the third neighbor of 
            ! element KEL which is the third neighbor of JEL, whereby KEL needs to 
            ! be found in the dynamic data structure first
            kel = rhadapt%p_IneighboursAtElement(3,jel)

            ! Ok, now we can proceed to its first neighbor
            iel1 = rhadapt%p_IneighboursAtElement(3,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3,kel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)

            ! Now, we can physically convert the three elements IEL,JEL and KEL
            ! into four similar quadrilaterals
            CALL convert_Quad3Tria(rhadapt,kel,iel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new element NEL0+1 has zero marker by construction but that of JEL must
            ! be nullified. The markers for the modified triangles also need to be adjusted.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),4)
            
          CASE(STATE_TRIA_OUTERINNER)
            ! Element IEL and JEL are the result of a Quad4Tria refinement, whereby
            ! element JEL is the inner triangles  and IEL is its left neighbor.
            kel = rhadapt%p_IneighboursAtElement(2,jel)
            lel = rhadapt%p_IneighboursAtElement(1,jel)

            ! Mark the edge of the element adjacent to IEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(1,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

            ! Mark the edge of the element adjacent to LEL for subdivision
            iel1 = rhadapt%p_IneighboursAtElement(3,lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(3,lel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

            ! Now, we can physically convert the four triangles into four quadrilaterals
            CALL convert_Quad4Tria(rhadapt,lel,kel,iel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! All four elements have to be converted from triangles to quadrilaterals.
            p_Imarker(iel) = ibset(0,0)
            p_Imarker(jel) = ibset(0,0)
            p_Imarker(kel) = ibset(0,0)
            p_Imarker(lel) = ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1,lel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,lel)
            IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel) = ibset(p_Imarker(lel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),4)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(1,kel)
            iel2 = rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel) = ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1 = rhadapt%p_IneighboursAtElement(4,iel)
            iel2 = rhadapt%p_ImidneighboursAtElement(4,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel) = ibset(p_Imarker(iel),4)

          CASE DEFAULT
            PRINT *, "redgreen_mark_refinement2D: Invalid state",istate,iel
            CALL sys_halt()
          END SELECT
          
        CASE DEFAULT
          PRINT *, "redgreen_mark_refinement2D: Invalid state",istate,iel
          CALL sys_halt()
        END SELECT
100     CONTINUE

        !------------------------------------------------------------
        ! Now, we can safely loop over adjacent cells of element IEL
        !------------------------------------------------------------
        DO ive=1,nve

          ! If element IEL is adjacent to two different elements along one edge
          ! then ignore this edge. This situation can arise, if the element IEL
          ! is adjacent to a green element JEL which has been converted previously,
          ! so that element IEL has two neighbors along the edge, i.e., IEL and one 
          ! new element number >NEL0
          IF (rhadapt%p_IneighboursAtElement(ive,iel) .NE. &
              rhadapt%p_ImidneighboursAtElement(ive,iel)) CYCLE
          
          ! If the edge shared by element IEL and its neighbor JEL has been 
          ! marked for refinement and JEL is not outside of the domain, i.e.,
          ! element IEL is located at the boundary, then proceed to element JEL.
          IF (btest(p_Imarker(iel),ive)) THEN
            
            ! Check if the current edge is located at the boundary, 
            ! than nothing needs to be done
            jel = rhadapt%p_IneighboursAtElement(ive,iel)
            IF (jel .EQ. 0) CYCLE
            
            ! Get number of vertices per element
            mve = get_NVE(rhadapt,jel)
            
            ! Are we triangular or quadrilateral element?
            SELECT CASE(mve)
            CASE(TRIA_NVETRI2D)
              p_Imarker(jel) = ibclr(p_Imarker(jel),0)

            CASE(TRIA_NVEQUAD2D)
              p_Imarker(jel) = ibset(p_Imarker(jel),0)

            CASE DEFAULT
              PRINT *, "redgreen_mark_refinement2D: Invalid number of vertices per element!"
              CALL sys_halt()
            END SELECT

            ! Now, we need to find the local position of element IEL in the
            ! adjacency list of the nieghboring element JEL
            DO jve=1,mve
              IF (rhadapt%p_IneighboursAtElement(jve,jel) == iel) EXIT
            END DO

            IF (jve > mve) THEN
              PRINT *, "redgreen_mark_refinement2D: Unable to find element!"
              CALL sys_halt()
            END IF

            ! If the edge is already marked for refinement then we can
            ! guarantee conformity for the edge shared by IEL and JEL.
            ! Otherwise, this edge needs to be marked "looking" from JEL
            IF (.NOT.btest(p_Imarker(jel),jve)) THEN
              p_Imarker(jel) = ibset(p_Imarker(jel),jve)
              isConform      = .FALSE.
              IF (p_Imodified(jel) .NE. imodifier) p_Imodified(jel)= -imodifier
            END IF

            ! Finally, "lock" all vertices of element JEL
            Kvert(1:mve) = rhadapt%p_IverticesAtElement(1:mve,jel)
            rhadapt%p_IvertexAge(Kvert(1:mve)) = -ABS(rhadapt%p_IvertexAge(Kvert(1:mve)))
          END IF
        END DO
      END DO
      
      ! Note that the number of elements (and vertices) may have changed
      ! due to conversion of green into red elements. Hence, adjust dimensions.
      rhadapt%NEL0 = rhadapt%NEL
      rhadapt%NVT0 = rhadapt%NVT

      !--------------------------------------------------------------
      ! We don't want to have elements with two and/or three divided edges, 
      ! aka, blue refinement for triangles and quadrilaterals, respectively.
      ! As a remedy, these elements are filtered and marked for red refinement.
      !--------------------------------------------------------------
      DO iel=1,rhadapt%NEL0
        IF (p_Imodified(iel)*imodifier == 0) CYCLE
        
        ! What type of element are we?
        SELECT CASE(p_Imarker(iel))
        CASE(MARK_REF_TRIA3TRIA_12,MARK_REF_TRIA3TRIA_13,MARK_REF_TRIA3TRIA_23)
          ! Blue refinement for triangles is not allowed.
          ! Hence, mark triangle for red refinement
          p_Imarker(iel)         =  MARK_REF_TRIA4TRIA
          p_Imodified(iel)       = -imodifier
          isConform              = .FALSE.
          rhadapt%increaseNVT    =  rhadapt%increaseNVT+1
          
        CASE(MARK_REF_QUADBLUE_123,MARK_REF_QUADBLUE_412,&
            MARK_REF_QUADBLUE_341,MARK_REF_QUADBLUE_234)
          ! Blue refinement for quadrilaterals is not allowed. 
          ! Hence, mark quadrilateral for red refinement
          p_Imarker(iel)         =  MARK_REF_QUAD4QUAD
          p_Imodified(iel)       = -imodifier
          isConform              = .FALSE.
          rhadapt%increaseNVT    =  rhadapt%increaseNVT+2

        CASE DEFAULT
          p_Imodified(iel) = (p_Imodified(iel)-imodifier)/2
        END SELECT
      END DO
      
      ! Reverse modifier
      imodifier = -imodifier
    END DO conformity

    ! As a last step, lock all vertices of those elements which are marked for refinement
    DO iel=1,SIZE(p_Imarker)
      
      ! What type of element are we?
      SELECT CASE(p_Imarker(iel))
      CASE(MARK_REF_TRIA2TRIA_1)
        i = rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
        
      CASE(MARK_REF_TRIA2TRIA_2)
        i = rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_TRIA2TRIA_3)
        i = rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_QUAD3TRIA_1)
        i = rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
        i = rhadapt%p_IverticesAtElement(4,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
        
      CASE(MARK_REF_QUAD3TRIA_2)
        i = rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
        i = rhadapt%p_IverticesAtElement(4,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_QUAD3TRIA_3)
        i = rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
        i = rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_QUAD3TRIA_4)
        i = rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
        i = rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i) = -ABS(rhadapt%p_IvertexAge(i))
      END SELECT
    END DO
      
    ! Free auxiliary storage
    CALL storage_free(h_Imodified)

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! For a given element IEL, mark the edge that connects IEL 
    ! to its neighbor JEL in the marker array at position JEL
    
    SUBROUTINE mark_edge(iel,jel)

      INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
      INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: jel

      ! local variables
      INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D)  :: Kvert
      INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE2D) :: Kadj
      INTEGER :: ive,nve

      ! Get number of vertices per element
      nve = get_NVE(rhadapt,jel)

      ! Get local data
      Kvert(1:nve) = rhadapt%p_IverticesAtElement(1:nve,jel)
      Kadj(1:nve)  = rhadapt%p_IneighboursAtElement(1:nve,jel)
      
      ! Find local position of element IEL in adjacency list of JEL
      SELECT CASE(nve)
      CASE(TRIA_NVETRI2D)
        ! Triangular elements
        DO ive=1,nve
          IF (Kadj(ive) .EQ. iel) THEN
            p_Imarker(jel) = ibclr(ibset(p_Imarker(jel),ive),0)
            isConform      =.FALSE.
            IF (p_Imodified(jel) .NE. imodifier) p_Imodified(jel) = -imodifier
            EXIT
          END IF
        END DO

      CASE(TRIA_NVEQUAD2D)
        ! Quadrilateral element
        DO ive=1,nve
          IF (Kadj(ive) .EQ. iel) THEN
            p_Imarker(jel) = ibset(ibset(p_Imarker(jel),ive),0)
            isConform      =.FALSE.
            IF (p_Imodified(jel) .NE. imodifier) p_Imodified(jel) = -imodifier
            EXIT
          END IF
        END DO

      CASE DEFAULT
        PRINT *, "mark_edge: Invalid number of vertices per element!"
        CALL sys_halt()
      END SELECT
    END SUBROUTINE mark_edge

    !**************************************************************
    ! For a given element IEL, check if the edge that connects IEL
    ! with its adjacent element JEL is marked

    FUNCTION ismarked_edge(iel,ielmid,jel) RESULT (bismarked)

      INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel,ielmid,jel
      LOGICAL :: bismarked

      ! local variables
      INTEGER :: ive,nve

      ! Are we at the boundary?
      IF (iel*ielmid .EQ.0) THEN
        bismarked=.FALSE.
        RETURN
      ELSEIF(iel .NE. ielmid) THEN
        bismarked=.TRUE.
        RETURN
      END IF
      
      ! Loop over all edges of element IEL and try to find neighbor JEL
      DO ive=1,get_NVE(rhadapt,iel)
        IF (rhadapt%p_IneighboursAtElement(ive,iel) .EQ. jel) THEN
          bismarked = btest(p_Imarker(iel),ive)
          RETURN
        END IF
      END DO

      PRINT *, "ismarked_egde: Unable to find common edge!"
      CALL sys_halt()
    END FUNCTION ismarked_edge
  END SUBROUTINE redgreen_mark_refinement2D
  
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
      PRINT *, "redgreen_refine: dynamic data structures are not generated &
          &or no marker for refinement is available!"
      CALL sys_halt()
    END IF
    
    ! Set pointers
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
        
    ! Perform red-green refinement
    DO iel=1,SIZE(p_Imarker)
      
      SELECT CASE(p_Imarker(iel))
      CASE(MARK_ASIS_TRIA,MARK_ASIS_QUAD)
        ! Do nothing for elements that should be kept 'as is'

      CASE(:MARK_CRS_GENERIC)
        ! Do nothing for element that have been marked for coarsening

      CASE(MARK_REF_TRIA4TRIA)
        ! Red refinement triangle
        CALL refine_Tria4Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)
        
      CASE(MARK_REF_QUAD4QUAD)
        ! Red refinement quadrilateral
        CALL refine_Quad4Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)
        
      CASE(MARK_REF_TRIA2TRIA_1,MARK_REF_TRIA2TRIA_2,MARK_REF_TRIA2TRIA_3)   
        ! Green refinement triangle
        CALL refine_Tria2Tria(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2
        
      CASE(MARK_REF_QUAD3TRIA_1,MARK_REF_QUAD3TRIA_2,&
           MARK_REF_QUAD3TRIA_3,MARK_REF_QUAD3TRIA_4)
        ! Green refinement quadrilateral
        CALL refine_Quad3Tria(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+3
        
      CASE(MARK_REF_QUAD2QUAD_13,MARK_REF_QUAD2QUAD_24)
        ! Green refinement quadrilateral
        CALL refine_Quad2Quad(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+2

      CASE(MARK_REF_QUAD4TRIA_12,MARK_REF_QUAD4TRIA_23,&
           MARK_REF_QUAD4TRIA_34,MARK_REF_QUAD4TRIA_14)
        ! Green refinement quadrilateral
        CALL refine_Quad4Tria(rhadapt,iel,p_Imarker(iel),rcollection,fcb_hadaptCallback)
        rhadapt%nGreenElements = rhadapt%nGreenElements+4

      CASE DEFAULT
        PRINT *, "redgreen_refine: Invalid refinement marker!"
        CALL sys_halt()
      END SELECT
    END DO

    ! Increase the number of refinement steps by one
    rhadapt%nRefinementSteps = rhadapt%nRefinementSteps+1
    
    ! Refinement has been performed.
    rhadapt%iSpec = IOR(rhadapt%ispec,HADAPT_REFINED)
    
    ! Hence, the markers are no longer valid
    rhadapt%iSpec = IAND(rhadapt%iSpec,NOT(HADAPT_MARKEDREFINE))
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

!<inputoutptut>
    ! adativity structure
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt
    
    ! OPTIONAL Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER,  DIMENSION(:), POINTER :: p_Imarker
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: iel,jel
    INTEGER(PREC_VERTEXIDX)  :: ivt,ivtReplace
    INTEGER :: ive

    ! Check if dynamic data structures are o.k. and if 
    ! cells are marked for coarsening
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA .OR.&
        IAND(rhadapt%iSpec,HADAPT_MARKEDCOARSEN).NE.HADAPT_MARKEDCOARSEN) THEN
      PRINT *, "redgreen_coarsen: dynamic data structures are not generated &
          &or no marker for coarsening is available!"
      CALL sys_halt()
    END IF
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
    
    ! Perform hierarchical red-green recoarsening
    DO iel=SIZE(p_Imarker),1,-1
      
      SELECT CASE(p_Imarker(iel))
      CASE(MARK_CRS_GENERIC:)
        ! Do nothing for elements ...
        ! - that should be kept 'as is'
        ! - that are only marked for generic recoarsening
        ! - that are marked for refinement.

      CASE(MARK_CRS_2TRIA1TRIA_LEFT)
        CALL coarsen_2Tria1Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)
        
      CASE(MARK_CRS_4TRIA1TRIA)
        CALL coarsen_4Tria1Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_4TRIA2TRIA_1)
        CALL coarsen_4Tria2Tria(rhadapt,iel,1,rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_4TRIA2TRIA_2)
        CALL coarsen_4Tria2Tria(rhadapt,iel,2,rcollection,fcb_hadaptCallback)

      CASE(MARK_CRS_4TRIA2TRIA_3)
        CALL coarsen_4Tria2Tria(rhadapt,iel,3,rcollection,fcb_hadaptCallback)

      CASE DEFAULT
        PRINT *, "redgreen_coarsen: Invalid recoarsening marker!",p_Imarker(iel)
        CALL sys_halt()
      END SELECT
    END DO

    PRINT *, "Elemental removal is done"
    
    ! Loop over all vertices present in the triangulation before refinement
    ! and check if they are free for vertex removal.
    DO ivt=rhadapt%NVT0,1,-1
      IF (rhadapt%p_IvertexAge(ivt) .GT. 0) THEN

        PRINT *, "Removing",ivt,"..."
        CALL remove_vertex2D(rhadapt,ivt,ivtReplace)

        ! If vertex has been replace, update vertex at element list
        IF (ivtReplace .NE. 0) THEN

          PRINT *, "... is replaced by",ivtReplace

          CALL arrlst_printArrayList(rhadapt%rElementsAtVertex,ivtReplace)

          ! Start with first element in list
          ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivtReplace,.TRUE.)
          DO WHILE(ipos .NE. ANULL)
            jel=rhadapt%rElementsAtVertex%IData(ipos)

            print *, "jel=",jel, rhadapt%p_IverticesAtElement(:,jel)

            ! Look for vertex ivtReplace in element JEL and replace it by IVT
            DO ive=1,get_NVE(rhadapt,jel)
              IF (rhadapt%p_IverticesAtElement(ive,jel) .EQ. ivtReplace) THEN
                rhadapt%p_IverticesAtElement(ive,jel) = ivt
                GOTO 99
              END IF
            END DO
            
            PRINT *, "NOT FOUND"
!            STOP

99          continue

            ! Proceed to next element
            ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivtReplace,.FALSE.)
          END DO

          ! Swap tables IVT and ivtReplace in arraylist and release table ivtReplace
          CALL arrlst_swapArrayList(rhadapt%rElementsAtVertex,ivt,ivtReplace)
          CALL arrlst_releaseArrayList(rhadapt%rElementsAtVertex,ivtReplace)

        ELSE

          PRINT *,"is last vertex"
          ! Release table IVT
          CALL arrlst_releaseArrayList(rhadapt%rElementsAtVertex,ivt)

        END IF
        
        PAUSE        
!        CALL fcb_removeVertex(ivt,ivtReplace)
      END IF
    END DO

    PAUSE

    ! Increase the number of recoarsening steps by one
    rhadapt%nCoarseningSteps = rhadapt%nCoarseningSteps+1

    ! Coarsening has been performed.
    rhadapt%iSpec = IOR(rhadapt%ispec,HADAPT_COARSENED)
    
    ! Hence, the markers are no longer valid
    rhadapt%iSpec = IAND(rhadapt%iSpec,NOT(HADAPT_MARKEDCOARSEN))
  END SUBROUTINE redgreen_coarsen

  ! ***************************************************************************

!<function>

  PURE FUNCTION redgreen_getState(rhadapt,iel) RESULT(istate)

!<description>
    ! This function encodes the state of an element iel.
    ! Note: The state of an element is solely determined from the 
    ! age information of its surrounding vertices.
!</description>

!<input>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(IN)        :: rhadapt

    ! Number of element for which state should be computed
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
!</input>

!<result>
    ! State of element
    INTEGER :: istate
!</result>
!</function>

    ! local variables
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D) :: Kvert
    
    ! Get local data
    Kvert = rhadapt%p_IverticesAtElement(1:4,iel)

    ! Are we triangular or quadrilateral element?
    IF (Kvert(4).EQ.0) THEN
      istate=redgreen_getStateTria(rhadapt%p_IvertexAge(Kvert(1:3)))
    ELSE
      istate=redgreen_getStateQuad(rhadapt%p_IvertexAge(Kvert(1:4)))
    END IF
  END FUNCTION redgreen_getState

  ! ***************************************************************************

!<function>

  PURE FUNCTION redgreen_getStateTria(IvertexAge) RESULT(istate)

!<description>
    ! This pure functions encodes the state of the given triangle
    ! into a unique positive integer value. If the state cannot
    ! be determined uniquely, then the position of the youngest
    ! vertex is encoded but with negative sign.
!</description>

!<input>
    ! Age of the three vertices
    INTEGER, DIMENSION(3), INTENT(IN) :: IvertexAge
!</input>

!<result>
    ! State of the triangle
    INTEGER :: istate
!</result>
!</function>

    ! local variables
    INTEGER :: ive,ipos

    ! Reset state
    istate = 0
    
    ! Check if triangle belongs to the initial triangulation. Then nothing
    ! needs to be done. Otherwise, perform additional checks
    IF (SUM(ABS(IvertexAge)).NE.0) THEN
      
      ! Check if any two vertices have the same age and mark that edge.
      ! In addition, determine the largest age and its position
      ipos=1
      DO ive=1,3
        IF (ABS(IvertexAge(ive)) .EQ. ABS(IvertexAge(MOD(ive,3)+1))) THEN
          ! Edge connects two nodes with the same age
          istate=ibset(istate,ive)
        ELSEIF (ABS(IvertexAge(ive)) .GT. ABS(IvertexAge(ipos))) THEN
          ! The age of the "next" node is larger, so than it might be the largest
          ipos=ive
        END IF
      END DO
      
      ! If ISTATE = 0, then there must be one youngest vertex which is returned
      ! but with negative sign. In this case the state of the triangle cannot be
      ! uniquely determined.
      ! Otherwise, either all nodes have the same age (1110 : inner red triangle)
      ! or exactly two nodes have the same age. In the latter case, we must check
      ! if the opposite vertex is the youngest one in the triple of nodes.
      SELECT CASE(istate)
      CASE(STATE_TRIA_ROOT)
        istate=-ibset(istate,ipos)

      CASE(STATE_TRIA_OUTERINNER1)
        IF (ipos .EQ. 3) istate=-ibset(0,ipos)

      CASE(STATE_TRIA_OUTERINNER)
        IF (ipos .EQ. 1) istate=-ibset(0,ipos)

      CASE(STATE_TRIA_OUTERINNER2)
        IF (ipos .EQ. 2) istate=-ibset(0,ipos)

      END SELECT
    END IF

    ! Mark state for triangle
    istate = ibclr(istate,0)
  END FUNCTION redgreen_getStateTria

  ! ***************************************************************************

!<function>

  PURE FUNCTION redgreen_getStateQuad(IvertexAge) RESULT(istate)

!<description>
    ! This pure functions encodes the state of the given quadrilateral
    ! into a unique integer value
!</description>

!<input>
    ! Age of the four vertices
    INTEGER, DIMENSION(4), INTENT(IN) :: IvertexAge
!</input>

!<result>
    ! State of the quadrilateral
    INTEGER :: istate
!</result>
!</function>

    ! local variables
    INTEGER :: ive

    ! Reset state
    istate = 0

    ! Check if quadrilateral belongs to the initial triangulation. Then nothing
    ! needs to be done. Otherwise, perform additional checks.
    IF (SUM(ABS(IvertexAge)).NE.0) THEN
      
      ! Check if any two vertices have the same age and mark that edge.
      ! After this procedure, ISTATE must be different from 0 and a unique
      ! state for the quadrilateral has been determined.
      DO ive=1,4
        IF (ABS(IvertexAge(ive)) .EQ. ABS(IvertexAge(MOD(ive,4)+1))) THEN
          ! Edge connects two nodes with the same age
          istate = ibset(istate,ive)
        END IF
      END DO
    END IF

    ! Mark state for quadrilateral
    istate = ibset(istate,0)
  END FUNCTION redgreen_getStateQuad

  ! ***************************************************************************

!<function>
  
  PURE FUNCTION redgreen_rotateState(istate,irotate) RESULT(inewstate)
    
!<description>
    ! This pure function "rotates" a given state
!</description>

!<input>
    ! Current state
    INTEGER, INTENT(IN) :: istate

    ! Positive rotation
    INTEGER, INTENT(IN) :: irotate
!</input>

!<result>
    ! New state
    INTEGER :: inewstate
!</result>
    !</function>
    
    ! Are we triangular or quadrilateral?
    IF (btest(istate,0)) THEN
      inewstate = ishft(istate,-1)
      inewstate = ishftc(inewstate,irotate,4)
      inewstate = ishft(inewstate,1)
      inewstate = ibset(inewstate,0)
    ELSE
      inewstate = ishft(istate,-1)
      inewstate = ishftc(inewstate,irotate,3)
      inewstate = ishft(inewstate,1)
      inewstate = ibclr(inewstate,0)
    END IF
  END FUNCTION redgreen_rotateState

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bisection_refine(rhadapt)

!<description>
    ! This subroutine performs longest edge bisection as proposed by M. Rivara.
!</description>

!<inputoutput>
    ! adaptive structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
        
!!$      ! Mark elements for refinement
!!$      SELECT CASE(ind)
!!$      CASE(NODE)
!!$         DO iel=1,mg%nel
!!$            IF (search(TT_elem,iel,ipred) == TFOUND) THEN
!!$               ipos  = TT_elem%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$               kvert = TT_elem%kb(ELivert1:ELivert3,ipos)
!!$               
!!$               IF (SUM(dind(kvert)) > 3*tol) THEN
!!$                  d(1)=seglength(TQ_dcorvg%dcorvg(1,kvert((/1,2/))),TQ_dcorvg%dcorvg(2,kvert((/1,2/))))
!!$                  d(2)=seglength(TQ_dcorvg%dcorvg(1,kvert((/2,3/))),TQ_dcorvg%dcorvg(2,kvert((/2,3/))))
!!$                  d(3)=seglength(TQ_dcorvg%dcorvg(1,kvert((/1,3/))),TQ_dcorvg%dcorvg(2,kvert((/1,3/))))
!!$                  
!!$                  kref=.FALSE.; kref(MAXLOC(d))=.TRUE.
!!$                  CALL refine_1to2(kref,iel)
!!$               END IF
!!$            END IF
!!$         END DO
!!$         
!!$      CASE(ELEMENT)
!!$         DO iel=1,mg%nel
!!$            IF (search(TT_elem,iel,ipred) == TFOUND) THEN
!!$               ipos  = TT_elem%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$               kvert = TT_elem%kb(ELivert1:ELivert3,ipos)
!!$               
!!$               IF (dind(iel) > tol) THEN
!!$             !     d(1)=seglength(TQ_dcorvg%dcorvg(1,kvert((/1,2/))),TQ_dcorvg%dcorvg(2,kvert((/1,2/))))
!!$             !     d(2)=seglength(TQ_dcorvg%dcorvg(1,kvert((/2,3/))),TQ_dcorvg%dcorvg(2,kvert((/2,3/))))
!!$             !     d(3)=seglength(TQ_dcorvg%dcorvg(1,kvert((/1,3/))),TQ_dcorvg%dcorvg(2,kvert((/1,3/))))
!!$                  
!!$                  IF (search(TT_node,kvert(1),ipred) /= TFOUND) CALL sys_halt()
!!$                  i1=TT_node%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$                  IF (search(TT_node,kvert(2),ipred) /= TFOUND) CALL sys_halt()
!!$                  i2=TT_node%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$                  IF (search(TT_node,kvert(3),ipred) /= TFOUND) CALL sys_halt()
!!$                  i3=TT_node%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$                  
!!$                  d(1)=ABS(TT_node%db(1,i1)-TT_node%db(1,i2))
!!$                  d(2)=ABS(TT_node%db(1,i2)-TT_node%db(1,i3))
!!$                  d(3)=ABS(TT_node%db(1,i3)-TT_node%db(1,i1))
!!$
!!$                  kref=.FALSE.; kref(MAXLOC(d))=.TRUE.
!!$                  CALL refine_1to2(kref,iel)
!!$               END IF
!!$            END IF
!!$         END DO
!!$      END SELECT
!!$
!!$      GOTO 1000
!!$      
!!$      DO
!!$         bmod=.FALSE.
!!$         mg%nel=TT_elem%na
!!$         ALLOCATE(r(mg%nel)); r=.FALSE.
!!$         DO iel=1,mg%nel
!!$            IF (search(TT_elem,iel,ipred) == TFOUND) THEN
!!$               ipos    = TT_elem%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$               kvert   = TT_elem%kb(ELivert1:ELivert3,ipos)
!!$               kadj    = TT_elem%kb(ELiadj1:ELiadj3,ipos)
!!$               kmidadj = TT_elem%kb(ELimidadj1:ELimidadj3,ipos)
!!$
!!$               IF (ANY(kadj/=kmidadj)) THEN
!!$                  r(iel)=.TRUE.; bmod=.TRUE.
!!$               END IF
!!$            END IF
!!$         END DO
!!$
!!$         IF (.NOT.bmod) EXIT
!!$
!!$         DO iel=1,TT_elem%na
!!$            IF (.NOT.r(iel)) CYCLE
!!$
!!$            IF (search(TT_elem,iel,ipred) == TFOUND) THEN
!!$               ipos    = TT_elem%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$               kvert   = TT_elem%kb(ELivert1:ELivert3,ipos)
!!$               kadj    = TT_elem%kb(ELiadj1:ELiadj3,ipos)
!!$               kmidadj = TT_elem%kb(ELimidadj1:ELimidadj3,ipos)
!!$
!!$    !           d(1)=seglength(TQ_dcorvg%dcorvg(1,kvert((/1,2/))),TQ_dcorvg%dcorvg(2,kvert((/1,2/))))
!!$    !           d(2)=seglength(TQ_dcorvg%dcorvg(1,kvert((/2,3/))),TQ_dcorvg%dcorvg(2,kvert((/2,3/))))
!!$    !           d(3)=seglength(TQ_dcorvg%dcorvg(1,kvert((/1,3/))),TQ_dcorvg%dcorvg(2,kvert((/1,3/))))
!!$               
!!$               IF (search(TT_node,kvert(1),ipred) /= TFOUND) CALL sys_halt()
!!$               i1=TT_node%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$               IF (search(TT_node,kvert(2),ipred) /= TFOUND) CALL sys_halt()
!!$               i2=TT_node%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$               IF (search(TT_node,kvert(3),ipred) /= TFOUND) CALL sys_halt()
!!$               i3=TT_node%kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
!!$
!!$               d(1)=ABS(TT_node%db(1,i1)-TT_node%db(1,i2))
!!$               d(2)=ABS(TT_node%db(1,i2)-TT_node%db(1,i3))
!!$               d(3)=ABS(TT_node%db(1,i3)-TT_node%db(1,i1))
!!$
!!$               ipos=ABS(SUM((kadj/=kmidadj)*(/1,2,3/)))
!!$
!!$               IF (ipos > 3) THEN
!!$                  kref=.TRUE.
!!$                  CALL refine_1to4(kref,iel)
!!$               ELSE
!!$                  IF (ABS(d(ipos)-MAXVAL(d))  <= 1e-3) THEN
!!$                     kref=kadj/=kmidadj
!!$                     CALL refine_1to2(kref,iel)
!!$                  ELSE
!!$                     kref=.FALSE.; kref(MAXLOC(d))=.TRUE.; kref=kref.OR.(kadj/=kmidadj)
!!$                     CALL refine_1to3(kref,iel)
!!$                  END IF
!!$               END IF
!!$            END IF
!!$         END DO
!!$         DEALLOCATE(r)
  END SUBROUTINE bisection_refine
  
END MODULE hadaptivity

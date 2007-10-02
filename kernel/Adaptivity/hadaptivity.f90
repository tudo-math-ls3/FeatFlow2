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
!#  6.) update_ElementNeighbors2D = update_ElemNeighb2D_1to2 /
!#                                  update_ElemNeighb2D_2to2
!#      -> update the list of neighboring elements  in 2D
!#
!#  7.) update_AllElementNeighbors2D
!#      -> update the lists of neighboring elements of ALL adjacent elements
!#
!#  8.) refine_Tria2Tria
!#      -> refine a triangle by subdivision into two triangles
!#
!#  9.) refine_Tria3Tria
!#      -> refine a triangle by subdivision into three triangles
!#
!# 10.) refine_Tria4Tria
!#      -> refine a triangle by subdivision into four triangles
!#
!# 11.) refine_Quad2Quad
!#      -> refine a quadrilateral by subdivision into two quadrilaterals
!#
!# 12.) refine_Quad3Tria
!#      -> refine a quadrilateral by subdivision into three triangles
!#
!# 13.) refine_Quad4Tria
!#      -> refine a quadrilateral by subdivision into four triangles
!#
!# 14.) refine_Quad4Quad
!#      -> refine a quadrilateral by subdivision into four quadrilaterals
!#
!# 15.) convert_Tria2Tria
!#      -> convert two neighboring triangles into four similar triangle
!#
!# 16.) convert_Quad2Quad
!#      -> convert two neighboring quadrilaterals into four similar quadrilaterals
!#
!# 17.) convert_Quad3Tria
!#      -> convert three neighboring triangles into four similar quadrilaterals
!#
!# 18.) convert_Quad4Tria
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
!# 22.) coarsen_4Quad1Quad
!#      -> coarsen four red quadrilaterals by combining them to the macro element
!#
!# 23.) coarsen_4Quad2Quad
!#      -> coarsen four red quadrilaterals by combining them to two green elements
!#
!# 24.) coarsen_4Quad3Tria
!#      -> coarsen four red quadrilaterals by combining them to three green elements
!#
!# 25.) coarsen_4Quad4Tria
!#      -> coarsen four red quadrilaterals by combining them to four green elements
!#
!# 26.) coarsen_2Quad1Quad
!#      -> coarsen two green quadrilaterals by combining them to the macro element
!#
!# 27.) coarsen_2Quad3Tria
!#      -> coarsen two green quadrilaterals by combining them to three green triangles
!#
!# 28.) coarsen_3Tria1Quad
!#      -> coarsen three green triangles by combining them to the macro element
!#
!# 29.) coarsen_4Tria1Quad
!#      -> coarsen four green triangles by combining them to the macro element
!#
!# 30.) coarsen_4Tria3Tria
!#      -> coarsen four green triangles by combining them to three green triangles
!#
!# 31.) mark_refinement2D
!#      -> mark elements for refinement in 2D
!#
!# 32.) redgreen_mark_coarsening2D
!#      -> mark elements for coarsening in 2D
!#
!# 33.) redgreen_mark_refinement2D
!#      -> mark elements for refinement due to Red-Green strategy in 2D
!#
!# 34.) redgreen_refine
!#      -> perform red-green refinement for marked elements
!#
!# 35.) redgreen_coarsen
!#      -> perform red-green coarsening for marked elements
!#
!# 36.) redgreen_getState
!#      -> return the state of an element
!#
!# 37.) redgreen_getStateTria
!#      -> return the state of a triangle
!#
!# 38.) redgreen_getStateQuad
!#      -> return the state of a quadrilateral
!#
!# 39.) redgreen_rotateState
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

  ! Operation identifier for coarsening: 2-tria : 1-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_2TRIA1TRIA   = 17

  ! Operation identifier for coarsening: 4-tria : 1-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA1TRIA   = 18

  ! Operation identifier for coarsening: 4-tria : 2-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA2TRIA1  = 19
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA2TRIA2  = 20
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA2TRIA3  = 21

  ! Operation identifier for coarsening: 4-quad : 1-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4QUAD1QUAD   = 22

  ! Operation identifier for coarsening: 4-quad : 2-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4QUAD2QUAD   = 23

  ! Operation identifier for coarsening: 4-quad : 3-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4QUAD3TRIA   = 24

  ! Operation identifier for coarsening: 4-quad : 4-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4QUAD4TRIA   = 25

  ! Operation identifier for coarsening: 2-quad : 1-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_2QUAD1QUAD   = 26

  ! Operation identifier for coarsening: 2-quad : 3-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_2QUAD3TRIA   = 27

  ! Operation identifier for coarsening: 3-tria : 1-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_3TRIA1QUAD   = 28

  ! Operation identifier for coarsening: 4-tria : 1-quad
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA1QUAD   = 29

  ! Operation identifier for coarsening: 4-tria : 3-tria
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA3TRIA2 = 30
  INTEGER, PARAMETER, PUBLIC :: HADAPT_OPR_CRS_4TRIA3TRIA3 = 31

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
  ! recoarsening into the macro quadrilateral together with its neighbors
  INTEGER, PARAMETER :: MARK_CRS_3TRIA1QUAD         = -6

  ! Mark (left) green triangle of a 1-tria : 2-tria refinement for 
  ! recoarsening into the macro triangle together with its right neighbor
  INTEGER, PARAMETER :: MARK_CRS_2TRIA1TRIA         = -7

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! recoarsening into the macro quadrilateral together with its three neighbors
  INTEGER, PARAMETER :: MARK_CRS_4TRIA1QUAD         = -8

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! conversion into three triangles keeping the left neighboring triangle
  INTEGER, PARAMETER :: MARK_CRS_4TRIA3TRIA_LEFT    = -9

  ! Mark "most inner" green triangle of a 1-quad : 4-tria refinement for
  ! conversion into three triangles keeping the right neighboring triangle
  INTEGER, PARAMETER :: MARK_CRS_4TRIA3TRIA_RIGHT   = -10

  ! Mark green quadrilateral of a 1-quad : 2-quad refinement for
  ! recoarsening into the macro quadrilateral together with its neighbor
  INTEGER, PARAMETER :: MARK_CRS_2QUAD1QUAD         = -11
  
  ! Mark green quadrilateral of a 1-quad : 2-quad refinement for
  ! conversion into three triangles keeping the second vertex
  INTEGER, PARAMETER :: MARK_CRS_2QUAD3TRIA         = -12

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for 
  ! recoarsening into the macro quadrilateral together with its neighbors
  INTEGER, PARAMETER :: MARK_CRS_4QUAD1QUAD         = -13

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for 
  ! recoarsening into two quadrilaterals
  INTEGER, PARAMETER :: MARK_CRS_4QUAD2QUAD         = -14

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for 
  ! recoarsening into three triangles
  INTEGER, PARAMETER :: MARK_CRS_4QUAD3TRIA         = -15

  ! Mark red quadrilateral of a 1-quad : 4-quad refinement for 
  ! recoarsening into four triangles
  INTEGER, PARAMETER :: MARK_CRS_4QUAD4TRIA         = -16
  
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
  INTEGER, PARAMETER :: STATE_QUAD_HALF1            = 5
  INTEGER, PARAMETER :: STATE_QUAD_HALF2            = 21
  

!</constantblock>


!<constantblock description="Constants for grid adaptation">
  
  ! Array position of the boundary
  INTEGER, PARAMETER :: BdrValue = 1
  
  ! Array position of the previous boundary vertex
  INTEGER, PARAMETER :: BdrPrev  = 1

  ! Array position of the next boundary vertex
  INTEGER, PARAMETER :: BdrNext  = 2

!</constantblock>


!<constantblock description="Duplication flags. Specifies which information is
!                            shared between adaptivity structures">

  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IMARKER                 = 2** 0
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IVERTEXAGE              = 2** 1
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_INODALPROPERTY          = 2** 2
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IVERTICESATELEMENT      = 2** 3
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_INEIGHATELEMENT         = 2** 4
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_IMIDNEIGHATELEMENT      = 2** 5
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_RVERTEXCOORDINATES      = 2** 6
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_RBOUNDARY               = 2** 7
  INTEGER(I32), PARAMETER :: HADAPT_SHARE_RELEMENTSATVERTEX       = 2** 8

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

    ! Duplication flag. Bitfield that indicates which information is
    ! shared with another adaptivity structure.
    ! When a bit is set to 1, the corresponding array is
    ! maintained by another adaptivity structure and must
    ! not be deleted by hadapt_releaseAdaptation. 
    ! When the bit is 0, the array is a real copy of another array 
    ! and must be deleted in hadapt_releaseAdaptation.
    INTEGER(I32) :: iduplicationFlag                 = 0

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
    TYPE(t_arraylist) :: rElementsAtVertex
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
    MODULE PROCEDURE update_ElemNeighb2D_1to2
    MODULE PROCEDURE update_ElemNeighb2D_2to2
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
    rhadapt%iSpec=HADAPT_HAS_PARAMETERS
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
      DO ibct=1,SIZE(rhadapt%rBoundary)
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
      rhadaptBackup%iRefinementStrategy  = rhadapt%iRefinementStrategy
      rhadaptBackup%iCoarseningStrategy  = rhadapt%iCoarseningStrategy
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
    rhadapt%iRefinementStrategy  = rhadaptBackup%iRefinementStrategy
    rhadapt%iCoarseningStrategy  = rhadaptBackup%iCoarseningStrategy
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
!</inputoutput
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
!</inputoutput
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
    CALL storage_getsize(rhadapt%h_IverticesAtElement,Isize)
    rhadapt%NELMAX=Isize(2)
    CALL storage_getsize(rhadapt%h_IneighboursAtElement,Isize)
    rhadapt%NELMAX=MIN(rhadapt%NELMAX,Isize(2))

    rhadapt%InelOfType0=rhadapt%InelOfType
    rhadapt%NVT0       =rhadapt%NVT
    rhadapt%NEL0       =rhadapt%NEL
    rhadapt%NVBD0      =rhadapt%NVBD
    rhadapt%increaseNVT=0
    
    ! What kind of grid refinement should be performed
    SELECT CASE(rhadapt%irefinementStrategy)
    CASE (HADAPT_NOREFINEMENT)   ! No grid refinement

      
    CASE (HADAPT_REDGREEN)       ! Red-green grid refinement

      ! Mark elements for refinement based on indicator function
      CALL mark_refinement2D(rhadapt,rindicator)

      ! Mark additional elements to restore conformity
      CALL redgreen_mark_refinement2D(rhadapt,rcollection,fcb_hadaptCallback)

      ! Mark element for recoarsening based on indicator function
      CALL redgreen_mark_coarsening2D(rhadapt,rindicator)

      ! Compute new dimensions
      nvt=rhadapt%NVT+rhadapt%increaseNVT
      nel=NumberOfElements(rhadapt)

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
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/nvt/); Ielements=(/0/)
        CALL fcb_hadaptCallback(rcollection,&
            HADAPT_OPR_ADJUSTVERTEXDIM,Ivertices,Ielements)
      END IF

      ! Perform refinement
      CALL redgreen_refine(rhadapt,rcollection,fcb_hadaptCallback)

      ! Perform coarsening
      CALL redgreen_coarsen(rhadapt,rcollection,fcb_hadaptCallback)

      ! Adjust nodal property array
      CALL storage_realloc('hadapt_performAdaptation',rhadapt%NVT,&
          rhadapt%h_InodalProperty,ST_NEWBLOCK_NOINIT,.TRUE.)
      

!!$    CASE (HADAPT_LONGESTEDGE)   ! Bisection of longest edge
!!$
!!$      ! Mark elements for refinement based on indicator function
!!$      CALL mark_refinement2D(rhadapt,rindicator)
!!$      
!!$      ! Perform refinement
!!$      CALL bisection_refine(rhadapt)

      
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
      
      CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
      
      ! Initialize number of elements by current number
      nel=rhadapt%NEL0
      
      ! Loop over all elements and check marker
      DO iel=1,rhadapt%NEL0
        
        SELECT CASE(p_Imarker(iel))
        CASE(MARK_REF_TRIA2TRIA_1,MARK_REF_TRIA2TRIA_2,MARK_REF_TRIA2TRIA_3,&
             MARK_REF_QUAD2QUAD_13,MARK_REF_QUAD2QUAD_24,MARK_CRS_2QUAD3TRIA)
          ! Interestingly enought, the coarsening of 2 quadrilaterals into three
          ! triangles which reduces the number of vertices by one also increases
          ! the number of elements by one which has to be taken into account here.
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
    CALL output_line('Strategy for grid refinement:                '//&
        TRIM(sys_siL(rhadapt%irefinementStrategy,3)))
    CALL output_line('Strategy for grid coarsening:                '//&
        TRIM(sys_siL(rhadapt%icoarseningStrategy,3)))
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
        nve=get_NVE(rhadapt,iel)
        
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
                TRIM(sys_siL(redgreen_getState(rhadapt,iel),4))//''','''//&
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
                TRIM(sys_siL(redgreen_getState(rhadapt,iel),4))//''','''//&
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
                TRIM(sys_siL(redgreen_getState(rhadapt,iel),4))//''','''//&
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
                TRIM(sys_siL(redgreen_getState(rhadapt,iel),4))//''','''//&
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
        nve=get_NVE(rhadapt,iel)
        
        ! For all new elements there si no marker available
        WRITE(iunit,FMT='(A)') '<polygon id="el'//TRIM(sys_siL(iel,9))//&
            '" fill="white" stroke="black" stroke-width="1"'

        ! Write data which is common to all elements
        SELECT CASE(nve)
        CASE(TRIA_NVETRI2D)
          WRITE(iunit,FMT='(A)') ' onmousemove="ShowElementInfo(evt,'''//&
              TRIM(sys_siL(iel,10))//''','''//&
              TRIM(sys_siL(redgreen_getState(rhadapt,iel),4))//''','''//&
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
              TRIM(sys_siL(redgreen_getState(rhadapt,iel),4))//''','''//&
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
      nve=get_NVE(rhadapt,iel)
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

    DO ivt=1,SIZE(rhadapt%p_IvertexAge)
      WRITE(UNIT=iunit,FMT=40) rhadapt%p_IvertexAge(ivt)
    END DO
    DO ivt=SIZE(rhadapt%p_IvertexAge)+1,rhadapt%NVT
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
      nve=get_NVE(rhadapt,iel)
      
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
          mve=get_NVE(rhadapt,jel)
          
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
          mve=get_NVE(rhadapt,jel)
          
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
          mve=get_NVE(rhadapt,jelmid)
          
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
      nve=get_NVE(rhadapt,iel)
      
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
        mve=get_NVE(rhadapt,jel)

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
        nve=get_NVE(rhadapt,iel)

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
        nve=get_NVE(rhadapt,iel)

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
    IF (rhadapt%ndim .EQ. NDIM2D) THEN
      
      ! Do we have quadrilaterals in the triangulation?
      IF (rhadapt%InelOfType(TRIA_NVEQUAD2D).EQ.0) THEN

        ! There are no quadrilaterals in the current triangulation.
        ! Hence, return TRIA_NVETRI2D by default.
        nve=TRIA_NVETRI2D
        
      ELSE
        
        ! There are quadrilaterals and possible also triangles in
        ! the current triangulatin. If the last entry of the
        ! vertices-at-element list is nonzero then TRIA_NVEQUAD2D vertices
        ! are present in the current element. Otherwise return TRIA_NVETRI2D.
        IF (rhadapt%p_IverticesAtElement(TRIA_NVEQUAD2D,iel).EQ.0) THEN
          nve=TRIA_NVETRI2D
        ELSE
          nve=TRIA_NVEQUAD2D
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
    REAL(DP), DIMENSION(1)      :: Ddata
    REAL(DP)                    :: x1,y1,x2,y2,dvbdp1,dvbdp2
    INTEGER(PREC_VERTEXIDX),  DIMENSION(3) :: Ivertices
    INTEGER(PREC_ELEMENTIDX), DIMENSION(1) :: Ielements
    INTEGER, DIMENSION(2)       :: Idata
    INTEGER(PREC_QTREEIDX)      :: inode
    INTEGER(PREC_TREEIDX)       :: ipred,ipos
    INTEGER                     :: ibct
    
    ! Get coordinates of edge vertices
    x1=qtree_getX(rhadapt%rVertexCoordinates2D,i1)
    y1=qtree_getY(rhadapt%rVertexCoordinates2D,i1)
    x2=qtree_getX(rhadapt%rVertexCoordinates2D,i2)
    y2=qtree_getY(rhadapt%rVertexCoordinates2D,i2)

    ! Compute coordinates of new vertex 
    Dcoord=0.5_DP*(/x1+x2,y1+y2/)

    ! Search for vertex coordinates in quadtree: 
    ! If the vertex already exists, e.g., it was added when the adjacent element
    ! was refined, then nothing needs to be done for this vertex
    IF (qtree_searchInQuadtree(rhadapt%rVertexCoordinates2D,&
        Dcoord,inode,ipos,i12) .EQ. QTREE_FOUND) RETURN
    
    ! Otherwise, update number of vertices
    rhadapt%NVT=rhadapt%NVT+1
    i12        =rhadapt%NVT

    ! Set age of vertex
    rhadapt%p_IvertexAge(i12)=&
        1+MAX(ABS(rhadapt%p_IvertexAge(i1)),ABS(rhadapt%p_IvertexAge(i2)))
    
    ! Set nodal property
    IF (e1 .EQ. 0) THEN
      rhadapt%p_InodalProperty(i12)=rhadapt%p_InodalProperty(i1)
    ELSE
      rhadapt%p_InodalProperty(i12)=0
    END IF
    
    ! Add new entry to vertex coordinates
    CALL qtree_insertIntoQuadtree(rhadapt%rVertexCoordinates2D,i12,Dcoord,inode)
    
    ! Are we at the boundary?
    IF (e1 .EQ. 0) THEN
      ! Increment number of boundary nodes
      rhadapt%NVBD=rhadapt%NVBD+1
      
      ! Get number of boundary component
      ibct=rhadapt%p_InodalProperty(i1)
      
      ! Get parameter values of the boundary nodes
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i1,ipred) .EQ. BTREE_NOT_FOUND) THEN
        CALL output_line('Unable to find first vertex in boudary data structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atEdgeMidpoint2D')
        CALL sys_halt()
      END IF
      ipos  =rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      dvbdp1=rhadapt%rBoundary(ibct)%p_DData(BdrValue,ipos)
      
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i2,ipred) .EQ. BTREE_NOT_FOUND) THEN
        CALL output_line('Unable to find second vertex in boudary data structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'add_vertex_atEdgeMidpoint2D')
        CALL sys_halt()
      END IF
      ipos  =rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      dvbdp2=rhadapt%rBoundary(ibct)%p_DData(BdrValue,ipos)
      
      ! If I2 is last(=first) node on boundary component IBCT round DVBDP2 to next integer
      IF (dvbdp2 .LE. dvbdp1) dvbdp2=CEILING(dvbdp1)
      
      ! Add new entry to boundary structure
      Idata=(/i1,i2/); Ddata=(/0.5_DP*(dvbdp1+dvbdp2)/)
      CALL btree_insertIntoTree(rhadapt%rBoundary(ibct),i12,Idata=Idata,Ddata=Ddata)
    END IF
      
    ! Optionally, invoke callback function
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i12,i1,i2/); Ielements=(/0/)
      CALL fcb_hadaptCallback(rcollection,&
          HADAPT_OPR_INSERTVERTEXEDGE,Ivertices,Ielements)
    END IF
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(1) :: Ielements
    INTEGER(PREC_VERTEXIDX), DIMENSION(5)  :: Ivertices
    REAL(DP) :: x1,y1,x2,y2,x3,y3,x4,y4,x21,y21,x31,y31,x24,y24,alpha
    INTEGER(PREC_QTREEIDX) :: inode,ipos
    
    ! Compute coordinates of new vertex
    x1=qtree_getX(rhadapt%rVertexCoordinates2D,i1)
    y1=qtree_getY(rhadapt%rVertexCoordinates2D,i1)
    x2=qtree_getX(rhadapt%rVertexCoordinates2D,i2)
    y2=qtree_getY(rhadapt%rVertexCoordinates2D,i2)
    x3=qtree_getX(rhadapt%rVertexCoordinates2D,i3)
    y3=qtree_getY(rhadapt%rVertexCoordinates2D,i3)
    x4=qtree_getX(rhadapt%rVertexCoordinates2D,i4)
    y4=qtree_getY(rhadapt%rVertexCoordinates2D,i4)

    x21=x2-x1; x31=x3-x1; x24=x2-x4
    y21=y2-y1; y31=y3-y1; y24=y2-y4
    alpha=(x21*y24-y21*x24)/(x31*y24-y31*x24)
    
    Dcoord=(/x1+alpha*x31,y1+alpha*y31/)

    ! Search for vertex coordinates in quadtree
    IF (qtree_searchInQuadtree(rhadapt%rVertexCoordinates2D,&
        Dcoord,inode,ipos,i5) .EQ. QTREE_NOT_FOUND) THEN
      
      ! Update number of vertices
      rhadapt%NVT=rhadapt%NVT+1
      i5         =rhadapt%NVT
      
      ! Set age of vertex
      rhadapt%p_IvertexAge(I5)=&
          1+MAX(ABS(rhadapt%p_IvertexAge(i1)),&
                ABS(rhadapt%p_IvertexAge(i2)),&
                ABS(rhadapt%p_IvertexAge(i3)),&
                ABS(rhadapt%p_IvertexAge(i4)))

      ! Set nodal property
      rhadapt%p_InodalProperty(i5)=0
      
      ! Add new entry to DCORVG
      CALL qtree_insertIntoQuadtree(rhadapt%rVertexCoordinates2D,i5,Dcoord,inode)
    END IF
    
    ! Optionally, invoke callback function
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i5,i1,i2,i3,i4/); Ielements=(/0/)
      CALL fcb_hadaptCallback(rcollection,&
          HADAPT_OPR_INSERTVERTEXCENTR,Ivertices,Ielements)
    END IF
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

    ! Remove vertex from coordinates and get number of replacement vertex
    IF (qtree_deleteFromQuadtree(rhadapt%rVertexCoordinates2D,&
        ivt,ivtReplace) .EQ. QTREE_NOT_FOUND) THEN
      CALL output_line('Unable to delete vertex coordinates!',&
          OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
      CALL sys_halt()
    END IF
    
    ! Decrease number of vertices by one
    rhadapt%NVT=rhadapt%NVT-1
    
    ! If IVT is a boundary node remove it from the boundary and
    ! connect its boundary neighbors with each other
    ibct=rhadapt%p_InodalProperty(ivt)
    
    ! Are we at the boundary?
    IF (ibct .NE. 0) THEN
      
      ! Find position of vertex IVT in boundary array
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),ivt,ipred) .EQ.&
          BTREE_NOT_FOUND) THEN
        CALL output_line('Unable to find vertex in boundary data structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        CALL sys_halt()
      END IF
      ipos=rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      
      ! Get the two boundary neighbors: I1 <- IVT -> I2
      i1=rhadapt%rBoundary(ibct)%p_IData(BdrPrev,ipos)
      i2=rhadapt%rBoundary(ibct)%p_IData(BdrNext,ipos)
      
      ! Connect the boundary neighbors with each other: I1 <=> I2
      ! First, set I2 as next neighboring of I1
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i1,ipred) .EQ.&
          BTREE_NOT_FOUND) THEN
        CALL output_line('Unable to find left neighboring vertex in boundary data structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        CALL sys_halt()
      END IF
      ipos=rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      rhadapt%rBoundary(ibct)%p_IData(BdrNext,ipos)=i2
      
      ! Second, set I1 as previous neighbor of I2
      IF (btree_searchInTree(rhadapt%rBoundary(ibct),i2,ipred) .EQ.&
          BTREE_NOT_FOUND) THEN
        CALL output_line('Unable to find right neighboring vertex in boundary data structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        CALL sys_halt()
      END IF
      ipos=rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
      rhadapt%rBoundary(ibct)%p_IData(BdrPrev,ipos)=i1
      
      ! And finally, delete IVT from the boundary
      IF (btree_deleteFromTree(rhadapt%rBoundary(ibct),ivt) .EQ.&
          BTREE_NOT_FOUND) THEN
        CALL output_line('Unable to delete vertex from boundary data structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
        CALL sys_halt()
      END IF
    END IF
    
    ! If IVT is not the last node then copy the data for vertex IVTREPLACE
    ! to IVT and prepare IVTREPLACE for elimination
    IF (ivt < ivtReplace) THEN
      
      ! If IVTREPLACE is a boundary node then remove IVTREPLACE from the boundary
      ! vector, insert IVT into the boundary vector instead, and connect the
      ! boundary neighbors of IVTREPLACE with IVT
      ibct=rhadapt%p_InodalProperty(ivtReplace)

      ! Are we at the boundary?
      IF (ibct .NE. 0) THEN

        IF (btree_searchInTree(rhadapt%rBoundary(ibct),ivtReplace,ipred) .EQ.&
            BTREE_NOT_FOUND) THEN
          CALL output_line('Unable to find replacement vertex in boundary data structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          CALL sys_halt()
        END IF
        ipos=rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        
        ! Insert IVT into the boundary vector
        CALL btree_insertIntoTree(rhadapt%rBoundary(ibct),ivt,&
            Idata=rhadapt%rBoundary(ibct)%p_IData(:,ipos),&
            Ddata=rhadapt%rBoundary(ibct)%p_DData(:,ipos))
        
        ! Get the two boundary neighbors: I1 <- IVTREPLACE -> I2
        i1=rhadapt%rBoundary(ibct)%p_IData(BdrPrev,ipos)
        i2=rhadapt%rBoundary(ibct)%p_IData(BdrNext,ipos)
        
        ! Connect the boundary neighbors with IVT: I1 <- IVT -> I2
        ! First, set IVT as next neighbor of I1
        IF (btree_searchInTree(rhadapt%rBoundary(ibct),i1,ipred) .EQ.&
            BTREE_NOT_FOUND) THEN
          CALL output_line('Unable to find left neighboring vertex in boundary data structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          CALL sys_halt()
        END IF
        ipos=rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        rhadapt%rBoundary(ibct)%p_IData(BdrNext,ipos)=ivt
        
        ! Second, set IVT as previous neighbor of I2
        IF (btree_searchInTree(rhadapt%rBoundary(ibct),i2,ipred) .EQ.&
            BTREE_NOT_FOUND) THEN
          CALL output_line('Unable to find right neighboring vertex in boundary data structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          CALL sys_halt()
        END IF
        ipos=rhadapt%rBoundary(ibct)%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
        rhadapt%rBoundary(ibct)%p_IData(BdrPrev,ipos)=ivt
        
        ! Finally, delete IVTREPLACE from the boundary
        IF (btree_deleteFromTree(rhadapt%rBoundary(ibct),ivtReplace) .EQ.&
            BTREE_NOT_FOUND) THEN
          CALL output_line('Unable to delete vertex from the boundary data structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'remove_vertex2D')
          CALL sys_halt()
        END IF
      END IF
      
      ! Copy data from node IVTREPLACE to node IVT
      rhadapt%p_InodalProperty(ivt)=rhadapt%p_InodalProperty(ivtReplace)
      rhadapt%p_IvertexAge(ivt)    =rhadapt%p_IvertexAge(ivtReplace)

      ! Clear data for node IVTREPLACE
      rhadapt%p_InodalProperty(ivtReplace)=0
      rhadapt%p_IvertexAge(ivtReplace)    =0
      
    ELSE

      ! IVT is the last vertex of the adaptivity structure
      ivtReplace=0

      ! Clear data for node IVT
      rhadapt%p_InodalProperty(ivt)=0
      rhadapt%p_IvertexAge(ivt)    =0
      
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
    rhadapt%p_IverticesAtElement(:,ipos)     =(/i1,i2,i3,0/)
    rhadapt%p_IneighboursAtElement(:,ipos)   =(/e1,e2,e3,0/)
    rhadapt%p_ImidneighboursAtElement(:,ipos)=(/e4,e5,e6,0/)    
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
    rhadapt%p_IverticesAtElement(:,ipos)     =(/i1,i2,i3,i4/)
    rhadapt%p_IneighboursAtElement(:,ipos)   =(/e1,e2,e3,e4/)
    rhadapt%p_ImidneighboursAtElement(:,ipos)=(/e5,e6,e7,e8/)
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
    rhadapt%NEL=rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)+1

    rhadapt%p_IverticesAtElement(:,rhadapt%NEL)     =(/i1,i2,i3,0/)
    rhadapt%p_IneighboursAtElement(:,rhadapt%NEL)   =(/e1,e2,e3,0/)
    rhadapt%p_ImidneighboursAtElement(:,rhadapt%NEL)=(/e4,e5,e6,0/)
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
    rhadapt%NEL=rhadapt%NEL+1
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)+1

    rhadapt%p_IverticesAtElement(:,rhadapt%NEL)     =(/i1,i2,i3,i4/)
    rhadapt%p_IneighboursAtElement(:,rhadapt%NEL)   =(/e1,e2,e3,e4/)
    rhadapt%p_ImidneighboursAtElement(:,rhadapt%NEL)=(/e5,e6,e7,e8/)
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
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_VERTEXIDX)    :: ivt
    INTEGER(PREC_ELEMENTIDX)   :: jel,jelmid
    INTEGER :: ive,jve
    LOGICAL :: bfound

    ! Which kind of element are we?
    SELECT CASE(get_NVE(rhadapt,iel))
    CASE(TRIA_NVETRI2D)
      rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)-1

    CASE(TRIA_NVEQUAD2D)
      rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)-1

    CASE DEFAULT
      CALL output_line('Invalid element type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'remove_element2D')
      CALL sys_halt()
    END SELECT

    ! Replace element by the last element and delete last element
    ielReplace=rhadapt%NEL

    ! Check if we are the last element
    IF (iel .NE. ielReplace) THEN
      
      ! Element is not the last one. Then the element that should be removed must
      ! have a smaller element number. If this is not the case, something is wrong.
      IF (iel > ielReplace) THEN
        CALL output_line('Number of replacement element must not be smaller than that of&
            & the removed elements!',OU_CLASS_ERROR,OU_MODE_STD,'remove_element2D')
        CALL sys_halt()
      END IF

      ! The element which formally was labeled ielReplace is now labeled IEL.
      ! This modification must be updated in the list of adjacent element
      ! neighbors of all surrounding elements. Moreover, the modified element
      ! number must be updated in the "elements-meeting-at-vertex" lists of the
      ! corner nodes of element IEL. Both operations are performed below.
      update: DO ive=1,get_NVE(rhadapt,ielReplace)
        
        ! Get vertex number of corner node
        ivt=rhadapt%p_IverticesAtElement(ive,ielReplace)

        ! Start with first element in "elements-meeting-at-vertex" list
        ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivt,.TRUE.)
        elements: DO WHILE(ipos .GT. ARRLST_NULL)
          
          ! Check if element number corresponds to the replaced element
          IF (rhadapt%rElementsAtVertex%p_IData(ipos) .EQ. ielReplace) THEN
            rhadapt%rElementsAtVertex%p_IData(ipos)=iel
            EXIT elements
          END IF
          
          ! Proceed to next element in list
          ipos=arrlst_getNextInArraylist(rhadapt%rElementsAtVertex,ivt,.FALSE.)
        END DO elements

                
        ! Get element number of element JEL and JELMID which 
        ! are (mid-)adjacent to element ielReplace
        jel    = rhadapt%p_IneighboursAtElement(ive,ielReplace)
        jelmid = rhadapt%p_ImidneighboursAtElement(ive,ielReplace)

        
        ! Do we have different neighbours along the first part and the
        ! second part of the common edge?
        IF (jel .NE. jelmid) THEN
          ! There are two elements sharing the common edge with ielReplace. 
          ! We have to find the position of element ielReplace in the lists
          ! of (mid-)adjacent elements for both element JEL and JELMID
          ! seperately. If we are at the boundary of ir element IEL to be 
          ! removed is adjacent to element ielReplace, the skip this edge.
          bfound=.FALSE.
          
          ! Find position of replacement element in adjacency list of 
          ! element JEL and update the entry to new element number IEL
          IF (jel .EQ. 0 .OR. jel .EQ. iel) THEN
            bfound=.TRUE.
          ELSE
            adjacent1: DO jve=1,get_NVE(rhadapt,jel)
              IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. ielReplace) THEN
                rhadapt%p_IneighboursAtElement(jve,jel)=iel
                rhadapt%p_ImidneighboursAtElement(jve,jel)=iel
                bfound=.TRUE.
                EXIT adjacent1
              END IF
            END DO adjacent1
          END IF
          
          ! Find position of replacement element in adjacentcy list of
          ! element JELMID and update the entry to new element number IEL
          IF (jelmid .EQ. 0 .OR. jelmid .EQ. iel) THEN
            bfound=bfound.AND..TRUE.
          ELSE
            adjacent2: DO jve=1,get_NVE(rhadapt,jelmid)
              IF (rhadapt%p_IneighboursAtElement(jve,jelmid) .EQ. ielReplace) THEN
                rhadapt%p_IneighboursAtElement(jve,jelmid)=iel
                rhadapt%p_ImidneighboursAtElement(jve,jelmid)=iel
                bfound=bfound.AND..TRUE.
                EXIT adjacent2
              END IF
            END DO adjacent2
          END IF

        ELSE
          ! There is only one element sharing the common edge with ielReplace.
          ! If we are at he boundary or if element IEL to be removed is
          ! adjacent to element ielReplace, then skip this edge.
          IF (jel .EQ. 0 .OR. jel .EQ. iel) CYCLE update

          ! We have to consider two possible situations. The neighbouring element
          ! JEL can be completely aligned with element ielReplace so that the
          ! element number must be updated in the list of adjacent and mid-
          ! adjacent elements. On the other hand, element JEL can share one
          ! half of the common edge with element ielReplace and the other half-
          ! edge with another element. In this case, the number ielReplace can
          ! only be found in either the list of adjacent or mid-adjacent elements.
          bfound=.FALSE.
          adjacent3: DO jve=1,get_NVE(rhadapt,jel)
            IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. ielReplace) THEN
              rhadapt%p_IneighboursAtElement(jve,jel)=iel
              bfound=bfound.OR..TRUE.
            END IF

            IF (rhadapt%p_ImidneighboursAtElement(jve,jel) .EQ. ielReplace) THEN
              rhadapt%p_ImidneighboursAtElement(jve,jel)=iel
              bfound=bfound.OR..TRUE.
            END IF
            
            IF (bfound) EXIT adjacent3
          END DO adjacent3
          
        END IF

        ! If the element could not be found, something is wrong
        IF (.NOT.bfound) THEN
          CALL output_line('Unable to update element neighbor!',&
              OU_CLASS_ERROR,OU_MODE_STD,'remove_element2D')         
          CALL sys_halt()
        END IF
        
      END DO update
      
      ! Copy data from element ielReplace to element IEL
      rhadapt%p_IverticesAtElement(:,iel)     =rhadapt%p_IverticesAtElement(:,ielReplace)
      rhadapt%p_IneighboursAtElement(:,iel)   =rhadapt%p_IneighboursAtElement(:,ielReplace)
      rhadapt%p_ImidneighboursAtElement(:,iel)=rhadapt%p_ImidneighboursAtElement(:,ielReplace)
      
    ELSE

      ! Element iel is the last element
      ielReplace=0
    END IF

    ! Decrease number of elements
    rhadapt%NEL=rhadapt%NEL-1
  END SUBROUTINE remove_element2D

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE update_ElemNeighb2D_1to2(rhadapt,jel,jelmid,iel0,iel,ielmid)

!<description>
    ! This subroutine updates the list of elements adjacent to another 
    ! element and the list of elements mid-adjacent to another element.
    !
    ! The situation is as follows:
    !
    !  +---------------------+            +---------------------+
    !  |                     |            |          .          |
    !  |                     |            |          .          |
    !  |        IEL0         |            |  IELMID  .    IEL   |
    !  |                     |            |          .          |
    !  |                     |            |          .          |
    !  +---------------------+    --->    +----------+----------+
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  |    JEL   .  JELMID  |            |   JEL    .  JELMID  |
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  +---------------------+            +---------------------+
    !
    ! "Sitting" on element IEL0 we want to update the element lists:
    !
    ! 1) If IEL0 is located at the boundary then nothing needs to be done.
    ! 2) If the adjacent element has not been subdivided, that is,
    !    JEL = JELMID then it suffices to update its adjacency
    !    list for the entry IEL0 adopting the values IEL and IELMID.
    ! 3) If the adjacent element has been subdivided, that is,
    !    JEL != JELMID then the adjacency list if each of these two 
    !    elements is updated by the value IEL and IELMID, respectively.
!</description>

!<input>
    ! Number of the neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: jel

    ! Number of the mid-neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: jelmid

    ! Number of the updated macro-element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel0

    ! Number of the new neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Number of the new mid-neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: ielmid
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER :: ive,nve
    LOGICAL :: bfound1,bfound2,bfound

    ! Do nothing for elements adjacent to the boundary
    IF (jel .EQ. 0 .OR. jelmid .EQ. 0) RETURN

    ! Check if adjacent and mid-adjacent elements are the same.
    IF (jel .EQ. jelmid) THEN

      ! Case 1: Adjacent element has not been subdivided.
      bfound1=.FALSE.; bfound2=.FALSE.

      ! What kind of element is neighboring element?
      nve=get_NVE(rhadapt,jel)
      
      ! Loop over all entries in the list of adjacent and/or mid-adjacent
      ! elements for element JEL and check if the value IEL0 is present. 
      ! The situation may occur that element IEL0 is only adjacent to one
      ! "half" of the edge of element JEL. This may be the case, if element
      ! IEL0 was a green triangle which has already been subdivided in the
      ! marking procedure, whiereby element JEL is marked for red refinement.
      ! In order to consider this situation we are looking for element IEL0
      ! both in the list of adjacent and mid-adjacent elements of JEL.
      ! It suffices if IEL0 is found in one of these lists.
      DO ive=1,nve
        IF (rhadapt%p_IneighboursAtElement(ive,jel) .EQ. iel0) THEN
          rhadapt%p_IneighboursAtElement(ive,jel)=iel
          bfound1=.TRUE.
        END IF

        IF (rhadapt%p_ImidneighboursAtElement(ive,jel) .EQ. iel0) THEN
          rhadapt%p_ImidneighboursAtElement(ive,jel)=ielmid
          bfound2=.TRUE.
        END IF
        
        ! Exit if IEL0 has been found in the either the adjacency or the 
        ! mid-adjacency list of element JEL.
        bfound=bfound1.OR.bfound2
        IF (bfound) EXIT
      END DO

    ELSE
      
      ! Case 2: Adjacent element has already been subdivided.
      bfound1=.FALSE.; bfound2=.FALSE.
      
      ! What kind of element is neighboring element 
      nve=get_NVE(rhadapt,jel)

      ! Loop over all entries in the list of adjacent elements for element
      ! JEL and check if the value IEL0 is present. 
      ! If this is the case,  then update the corrseponding entries in the 
      ! lists of (mid-)adjacent element neighbors.
      DO ive=1,nve
        IF (rhadapt%p_IneighboursAtElement(ive,jel) .EQ. iel0) THEN
          rhadapt%p_IneighboursAtElement(ive,jel)   =ielmid
          rhadapt%p_ImidneighboursAtElement(ive,jel)=ielmid
          bfound1=.TRUE.
          EXIT
        END IF
      END DO
      
      ! What kind of element is neighboring element 
      nve=get_NVE(rhadapt,jelmid)
      
      ! Loop over all entries in the list of adjacent elements for element
      ! JELMID and check if the value IEL0 is present.
      ! If this is the case, then update the corrseponding entries in the 
      ! lists of (mid-)adjacent element neighbors.
      DO ive=1,nve
        IF (rhadapt%p_IneighboursAtElement(ive,jelmid) .EQ. iel0) THEN
          rhadapt%p_IneighboursAtElement(ive,jelmid)   =iel
          rhadapt%p_ImidneighboursAtElement(ive,jelmid)=iel
          bfound2=.TRUE.
          EXIT
        END IF
      END DO

      ! Check success of both searches
      bfound=bfound1.AND.bfound2

    END IF

    IF (.NOT.bfound) THEN
      CALL output_line('Inconsistent adjacency lists!',&
          OU_CLASS_ERROR,OU_MODE_STD,'update_ElemNeighb2D_1to2')
      CALL sys_halt()
    END IF
  END SUBROUTINE update_ElemNeighb2D_1to2

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE update_ElemNeighb2D_2to2(rhadapt,jel,jelmid,iel0,ielmid0,iel,ielmid)

!<description>
    ! This subroutine updates the list of elements adjacent to another 
    ! element and the list of elements mid-adjacent to another element.
    !
    ! The situation is as follows:
    !
    !  +---------------------+            +---------------------+
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  | IELMID0  .   IEL0   |            | IELMID   .   IEL    |
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  +---------------------+    --->    +----------+----------+
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  |    JEL   .  JELMID  |            |    JEL   .   JELMID |
    !  |          .          |            |          .          |
    !  |          .          |            |          .          |
    !  +---------------------+            +---------------------+
    !
    ! If IEL0 and IELMID0 are the same, then thins subroutine is identical
    ! to subroutine update_EkemNeighb2D_1to2 which is called in this case
!</description>

!<input>
    ! Number of the neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: jel

    ! Number of the mid-neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: jelmid

    ! Number of the updated macro-element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel0

    ! Number of the updated macro-element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: ielmid0

    ! Number of the new neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Number of the new mid-neighboring element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: ielmid
!</input>

!<inputoutput>
    ! Adaptive data structure
    TYPE(t_hadapt), INTENT(INOUT) :: rhadapt
!</inputoutput>
!</subroutine>
    
    ! local variables
    INTEGER :: ive,nve
    LOGICAL :: bfound1,bfound2,bfound
    
    ! Check if IEL0 and IELMID0 are the same. 
    ! In this case call the corresponding subroutine
    IF (iel0 .EQ. ielmid0) THEN
      CALL update_ElemNeighb2D_1to2(rhadapt,jel,jelmid,iel0,iel,ielmid)
      RETURN
    END IF

    ! Do nothing for elements adjacent to the boundary
    IF (jel .EQ. 0 .OR. jelmid .EQ. 0) RETURN
    
    ! Check if adjacent and mid-adjacent elements are the same.
    IF (jel .EQ. jelmid) THEN

      ! Case 1: Adjacent element has not been subdivided, that is, the 
      !         current edge contains a hanging node for the adjacent element.
      bfound=.FALSE.

      ! What kind of element is neighboring element?
      nve=get_NVE(rhadapt,jel)
      
      ! Loop over all entries in the list of adjacent elements for element 
      ! JEL and check if the value IEL0 or IELMID0 is present. 
      ! If this is the case, then update the corrseponding entries in the 
      ! lists of (mid-)adjacent element neighbors.
      DO ive=1,nve
        IF (rhadapt%p_IneighboursAtElement(ive,jel)    .EQ. iel0 .AND.&
            rhadapt%p_ImidneighboursAtElement(ive,jel) .EQ. ielmid0) THEN
          rhadapt%p_IneighboursAtElement(ive,jel)   =iel
          rhadapt%p_ImidneighboursAtElement(ive,jel)=ielmid
          bfound=.TRUE.
          EXIT
        END IF
      END DO
                
    ELSE
      
      ! Case 2: Adjacent element has already been subdivided.
      bfound1=.FALSE.; bfound2=.FALSE.
      
      ! What kind of element is neighboring element 
      nve=get_NVE(rhadapt,jel)

      ! Loop over all entries in the list of adjacent elements for element
      ! JEL and check if the value IELMID0 is present. 
      ! If this is the case,  then update the corrseponding entries in the 
      ! lists of (mid-)adjacent element neighbors.
      DO ive=1,nve
        IF (rhadapt%p_IneighboursAtElement(ive,jel) .EQ. ielmid0 .AND.&
            rhadapt%p_ImidneighboursAtElement(ive,jel) .EQ. ielmid0) THEN
          rhadapt%p_IneighboursAtElement(ive,jel)   =ielmid
          rhadapt%p_ImidneighboursAtElement(ive,jel)=ielmid
          bfound1=.TRUE.
          EXIT
        END IF
      END DO
      
      ! What kind of element is neighboring element 
      nve=get_NVE(rhadapt,jelmid)
      
      ! Loop over all entries in the list of adjacent elements for element
      ! JELMID and check if the value IEL0 is present.
      ! If this is the case, then update the corrseponding entries in the 
      ! lists of (mid-)adjacent element neighbors.
      DO ive=1,nve
        IF (rhadapt%p_IneighboursAtElement(ive,jelmid) .EQ. iel0 .AND.&
            rhadapt%p_ImidneighboursAtElement(ive,jelmid) .EQ. iel0) THEN
          rhadapt%p_IneighboursAtElement(ive,jelmid)   =iel
          rhadapt%p_ImidneighboursAtElement(ive,jelmid)=iel
          bfound2=.TRUE.
          EXIT
        END IF
      END DO

      bfound=bfound1.AND.bfound2
    END IF

    IF (.NOT.bfound) THEN
      CALL output_line('Inconsistent adjacency lists!',&
          OU_CLASS_ERROR,OU_MODE_STD,'update_ElemNeighb2D_2to2')
      CALL sys_halt()
    END IF
  END SUBROUTINE update_ElemNeighb2D_2to2

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE update_AllElementNeighbors2D(rhadapt,iel0,iel)

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
    LOGICAL :: bfound

    ! Check if the old element is still present in the triangulation
    IF (iel0 > rhadapt%NEL) THEN
      RETURN
    END IF

    ! Loop over adjacent elements
    adjacent: DO ive=1,get_NVE(rhadapt,iel0)
      ! Get number of adjacent element
      jel=rhadapt%p_IneighboursAtElement(ive,iel0)

      ! Are we at the boundary?
      IF (jel .EQ. 0) CYCLE adjacent

      ! Initialize indicator
      bfound=.FALSE.

      ! Find position of element IEL0 in adjacent element JEL
      DO jve=1,get_NVE(rhadapt,jel)
        IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. iel0) THEN
          rhadapt%p_IneighboursAtElement(jve,jel)=iel
          bfound=.TRUE.
          EXIT
        END IF
      END DO

      ! Find position of element IEL0 in mid-adjacent element JEL
      DO jve=1,get_NVE(rhadapt,jel)
        IF (rhadapt%p_ImidneighboursAtElement(jve,jel) .EQ. iel0) THEN
          rhadapt%p_ImidneighboursAtElement(jve,jel)=iel
          bfound=(bfound .AND. .TRUE.)
          EXIT
        END IF
      END DO

      ! If the old element number was not found in adjacent element JEL
      ! then something went wrong and we should not proceed.
      IF (.NOT.bfound) THEN
        CALL output_line('Inconsistent adjacency lists!',&
            OU_CLASS_ERROR,OU_MODE_STD,'update_AllElementNeighbors2D')
        CALL sys_halt()
      END IF
    END DO adjacent
  END SUBROUTINE update_AllElementNeighbors2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Tria2Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into two triangular 
    ! elements by subdividing one edge. The local number of the edge that is
    ! subdivided can be uniquely determined from the element marker.
    !
    ! In the illustration, i1,i2,i3 and i4 denote the vertices whereby i4 
    ! is the new vertex. Moreover, (e1)-(e6) stand for the element numbers
    ! which are adjecent and mid-adjacent to the current element.
    ! The new element is assigned the total number of elements currently present
    ! in the triangulation increased by one.
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(6) :: Ielements
    INTEGER(PREC_VERTEXIDX), DIMENSION(4)  :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4
    INTEGER                    :: loc1,loc2,loc3
    
    ! Determine the local position of edge to 
    ! be subdivided from the element marker
    SELECT CASE(imarker)
    CASE(MARK_REF_TRIA2TRIA_1)
      loc1=1; loc2=2; loc3=3

    CASE(MARK_REF_TRIA2TRIA_2)
      loc1=2; loc2=3; loc3=1

    CASE(MARK_REF_TRIA2TRIA_3)
      loc1=3; loc2=1; loc3=2

    CASE DEFAULT
      CALL output_line('Invalid element marker',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria2Tria')
      CALL sys_halt()
    END SELECT
    
    ! Store vertex- and element-values of the current element
    i1=rhadapt%p_IverticesAtElement(loc1,iel)
    i2=rhadapt%p_IverticesAtElement(loc2,iel)
    i3=rhadapt%p_IverticesAtElement(loc3,iel)

    e1=rhadapt%p_IneighboursAtElement(loc1,iel)
    e2=rhadapt%p_IneighboursAtElement(loc2,iel)
    e3=rhadapt%p_IneighboursAtElement(loc3,iel)
    
    e4=rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e5=rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e6=rhadapt%p_ImidneighboursAtElement(loc3,iel)
    
    ! Store total number of elements before refinement
    nel0=rhadapt%NEL

    
    ! Add one new vertex I4 at the midpoint of edge (I1,I2).
    CALL add_vertex2D(rhadapt,i1,i2,e1,i4,rcollection,fcb_hadaptCallback)
    
    ! Replace element IEL and add one new element numbered NEL0+1
    CALL replace_element2D(rhadapt,iel,i1,i4,i3,e1,nel0+1,e3,e1,nel0+1,e6)
    CALL add_element2D(rhadapt,i2,i3,i4,e2,iel,e4,e5,iel,e4)

    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,nel0+1,nel0+1)

    
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria2Tria')
      CALL sys_halt()
    END IF

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel,   ipos)


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4/); Ielements=(/e1,e2,e3,e4,e5,e6/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA2TRIA,&
         Ivertices,Ielements)
    END IF
  END SUBROUTINE refine_Tria2Tria

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE refine_Tria3Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into three triangular 
    ! elements by subdividing the longest edge and connecting the new vertex 
    ! with the opposite midpoint. The local numbers of the two edges that are
    ! subdivided can be uniquely determined from the element marker.
    !
    ! In the illustration, i1,i2,i3,i4 and i5 denote the vertices whereby i4 and
    ! i5 are the new vertices. Moreover, (e1)-(e6) stand for the element numbers
    ! which are adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently present
    ! in the triangulation increased by one and two, respectively.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           / \                       /|\
    !     (e3) /   \ (e5)           (e3) / |*\ (e5)
    !         /     \                   /  |ne\
    !        /       \          ->     /   |l+2+i5
    !       /   iel   \               /    |  / \
    ! (e6) /           \ (e2)   (e6) / iel | /nel\ (e2)
    !     /*            \           /*     |/ +1 *\
    !    +---------------+         +-------+-------+
    !   i1 (e1)     (e4) i2       i1 (e1) i4  (e4) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(5) :: Ielements
    INTEGER(PREC_VERTEXIDX), DIMENSION(5)  :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5
    INTEGER                    :: loc1,loc2,loc3
    REAL(DP)                   :: dlen12,dlen23,x,y
    
    ! Find corresponding edges according to convention to be subdivided
    SELECT CASE(imarker)
    CASE(MARK_REF_TRIA3TRIA_12)
      loc1=1; loc2=2; loc3=3

    CASE(MARK_REF_TRIA3TRIA_23)
      loc1=2; loc2=3; loc3=1

    CASE(MARK_REF_TRIA3TRIA_13)
      loc1=3; loc2=1; loc3=2

    CASE DEFAULT
      CALL output_line('Invalid element marker!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria3Tria')
      CALL sys_halt()
    END SELECT

    ! Store vertex- and element-values of the current element
    i1=rhadapt%p_IverticesAtElement(loc1,iel)
    i2=rhadapt%p_IverticesAtElement(loc2,iel)
    i3=rhadapt%p_IverticesAtElement(loc3,iel)

    e1=rhadapt%p_IneighboursAtElement(loc1,iel)
    e2=rhadapt%p_IneighboursAtElement(loc2,iel)
    e3=rhadapt%p_IneighboursAtElement(loc3,iel)

    e4=rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e5=rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e6=rhadapt%p_ImidneighboursAtElement(loc3,iel)
    
    ! Store total number of elements before refinement
    nel0=rhadapt%NEL
    

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
      
      ! Replace element IEL and add two new elements numbered NEL0+1 and NEL0+2
      CALL replace_element2D(rhadapt,iel,i1,i4,i3,e1,nel0+2,e3,e1,nel0+2,e6)
      CALL add_element2D(rhadapt,i2,i5,i4,e2,nel0+2,e4,e2,nel0+2,e4)
      CALL add_element2D(rhadapt,i3,i4,i5,iel,nel0+1,e5,iel,nel0+1,e5)
      

      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+1,iel)
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,nel0+2,nel0+1)

            
      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria3Tria')
        CALL sys_halt()
      END IF
      
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel,   ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel,   ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)


      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/i1,i2,i3,i4,i5/); Ielements=(/e1,e2,e3,e4,e5/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA3TRIA12,&
            Ivertices,Ielements)
      END IF
      
    ELSE
      
      ! 2nd CASE: longest edge is (I2,I3)
      
      ! Replace element IEL and add two new elements numbered NEL0+1 and NEL0+2
      CALL replace_element2D(rhadapt,iel,i1,i5,i3,nel0+2,e5,e3,nel0+2,e5,e6)
      CALL add_element2D(rhadapt,i2,i5,i4,e2,nel0+2,e4,e2,nel0+2,e4)
      CALL add_element2D(rhadapt,i1,i4,i5,e1,nel0+1,iel,e1,nel0+1,iel)
      
      
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+1,nel0+2)
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,iel,nel0+1)
      
      
      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria3Tria')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,nel0+2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,   ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)


      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/i1,i2,i3,i4,i5/); Ielements=(/e1,e2,e3,e4,e5/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA3TRIA23,&
            Ivertices,Ielements)
      END IF
    END IF
  END SUBROUTINE refine_Tria3Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Tria4Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one triangular element into four similar 
    ! triangular elements by connecting the three edge midpoints.
    !
    ! In the illustration, i1-i6 denote the vertices whereby i4-i6 are the 
    ! new vertices. Moreover, (e1)-(e6) stand for the element numbers which
    ! are adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently present
    ! in the triangulation increased by one, two and three, respectively.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           / \                       /*\
    !     (e3) /   \ (e5)           (e3) /nel\ (e5)
    !         /     \                   / +3  \
    !        /       \          ->  i6 +-------+ i5
    !       /   iel   \               / \ iel / \
    ! (e6) /           \ (e2)   (e6) /nel\   /nel\ (e2)
    !     /*            \           /* +1 \*/ +2 *\
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(6) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(6) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6

    ! Store vertex- and element-values of the current element
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i2=rhadapt%p_IverticesAtElement(2,iel)
    i3=rhadapt%p_IverticesAtElement(3,iel)

    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e2=rhadapt%p_IneighboursAtElement(2,iel)
    e3=rhadapt%p_IneighboursAtElement(3,iel)

    e4=rhadapt%p_ImidneighboursAtElement(1,iel)
    e5=rhadapt%p_ImidneighboursAtElement(2,iel)
    e6=rhadapt%p_ImidneighboursAtElement(3,iel)
        
    ! Store total number of elements before refinement
    nel0=rhadapt%NEL


    ! Add three new vertices I4,I5 and I6 at the midpoint of edges (I1,I2),
    ! (I2,I3) and (I1,I3), respectively. 
    CALL add_vertex2D(rhadapt,i1,i2,e1,i4,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i2,i3,e2,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i1,e3,i6,rcollection,fcb_hadaptCallback)
    

    ! Replace element IEL and add three new elements NEL0+1, NEL0+2 and NEL0+3
    CALL replace_element2D(rhadapt,iel,i4,i5,i6,nel0+2,nel0+3,nel0+1,nel0+2,nel0+3,nel0+1)
    CALL add_element2D(rhadapt,i1,i4,i6,e1,iel,e6,e1,iel,e6)
    CALL add_element2D(rhadapt,i2,i5,i4,e2,iel,e4,e2,iel,e4)
    CALL add_element2D(rhadapt,i3,i6,i5,e3,iel,e5,e3,iel,e5)

    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel,nel0+2,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel,nel0+3,nel0+2)
    CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel,nel0+1,nel0+3)

    
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria4Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria4Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Tria4Tria')
      CALL sys_halt()
    END IF
   
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel   ,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel   ,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel   ,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+3,ipos)


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_TRIA4TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE refine_Tria4Tria

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE refine_Quad2Quad(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element into two quadrilateral
    ! elements. The local number of the edge that is subdivided can be uniquely
    ! determined from the element marker.
    !
    ! In the illustration, i1-i6 denote the vertices whereby i5 and i6
    ! are the new vertices. Moreover, (e1)-(e8) stand for the element numbers
    ! which are adjecent and mid-adjacent to the current element.
    ! The new element is assigned the total number of elements currently present
    ! in the triangulation increased by one.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i3          i4 (e7)  i6 (e3) i3
    !      +---------------+            +-------+-------+
    !      |               |            |       |      *|
    ! (e4) |               | (e6)  (e4) |       |       | (e6)
    !      |               |            |       |       |
    !      |     iel       |     ->     |  iel  | nel+1 |
    !      |               |            |       |       |
    ! (e8) |               | (e2)  (e8) |       |       | (e2)
    !      |*              |            |*      |       |
    !      +---------------+            +-------+-------+
    !     i1 (e1)     (e5) i2          i1 (e1)  i5 (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(6) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6
    INTEGER                    :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    SELECT CASE(imarker)
    CASE(MARK_REF_QUAD2QUAD_13)
      loc1=1; loc2=2; loc3=3; loc4=4

    CASE(MARK_REF_QUAD2QUAD_24)
      loc1=2; loc2=3; loc3=4; loc4=1

    CASE DEFAULT
      CALL output_line('Invalid element marker!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad2Quad')
      CALL sys_halt()
    END SELECT
    
    ! Store vertex- and element-values of the current element
    i1=rhadapt%p_IverticesAtElement(loc1,iel)
    i2=rhadapt%p_IverticesAtElement(loc2,iel)
    i3=rhadapt%p_IverticesAtElement(loc3,iel)
    i4=rhadapt%p_IverticesAtElement(loc4,iel)

    e1=rhadapt%p_IneighboursAtElement(loc1,iel)
    e2=rhadapt%p_IneighboursAtElement(loc2,iel)
    e3=rhadapt%p_IneighboursAtElement(loc3,iel)
    e4=rhadapt%p_IneighboursAtElement(loc4,iel)

    e5=rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e6=rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e7=rhadapt%p_ImidneighboursAtElement(loc3,iel)
    e8=rhadapt%p_ImidneighboursAtElement(loc4,iel)
    
    ! Store total number of elements before refinement
    nel0=rhadapt%NEL
    
    ! Add two new vertices I5 and I6 at the midpoint of edges (I1,I2) and
    ! (I3,I4), respectively.
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i4,e3,i6,rcollection,fcb_hadaptCallback)
    

    ! Replace element IEL and add one new element NEL0+1
    CALL replace_element2D(rhadapt,iel,i1,i5,i6,i4,e1,nel0+1,e7,e4,e1,nel0+1,e7,e8)
    CALL add_element2D(rhadapt,i3,i6,i5,i2,e3,iel,e5,e2,e3,iel,e5,e6)
    

    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel,nel0+1,iel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel,nel0+1,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel,nel0+1)

    
    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad2Quad')
      CALL sys_halt()
    END IF

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD2QUAD,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE refine_Quad2Quad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Quad3Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element into three triangular
    ! elements. The local number of the edge that is
    ! subdivided can be uniquely determined from the element marker.
    !
    ! In the illustration, i1-i5 denote the vertices whereby i5 is the new 
    ! vertex. Moreover, (e1)-(e8) stand for the element numbers which are 
    ! adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently 
    ! present in the triangulation increased by one and two, respectively.
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
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(5) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5
    INTEGER                    :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    SELECT CASE(imarker)
    CASE(MARK_REF_QUAD3TRIA_1)
      loc1=1; loc2=2; loc3=3; loc4=4

    CASE(MARK_REF_QUAD3TRIA_2)
      loc1=2; loc2=3; loc3=4; loc4=1

    CASE(MARK_REF_QUAD3TRIA_3)
      loc1=3; loc2=4; loc3=1; loc4=2

    CASE(MARK_REF_QUAD3TRIA_4)
      loc1=4; loc2=1; loc3=2; loc4=3

    CASE DEFAULT
      CALL output_line('Invalid marker element!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad3Tria')
      CALL sys_halt()
    END SELECT
    

    ! Store vertex- and element-values of the current element
    i1=rhadapt%p_IverticesAtElement(loc1,iel)
    i2=rhadapt%p_IverticesAtElement(loc2,iel)
    i3=rhadapt%p_IverticesAtElement(loc3,iel)
    i4=rhadapt%p_IverticesAtElement(loc4,iel)


    e1=rhadapt%p_IneighboursAtElement(loc1,iel)
    e2=rhadapt%p_IneighboursAtElement(loc2,iel)
    e3=rhadapt%p_IneighboursAtElement(loc3,iel)
    e4=rhadapt%p_IneighboursAtElement(loc4,iel)

    e5=rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e6=rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e7=rhadapt%p_ImidneighboursAtElement(loc3,iel)
    e8=rhadapt%p_ImidneighboursAtElement(loc4,iel)

    ! Store total number of elements before refinement
    nel0=rhadapt%NEL
    

    ! Add one new vertex I5 it the midpoint of edge (I1,I2)
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)

    
    ! Replace element IEL and add two new elements NEL0+1 and NEL0+2
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
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad3Tria')
      CALL sys_halt()
    END IF   
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad3Tria')
      CALL sys_halt()
    END IF

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+2,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)+1
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)-1
    

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD3TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE refine_Quad3Tria
  
! ***************************************************************************

!<subroutine>

  SUBROUTINE refine_Quad4Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element into four triangular
    ! elements. The local numbers of the edges that are subdivided can be 
    ! uniquely determined from the element marker.
    !
    ! In the illustration, i1-i6 denote the vertices whereby i5 and i6 are the 
    ! new vertices. Moreover, (e1)-(e8) stand for the element numbers which are
    ! adjecent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently 
    ! present in the triangulation increased by one, two an three, respectively.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7)     (e3) i3
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(6) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6
    INTEGER                    :: loc1,loc2,loc3,loc4
    
    ! Find local position of edge to be subdivided
    SELECT CASE(imarker)
    CASE(MARK_REF_QUAD4TRIA_12)
      loc1=1; loc2=2; loc3=3; loc4=4

    CASE(MARK_REF_QUAD4TRIA_23)
      loc1=2; loc2=3; loc3=4; loc4=1

    CASE(MARK_REF_QUAD4TRIA_34)
      loc1=3; loc2=4; loc3=1; loc4=2

    CASE(MARK_REF_QUAD4TRIA_14)
      loc1=4; loc2=1; loc3=2; loc4=3

    CASE DEFAULT
      CALL output_line('Invalid element marker!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Tria')
      CALL sys_halt()
    END SELECT

    ! Store vertex- and element-values of the current element
    i1=rhadapt%p_IverticesAtElement(loc1,iel)
    i2=rhadapt%p_IverticesAtElement(loc2,iel)
    i3=rhadapt%p_IverticesAtElement(loc3,iel)
    i4=rhadapt%p_IverticesAtElement(loc4,iel)

    e1=rhadapt%p_IneighboursAtElement(loc1,iel)
    e2=rhadapt%p_IneighboursAtElement(loc2,iel)
    e3=rhadapt%p_IneighboursAtElement(loc3,iel)
    e4=rhadapt%p_IneighboursAtElement(loc4,iel)

    e5=rhadapt%p_ImidneighboursAtElement(loc1,iel)
    e6=rhadapt%p_ImidneighboursAtElement(loc2,iel)
    e7=rhadapt%p_ImidneighboursAtElement(loc3,iel)
    e8=rhadapt%p_ImidneighboursAtElement(loc4,iel)
    
    ! Store total number of elements before refinement
    nel0=rhadapt%NEL
    

    ! Add new vertices I5 and I6 at the midpoint of edges (I1,I2) and
    ! (I2,I3), respectively.
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i2,i3,e2,i6,rcollection,fcb_hadaptCallback)

    
    ! Replace element IEL and add three new elements NEL0+1,NEL0+2 and NEL0+3
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
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Tria')
      CALL sys_halt()
    END IF
    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+3,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)+1
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)-1


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD4TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE refine_Quad4Tria
    
    ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE refine_Quad4Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine subdivides one quadrilateral element four similar 
    ! quadrilateral elements.
    !
    ! In the illustration, i1-i9 denote the vertices whereby i5-i9 are the new 
    ! vertices. Moreover, (e1)-(e8) stand for the element numbers which are 
    ! adjacent and mid-adjacent to the current element.
    ! The new elements are assigned the total number of elements currently present
    ! in the triangulation increased by one, two and three, respectively.
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Store vertex- and element-values of the current element
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i2=rhadapt%p_IverticesAtElement(2,iel)
    i3=rhadapt%p_IverticesAtElement(3,iel)
    i4=rhadapt%p_IverticesAtElement(4,iel)

    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e2=rhadapt%p_IneighboursAtElement(2,iel)
    e3=rhadapt%p_IneighboursAtElement(3,iel)
    e4=rhadapt%p_IneighboursAtElement(4,iel)

    e5=rhadapt%p_ImidneighboursAtElement(1,iel)
    e6=rhadapt%p_ImidneighboursAtElement(2,iel)
    e7=rhadapt%p_ImidneighboursAtElement(3,iel)
    e8=rhadapt%p_ImidneighboursAtElement(4,iel)

    ! Store total number of elements before refinement
    nel0=rhadapt%NEL
    
    
    ! Add five new vertices I5,I6,I7,I8 and I9 at the midpoint of edges
    ! (I1,I2), (I2,I3), (I3,I4) and (I1,I4) and at the center of element IEL
    CALL add_vertex2D(rhadapt,i1,i2,e1,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i2,i3,e2,i6,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i4,e3,i7,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i4,i1,e4,i8,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i1,i2,i3,i4,i9,rcollection,fcb_hadaptCallback)
    

    ! Replace element IEL and add three new elements NEL0+1,NEL0+2 and NEL0+3
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
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'refine_Quad4Quad')
      CALL sys_halt()
    END IF
    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel,   ipos)
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
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_REF_QUAD4QUAD,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE refine_Quad4Quad
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE convert_Tria2Tria(rhadapt,iel,jel,rcollection,fcb_hadaptCallback)

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
    !        /   |   \        ->      i6+-------+i5
    !       /    |    \                / \nel+2/ \
    ! (e6) / iel | jel \ (e2)    (e6) /   \   /   \ (e2)
    !     /*     |     *\            /* iel\*/jel *\
    !    +-------+-------+          +-------+-------+
    !   i1(e1,e7)i4(e4,e8)i2       i1(e1,e7)i4(e4,e8)i2
    !
!</description>

!<input>
    ! Number of first (left) element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Number of second (right) element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: jel

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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(6) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6,e7,e8

    ! Find local positions of elements IEL and JEL
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i4=rhadapt%p_IverticesAtElement(2,iel)
    i3=rhadapt%p_IverticesAtElement(3,iel)
    i2=rhadapt%p_IverticesAtElement(1,jel)
    
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e3=rhadapt%p_IneighboursAtElement(3,iel)
    e2=rhadapt%p_IneighboursAtElement(1,jel)
    e4=rhadapt%p_IneighboursAtElement(3,jel)
    
    e7=rhadapt%p_ImidneighboursAtElement(1,iel)
    e6=rhadapt%p_ImidneighboursAtElement(3,iel)
    e5=rhadapt%p_ImidneighboursAtElement(1,jel)
    e8=rhadapt%p_ImidneighboursAtElement(3,jel)
    
    ! Store total number of elements before conversion
    nel0=rhadapt%NEL

    
    ! Add two new vertices I5 and I6 at the midpoint of edges (I2,I3)
    ! and (I1,I3), respectively.
    CALL add_vertex2D(rhadapt,i2,i3,e2,i5,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i3,i1,e3,i6,rcollection,fcb_hadaptCallback)


    ! Replace elements IEL and JEL and add two new elements NEL0+1 and NEL0+2
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
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Tria2Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Tria2Tria')
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


    ! "Lock" all vertices connected to the four triangles
    rhadapt%p_IvertexAge(i1)=-ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2)=-ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3)=-ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4)=-ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5)=-ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6)=-ABS(rhadapt%p_IvertexAge(i6))
    

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_TRIA2TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE convert_Tria2Tria

! ***************************************************************************

!<subroutine>

  SUBROUTINE convert_Quad2Quad(rhadapt,iel,jel,rcollection,fcb_hadaptCallback)

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
    !     i4(e7,e3)i7(f5,f1)i3         i4(e7,e3)i7(f5,f1)i3
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
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Number of second element
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: jel
    
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e3,e4,e5,e7,e8,f1,f3,f4,f5,f7,f8
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Find local positions of elements IEL and JEL
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i7=rhadapt%p_IverticesAtElement(3,iel)
    i4=rhadapt%p_IverticesAtElement(4,iel)
    i3=rhadapt%p_IverticesAtElement(1,jel)
    i2=rhadapt%p_IverticesAtElement(4,jel)

    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e3=rhadapt%p_IneighboursAtElement(3,iel)
    e4=rhadapt%p_IneighboursAtElement(4,iel)
    f1=rhadapt%p_IneighboursAtElement(1,jel)   
    f3=rhadapt%p_IneighboursAtElement(3,jel)
    f4=rhadapt%p_IneighboursAtElement(4,jel)
    
    e5=rhadapt%p_ImidneighboursAtElement(1,iel)
    e7=rhadapt%p_ImidneighboursAtElement(3,iel)
    e8=rhadapt%p_ImidneighboursAtElement(4,iel)
    f5=rhadapt%p_ImidneighboursAtElement(1,jel)
    f7=rhadapt%p_ImidneighboursAtElement(3,jel)
    f8=rhadapt%p_ImidneighboursAtElement(4,jel)

    ! Store total number of elements before conversion
    nel0=rhadapt%NEL

    
    ! Add two new vertices I6, I8, and I9 at the midpoint of edges (I2,I3),
    ! (I1,I4) and (I1,I2), respectively.
    CALL add_vertex2D(rhadapt,i2,i3,f4,i6,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i4,i1,e4,i8,rcollection,fcb_hadaptCallback)
    CALL add_vertex2D(rhadapt,i1,i2,i3,i4,i9,rcollection,fcb_hadaptCallback)


    ! Replace element IEL and JEL and add two new elements NEL0+1 and NEL0+2
    CALL replace_element2D(rhadapt,iel,i1,i5,i9,i8,e1,jel,nel0+2,e8,e5,jel,nel0+2,e8)
    CALL replace_element2D(rhadapt,jel,i2,i6,i9,i5,f4,nel0+1,iel,f3,f4,nel0+1,iel,f7)
    CALL add_element2D(rhadapt,i3,i7,i9,i6,f1,nel0+2,jel,f8,f5,nel0+2,jel,f8)
    CALL add_element2D(rhadapt,i4,i8,i9,i7,e4,iel,nel0+1,e3,e4,iel,nel0+1,e7)


    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,f4,f8,jel,nel0+1,jel)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel,nel0+2)
    CALL update_ElementNeighbors2D(rhadapt,f1,f5,jel,nel0+1,nel0+1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,nel0+2,nel0+2)


    ! Update list of elements meeting at vertices
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,jel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad2Quad')
      CALL sys_halt()
    END IF

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,nel0+2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,jel,   ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+2,ipos)


    ! "Lock" all vertices connected to the four quadrilaterals
    rhadapt%p_IvertexAge(i1)=-ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2)=-ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3)=-ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4)=-ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5)=-ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6)=-ABS(rhadapt%p_IvertexAge(i6))
    rhadapt%p_IvertexAge(i7)=-ABS(rhadapt%p_IvertexAge(i7))
    rhadapt%p_IvertexAge(i8)=-ABS(rhadapt%p_IvertexAge(i8))
    rhadapt%p_IvertexAge(i9)=-ABS(rhadapt%p_IvertexAge(i9))


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,f4,f1,e4,f3,f8,e3,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_QUAD2QUAD,&
          Ivertices,Ielements)
    END IF
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
    TYPE(t_hadapt), INTENT(INOUT)               :: rhadapt

    ! OPTIONAL: Collection
    TYPE(t_collection), INTENT(INOUT), OPTIONAL :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: nel0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Get local data from elements IEL1, IEL2 and IEL3
    i1=rhadapt%p_IverticesAtElement(1,iel1)
    i5=rhadapt%p_IverticesAtElement(2,iel1)
    i4=rhadapt%p_IverticesAtElement(3,iel1)
    i2=rhadapt%p_IverticesAtElement(1,iel2)
    i3=rhadapt%p_IverticesAtElement(2,iel2)

    e1=rhadapt%p_IneighboursAtElement(1,iel1)
    e4=rhadapt%p_IneighboursAtElement(3,iel1)
    e2=rhadapt%p_IneighboursAtElement(1,iel2)
    e5=rhadapt%p_IneighboursAtElement(3,iel2)
    e3=rhadapt%p_IneighboursAtElement(2,iel3)
    
    e9 =rhadapt%p_ImidneighboursAtElement(1,iel1)
    e8 =rhadapt%p_ImidneighboursAtElement(3,iel1)
    e6 =rhadapt%p_ImidneighboursAtElement(1,iel2)
    e10=rhadapt%p_ImidneighboursAtElement(3,iel2)
    e7 =rhadapt%p_ImidneighboursAtElement(2,iel3)
    
    ! Store total number of elements before conversion
    nel0=rhadapt%NEL

    
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
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad3Tria')
      CALL sys_halt()
    END IF
    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel2,  ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,iel3,  ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,iel3,  ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel1,  ipos)    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel1,  ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel2,  ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,iel3,  ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i9,nel0+1,ipos)


    ! Finally, adjust numbers of triangles/quadrilaterals
    rhadapt%InelOfType(TRIA_NVETRI2D) =rhadapt%InelOfType(TRIA_NVETRI2D)-3
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)+3


    ! "Lock" all vertices connected to the four quadrilaterals
    rhadapt%p_IvertexAge(i1)=-ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2)=-ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3)=-ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4)=-ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5)=-ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6)=-ABS(rhadapt%p_IvertexAge(i6))
    rhadapt%p_IvertexAge(i7)=-ABS(rhadapt%p_IvertexAge(i7))
    rhadapt%p_IvertexAge(i8)=-ABS(rhadapt%p_IvertexAge(i8))
    rhadapt%p_IvertexAge(i9)=-ABS(rhadapt%p_IvertexAge(i9))
    

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_QUAD3TRIA,&
          Ivertices,Ielements)
    END IF
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Get local data from elements IEL1, IEL2 and IEL3
    i1=rhadapt%p_IverticesAtElement(1,iel1)
    i5=rhadapt%p_IverticesAtElement(2,iel1)
    i4=rhadapt%p_IverticesAtElement(3,iel1)
    i2=rhadapt%p_IverticesAtElement(1,iel2)
    i6=rhadapt%p_IverticesAtElement(2,iel2)
    i3=rhadapt%p_IverticesAtElement(1,iel3)

    e1=rhadapt%p_IneighboursAtElement(1,iel1)
    e4=rhadapt%p_IneighboursAtElement(3,iel1)
    e2=rhadapt%p_IneighboursAtElement(1,iel2)
    e5=rhadapt%p_IneighboursAtElement(3,iel2)
    e3=rhadapt%p_IneighboursAtElement(1,iel3)
    e6=rhadapt%p_IneighboursAtElement(3,iel3)

    e9 =rhadapt%p_ImidneighboursAtElement(1,iel1)
    e8 =rhadapt%p_ImidneighboursAtElement(3,iel1)
    e11=rhadapt%p_ImidneighboursAtElement(1,iel2)
    e10=rhadapt%p_ImidneighboursAtElement(3,iel2)
    e7 =rhadapt%p_ImidneighboursAtElement(1,iel3)
    e12=rhadapt%p_ImidneighboursAtElement(3,iel3)
    

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
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel4).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'convert_Quad4Tria')
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
    rhadapt%InelOfType(TRIA_NVETRI2D) =rhadapt%InelOfType(TRIA_NVETRI2D)-4
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)+4


    ! "Lock" all vertices of the four quadrilaterals
    rhadapt%p_IvertexAge(i1)=-ABS(rhadapt%p_IvertexAge(i1))
    rhadapt%p_IvertexAge(i2)=-ABS(rhadapt%p_IvertexAge(i2))
    rhadapt%p_IvertexAge(i3)=-ABS(rhadapt%p_IvertexAge(i3))
    rhadapt%p_IvertexAge(i4)=-ABS(rhadapt%p_IvertexAge(i4))
    rhadapt%p_IvertexAge(i5)=-ABS(rhadapt%p_IvertexAge(i5))
    rhadapt%p_IvertexAge(i6)=-ABS(rhadapt%p_IvertexAge(i6))
    rhadapt%p_IvertexAge(i7)=-ABS(rhadapt%p_IvertexAge(i7))
    rhadapt%p_IvertexAge(i8)=-ABS(rhadapt%p_IvertexAge(i8))
    rhadapt%p_IvertexAge(i9)=-ABS(rhadapt%p_IvertexAge(i9))


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CVT_QUAD4TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE convert_Quad4Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_2Tria1Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines two triangles resulting from a 1-tria : 2-tria
    ! green refinement into the original macro triangle. Note that element IEL 
    ! must(!) be the left triangle which is combined with its right neighbour.
    ! This is not checked in the routine since this check is already performed
    ! in the marking routine.
    !
    ! Due to the numbering convention used in the refinement procedure
    ! the "second" element neighbor is the other green triangle that is removed.
    ! In order to restore the original orientation the outcome of the resulting
    ! triangle is considered. If the macro element belongs to the inital
    ! triangulation of if it is an inner red triangle, an arbitrary orientation
    ! is adopted. If the macro element is an outer red triangle, then the states
    ! 2 and 8 are trapped and converted into state 4.
    !
    !     initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /|\                       / \
    !     (e3) / | \ (e5)           (e3) /   \ (e5)
    !         /  |  \                   /     \
    !        /   |   \        ->       /       \
    !       /    |    \               /   jel   \
    ! (e6) / iel | iel1\ (e2)   (e6) /           \ (e2)
    !     /*     |     *\           /*            \
    !    +-------+-------+         +---------------+
    !    i1 (e1) i4 (e4) i2        i1 (e1)    (e4) i2
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Callback routines
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(4) :: Ivertices
    INTEGER(PREC_VERTEXIDX),  DIMENSION(3) :: ImacroVertices
    INTEGER, DIMENSION(TRIA_NVETRI2D)      :: IvertexAge
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,e1,e2,e3,e4,e5,e6,jel,ielRemove,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4
    INTEGER :: istate
    
          
    ! Get right-adjacent green elements
    iel1=rhadapt%p_IneighboursAtElement(2,iel)
    
    ! Determine element with smaller element number
    IF (iel < iel1) THEN
      jel=iel; ielRemove=iel1
    ELSE
      jel=iel1; ielRemove=iel
    END IF
    
    ! Store vertex- and element values of the two elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i4=rhadapt%p_IverticesAtElement(2,iel)
    i3=rhadapt%p_IverticesAtElement(3,iel)
    i2=rhadapt%p_IverticesAtElement(1,iel1)
    
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e3=rhadapt%p_IneighboursAtElement(3,iel)
    e2=rhadapt%p_IneighboursAtElement(1,iel1)
    e4=rhadapt%p_IneighboursAtElement(3,iel1)
    
    e5=rhadapt%p_ImidneighboursAtElement(1,iel1)
    e6=rhadapt%p_ImidneighboursAtElement(3,iel)

    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel1,iel,jel,jel)
    IF (iel < iel1) THEN
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel1,jel,jel)
    ELSE
      CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel,jel,jel)
    END IF

    
    ! The resulting triangle will possess one of the states STATE_TRIA_OUTERINNERx, whereby
    ! x can be blank, 1 or 2. Due to our refinement convention, the two states x=1 and x=2
    ! should not appear, that is, local numbering of the resulting triangle starts at the
    ! vertex which is opposite to the inner red triangle. To this end, we check the state of
    ! the provisional triangle (I1,I2,I3) and transform the orientation accordingly.
    ImacroVertices = (/i1,i2,i3/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate = redgreen_getstateTria(IvertexAge)
    
    SELECT CASE(istate)
    CASE(STATE_TRIA_OUTERINNER,&
        STATE_TRIA_ROOT,STATE_TRIA_REDINNER)
      ! Update element JEL = (I1,I2,I3)
      CALL replace_element2D(rhadapt,jel,i1,i2,i3,e1,e2,e3,e4,e5,e6)
      
    CASE(STATE_TRIA_OUTERINNER2)
      ! Update element JEL = (I2,I3,I1)
      CALL replace_element2D(rhadapt,jel,i2,i3,i1,e2,e3,e1,e5,e6,e4)
      
    CASE(STATE_TRIA_OUTERINNER1)
      ! Update element JEL = (I3,I1,I2)
      CALL replace_element2D(rhadapt,jel,i3,i1,i2,e3,e1,e2,e6,e4,e5)
      
    CASE DEFAULT
      CALL output_line('Invalid state of resulting triangle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
      CALL sys_halt()
    END SELECT
    

    ! Delete element IEL or IEL1 depending on which element has smaller
    ! element number, that is, is not equal to JEL
    CALL remove_element2D(rhadapt,ielRemove,ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,ielRemove)
    
    
    ! Update list of elements meeting at vertices
    IF (ielRemove .EQ. iel1) THEN
      
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,ielRemove).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,ielRemove).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel,ipos)
      
    ELSE
      
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,ielRemove).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,ielRemove).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Tria1Tria')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel,ipos)
      
    END IF
    
    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4/); Ielements=(/e1,e2,e3,e4/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_2TRIA1TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_2Tria1Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Tria1Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-tria : 4-tria refinement into the original macro triangle.
    ! By definition, iel is the number of the inner red triangle.
    !
    !    initial triangle           subdivided triangle
    !
    !            i3                        i3
    !            +                         +
    !           /*\                       / \
    !     (e3) /   \ (e5)           (e3) /   \ (e5)
    !         / iel3\                   /     \
    !      i6+-------+i5      ->       /       \
    !       / \ iel / \               /   jel   \
    ! (e6) /   \   /   \ (e2)   (e6) /           \ (e2)
    !     /*iel1\*/iel2*\           /             \
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(6) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IsortedElements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(6) :: Ivertices
    INTEGER(PREC_VERTEXIDX),  DIMENSION(3) :: ImacroVertices
    INTEGER, DIMENSION(TRIA_NVETRI2D)      :: IvertexAge
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,jel,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6
    INTEGER                    :: istate
    
    ! Retrieve patch of elements
    iel2=rhadapt%p_IneighboursAtElement(1,iel)
    iel3=rhadapt%p_IneighboursAtElement(2,iel)
    iel1=rhadapt%p_IneighboursAtElement(3,iel)


    ! Store vertex- and element-values of the three neighboring elements
    i4=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i6=rhadapt%p_IverticesAtElement(3,iel)
    i1=rhadapt%p_IverticesAtElement(1,iel1)
    i2=rhadapt%p_IverticesAtElement(1,iel2)
    i3=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel1)
    e6=rhadapt%p_IneighboursAtElement(3,iel1)
    e2=rhadapt%p_IneighboursAtElement(1,iel2)
    e4=rhadapt%p_IneighboursAtElement(3,iel2)
    e3=rhadapt%p_IneighboursAtElement(1,iel3)
    e5=rhadapt%p_IneighboursAtElement(3,iel3)


    ! Sort the four elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements=(/iel,iel1,iel2,iel3/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel=IsortedElements(1)


    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel2,iel1,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel3,iel2,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel1,iel3,jel,jel)

    
    ! The resulting triangle will posses one of the states STATE_TRIA_REDINNER or
    ! STATE_TRIA_OUTERINNERx, whereby x can be blank, 1 or 2. Due to our refinement 
    ! convention, the two states x=1 and x=2 should not appear, that is, local
    ! numbering of the resulting triangle starts at the vertex which is opposite 
    ! to the inner red triangle. To this end, we check the state of the provisional
    ! triangle (I1,I2,I3) and transform the orientation accordingly.
    ImacroVertices = (/i1,i2,i3/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate = redgreen_getstateTria(IvertexAge)
    
    SELECT CASE(istate)
    CASE(STATE_TRIA_OUTERINNER,&
         STATE_TRIA_ROOT,STATE_TRIA_REDINNER)
      ! Update element JEL = (I1,I2,I3)
      CALL replace_element2D(rhadapt,jel,i1,i2,i3,e1,e2,e3,e4,e5,e6)
      
    CASE(STATE_TRIA_OUTERINNER1)
      ! Update element JEL = (I3,I1,I2)
      CALL replace_element2D(rhadapt,jel,i3,i1,i2,e3,e1,e2,e6,e4,e5)

    CASE(STATE_TRIA_OUTERINNER2)
      ! Update element JEL = (I2,I3,I1)
      CALL replace_element2D(rhadapt,jel,i2,i3,i1,e2,e3,e1,e5,e6,e4)

    CASE DEFAULT
      CALL output_line('Invalid state of resulting triangle!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      CALL sys_halt()
    END SELECT


    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
    ! element corresponds to element with minimum number JEL
    CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))

    CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))

    CALL remove_element2D(rhadapt,IsortedElements(2),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(2))


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Tria')
      CALL sys_halt()
    END IF

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel,ipos)


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA1TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_4Tria1Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Tria2Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-tria : 4-tria refinement into two green triangle.
    ! The local position (node 1,2,3 of the interior element) of the midpoint 
    ! vertex that should be kept is identified by the marker imarker. 
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
    !         /iel3 \                   /  |  \
    !     i6 +-------+ i5     ->       /   |   \
    !       / \ iel / \               /    |    \
    ! (e6) /   \   /   \ (e2)   (e6) /     |     \ (e2)
    !     /*iel1\*/iel2*\           /* jel1|jel2 *\
    !    +-------+-------+         +-------+-------+
    !    i1 (e1) i4 (e4) i2        i1 (e1) i4 (e4)  i2
    !
!</description>

!<input>
    ! Element number of the inner red triangle
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Identifiert for element marker
    INTEGER, INTENT(IN)                  :: imarker

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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(6) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IsortedElements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(6) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX) :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,jel1,jel2,ielReplace
    INTEGER(PREC_VERTEXIDX)  :: i1,i2,i3,i4,i5,i6
    
    ! Retrieve patch of elements
    iel2=rhadapt%p_IneighboursAtElement(1,iel)
    iel3=rhadapt%p_IneighboursAtElement(2,iel)
    iel1=rhadapt%p_IneighboursAtElement(3,iel)

    ! Store vertex- and element-values of the three neighboring elements
    i4=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i6=rhadapt%p_IverticesAtElement(3,iel)
    i1=rhadapt%p_IverticesAtElement(1,iel1)
    i2=rhadapt%p_IverticesAtElement(1,iel2)
    i3=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel1)
    e6=rhadapt%p_IneighboursAtElement(3,iel1)
    e2=rhadapt%p_IneighboursAtElement(1,iel2)
    e4=rhadapt%p_IneighboursAtElement(3,iel2)
    e3=rhadapt%p_IneighboursAtElement(1,iel3)
    e5=rhadapt%p_IneighboursAtElement(3,iel3)


    ! Sort the four elements according to their number and
    ! determine the two elements withthe smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel1=IsortedElements(1)
    jel2=IsortedElements(2)


    ! Which midpoint vertex should be kept?
    SELECT CASE(imarker)
    CASE(MARK_CRS_4TRIA2TRIA_1)     
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel2,iel1,jel2,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel3,iel2,jel2,jel2)
      CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel1,iel3,jel1,jel1)


      ! Update elements JEL1 and JEL2
      CALL replace_element2D(rhadapt,jel1,i1,i4,i3,e1,jel2,e3,e1,jel2,e6)
      CALL replace_element2D(rhadapt,jel2,i2,i3,i4,e2,jel1,e4,e5,jel1,e4)
      

      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
      ! elements correspond to the two elements with smallest numbers
      CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))

      CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))

      
      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF

      ! Note, this can be improved by checking against JEL1 and JEL2
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel2,ipos)


      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA2TRIA1,&
            Ivertices,Ielements)
      END IF


    CASE(MARK_CRS_4TRIA2TRIA_2)
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel2,iel1,jel1,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel3,iel2,jel2,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel1,iel3,jel2,jel2)


      ! Update elements JEL1 and JEL2
      CALL replace_element2D(rhadapt,jel1,i2,i5,i1,e2,jel2,e1,e2,jel2,e4)
      CALL replace_element2D(rhadapt,jel2,i3,i1,i5,e3,jel1,e5,e6,jel1,e5)


      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
      ! elements correspond to the two elements with smallest numbers
      CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))

      CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))


      ! Update list of elements meeting at vertices
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF

      ! Note, this can be improved by checking against JEL1 and JEL2
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel2,ipos)


      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA2TRIA2,&
            Ivertices,Ielements)
      END IF


    CASE(MARK_CRS_4TRIA2TRIA_3)
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e4,iel2,iel1,jel2,jel2)
      CALL update_ElementNeighbors2D(rhadapt,e2,e5,iel3,iel2,jel1,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e3,e6,iel1,iel3,jel2,jel1)
      

      ! Update elements JEL1 and JEL2
      CALL replace_element2D(rhadapt,jel1,i3,i6,i2,e3,jel2,e2,e3,jel2,e5)
      CALL replace_element2D(rhadapt,jel2,i1,i2,i6,e1,jel1,e6,e4,jel1,e6)


      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
      ! elements correspond to the two elements with smallest numbers
      CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))

      CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))


      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF

      ! Note, this can be improved by checking against JEL1 and JEL2
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel2,ipos)
      

      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/i1,i2,i3,i4,i5,i6/); Ielements=(/e1,e2,e3,e4,e5,e6/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA2TRIA3,&
            Ivertices,Ielements)
      END IF


    CASE DEFAULT
      CALL output_line('Invalid position of midpoint vertex!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria2Tria')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE coarsen_4Tria2Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Quad1Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into the original macro quadrilateral.
    ! The remaining element JEL is labeled by the smallest element number
    ! of the four neighboring elements. Due to the numbering convention
    ! all elements can be determined by visiting neighbors along the
    ! second edge starting at the given element IEL. The resulting element
    ! JEL is oriented such that local numbering starts at the vertex of the
    ! macro element, that is, such that the first vertex is the oldest one.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)      (e3) i3          i4 (e7) i7  (e3) i3
    !      +-------+-------+             +---------------+
    !      |*      |      *|             |               |
    ! (e4) | iel3  |  iel2 | (e6)   (e4) |               | (e6)
    !      |       |       |             |               |
    !    i8+-------+-------+i6    ->     |      jel      |
    !      |     i9|       |             |               |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |               | (e2)
    !      |*      |      *|             |*              |
    !      +-------+-------+             +---------------+
    !     i1 (e1) i5  (e5) i2           i1 (e1)     (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IsortedElements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_VERTEXIDX),  DIMENSION(TRIA_NVEQUAD2D) :: ImacroVertices
    INTEGER, DIMENSION(TRIA_NVEQUAD2D)                  :: IvertexAge
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8,jel,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    INTEGER                    :: istate

    ! Retrieve patch of elements
    iel1=rhadapt%p_IneighboursAtElement(2,iel)
    iel2=rhadapt%p_IneighboursAtElement(2,iel1)
    iel3=rhadapt%p_IneighboursAtElement(2,iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i9=rhadapt%p_IverticesAtElement(3,iel)
    i8=rhadapt%p_IverticesAtElement(4,iel)
    i2=rhadapt%p_IverticesAtElement(1,iel1)
    i6=rhadapt%p_IverticesAtElement(2,iel1)
    i3=rhadapt%p_IverticesAtElement(1,iel2)
    i7=rhadapt%p_IverticesAtElement(2,iel2)
    i4=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e8=rhadapt%p_IneighboursAtElement(4,iel)   
    e2=rhadapt%p_IneighboursAtElement(1,iel1)
    e5=rhadapt%p_IneighboursAtElement(4,iel1)
    e3=rhadapt%p_IneighboursAtElement(1,iel2)
    e6=rhadapt%p_IneighboursAtElement(4,iel2)
    e4=rhadapt%p_IneighboursAtElement(1,iel3)
    e7=rhadapt%p_IneighboursAtElement(4,iel3)


    ! Sort the four elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements=(/iel,iel1,iel2,iel3/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel=IsortedElements(1)


    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel3,jel,jel)


    ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
    ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
    ! should not appear, that is, local numbering of the resulting quadrilateral
    ! starts at the oldest vertex. to this end, we check the state of the 
    ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
    ImacroVertices = (/i1,i2,i3,i4/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate = redgreen_getstateQuad(IvertexAge)

    SELECT CASE(istate)
    CASE(STATE_QUAD_ROOT,STATE_QUAD_RED4)
      ! Update element JEL = (I1,I2,I3,I4)
      CALL replace_element2D(rhadapt,jel,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)

    CASE(STATE_QUAD_RED1)
      ! Update element JEL = (I2,I3,I4,I1)
      CALL replace_element2D(rhadapt,jel,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)

    CASE(STATE_QUAD_RED2)
      ! Update element JEL = (I3,I4,I1,I2)
      CALL replace_element2D(rhadapt,jel,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)

    CASE(STATE_QUAD_RED3)
      ! Update element JEL = (I4,I1,I2,I3)
      CALL replace_element2D(rhadapt,jel,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
      
    CASE DEFAULT
      CALL output_line('Invalid state of resulting quadrilateral!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      CALL sys_halt()
    END SELECT


    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
    ! element corresponds to element with minimum number JEL
    CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))

    CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))

    CALL remove_element2D(rhadapt,IsortedElements(2),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(2))


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad1Quad')
      CALL sys_halt()
    END IF

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel,ipos)


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4QUAD1QUAD,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_4Quad1Quad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Quad2Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into two green quadrilaterals.
    ! The remaining elements JEL1 and JEL2 are labeled by the  two smallest
    ! element numbers of the four neighboring elements. Due to the numbering
    ! convention all elements can be determined by visiting neighbors along 
    ! the second edge starting at the given element IEL. The resulting elements
    ! JEL1 and JEL2 are oriented such that local numbering starts at the vertices
    ! of the macro element, that is, such that the first vertex is the oldest one.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)  i7  (e3) i3
    !      +-------+-------+             +-------+-------+
    !      |*      |      *|             |       |      *|
    ! (e4) | iel3  |  iel2 | (e6)   (e4) |       |       | (e6)
    !      |       |       |             |       |       |
    !    i8+-------+-------+i6    ->     | jel1  | jel2  |
    !      |     i9|       |             |       |       |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |       |       | (e2)
    !      |*      |      *|             |*      |       |
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IsortedElements    
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8,jel1,jel2,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Retrieve patch of elements
    iel1=rhadapt%p_IneighboursAtElement(2,iel)
    iel2=rhadapt%p_IneighboursAtElement(2,iel1)
    iel3=rhadapt%p_IneighboursAtElement(2,iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i9=rhadapt%p_IverticesAtElement(3,iel)
    i8=rhadapt%p_IverticesAtElement(4,iel)
    i2=rhadapt%p_IverticesAtElement(1,iel1)
    i6=rhadapt%p_IverticesAtElement(2,iel1)
    i3=rhadapt%p_IverticesAtElement(1,iel2)
    i7=rhadapt%p_IverticesAtElement(2,iel2)
    i4=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e8=rhadapt%p_IneighboursAtElement(4,iel)   
    e2=rhadapt%p_IneighboursAtElement(1,iel1)
    e5=rhadapt%p_IneighboursAtElement(4,iel1)
    e3=rhadapt%p_IneighboursAtElement(1,iel2)
    e6=rhadapt%p_IneighboursAtElement(4,iel2)
    e4=rhadapt%p_IneighboursAtElement(1,iel3)
    e7=rhadapt%p_IneighboursAtElement(4,iel3)


    ! Sort the four elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements=(/iel,iel1,iel2,iel3/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel1=IsortedElements(1)
    jel2=IsortedElements(2)


    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,jel2,jel1)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel2,jel2)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel1,jel2)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel3,jel1,jel1)


    ! Update elements JEL1 = (I1,I5,I7,I4) and JEL2=(I3,I7,I5,I2)
    CALL replace_element2D(rhadapt,jel1,i1,i5,i7,i4,e1,jel2,e7,e4,e1,jel2,e7,e8)
    CALL replace_element2D(rhadapt,jel2,i3,i7,i5,i2,e3,jel1,e5,e2,e3,jel1,e5,e6)

    
    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
    ! elements corresponds to elements with minimum numbers JEL1 and JEL2
    CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))
    
    CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad2Quad')
      CALL sys_halt()
    END IF
   
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel2,ipos)


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4QUAD2QUAD,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_4Quad2Quad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Quad3Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into three green triangles.
    ! The elements remaining JEL1,JEL2, and JEL3 are the ones with
    ! the smallest element number, that is, only the largest element is
    ! removed. The numbering convention follows the strategy used for
    ! refinement, that is, the local numbering of the two outer triangles
    ! starts at the vertices of the macro element and the local numbering
    ! numbering of the interior triangle starts at the edge midpoint.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*      |      *|             |\             /|
    ! (e4) | iel3  |  iel2 | (e6)   (e4) | \           / | (e6)
    !      |       |       |             |  \         /  |
    !    i8+-------+-------+i6    ->     |   \  jel3 /   |
    !      |     i9|       |             |    \     /    |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |jel1 \   / jel2| (e2)
    !      |*      |      *|             |*     \*/     *|
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IsortedElements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_ELEMENTIDX)   :: jel1,jel2,jel3,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Retrieve patch of elements
    iel1=rhadapt%p_IneighboursAtElement(2,iel)
    iel2=rhadapt%p_IneighboursAtElement(2,iel1)
    iel3=rhadapt%p_IneighboursAtElement(2,iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i9=rhadapt%p_IverticesAtElement(3,iel)
    i8=rhadapt%p_IverticesAtElement(4,iel)
    i2=rhadapt%p_IverticesAtElement(1,iel1)
    i6=rhadapt%p_IverticesAtElement(2,iel1)
    i3=rhadapt%p_IverticesAtElement(1,iel2)
    i7=rhadapt%p_IverticesAtElement(2,iel2)
    i4=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e8=rhadapt%p_IneighboursAtElement(4,iel)   
    e2=rhadapt%p_IneighboursAtElement(1,iel1)
    e5=rhadapt%p_IneighboursAtElement(4,iel1)
    e3=rhadapt%p_IneighboursAtElement(1,iel2)
    e6=rhadapt%p_IneighboursAtElement(4,iel2)
    e4=rhadapt%p_IneighboursAtElement(1,iel3)
    e7=rhadapt%p_IneighboursAtElement(4,iel3)


    ! Sort the four elements according to their number and
    ! determine the elements with the smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel1=IsortedElements(1)
    jel2=IsortedElements(2)
    jel3=IsortedElements(3)


    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,jel2,jel1)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel2,jel2)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel3,jel3)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel3,jel1,jel1)


    ! Update elements JEL1 = (I1,I5,I4), JEL2 = (I2,I3,I5), and JEL3 = (I5,I3,I4)
    CALL replace_element2D(rhadapt,jel1,i1,i5,i4,e1,jel3,e4,e1,jel3,e8)
    CALL replace_element2D(rhadapt,jel2,i2,i3,i5,e2,jel3,e5,e6,jel3,e5)
    CALL replace_element2D(rhadapt,jel3,i5,i3,i4,jel2,e3,jel1,jel2,e7,jel1)


    ! Delete the element with the largest element number
    CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i5,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Quad3Tria')
      CALL sys_halt()
    END IF

    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,jel3,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)+3
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)-3
    
    
    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4QUAD3TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_4Quad3Tria

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Quad4Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four quadrilaterals resulting from a
    ! 1-quad : 4-quad refinement into four green triangles. At first glance
    ! this "coarsening" does not make sense since the number of elements
    ! is not reduced. However, the number of vertices is decreased by one.
    ! Moreover, the removal of this vertex may lead to further coarsening
    ! starting at the elements which meet at the deleted vertex.
    ! The numbering convention follows the strategy used for
    ! refinement, that is, the local numbering of the two outer triangles
    ! starts at the vertices of the macro element and the local numbering
    ! numbering of the interior triangle starts at the edge midpoint.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*      |      *|             |*          ---/|
    ! (e4) | iel3  |  iel2 | (e6)   (e4) |iel3   ---/ */ | (e6)
    !      |       |       |             |   ---/     /  |
    !    i8+-------+-------+i6    ->   i8+--/  iel2  /   |
    !      |     i9|       |             |-\        /    |
    ! (e8) | iel   |  iel1 | (e2)   (e8) |  --\    /     | (e2)
    !      |*      |      *|             |*iel -\ / iel1*|
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(9) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i6,i7,i8,i9
    
    ! Retrieve patch of elements
    iel1=rhadapt%p_IneighboursAtElement(2,iel)
    iel2=rhadapt%p_IneighboursAtElement(2,iel1)
    iel3=rhadapt%p_IneighboursAtElement(2,iel2)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i9=rhadapt%p_IverticesAtElement(3,iel)
    i8=rhadapt%p_IverticesAtElement(4,iel)
    i2=rhadapt%p_IverticesAtElement(1,iel1)
    i6=rhadapt%p_IverticesAtElement(2,iel1)
    i3=rhadapt%p_IverticesAtElement(1,iel2)
    i7=rhadapt%p_IverticesAtElement(2,iel2)
    i4=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e8=rhadapt%p_IneighboursAtElement(4,iel)   
    e2=rhadapt%p_IneighboursAtElement(1,iel1)
    e5=rhadapt%p_IneighboursAtElement(4,iel1)
    e3=rhadapt%p_IneighboursAtElement(1,iel2)
    e6=rhadapt%p_IneighboursAtElement(4,iel2)
    e4=rhadapt%p_IneighboursAtElement(1,iel3)
    e7=rhadapt%p_IneighboursAtElement(4,iel3)


    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,iel1,iel1)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,iel3,iel3)

    
    ! Update the four element
    CALL replace_element2D(rhadapt,iel,i1,i5,i8,e1,iel2,e8,e1,iel2,e8)
    CALL replace_element2D(rhadapt,iel1,i2,i3,i5,e2,iel2,e5,e6,iel2,e5)
    CALL replace_element2D(rhadapt,iel2,i3,i8,i5,iel3,iel,iel1,iel3,iel,iel1)
    CALL replace_element2D(rhadapt,iel3,i4,i8,i3,e4,iel2,e3,e4,iel2,e7)


    ! Update list of elements meeting at vertices. Note that we only have to
    ! add some elements to the vertices since all four elements are already
    ! "attached" to one of the four corner nodes
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel3,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,iel2,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i8,iel2,ipos)
    
    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)+4
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)-4
    
    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,i6,i7,i8,i9/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4QUAD4TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_4Quad4Tria
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_2Quad1Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines two quadrilaterals resulting from a
    ! 1-quad : 2-quad refinement into the macro quadrilateral.
    ! The remaining element JEL is labeled by the smallest element number
    ! of the four neighboring elements. Due to the numbering convention
    ! all elements can be determined by visiting neighbors along the
    ! second edge starting at the given element IEL. The resulting element
    ! JEL is oriented such that local numbering starts at the vertex of the
    ! macro element, that is, such that the first vertex is the oldest one.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |       |      *|             |               |
    ! (e4) |       |       | (e6)   (e4) |               | (e6)
    !      |       |       |             |               |
    !      |  iel  |  iel1 |      ->     |      jel      |
    !      |       |       |             |               |
    ! (e8) |       |       | (e2)   (e8) |               | (e2)
    !      |*      |      *|             |*              |
    !      +-------+-------+             +---------------+
    !     i1 (e1) i5  (e5) i2           i1 (e1)     (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(8) :: Ivertices
    INTEGER(PREC_VERTEXIDX),  DIMENSION(TRIA_NVEQUAD2D) :: ImacroVertices
    INTEGER, DIMENSION(TRIA_NVEQUAD2D)                  :: IvertexAge
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,e1,e2,e3,e4,e5,e6,e7,e8,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i7
    INTEGER                    :: istate
    
    ! Retrieve neighboring element
    iel1=rhadapt%p_IneighboursAtElement(2,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i7=rhadapt%p_IverticesAtElement(3,iel)
    i4=rhadapt%p_IverticesAtElement(4,iel)
    i3=rhadapt%p_IverticesAtElement(1,iel1)
    i2=rhadapt%p_IverticesAtElement(4,iel1)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e7=rhadapt%p_IneighboursAtElement(3,iel)
    e4=rhadapt%p_IneighboursAtElement(4,iel)    
    e3=rhadapt%p_IneighboursAtElement(1,iel1)
    e5=rhadapt%p_IneighboursAtElement(3,iel1)
    e2=rhadapt%p_IneighboursAtElement(4,iel1)

    e8=rhadapt%p_ImidneighboursAtElement(4,iel)
    e6=rhadapt%p_ImidneighboursAtElement(4,iel1)


    ! Which is the smaller element?
    IF (iel < iel1) THEN
      
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,iel,iel)
      CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel1,iel,iel)
      CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel1,iel,iel)


      ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
      ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
      ! should not appear, that is, local numbering of the resulting quadrilateral
      ! starts at the oldest vertex. to this end, we check the state of the 
      ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
      ImacroVertices = (/i1,i2,i3,i4/)
      IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
      istate = redgreen_getstateQuad(IvertexAge)
      
      SELECT CASE(istate)
      CASE(STATE_QUAD_ROOT,STATE_QUAD_RED4)
        ! Update element IEL = (I1,I2,I3,I4)
        CALL replace_element2D(rhadapt,iel,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)
        
      CASE(STATE_QUAD_RED1)
        ! Update element IEL = (I2,I3,I4,I1)
        CALL replace_element2D(rhadapt,iel,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)

      CASE(STATE_QUAD_RED2)
        ! Update element IEL = (I3,I4,I1,I2)
        CALL replace_element2D(rhadapt,iel,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)

      CASE(STATE_QUAD_RED3)
        ! Update element IEL = (I4,I1,I2,I3)
        CALL replace_element2D(rhadapt,iel,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
        
      CASE DEFAULT
        CALL output_line('Invalid state of resulting quadrilateral!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        CALL sys_halt()
      END SELECT


      ! Delete element IEL1
      CALL remove_element2D(rhadapt,iel1,ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,iel1)


      ! Update list of elements meeting at vertices.
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,iel,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,iel,ipos)

    ELSE

      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,iel,iel1,iel1)
      CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel1,iel1,iel1)
      CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel,iel1,iel1)
      
      ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
      ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
      ! should not appear, that is, local numbering of the resulting quadrilateral
      ! starts at the oldest vertex. to this end, we check the state of the 
      ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
      ImacroVertices = (/i1,i2,i3,i4/)
      IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
      istate = redgreen_getstateQuad(IvertexAge)
      
      SELECT CASE(istate)
      CASE(STATE_QUAD_ROOT,STATE_QUAD_RED4)
        ! Update element IEL1 = (I1,I2,I3,I4)
        CALL replace_element2D(rhadapt,iel1,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)
        
      CASE(STATE_QUAD_RED1)
        ! Update element IEL1 = (I2,I3,I4,I1)
        CALL replace_element2D(rhadapt,iel1,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)

      CASE(STATE_QUAD_RED2)
        ! Update element IEL1 = (I3,I4,I1,I2)
        CALL replace_element2D(rhadapt,iel1,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)

      CASE(STATE_QUAD_RED3)
        ! Update element IEL1 = (I4,I1,I2,I3)
        CALL replace_element2D(rhadapt,iel1,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
        
      CASE DEFAULT
        CALL output_line('Invalid state of resulting quadrilateral!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        CALL sys_halt()
      END SELECT


      ! Delete element IEL
      CALL remove_element2D(rhadapt,iel,ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,iel)
      

      ! Update list of elements meeting at vertices.
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_2Quad1Quad')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,iel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,iel1,ipos)

    END IF


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,0,i7,0/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_2QUAD1QUAD,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_2Quad1Quad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_2Quad3Tria(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines two quadrilaterals resulting from a
    ! 1-quad : 2-quad refinement into three green triangles. At first glance
    ! this "coarsening" does not make sense since the number of elements
    ! is even increased. However, the number of vertices is decreased by one.
    ! Moreover, the removal of this vertex may lead to further coarsening
    ! starting at the elements which meet at the deleted vertex.
    ! The numbering convention follows the strategy used for
    ! refinement, that is, the local numbering of the two outer triangles
    ! starts at the vertices of the macro element and the local numbering
    ! numbering of the interior triangle starts at the edge midpoint.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7) i7  (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |       |      *|             |\             /|
    ! (e4) |       |       | (e6)   (e4) | \           / | (e6)
    !      |       |       |             |  \  nel+1  /  |
    !      |  iel  |  iel1 |      ->     |   \       /   |
    !      |       |       |             |    \     /    |
    ! (e8) |       |       | (e2)   (e8) | iel \   / iel1| (e2)
    !      |*      |      *|             |*     \*/     *|
    !      +-------+-------+             +-------+-------+
    !     i1 (e1) i5  (e5) i2           i1 (e1) i5  (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(8) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,nel0,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5,i7

    ! Retrieve neighboring element
    iel1=rhadapt%p_IneighboursAtElement(2,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i5=rhadapt%p_IverticesAtElement(2,iel)
    i7=rhadapt%p_IverticesAtElement(3,iel)
    i4=rhadapt%p_IverticesAtElement(4,iel)
    i3=rhadapt%p_IverticesAtElement(1,iel1)
    i2=rhadapt%p_IverticesAtElement(4,iel1)

    ! Store values of the elements adjacent to the resulting macro element
    e1=rhadapt%p_IneighboursAtElement(1,iel)
    e7=rhadapt%p_IneighboursAtElement(3,iel)
    e4=rhadapt%p_IneighboursAtElement(4,iel)
    e3=rhadapt%p_IneighboursAtElement(1,iel1)
    e5=rhadapt%p_IneighboursAtElement(3,iel1)
    e2=rhadapt%p_IneighboursAtElement(4,iel1)

    e8=rhadapt%p_ImidneighboursAtElement(4,iel)
    e6=rhadapt%p_ImidneighboursAtElement(4,iel1)

    ! Store total number of elements before refinement
    nel0=rhadapt%NEL
    

    ! Replace elements IEL and IEL1 and add one new element JEL
    CALL replace_element2D(rhadapt,iel,i1,i5,i4,e1,nel0+1,e4,e1,nel0+1,e8)
    CALL replace_element2D(rhadapt,iel1,i2,i3,i5,e2,nel0+1,e5,e6,nel0+1,e5)
    CALL add_element2D(rhadapt,i5,i3,i4,iel1,e3,iel,iel1,e7,iel)


    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,iel1,nel0+1,nel0+1)


    ! Update list of elements meeting at vertices.
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,nel0+1,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i5,nel0+1,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)+2
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)-2

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5,0,i7,0/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_2QUAD3TRIA,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_2Quad3Tria
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_3Tria1Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines three triangles resulting from a
    ! 1-quad : 3-tria refinement into the macro quadrilateral.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)     (e3) i3          i4 (e7)      (e3) i3
    !      +---------------+             +---------------+
    !      |\             /|             |               |
    ! (e4) | \           / | (e6)   (e4) |               | (e6)
    !      |  \   iel   /  |             |               |
    !      |   \       /   |      ->     |      jel      |
    !      |    \     /    |             |               |
    ! (e8) | iel1\   /iel2 | (e2)   (e8) |               | (e2)
    !      |*     \*/     *|             |*              |
    !      +-------+-------+             +---------------+
    !     i1 (e1) i5  (e5) i2           i1 (e1)     (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(3) :: IsortedElements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(5) :: Ivertices
    INTEGER(PREC_VERTEXIDX),  DIMENSION(TRIA_NVEQUAD2D) :: ImacroVertices
    INTEGER, DIMENSION(TRIA_NVEQUAD2D)                  :: IvertexAge
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,e1,e2,e3,e4,e5,e6,e7,e8,jel,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i5
    INTEGER                    :: istate

    ! Retrieve neighboring element
    iel2=rhadapt%p_IneighboursAtElement(1,iel)
    iel1=rhadapt%p_IneighboursAtElement(3,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i5=rhadapt%p_IverticesAtElement(1,iel)
    i3=rhadapt%p_IverticesAtElement(2,iel)
    i4=rhadapt%p_IverticesAtElement(3,iel)
    i1=rhadapt%p_IverticesAtElement(1,iel1)
    i2=rhadapt%p_IverticesAtElement(1,iel2)

    ! Store values of the elements adjacent to the resulting macro element
    e3=rhadapt%p_IneighboursAtElement(2,iel)
    e1=rhadapt%p_IneighboursAtElement(1,iel1)
    e4=rhadapt%p_IneighboursAtElement(3,iel1)
    e2=rhadapt%p_IneighboursAtElement(1,iel2)
    e5=rhadapt%p_IneighboursAtElement(3,iel2)

    e7=rhadapt%p_ImidneighboursAtElement(2,iel)    
    e8=rhadapt%p_ImidneighboursAtElement(3,iel1)
    e6=rhadapt%p_ImidneighboursAtElement(1,iel2)


    ! Sort the three elements according to their number and
    ! determine the element with the smallest element number
    IsortedElements=(/iel,iel1,iel2/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel=IsortedElements(1)

    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel2,iel1,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel1,jel,jel)

    ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
    ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
    ! should not appear, that is, local numbering of the resulting quadrilateral
    ! starts at the oldest vertex. to this end, we check the state of the 
    ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
    ImacroVertices = (/i1,i2,i3,i4/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate = redgreen_getstateQuad(IvertexAge)

    SELECT CASE(istate)
    CASE(STATE_QUAD_ROOT,STATE_QUAD_RED4)
      ! Update element JEL = (I1,I2,I3,I4)
      CALL replace_element2D(rhadapt,jel,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)
      
    CASE(STATE_QUAD_RED1)
      ! Update element JEL = (I2,I3,I4,I1)
      CALL replace_element2D(rhadapt,jel,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)
      
    CASE(STATE_QUAD_RED2)
      ! Update element JEL = (I3,I4,I1,I2)
      CALL replace_element2D(rhadapt,jel,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)
      
    CASE(STATE_QUAD_RED3)
      ! Update element JEL = (I4,I1,I2,I3)
      CALL replace_element2D(rhadapt,jel,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
      
    CASE DEFAULT
      CALL output_line('Invalid state of resulting quadrilateral!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      CALL sys_halt()
    END SELECT


    ! Delete the two elements with the largest element numbers
    CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))

    CALL remove_element2D(rhadapt,IsortedElements(2),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(2))


    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_3Tria1Quad')
      CALL sys_halt()
    END IF
    
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel,ipos)      
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel,ipos)


    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)-1
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)+1


    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,i5/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_3TRIA1QUAD,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_3Tria1Quad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Tria1Quad(rhadapt,iel,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-quad : 4-tria refinement into the macro quadrilateral.
    ! The remaining element JEL is the ones with the smallest element number.
    ! The numbering convention follows the strategy used for refinement.
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)  i7 (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*     / \-iel2*|             |               |
    ! (e4) |iel3 /    \--  | (e6)   (e4) |               | (e6)
    !      |    /        \-|             |               |
    !      |   /  iel    --+i6    ->     |      jel      |
    !      |  /      ---/  |             |               |
    ! (e8) | /*  ---/      | (e2)   (e8) |               | (e2)
    !      |/---/    iel1 *|             |*              |
    !      +---------------+             +---------------+
    !     i1 (e1)     (e5) i2           i1 (e1)     (e5) i2
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IsortedElements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(8) :: Ivertices
    INTEGER(PREC_VERTEXIDX),  DIMENSION(TRIA_NVEQUAD2D) :: ImacroVertices
    INTEGER, DIMENSION(TRIA_NVEQUAD2D)                  :: IvertexAge
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8,jel,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i6,i7
    INTEGER                    :: istate

    ! Retrieve patch of elements
    iel1=rhadapt%p_IneighboursAtElement(1,iel)
    iel2=rhadapt%p_IneighboursAtElement(2,iel)
    iel3=rhadapt%p_IneighboursAtElement(3,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i6=rhadapt%p_IverticesAtElement(2,iel)
    i7=rhadapt%p_IverticesAtElement(3,iel)
    i2=rhadapt%p_IverticesAtElement(1,iel1)
    i3=rhadapt%p_IverticesAtElement(1,iel2)
    i4=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e2=rhadapt%p_IneighboursAtElement(1,iel1)  
    e1=rhadapt%p_IneighboursAtElement(3,iel1)
    e3=rhadapt%p_IneighboursAtElement(1,iel2)
    e6=rhadapt%p_IneighboursAtElement(3,iel2)
    e4=rhadapt%p_IneighboursAtElement(1,iel3)
    e7=rhadapt%p_IneighboursAtElement(3,iel3)

    e5=rhadapt%p_ImidneighboursAtElement(3,iel1)
    e8=rhadapt%p_ImidneighboursAtElement(1,iel3)

    
    ! Sort the four elements according to their number and
    ! determine the elements with the smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel=IsortedElements(1)

    
    ! Update list of neighboring elements
    CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel,jel)
    CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel3,jel,jel)


    ! The resulting quadrilateral will posses one of the states STATES_QUAD_REDx,
    ! whereby x can be 1,2,3 and 4. Die to our refinement convention, x=1,2,3
    ! should not appear, that is, local numbering of the resulting quadrilateral
    ! starts at the oldest vertex. to this end, we check the state of the 
    ! provisional quadrilateral (I1,I2,I3,I4) and transform the orientation.
    ImacroVertices = (/i1,i2,i3,i4/)
    IvertexAge = rhadapt%p_IvertexAge(ImacroVertices)
    istate = redgreen_getstateQuad(IvertexAge)

    SELECT CASE(istate)
    CASE(STATE_QUAD_ROOT,STATE_QUAD_RED4)
      ! Update element JEL = (I1,I2,I3,I4)
      CALL replace_element2D(rhadapt,jel,i1,i2,i3,i4,e1,e2,e3,e4,e5,e6,e7,e8)

    CASE(STATE_QUAD_RED1)
      ! Update element JEL = (I2,I3,I4,I1)
      CALL replace_element2D(rhadapt,jel,i2,i3,i4,i1,e2,e3,e4,e1,e6,e7,e8,e5)

    CASE(STATE_QUAD_RED2)
      ! Update element JEL = (I3,I4,I1,I2)
      CALL replace_element2D(rhadapt,jel,i3,i4,i1,i2,e3,e4,e1,e2,e7,e8,e5,e6)

    CASE(STATE_QUAD_RED3)
      ! Update element JEL = (I4,I1,I2,I3)
      CALL replace_element2D(rhadapt,jel,i4,i1,i2,i3,e4,e1,e2,e3,e8,e5,e6,e7)
      
    CASE DEFAULT
      CALL output_line('Invalid state of resulting quadrilateral!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      CALL sys_halt()
    END SELECT

    
    ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
    ! element corresponds to element with minimum number JEL
    CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))

    CALL remove_element2D(rhadapt,IsortedElements(3),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(3))
    
    CALL remove_element2D(rhadapt,IsortedElements(2),ielReplace)
    IF (ielReplace.NE.0)&
        CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(2))

    
    ! Update list of elements meeting at vertices.
    ! Note that all elements are removed in the first step. Afterwards,
    ! element JEL is appended to the list of elements meeting at each vertex
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      CALL sys_halt()
    END IF
    IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
        ARRAYLIST_NOT_FOUND) THEN
      CALL output_line('Unable to delete element from vertex list!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria1Quad')
      CALL sys_halt()
    END IF
        
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel,ipos)
    CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel,ipos)
    
    
    ! Adjust number of elements
    rhadapt%InelOfType(TRIA_NVETRI2D)=rhadapt%InelOfType(TRIA_NVETRI2D)-1
    rhadapt%InelOfType(TRIA_NVEQUAD2D)=rhadapt%InelOfType(TRIA_NVEQUAD2D)+1
    

    ! Optionally, invoke callback routine
    IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
      Ivertices=(/i1,i2,i3,i4,0,i6,i7,0/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
      CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA1QUAD,&
          Ivertices,Ielements)
    END IF
  END SUBROUTINE coarsen_4Tria1Quad

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coarsen_4Tria3Tria(rhadapt,iel,imarker,rcollection,fcb_hadaptCallback)

!<description>
    ! This subroutine combines four triangles resulting from a
    ! 1-quad : 4-tria refinement into three green triangles. The remaining 
    ! elements JEL1, JEL2 and JEL3 are the ones with the smallest element number.
    ! The numbering convention follows the strategy used for refinement. 
    !
    !    initial quadrilateral      subdivided quadrilateral
    !
    !     i4 (e7)  i7 (e3) i3          i4 (e7)      (e3) i3
    !      +-------+-------+             +---------------+
    !      |*     / \-iel2*|             |---\          *|
    ! (e4) |iel3 /    \--  | (e6)   (e4) |    ---\ jel2  | (e6)
    !      |    /        \-|             |        ---\   |
    !      |   /  iel    --+i6    ->     |  jel3     *-- +i6
    !      |  /      ---/  |             |        ---/   |
    ! (e8) | /*  ---/      | (e2)   (e8) |    ---/       | (e2)
    !      |/---/    iel1 *|             |---/     jel1 *|
    !      +---------------+             +---------------+
    !     i1 (e1)     (e5) i2           i1 (e1)     (e5) i2
!</description>

!<input>
    ! Number of element to be refined
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel

    ! Identifiert for element marker
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
    INTEGER(PREC_ELEMENTIDX), DIMENSION(8) :: Ielements
    INTEGER(PREC_ELEMENTIDX), DIMENSION(4) :: IsortedElements
    INTEGER(PREC_VERTEXIDX),  DIMENSION(8) :: Ivertices
    INTEGER(PREC_ARRAYLISTIDX) :: ipos
    INTEGER(PREC_ELEMENTIDX)   :: iel1,iel2,iel3,e1,e2,e3,e4,e5,e6,e7,e8
    INTEGER(PREC_ELEMENTIDX)   :: jel1,jel2,jel3,ielReplace
    INTEGER(PREC_VERTEXIDX)    :: i1,i2,i3,i4,i6,i7

    ! Retrieve patch of elements
    iel1=rhadapt%p_IneighboursAtElement(1,iel)
    iel2=rhadapt%p_IneighboursAtElement(2,iel)
    iel3=rhadapt%p_IneighboursAtElement(3,iel)

    ! Store vertex- and element-values of the four neighboring elements
    i1=rhadapt%p_IverticesAtElement(1,iel)
    i6=rhadapt%p_IverticesAtElement(2,iel)
    i7=rhadapt%p_IverticesAtElement(3,iel)
    i2=rhadapt%p_IverticesAtElement(1,iel1)
    i3=rhadapt%p_IverticesAtElement(1,iel2)
    i4=rhadapt%p_IverticesAtElement(1,iel3)

    ! Store values of the elements adjacent to the resulting macro element
    e2=rhadapt%p_IneighboursAtElement(1,iel1)  
    e1=rhadapt%p_IneighboursAtElement(3,iel1)
    e3=rhadapt%p_IneighboursAtElement(1,iel2)
    e6=rhadapt%p_IneighboursAtElement(3,iel2)
    e4=rhadapt%p_IneighboursAtElement(1,iel3)
    e7=rhadapt%p_IneighboursAtElement(3,iel3)

    e5=rhadapt%p_ImidneighboursAtElement(3,iel1)
    e8=rhadapt%p_ImidneighboursAtElement(1,iel3)


    ! Sort the four elements according to their number and
    ! determine the elements with the smallest element numbers
    IsortedElements=(/iel,iel1,iel2,iel3/)
    CALL sort_I32(IsortedElements,SORT_INSERT)
    jel1=IsortedElements(1)
    jel2=IsortedElements(2)
    jel3=IsortedElements(3)
    
    ! Which midpoint vertex should be kept?
    SELECT CASE(imarker)
    CASE(MARK_CRS_4TRIA3TRIA_RIGHT)
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,jel1,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel2,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel2,jel2)
      CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel3,jel3,jel3)


      ! Update elements JEL1, JEL2 and JEL3
      CALL replace_element2D(rhadapt,jel1,i2,i6,i1,e2,jel3,e1,e2,jel3,e5)
      CALL replace_element2D(rhadapt,jel2,i3,i4,i6,e3,jel3,e6,e7,jel3,e6)
      CALL replace_element2D(rhadapt,jel3,i6,i4,i1,jel2,e4,jel1,jel2,e8,jel1)

      
      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
      ! elements correspond to the two elements with smallest numbers
      CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))


      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i6,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel3,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel3,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i6,jel3,ipos)

      
      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/i1,i2,i3,i4,0,i6,i7,0/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA3TRIA2,&
            Ivertices,Ielements)
      END IF
      

    CASE(MARK_CRS_4TRIA3TRIA_LEFT)
      ! Update list of neighboring elements
      CALL update_ElementNeighbors2D(rhadapt,e1,e5,iel1,jel3,jel3)
      CALL update_ElementNeighbors2D(rhadapt,e2,e6,iel2,iel1,jel1,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e3,e7,iel3,iel2,jel2,jel1)
      CALL update_ElementNeighbors2D(rhadapt,e4,e8,iel3,jel2,jel2)

      
      ! Update elements JEL1, JEL2 and JEL3
      CALL replace_element2D(rhadapt,jel1,i3,i7,i2,e3,jel3,e2,e3,jel3,e6)
      CALL replace_element2D(rhadapt,jel2,i4,i1,i7,e4,jel3,e7,e8,jel3,e7)
      CALL replace_element2D(rhadapt,jel3,i7,i1,i2,jel2,e1,jel1,jel2,e5,jel1)

      
      ! Delete elements IEL, IEL1, IEL2 and IEL3 depending on which 
      ! elements correspond to the two elements with smallest numbers
      CALL remove_element2D(rhadapt,IsortedElements(4),ielReplace)
      IF (ielReplace.NE.0)&
          CALL update_AllElementNeighbors2D(rhadapt,ielReplace,IsortedElements(4))


      ! Update list of elements meeting at vertices.
      ! Note that all elements are removed in the first step. Afterwards,
      ! element JEL is appended to the list of elements meeting at each vertex
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i1,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i2,iel1).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i3,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i4,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel2).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF
      IF (arrlst_deleteFromArraylist(rhadapt%relementsAtVertex,i7,iel3).EQ.&
          ARRAYLIST_NOT_FOUND) THEN
        CALL output_line('Unable to delete element from vertex list!',&
            OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
        CALL sys_halt()
      END IF

      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i1,jel3,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i2,jel3,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i3,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i4,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel1,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel2,ipos)
      CALL arrlst_appendToArraylist(rhadapt%relementsAtVertex,i7,jel3,ipos)

      ! Optionally, invoke callback routine
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/i1,i2,i3,i4,0,i6,i7,0/); Ielements=(/e1,e2,e3,e4,e5,e6,e7,e8/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_CRS_4TRIA3TRIA3,&
            Ivertices,Ielements)
      END IF

      
    CASE DEFAULT
      CALL output_line('Invalid position of midpoint vertex!',&
          OU_CLASS_ERROR,OU_MODE_STD,'coarsen_4Tria3Tria')
      CALL sys_halt()
    END SELECT
  END SUBROUTINE coarsen_4Tria3Tria

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
    INTEGER, DIMENSION(:),  POINTER    :: p_Imarker
    INTEGER(PREC_VERTEXIDX),  DIMENSION(TRIA_MAXNVE2D) :: p_IverticesAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE2D) :: p_IneighboursAtElement
    INTEGER(PREC_VERTEXIDX)  :: ivt
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER :: ive,nve

    ! Check if dynamic data structures are generated and contain data
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA) THEN
      CALL output_line('Dynamic data structures are not generated!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mark_refinement2D')
      CALL sys_halt()
    END IF

    ! Initialize marker structure for NEL0 elements
    IF (rhadapt%h_Imarker.NE.ST_NOHANDLE) &
        CALL storage_free(rhadapt%h_Imarker)
    CALL storage_new('mark_refinement2D','Imarker',rhadapt%NEL0,&
        ST_INT,rhadapt%h_Imarker,ST_NEWBLOCK_ZERO)
       
    ! Set pointers
    CALL lsyssc_getbase_double(rindicator,p_Dindicator)
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)

    ! Set state of all vertices to "free". Note that vertices of the
    ! initial triangulation are always "locked", i.e. have no positive age.
    DO ivt=1,SIZE(rhadapt%p_IvertexAge)
      rhadapt%p_IvertexAge(ivt)=ABS(rhadapt%p_IvertexAge(ivt))
    END DO

    ! Loop over all elements and mark those for which
    ! the indicator is greater than the prescribed treshold
    mark: DO iel=1,SIZE(p_Dindicator)

      IF (p_Dindicator(iel) .GT. rhadapt%drefinementTolerance) THEN
        
        ! Get number of vertices per elements
        nve=get_NVE(rhadapt,iel)

        ! Mark element for refinement
        p_IverticesAtElement(1:nve)=rhadapt%p_IverticesAtElement(1:nve,iel)
        p_IneighboursAtElement(1:nve)=rhadapt%p_IneighboursAtElement(1:nve,iel)

        ! An element can only be refined, if all of its vertices do
        ! not exceed the number of admissible subdivision steps.
        ! So, check the element type and loop over the 3 or 4 vertices.
        SELECT CASE(nve)
        CASE(TRIA_NVETRI2D)
          ! If triangle has reached maximum number of refinement levels,
          ! then enforce no further refinement of this element
          IF (ANY(ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(1:TRIA_NVETRI2D))).EQ.&
              rhadapt%NSUBDIVIDEMAX)) THEN
            p_Imarker(iel)=MARK_ASIS
            
            ! According to the indicator, this element should be refined. Since the 
            ! maximum admissible refinement level has been reached no refinement 
            ! was performed. At the same time, all vertices of the element should
            ! be "locked" to prevent this element from coarsening
            DO ive=1,TRIA_NVETRI2D
              rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                  -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive))) 
            END DO
            
            CYCLE mark
          END IF
          
          ! Otherwise, we can mark the triangle for refinement
          p_Imarker(iel)=MARK_REF_TRIA4TRIA
          
          ! Moreover, we have to "lock" its vertices from recoarsening
          DO ive=1,TRIA_NVETRI2D
            rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
          END DO
          
          ! Update number of new vertices. In principle, we can increase the number of 
          ! new vertices by 3, i.e., one for each each. However, some care must be taken
          ! if the edge belongs to some adjacent element which has been marked for
          ! refinement previously. Hence, a new vertex is only added, if the edge
          ! is connected to the boundary or if the adjacent element has not been marked.
          DO ive=1,TRIA_NVETRI2D
            IF (p_IneighboursAtElement(ive).EQ.0) THEN
              rhadapt%increaseNVT=rhadapt%increaseNVT+1
            ELSEIF((p_Imarker(p_IneighboursAtElement(ive)).NE.MARK_REF_TRIA4TRIA) .AND.&
                   (p_Imarker(p_IneighboursAtElement(ive)).NE.MARK_REF_QUAD4QUAD)) THEN
              rhadapt%increaseNVT=rhadapt%increaseNVT+1
            END IF
          END DO

        CASE(TRIA_NVEQUAD2D)
          ! If quadrilateral has reached maximum number of refinement levels,
          ! then enforce no further refinement of this element
          IF (ANY(ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(1:TRIA_NVEQUAD2D))).EQ.&
              rhadapt%NSUBDIVIDEMAX)) THEN
            p_Imarker(iel)=MARK_ASIS

            ! According to the indicator, this element should be refined. Since the 
            ! maximum admissible refinement level has been reached no refinement 
            ! was performed. At the same time, all vertices of the element should
            ! be "locked" to prevent this element from coarsening
            DO ive=1,TRIA_NVEQUAD2D
              rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                  -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive))) 
            END DO

            CYCLE mark
          END IF
          
          ! Otherwise, we can mark the quadrilateral for refinement
          p_Imarker(iel)=MARK_REF_QUAD4QUAD
          
          ! Moreover, we have to "lock" its vertices from recoarsening
          DO ive=1,TRIA_NVEQUAD2D
            rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
          END DO
          
          ! Update number of new vertices
          DO ive=1,TRIA_NVEQUAD2D
            IF (p_IneighboursAtElement(ive).EQ.0) THEN
              rhadapt%increaseNVT=rhadapt%increaseNVT+1
            ELSEIF((p_Imarker(p_IneighboursAtElement(ive)).NE.MARK_REF_TRIA4TRIA) .AND.&
                   (p_Imarker(p_IneighboursAtElement(ive)).NE.MARK_REF_QUAD4QUAD)) THEN
              rhadapt%increaseNVT=rhadapt%increaseNVT+1
            END IF
          END DO

!!$          ! And don't forget the new vertex in the interior of the quadrilateral
!!$          rhadapt%increaseNVT=rhadapt%increaseNVT+1

        CASE DEFAULT
          CALL output_line('Invalid element type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mark_refinement2D')
          CALL sys_halt()
        END SELECT
        
      ELSE
        ! Unmark element for refinement
        p_Imarker(iel)=MARK_ASIS
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
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D)  :: p_IverticesAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE2D) :: p_ImacroElement
    INTEGER, DIMENSION(TRIA_MAXNVE)                    :: IvertexAge
    INTEGER(PREC_VERTEXIDX)  :: i,j
    INTEGER(PREC_ELEMENTIDX) :: iel,jel,kel
    INTEGER :: ive,nve,istate,jstate,ivertexLock
    LOGICAL :: isModified

    ! Check if dynamic data structures are generated and contain data
    IF (IAND(rhadapt%iSpec,HADAPT_HAS_DYNAMICDATA).NE.HADAPT_HAS_DYNAMICDATA .OR.&
        IAND(rhadapt%iSpec,HADAPT_MARKEDREFINE).NE.HADAPT_MARKEDREFINE) THEN
      CALL output_line('Dynamic data structures are not &
          & generated or no marker for grid refinement exists!',&
          OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
      CALL sys_halt()
    END IF
    
    ! Set pointers
    IF (rhadapt%h_Imarker .EQ. ST_NOHANDLE) THEN
      CALL output_line('Marker array is not available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
      CALL sys_halt()
    END IF
    CALL storage_getbase_int(rhadapt%h_Imarker,p_Imarker)
    CALL lsyssc_getbase_double(rindicator,p_Dindicator)

    ! All nodes of the initial triangulation have age-0 and, hence, will
    ! never be deleted. Loop over all elements and "lock" nodes which 
    ! should not be removed from the triangulation due to the indicator.
    DO iel=1,SIZE(p_Dindicator)


      ! Check if the current element is marked for refinement. Then
      ! all vertices are/should be locked by the refinement procedure.
      SELECT CASE(p_Imarker(iel))
      CASE(MARK_REF_QUAD3TRIA_1,MARK_REF_QUAD3TRIA_2,&
           MARK_REF_QUAD3TRIA_3,MARK_REF_QUAD3TRIA_4,&
           MARK_REF_QUAD4TRIA_12,MARK_REF_QUAD4TRIA_23,&
           MARK_REF_QUAD4TRIA_34,MARK_REF_QUAD4TRIA_14,&
           MARK_REF_QUAD2QUAD_13,MARK_REF_QUAD2QUAD_24)
        ! The current element is a quadrilateral that is marked for green 
        ! refinement. In order to increase the performance of "phase 3" only
        ! those elements are considered which are marked for coarsening.
        ! Hence, we have to explicitely "lock" all vertices of quadrilaterals
        ! which are marked or refinement "by hand"
        DO ive=1,TRIA_NVEQUAD2D
          i=rhadapt%p_IverticesAtElement(ive,iel)
          rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        END DO
        CYCLE

      CASE(MARK_ASIS_TRIA,MARK_ASIS_QUAD)
        ! The current element is a candidate for coarsening
        
      CASE DEFAULT
        ! The current element is markerd for refinement
        CYCLE
      END SELECT


      ! Get number of vertices per element
      nve=get_NVE(rhadapt,iel)

      ! Get local data for element iel
      p_IverticesAtElement(1:nve)=rhadapt%p_IverticesAtElement(1:nve,iel)

      ! Get state of current element
      SELECT CASE(nve)
      CASE(TRIA_NVETRI2D)
        IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
            p_IverticesAtElement(1:TRIA_NVETRI2D))
        istate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))

      CASE(TRIA_NVEQUAD2D)
        IvertexAge(1:TRIA_NVEQUAD2D) = rhadapt%p_IvertexAge(&
            p_IverticesAtElement(1:TRIA_NVEQUAD2D))
        istate=redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))

      CASE DEFAULT
        CALL output_line('Invalid number of vertices per element!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
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
        ! Element IEL is red sub-element of a Quad4Quad refinement.
        ! Due to our orientation convention these cases should not appear.
        CALL output_line('This state should not appear!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
        CALL sys_halt()

      CASE(STATE_QUAD_RED4)
        ! Element IEL is red sub-element of a Quad4Quad refinement.

        ! Should the element be coarsened?
        IF (p_Dindicator(iel) .GE. rhadapt%dcoarseningTolerance) THEN

          ! If this is not the case, then "lock" all of its four nodes
          DO ive=1,TRIA_NVEQUAD2D
            i=rhadapt%p_IverticesAtElement(ive,iel)
            rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
          END DO
        ELSE

          ! Otherwise, "lock" only the node from the macro element.
          ! Due to the refinement strategy, this is the first vertex.
          i=rhadapt%p_IverticesAtElement(1,iel)
          rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

          ! Provisionally, mark element for generic coarsening
          p_Imarker(iel)=MARK_CRS_GENERIC
        END IF

      CASE(STATE_TRIA_REDINNER)
        ! Element IEL is inner red triangle of a Tria4Tria refinement.
        
        ! Should the element be coarsened?
        IF (p_Dindicator(iel) .GE. rhadapt%dcoarseningTolerance) THEN
          
          ! If this is not the case, then "lock" all of its three nodes
          DO ive=1,TRIA_NVETRI2D
            i=rhadapt%p_IverticesAtElement(ive,iel)
            rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
          END DO

        ELSE
          ! Otherwise, mark element for generic coarsening
          p_Imarker(iel)=MARK_CRS_GENERIC
        END IF

      CASE(STATE_TRIA_OUTERINNER)
        ! Element IEL is either an outer red triangle of a 1-tria : 4-tria re-
        ! finement or one of the two inner triangles of a 1-quad : 4-tria refinement
        
        ! Are we outer red triangle that was not created during the 
        ! red-green marking routine as a result of element conversion?
        jel   =rhadapt%p_IneighboursAtElement(2,iel)
        IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
                 rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,jel))
        jstate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
        
        IF (jstate .EQ. STATE_TRIA_REDINNER) THEN
          
          ! Should the element be coarsened?
          IF (p_Dindicator(iel) .GE. rhadapt%dcoarseningTolerance) THEN

            ! If this is not the case, then "lock" all of its three nodes
            DO ive=1,TRIA_NVETRI2D
              i=rhadapt%p_IverticesAtElement(ive,iel)
              rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
            END DO
          ELSE
            
            ! Otherwise, "lock" only the node from the macro element.
            ! Due to the refinement strategy, this is the first vertex.
            i=rhadapt%p_IverticesAtElement(1,iel)
            rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

            ! Provisionally, mark element for generic coarsening
            p_Imarker(iel)=MARK_CRS_GENERIC
          END IF
        ELSE
          
          ! We are an inner green triangle resulting from a 1-quad : 4-tria
          ! refinement. By construction, the first vertex belongs to the macro
          ! element and, consequently, has to be "locked" from removal
          i=rhadapt%p_IverticesAtElement(1,iel)
          rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

          ! Provisionally, mark element for generic coarsening
          p_Imarker(iel)=MARK_CRS_GENERIC
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
        i=rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

        i=rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

        ! Provisionally, mark element for generic coarsening
        p_Imarker(iel)=MARK_CRS_GENERIC
        
      CASE(STATE_TRIA_GREENOUTER_LEFT)
        ! Element IEL is left green triangle resulting from a 1-tria : 2-tria
        ! refinement. Here, the second vertex is the one which was last insrted.
        ! Hence, the first vertex is older than the second one by construction.
        i=rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

        i=rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

        ! Provisionally, mark element for generic coarsening
        p_Imarker(iel)=MARK_CRS_GENERIC

      CASE(STATE_TRIA_GREENINNER)
        ! Element IEL is the inner green triangle resulting from a 1-quad :
        ! 3-tria refinement. By construction, the first vertex is the newly
        ! introduced one which is the youngest vertex of the element. Hence,
        ! "lock" both the second and the third vertex which belong to the
        ! original macro element which was a quadrilateral
        i=rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        
        j=rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(j)=-ABS(rhadapt%p_IvertexAge(j))

        ! Provisionally, mark element for generic coarsening
        p_Imarker(iel)=MARK_CRS_GENERIC

      CASE(STATE_QUAD_HALF1,STATE_QUAD_HALF2)
        ! Element IEL is a green sub-element of a Quad2Quad refinement.
        
        ! Should the element be coarsened?
        IF (p_Dindicator(iel) .GE. rhadapt%dcoarseningTolerance) THEN

          ! If this is not the case, then "lock" all of its four nodes
          DO ive=1,TRIA_NVEQUAD2D
            i=rhadapt%p_IverticesAtElement(ive,iel)
            rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
          END DO
        ELSE
          
          ! Otherwise, "lock" only the node from the macro element.
          ! Due to the refinement strategy, this is the first and fourth vertex.
          i=rhadapt%p_IverticesAtElement(1,iel)
          rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

          i=rhadapt%p_IverticesAtElement(4,iel)
          rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

          ! Provisionally, mark element for generic coarsening
          p_Imarker(iel)=MARK_CRS_GENERIC
        END IF

      CASE DEFAULT
        CALL output_line('Invalid element state!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
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

    isModified=.TRUE.
    DO WHILE(isModified)
      isModified=.FALSE.

      ! Loop over all elements (from the initial triangulation + those which
      ! have been created during the conversion of green into red elements)
      DO iel=1,rhadapt%NEL

        ! Only process those elements which are candidates for removal
        IF (p_Imarker(iel) .NE. MARK_CRS_GENERIC) CYCLE

        ! Get number of vertices per element
        nve=get_NVE(rhadapt,iel)
        
        ! Get local data for element iel
        p_IverticesAtElement(1:nve)=rhadapt%p_IverticesAtElement(1:nve,iel)

        ! Check if all vertices of the element are locked then delete element
        ! from the list of removable elements. This could also be done in the
        ! third phase when the locked vertices are translated into coarsening
        ! rules. However, this check is relatively cheap as compared to the
        ! computation of the exact element state for some elements. Heuristically,
        ! a large number of elements can be filtered out by this check in advance
        ! so that no detailed investigation of their vertices is required.
        IF (ALL(rhadapt%p_IvertexAge(p_IverticesAtElement(1:nve)).LE.0)) THEN
          p_Imarker(iel)=MERGE(MARK_ASIS_TRIA,MARK_ASIS_QUAD,nve .EQ. 3)
          CYCLE
        END IF
        
        ! Get state of current element
        SELECT CASE(nve)
        CASE(TRIA_NVETRI2D)
          IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
              p_IverticesAtElement(1:TRIA_NVETRI2D))
          istate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))

        CASE(TRIA_NVEQUAD2D)
          IvertexAge(1:TRIA_NVEQUAD2D) = rhadapt%p_IvertexAge(&
              p_IverticesAtElement(1:TRIA_NVEQUAD2D))
          istate=redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))

        CASE DEFAULT
          CALL output_line('Invalid number of vertices per element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
          CALL sys_halt()
        END SELECT
        
        ! Do we have to "lock" some vertices?
        SELECT CASE(istate)

        CASE(STATE_QUAD_RED4)
          ! Element IEL is one of four red quadrilaterals of a 1-quad : 4-quad refinement.

          ! If the interior vertex is locked then lock all vertices of the element
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN

            ! Lock all other vertices. Note that the first vertiex and the third
            ! one are already locked, so that we have to consider vertices 2 and 4
            DO ive=2,TRIA_NVEQUAD2D,2
              rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                  -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
            END DO

            ! Delete element from list of removable elements
            p_Imarker(iel)=MARK_ASIS_QUAD
            
            ! We modified some vertex in this iteration
            isModified=.TRUE.

          ! If three midpoint vertices of adjacent quadrilaterals are locked
          ! then all vertices of the four red sub-quadrilaterals must be locked
          ! in order to prevent the generation of blue quadrilaterals
          ELSEIF(rhadapt%p_IvertexAge(p_IverticesAtElement(2)) .LE. 0 .AND.&
                 rhadapt%p_IvertexAge(p_IverticesAtElement(4)) .LE. 0) THEN
            
            ! Check the midpoint of the counterclockwise neighboring element
            jel=rhadapt%p_IneighboursAtElement(2,iel)
            i  =rhadapt%p_IverticesAtElement(2,jel)
            IF (rhadapt%p_IvertexAge(i) .LE. 0) THEN
              
              ! Lock all vertices of element IEL. Due to the fact that the interior
              ! vertex is locked, the nodes of all other quadrilaterals will be
              ! locked in the next loop of phase 2.
              DO ive=1,TRIA_NVEQUAD2D
                rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                    -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
              END DO

              ! Delete element from list of removable elements
              p_Imarker(iel)=MARK_ASIS_QUAD
              
              ! We modified some vertex in this iteration
              isModified=.TRUE.
              CYCLE
            END IF

            ! Check the midpoint of the clockwise neighboring element
            jel=rhadapt%p_IneighboursAtElement(3,iel)
            i  =rhadapt%p_IverticesAtElement(4,jel)
            IF (rhadapt%p_IvertexAge(i) .LE. 0) THEN
              
              ! Lock all vertices of element IEL. Due to the fact that the interior
              ! vertex is locked, the nodes of all other quadrilaterals will be
              ! locked in the next loop of phase 2.
              DO ive=1,TRIA_NVEQUAD2D
                rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                    -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
              END DO

              ! Delete element from list of removable elements
              p_Imarker(iel)=MARK_ASIS_QUAD
              
              ! We modified some vertex in this iteration
              isModified=.TRUE.
              CYCLE
            END IF
          END IF

        CASE(STATE_TRIA_REDINNER)
          ! Element IEL is inner red triangle of a 1-tria : 4-tria refinement.
          
          ! Determine number of locked vertices of the inner triangle.
          ivertexLock=0
          DO ive=1,TRIA_NVETRI2D
            IF (rhadapt%p_IvertexAge(p_IverticesAtElement(ive)) .LE. 0)&
                ivertexLock=ivertexLock+1
          END DO
          
          ! How many vertices are locked?
          SELECT CASE(ivertexLock)
          CASE(2)
            ! If exactly two vertices are locked, then "lock" the third one, too.
            DO ive=1,TRIA_NVETRI2D
              rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                  -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
            END DO
            
            ! Delete element from list of removable elements
            p_Imarker(iel)=MARK_ASIS
            
            ! We modified some vertex in this iteration
            isModified=.TRUE.
            
          CASE(3)
            ! If exactly three vertices are locked, then delete element from 
            ! list of removable elements
            p_Imarker(iel)=MARK_ASIS
          END SELECT

        CASE(STATE_TRIA_GREENOUTER_RIGHT)
          ! Element IEL is right green triangle of a 1-tria : 2-tria refinement.
          
          ! If the third vertex is locked, then "lock" all other vertices, too.
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
            
            ! Lock all other vertices
            DO ive=1,2
              rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                  -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
            END DO
            
            ! Delete element from list of removable elements              
            p_Imarker(iel)=MARK_ASIS
            
            ! We modified some vertex in this iteration
            isModified=.TRUE.
          END IF
          
        CASE(STATE_TRIA_GREENOUTER_LEFT)
          ! Element IEL is left green triangle of a 1-tria : 2-tria refinement.
          
          ! If the second vertex is locked, then "lock" all other vertices, too.
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(2)) .LE. 0) THEN
            
            ! Lock all other vertices
            DO ive=1,TRIA_NVETRI2D,2
              rhadapt%p_IvertexAge(p_IverticesAtElement(ive))=&
                  -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(ive)))
            END DO
            
            ! Delete element from list of removable elements              
            p_Imarker(iel)=MARK_ASIS
            
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
      nve=get_NVE(rhadapt,iel)

      ! Get local data for element IEL
      p_IverticesAtElement(1:nve)=rhadapt%p_IverticesAtElement(1:nve,iel)
      
      ! Get state of current element
      SELECT CASE(nve)
      CASE(TRIA_NVETRI2D)
        IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
            p_IverticesAtElement(1:TRIA_NVETRI2D))
        istate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))

      CASE(TRIA_NVEQUAD2D)
        IvertexAge(1:TRIA_NVEQUAD2D) = rhadapt%p_IvertexAge(&
            p_IverticesAtElement(1:TRIA_NVEQUAD2D))
        istate=redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))

      CASE DEFAULT
        CALL output_line('Invalid number of vertices per element!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
        CALL sys_halt()
      END SELECT
      
      ! What state are we?
      SELECT CASE(istate)
      CASE(STATE_TRIA_OUTERINNER,STATE_TRIA_GREENINNER)
        ! Dummy CASE to prevent immediate stop in the default branch.

      CASE(STATE_TRIA_REDINNER)
        ! If all vertices of the element are free, then the inner red triangle
        ! can be combined with its three neighboring red triangles so as to
        ! recover the original macro triangle. If one vertex of the inner
        ! triangle is locked, then the inner red tiangle together with its three
        ! neighboring elements can be convertec into two green triangles.
        IF (rhadapt%p_IvertexAge(p_IverticesAtElement(1)) .LE. 0) THEN
          p_Imarker(iel)=MARK_CRS_4TRIA2TRIA_1

        ELSEIF(rhadapt%p_IvertexAge(p_IverticesAtElement(2)) .LE. 0) THEN
          p_Imarker(iel)=MARK_CRS_4TRIA2TRIA_2
          
        ELSEIF(rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
          p_Imarker(iel)=MARK_CRS_4TRIA2TRIA_3
          
        ELSE
          p_Imarker(iel)=MARK_CRS_4TRIA1TRIA
        END IF
        
      CASE(STATE_TRIA_GREENOUTER_LEFT)
        ! We have to considered several situations depending on the state
        ! of the adjacent element that shares the second edge.
        jel=rhadapt%p_IneighboursAtElement(2,iel)
        p_IverticesAtElement(1:TRIA_NVETRI2D) = &
            rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,jel)
        IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
            p_IverticesAtElement(1:TRIA_NVETRI2D))
        jstate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
        
        ! What is the state of the "second-edge" neighbor?
        SELECT CASE(jstate)
        CASE(STATE_TRIA_GREENOUTER_RIGHT)
          ! Elements IEL and JEL represent the two green triangles which result from
          ! a 1-tria : 2-tria refinement. Since the youngest vertex is not locked
          ! we can mark the left triangle for coarsening and leave the right on as is.
          p_Imarker(iel)=MARK_CRS_2TRIA1TRIA
          p_Imarker(jel)=MARK_ASIS          
          
        CASE(STATE_TRIA_GREENINNER)
          ! If all vertices of element JEL are locked, then we can delete the
          ! element from the list of removable elements. Recall that both 
          ! vertices of the macro element are locked by definition so that it
          ! suffices to check the remaining (first) vertex.
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(1)) .LE. 0) THEN
            ! Delete element from list of removable elements              
            p_Imarker(iel)=MARK_ASIS
            p_Imarker(jel)=MARK_ASIS
            kel=rhadapt%p_IneighboursAtElement(1,jel)
            p_Imarker(kel)=MARK_ASIS
          ELSE
            ! Mark element for recoarsening together with its two neighbors
            p_Imarker(iel)=MARK_ASIS
            p_Imarker(jel)=MARK_CRS_3TRIA1QUAD
            kel=rhadapt%p_IneighboursAtElement(1,jel)
            p_Imarker(kel)=MARK_ASIS
          END IF

        CASE(STATE_TRIA_OUTERINNER)
          ! If all vertices of element JEL are locked, then we can delete the
          ! element from the list of removable elements. Recall that the vertex
          ! of the macro element is locked by definition so that it suffices to
          ! check the remaining second and third vertex individually.
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(2)) .LE. 0) THEN
            IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
              ! Delete element from list of removable elements              
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(3,jel)
              p_Imarker(kel)=MARK_ASIS
            ELSE
              ! Mark element for recoarsening together with its three neighbors,
              ! whereby the green outer triangle to the right is preserved.
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_CRS_4TRIA3TRIA_RIGHT
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(3,jel)
              p_Imarker(kel)=MARK_ASIS
            END IF
          ELSE
            IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
              ! Mark element for recoarsening together with its three neighbors,
              ! whereby the green outer triangle to the left is preserved.
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_CRS_4TRIA3TRIA_LEFT
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(3,jel)
              p_Imarker(kel)=MARK_ASIS
            ELSE
              ! Mark element for recoarsening together with its three neighbors
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_CRS_4TRIA1QUAD
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(3,jel)
              p_Imarker(kel)=MARK_ASIS
            END IF
          END IF
            
        CASE DEFAULT
          CALL output_line('Invalid element state!',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
          CALL sys_halt()
        END SELECT
        
      CASE(STATE_TRIA_GREENOUTER_RIGHT)
        ! This case is skew-symmetric to STATE_TRIA_GREENOUTER_LEFT so that 
        ! typically nothing needs to be done here. However, there is one exception
        ! to this rule. If element IEL is the right outer green triangle resulting
        ! from a 1-quad : 4-tria refinement and all vertices of the left outer 
        ! green triangle are locked, then the above CASE will never be reached.
        jel=rhadapt%p_IneighboursAtElement(2,iel)
        p_IverticesAtElement(1:TRIA_NVETRI2D) = &
            rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,jel)
        IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
            p_IverticesAtElement(1:TRIA_NVETRI2D))
        jstate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
        
        ! What is the state of the "second-edge" neighbor?
        SELECT CASE(jstate)
        CASE(STATE_TRIA_GREENOUTER_LEFT)
          ! Elements IEL and JEL represent the two green triangles which result from
          ! a 1-tria : 2-tria refinement. Since the youngest vertex is not locked
          ! we can mark the left triangle for coarsening and leave the right on as is.
          p_Imarker(iel)=MARK_ASIS
          p_Imarker(jel)=MARK_CRS_2TRIA1TRIA

        CASE(STATE_TRIA_GREENINNER)
          ! If all vertices of element JEL are locked, then we can delete the
          ! element from the list of removable elements. Recall that both 
          ! vertices of the macro element are locked by definition so that it
          ! suffices to check the remaining (first) vertex.
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(1)) .LE. 0) THEN
            ! Delete element from list of removable elements              
            p_Imarker(iel)=MARK_ASIS
            p_Imarker(jel)=MARK_ASIS
            kel=rhadapt%p_IneighboursAtElement(3,jel)
            p_Imarker(kel)=MARK_ASIS
          ELSE
            ! Mark element for recoarsening together with its two neighbors
            p_Imarker(iel)=MARK_ASIS
            p_Imarker(jel)=MARK_CRS_3TRIA1QUAD
            kel=rhadapt%p_IneighboursAtElement(3,jel)
            p_Imarker(kel)=MARK_ASIS
          END IF

        CASE(STATE_TRIA_OUTERINNER)
          ! If all vertices of element JEL are locked, then we can delete the
          ! element from the list of removable elements. Recall that the vertex
          ! of the macro element is locked by definition so that it suffices to
          ! check the remaining second and third vertex individually.
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(2)) .LE. 0) THEN
            IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
              ! Delete element from list of removable elements              
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(1,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
            ELSE
              ! Mark element for recoarsening together with its three neighbors,
              ! whereby the green outer triangle to the right is preserved.
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_CRS_4TRIA3TRIA_RIGHT
              kel=rhadapt%p_IneighboursAtElement(1,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
            END IF
          ELSE
            IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
              ! Mark element for recoarsening together with its three neighbors,
              ! whereby the green outer triangle to the left is preserved.
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_CRS_4TRIA3TRIA_LEFT
              kel=rhadapt%p_IneighboursAtElement(1,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
            ELSE
              ! Mark element for recoarsening together with its three neighbors
              p_Imarker(iel)=MARK_ASIS
              p_Imarker(jel)=MARK_CRS_4TRIA1QUAD
              kel=rhadapt%p_IneighboursAtElement(1,jel)
              p_Imarker(kel)=MARK_ASIS
              kel=rhadapt%p_IneighboursAtElement(2,jel)
              p_Imarker(kel)=MARK_ASIS
            END IF
          END IF
          
        CASE DEFAULT
          CALL output_line('Invalid element state!',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
          CALL sys_halt()
        END SELECT

      CASE(STATE_QUAD_HALF1,STATE_QUAD_HALF2)
        ! If the second vertex of element IEL is locked, then remove the element
        ! from the list of removable elements. However, it is still possible that
        ! the second vertex of the adjacent element is not locked so that "coarsening"
        ! into three triangles is initiated from that element. If the second vertex
        ! is not locked, then we check if the third vertex is not locked, either.
        ! in this case, we can mark the element for coarsening into the macro element
        ! and remove the neighboring element from the list of removable elements.
        jel=rhadapt%p_IneighboursAtElement(2,iel)
        
        ! Check if the second vertex is locked
        IF (rhadapt%p_IvertexAge(p_IverticesAtElement(2)) .LE. 0) THEN
          ! Check if the third vertex is also locked
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
            ! Delete both elements from list of removable elements
            p_Imarker(iel)=MARK_ASIS
            p_Imarker(jel)=MARK_ASIS
          ELSE
            ! Mark element IEL for coarsening into three triangles, 
            ! whereby the second vertex of element IEL is preserved
            p_Imarker(iel)=MARK_CRS_2QUAD3TRIA
            p_Imarker(jel)=MARK_ASIS
          END IF
        ELSE
          IF (rhadapt%p_IvertexAge(p_IverticesAtElement(3)) .LE. 0) THEN
            ! Mark element JEL for coarsening into three triangles,
            ! whereby the second vertex of element JEL is preserved
            p_Imarker(iel)=MARK_ASIS
            p_Imarker(jel)=MARK_CRS_2QUAD3TRIA
          ELSE
            ! Mark element IEL for recoarsening into the macro element
            p_Imarker(iel)=MARK_CRS_2QUAD1QUAD
            p_Imarker(jel)=MARK_ASIS
          END IF
        END IF
                
      CASE(STATE_QUAD_RED4)
        ! This is one of four quadrilaterals which result from a 1-quad : 4-quad
        ! refinement. Here, the situation is slightly more difficult because 
        ! multiple rotationally-symmetric situations have to be considered. The
        ! four elements forming the macro element can be collected by visiting the
        ! neighboring element along the second edge starting at element IEL. 

        ! As a first step, determine the number of locked midpoints. This is done
        ! by visiting all four elements and checking the second vertex individually.
        ! Note that we do not only count the number of locked vertices but also
        ! remember its position by setting or clearing the corresponding bit of
        ! the integer ivertexLocked. In the same loop, we determine the numbers 
        ! of the four elements of the macro element.
        p_ImacroElement(1)=iel
        ivertexLock=MERGE(2,0,&
            rhadapt%p_IvertexAge(rhadapt%p_IverticesAtElement(2,iel)) .LE. 0)
        
        DO ive=2,TRIA_NVEQUAD2D
          p_ImacroElement(ive)=rhadapt%p_IneighboursAtElement(2,p_ImacroElement(ive-1))
          IF (rhadapt%p_IvertexAge(rhadapt%p_IverticesAtElement(2,&
              p_ImacroElement(ive))) .LE. 0) ivertexLock=ibset(ivertexLock,ive)
        END DO
        
        ! How many vertices are locked?
        SELECT CASE(ivertexLock)
        CASE(0)
          ! Mark element IEL for recoarsening into the macro element and delete
          ! all remaining elements from the list of removable elements.
          p_Imarker(iel)=MARK_CRS_4QUAD1QUAD
          DO ive=2,TRIA_NVEQUAD2D
            p_Imarker(p_ImacroElement(ive))=MARK_ASIS
          END DO

        CASE(2)
          ! There is one vertex locked which is the second vertex of the first 
          ! element. All other elements are deleted from the list of removable elements
          p_Imarker(p_ImacroElement(1))=MARK_CRS_4QUAD3TRIA
          p_Imarker(p_ImacroElement(2))=MARK_ASIS
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(4)
          ! There is one vertex locked which is the second vertex of the second
          ! element. All other elements are deleted from the list of removable elements
          p_Imarker(p_ImacroElement(1))=MARK_ASIS
          p_Imarker(p_ImacroElement(2))=MARK_CRS_4QUAD3TRIA
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(8)
          ! There is one vertex locked which is the second vertex of the third
          ! element. All other elements are deleted from the list of removable elements
          p_Imarker(p_ImacroElement(1))=MARK_ASIS
          p_Imarker(p_ImacroElement(2))=MARK_ASIS
          p_Imarker(p_ImacroElement(3))=MARK_CRS_4QUAD3TRIA
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(16)
          ! There is one vertex locked which is the second vertex of the fourth
          ! element. All other elements are deleted from the list of removable elements
          p_Imarker(p_ImacroElement(1))=MARK_ASIS
          p_Imarker(p_ImacroElement(2))=MARK_ASIS
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_CRS_4QUAD3TRIA

        CASE(10)
          ! There are two vertices locked which are the second vertices of the
          ! first and third elements. Mark the first element for recoarsening.
          p_Imarker(p_ImacroElement(1))=MARK_CRS_4QUAD2QUAD
          p_Imarker(p_ImacroElement(2))=MARK_ASIS
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(20)
          ! There are two vertices locked which are the second vertices of the
          ! second and fourth elements. Mark the second element for recoarsening.
          p_Imarker(p_ImacroElement(1))=MARK_ASIS
          p_Imarker(p_ImacroElement(2))=MARK_CRS_4QUAD2QUAD
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(18)
          ! There are two vertices locked which are the second and fourth vertices
          ! of the firth elements. Mark the firth element for recoarsening.
          p_Imarker(p_ImacroElement(1))=MARK_CRS_4QUAD4TRIA
          p_Imarker(p_ImacroElement(2))=MARK_ASIS
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(6)
          ! There are two vertices locked which are the second and fourth vertices
          ! of the second elements. Mark the second element for recoarsening.
          p_Imarker(p_ImacroElement(1))=MARK_ASIS
          p_Imarker(p_ImacroElement(2))=MARK_CRS_4QUAD4TRIA
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(12)
          ! There are two vertices locked which are the second and fourth vertices
          ! of the third elements. Mark the third element for recoarsening.
          p_Imarker(p_ImacroElement(1))=MARK_ASIS
          p_Imarker(p_ImacroElement(2))=MARK_ASIS
          p_Imarker(p_ImacroElement(3))=MARK_CRS_4QUAD4TRIA
          p_Imarker(p_ImacroElement(4))=MARK_ASIS

        CASE(24)
          ! There are two vertices locked which are the second and fourth vertices
          ! of the fourth elements. Mark the fourth element for recoarsening.
          p_Imarker(p_ImacroElement(1))=MARK_ASIS
          p_Imarker(p_ImacroElement(2))=MARK_ASIS
          p_Imarker(p_ImacroElement(3))=MARK_ASIS
          p_Imarker(p_ImacroElement(4))=MARK_CRS_4QUAD4TRIA

        CASE DEFAULT
          ! Delete all four elements from list of removable elements
          DO ive=1,TRIA_NVEQUAD2D
            p_Imarker(p_ImacroElement(ive))=MARK_ASIS
          END DO
        END SELECT
        
      CASE DEFAULT
        CALL output_line('Invalid number of locked vertices!',&
            OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_coarsening2D')
        CALL sys_halt()
      END SELECT
    END DO

    ! Set specifier to "marked for coarsening"
    rhadapt%iSpec=IOR(rhadapt%iSpec,HADAPT_MARKEDCOARSEN)
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
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D) :: p_IverticesAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(1) :: Ielements
    INTEGER(PREC_VERTEXIDX), DIMENSION(1)  :: Ivertices
    INTEGER(PREC_VERTEXIDX)  :: i,nvt
    INTEGER(PREC_ELEMENTIDX) :: nel,iel,jel,kel,lel,iel1,iel2
    INTEGER :: ive,jve,nve,mve,istate,jstate,kstate
    INTEGER :: h_Imodified,imodifier
    LOGICAL :: isConform
    INTEGER, DIMENSION(TRIA_MAXNVE) :: IvertexAge
    
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
      nvt=CEILING(rhadapt%NVT0+1.5*rhadapt%nGreenElements)
      nel=rhadapt%NEL0+2*rhadapt%nGreenElements

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
      IF (PRESENT(fcb_hadaptCallback).AND.PRESENT(rcollection)) THEN
        Ivertices=(/nvt/); Ielements=(/0/)
        CALL fcb_hadaptCallback(rcollection,HADAPT_OPR_ADJUSTVERTEXDIM,&
            Ivertices,Ielements)
      END IF

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
    ! p_Imodified(iel) .EQ. imodifier. In order to pre-select elements for investigation
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
      IF (p_Imarker(iel) .NE. MARK_ASIS) p_Imodified(iel)=imodifier
    END DO

    isConform=.FALSE.
    conformity: DO

      ! If conformity is guaranteed for all cells, then exit
      IF (isConform) EXIT conformity
      isConform=.TRUE.
      
      ! Otherwise, loop over all elements present in the initial grid (IEL <= NEL0)
      ! which are modified and mark their neighbors for further refinement 
      ! if conformity is violated for some edge     
      DO iel=1,rhadapt%NEL0
        
        ! Skip those element which have not been modified
        IF (p_Imodified(iel) .NE. imodifier) CYCLE

        ! Get number of vertices per element
        nve=get_NVE(rhadapt,iel)
        
        ! Get local data for element iel
        p_IverticesAtElement(1:nve)=rhadapt%p_IverticesAtElement(1:nve,iel)
      
        ! Are we triangular or quadrilateral element?
        SELECT CASE(nve)
        CASE(TRIA_NVETRI2D)
          p_Imarker(iel)=ibclr(p_Imarker(iel),0)
          IvertexAge(1:TRIA_NVETRI2D) = &
              rhadapt%p_IvertexAge(p_IverticesAtElement(1:TRIA_NVETRI2D))
          istate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))

        CASE(TRIA_NVEQUAD2D)
          p_Imarker(iel)=ibset(p_Imarker(iel),0)
          IvertexAge(1:TRIA_NVEQUAD2D) = &
              rhadapt%p_IvertexAge(p_IverticesAtElement(1:TRIA_NVEQUAD2D))
          istate=redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))

        CASE DEFAULT
          CALL output_line('Invalid number of vertices per element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_refinement2D')
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
          CALL output_line('These states must not occur!',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_refinement2D')
          CALL sys_halt()
          
        CASE(STATE_TRIA_OUTERINNER)
          ! The triangle can either be an outer red element resulting from a Tria4Tria 
          ! refinement or one of the two opposite triangles resulting from Quad4Tria 
          ! refinement. This can be easily checked. If the opposite element is not 
          ! an inner red triangle (14), then the element IEL and its opposite neighbor 
          ! make up the inner "diamond" resulting from a Quad4Tria refinement. 
          jel=rhadapt%p_IneighboursAtElement(2,iel)
       
          ! First, we need to check a special case. If the second edge of element IEL
          ! has two different neighbors, then the original neighbor along this edge
          ! was a green triangle that has been converted into an (outer) red one. In
          ! this case, the current element IEL does not need further modifications.
          IF (jel .NE. rhadapt%p_ImidneighboursAtElement(2,iel)) GOTO 100
          
          ! Otherwise, determine the state of the edge neighbor JEL
          jstate=redgreen_getState(rhadapt,jel)
          
          IF (jstate .EQ. STATE_TRIA_OUTERINNER) THEN
            ! We know that element IEL and JEL make up the inner "diamond" resulting
            ! from a Quad4Tria refinement. It remains to find the two green triangle
            ! which make up the outer part of the macro element. Since we do not know
            ! if they are adjacent to element IEL or JEL we have to perform additional
            ! checks: The element along the first edge of the inner triangle must
            ! have state STATE_TRIA_GREENOUTER_LEFT.
            kel=rhadapt%p_IneighboursAtElement(1,iel)
            kstate=redgreen_getState(rhadapt,kel)

            IF (kstate .NE. STATE_TRIA_GREENOUTER_LEFT .OR.&
                rhadapt%p_IneighboursAtElement(2,kel).NE.iel) THEN
              ! At this stage we can be sure that element IEL is not (!) the inner
              ! triangle, hence, it must be element JEL
              
              ! To begin with, we need to find the two missing triangles KEL and LEL
              ! which make up the original quadrilateral
              kel=rhadapt%p_IneighboursAtElement(1,jel)
              lel=rhadapt%p_IneighboursAtElement(3,jel)

              ! Mark the edge of the element adjacent to KEL for subdivision
              iel1=rhadapt%p_IneighboursAtElement(3,kel)
              iel2=rhadapt%p_ImidneighboursAtElement(3,kel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)

              ! Mark the edge of the element adjacent to LEL for subdivision
              iel1=rhadapt%p_IneighboursAtElement(1,lel)
              iel2=rhadapt%p_ImidneighboursAtElement(1,lel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

              ! Now, we can physically convert the four triangles into four quadrilaterals
              CALL convert_Quad4Tria(rhadapt,kel,iel,lel,jel,rcollection,fcb_hadaptCallback)
              isConform=.FALSE.

              ! All four elements have to be converted from triangles to quadrilaterals.
              p_Imarker(jel)=ibset(0,0)

              ! The second and third edge of KEL must be unmarked. Moreover, the first edge is
              ! marked for refinement if and only if it is also marked from the adjacent element.
              iel1=rhadapt%p_IneighboursAtElement(1,kel)
              iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
              IF (ismarked_edge(iel1,iel2,kel)) THEN
                p_Imarker(kel)=ibset(ibset(0,1),0)
              ELSE
                p_Imarker(kel)=ibset(0,0)
              END IF

              ! The second and third edge of IEL must be unmarked.
              p_Imarker(iel)=ibset(0,0)

              ! The fourth edge is only marked if it is also marked from the adjacent element
              iel1=rhadapt%p_IneighboursAtElement(4,iel)
              iel2=rhadapt%p_ImidneighboursAtElement(4,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),4)
              
              ! The first edge is only marked if it is also marked from the adjacent element              
              iel1=rhadapt%p_IneighboursAtElement(1,iel)
              iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
              IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),1)
                      
              ! The first, second and third edge of LEL must be unmarked.
              p_Imarker(lel)=ibset(0,0)
              
              ! The fourth edge is only marked if it is also marked from the adjacent element              
              iel1=rhadapt%p_IneighboursAtElement(4,lel)
              iel2=rhadapt%p_ImidneighboursAtElement(4,lel)
              IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel)=ibset(p_Imarker(lel),4)

            ELSE

              ! Otherwise, element IEL is (!) the inner triangle.
              
              ! To begin with, we need to find the two missing triangles KEL and LEL
              ! which make up the original quadrilateral
              kel=rhadapt%p_IneighboursAtElement(1,iel)
              lel=rhadapt%p_IneighboursAtElement(3,iel)

              ! Mark the edge of the element adjacent to KEL for subdivision
              iel1=rhadapt%p_IneighboursAtElement(3,kel)
              iel2=rhadapt%p_ImidneighboursAtElement(3,kel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)
              
              ! Mark the edge of the element adjacent to LEL for subdivision
              iel1=rhadapt%p_IneighboursAtElement(1,lel)
              iel2=rhadapt%p_ImidneighboursAtElement(1,lel)
              IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

              ! Now, we can physically convert the four triangles into four quadrilaterals
              CALL convert_Quad4Tria(rhadapt,kel,jel,lel,iel,rcollection,fcb_hadaptCallback)
              isConform=.FALSE.

              ! All four elements have to be converted from triangles to quadrilaterals.
              p_Imarker(iel)=ibset(0,0)

              ! The second and third edge of KEL must be unmarked. Moreover, the first edge is
              ! marked for refinement if and only if it is also marked from the adjacent element
              iel1=rhadapt%p_IneighboursAtElement(1,kel)
              iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
              IF (ismarked_edge(iel1,iel2,kel)) THEN
                p_Imarker(kel)=ibset(ibset(0,1),0)
              ELSE
                p_Imarker(kel)=ibset(0,0)
              END IF

              ! The second and third edge of JEL must be unmarked. 
              p_Imarker(jel)=ibset(0,0)

              ! The fourth edge is only marked if it is also marked from the adjacent element
              iel1=rhadapt%p_IneighboursAtElement(4,jel)
              iel2=rhadapt%p_ImidneighboursAtElement(4,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel)=ibset(p_Imarker(jel),4)
              
              ! The first edge is only marked if it is also marked from the adjacent element
              iel1=rhadapt%p_IneighboursAtElement(1,jel)
              iel2=rhadapt%p_ImidneighboursAtElement(1,jel)
              IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel)=ibset(p_Imarker(jel),1)
              
              ! The first, second and third edge of LEL must be unmarked.
              p_Imarker(lel)=ibset(0,0)

              ! The fourth edge is only marked if it is also marked from the adjacent element
              iel1=rhadapt%p_IneighboursAtElement(4,lel)
              iel2=rhadapt%p_ImidneighboursAtElement(4,lel)
              IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel)=ibset(p_Imarker(lel),4)
              
            END IF
          END IF
                    
        CASE(STATE_QUAD_HALF1,STATE_QUAD_HALF2)
          ! Element is quadrilateral that results from a Quad2Quad refinement.
          ! Due to our orientation convention the other "halved" element is 
          ! adjacent to the second edge. 
          jel=rhadapt%p_IneighboursAtElement(2,iel)

          ! Mark the fourth edge of elements IEL for subdivision
          iel1=rhadapt%p_IneighboursAtElement(4,iel)
          iel2=rhadapt%p_ImidneighboursAtElement(4,iel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

          ! Mark the fourth edge of elements JEL for subdivision
          iel1=rhadapt%p_IneighboursAtElement(4,jel)
          iel2=rhadapt%p_ImidneighboursAtElement(4,jel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
          
          ! Now, we can physically convert the two quadrilaterals into four quadrilaterals
          CALL convert_Quad2Quad(rhadapt,iel,jel,rcollection,fcb_hadaptCallback)
          isConform=.FALSE.

          ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
          kel= rhadapt%p_IneighboursAtElement(2,jel)
          lel= rhadapt%p_IneighboursAtElement(3,iel)
          
          ! As a first step, clear all four elements and mark them as quadrilaterals.
          p_Imarker(iel)=ibset(0,0)
          p_Imarker(jel)=ibset(0,0)
          p_Imarker(kel)=ibset(0,0)
          p_Imarker(lel)=ibset(0,0)
          
          ! In addition, we have to transfer some markers from element IEL and JEL         
          ! to the new elements NEL0+1, NEL0+2.           
          iel1=rhadapt%p_IneighboursAtElement(1,iel)
          iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
          IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),1)
          
          iel1=rhadapt%p_IneighboursAtElement(4,jel)
          iel2=rhadapt%p_ImidneighboursAtElement(4,jel)
          IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel)=ibset(p_Imarker(jel),4)
          
          iel1=rhadapt%p_IneighboursAtElement(1,kel)
          iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
          IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),1)
          
          iel1=rhadapt%p_IneighboursAtElement(4,lel)
          iel2=rhadapt%p_ImidneighboursAtElement(4,lel)
          IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel)=ibset(p_Imarker(lel),4)
         
        CASE(STATE_TRIA_GREENINNER)
          ! We are processing a green triangle. Due to our refinement convention, element IEL
          ! can only be the inner triangle resulting from a Quad3Tria refinement.
          jel=rhadapt%p_IneighboursAtElement(3,iel)
          kel=rhadapt%p_IneighboursAtElement(1,iel)

          ! To begin with, we need to mark the edges of the adjacent elements for subdivision.
          iel1=rhadapt%p_IneighboursAtElement(2,iel)
          iel2=rhadapt%p_ImidneighboursAtElement(2,iel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
          
          ! Mark the edge of the element adjacent to JEL for subdivision
          iel1=rhadapt%p_IneighboursAtElement(3,jel)
          iel2=rhadapt%p_ImidneighboursAtElement(3,jel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)

          ! Mark the edge of the element adjacent to KEL for subdivision
          iel1=rhadapt%p_IneighboursAtElement(1,kel)
          iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
          IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)

          ! Now, we can physically convert the three elements IEL,JEL and KEL
          ! into four similar quadrilaterals
          CALL convert_Quad3Tria(rhadapt,jel,kel,iel,rcollection,fcb_hadaptCallback)
          isConform=.FALSE.
          
          ! The new element NEL0+1 has zero marker by construction but that of IEL must
          ! be nullified. The markers for the modified triangles also need to be adjusted.
          p_Imarker(iel)=ibset(0,0)
          p_Imarker(jel)=ibset(0,0)
          p_Imarker(kel)=ibset(0,0)

          ! The first edge is only marked if it is also marked from the adjacent element              
          iel1=rhadapt%p_IneighboursAtElement(1,jel)
          iel2=rhadapt%p_ImidneighboursAtElement(1,jel)
          IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel)=ibset(p_Imarker(jel),1)
          
          ! The fourth edge is only marked if it is also marked from the adjacent element              
          iel1=rhadapt%p_IneighboursAtElement(4,kel)
          iel2=rhadapt%p_ImidneighboursAtElement(4,kel)
          IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),4)
          
        CASE(STATE_TRIA_GREENOUTER_LEFT)
          ! We are processing a green triangle. Here, we have to consider several cases.
          ! First, let us find out the state of the neighboring element JEL.
          jel   =rhadapt%p_IneighboursAtElement(2,iel)
          jstate=redgreen_getState(rhadapt,jel)

          ! What state is element JEL
          SELECT CASE(jstate)
          CASE(STATE_TRIA_GREENOUTER_RIGHT)
            ! Element IEL and JEL are the result of a Tria2Tria refinement, whereby
            ! triangle IEL is located left to element JEL. We can safely convert both
            ! elements into one and perform regular refinement afterwards.

            ! To begin with, we need to mark the edges of the elements adjacent 
            ! to IEL and JEL for subdivision. 
            iel1=rhadapt%p_IneighboursAtElement(3,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
            
            ! The same procedure must be applied to the neighbor of element JEL
            iel1=rhadapt%p_IneighboursAtElement(1,jel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
            
            ! Now, we can physically convert the two elements IEL and JEL into four similar triangles
            CALL convert_Tria2Tria(rhadapt,iel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
            ! The markers for the modified elements IEL and JEL need to be adjusted.
            p_Imarker(iel)=0
            p_Imarker(jel)=0

            iel1=rhadapt%p_IneighboursAtElement(1,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),1)
            
            iel1=rhadapt%p_IneighboursAtElement(3,jel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,jel)
            IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel)=ibset(p_Imarker(jel),3)

          
          CASE(STATE_TRIA_GREENINNER)
            ! Element IEL and JEL are the result of a Quad3Tria refinement, whereby
            ! element JEL is the inner triangle and IEL is its left neighbor.

            ! To begin with, we need to mark the edges of the elements adjacent to IEL,
            ! JEL and KEL for subdivision. Let us start with the third neighbor of IEL.
            iel1=rhadapt%p_IneighboursAtElement(3,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

            ! The same procedure must be applied to the second neighbor of element JEL
            iel1=rhadapt%p_IneighboursAtElement(2,jel)
            iel2=rhadapt%p_ImidneighboursAtElement(2,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
            
            ! And again, the same procedure must be applied to the first neighbor of 
            ! element KEL which is the first neighbor of JEL, whereby KEL needs to 
            ! be found in the dynamic data structure first
            kel=rhadapt%p_IneighboursAtElement(1,jel)
            
            ! Ok, now we can proceed to its first neighbor
            iel1=rhadapt%p_IneighboursAtElement(1,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)
            
            ! Now, we can physically convert the three elements IEL,JEL and KEL
            ! into four similar quadrilaterals
            CALL convert_Quad3Tria(rhadapt,iel,kel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new element NEL0+1 has zero marker by construction but that of JEL must 
            ! be nullified. The markers for the modified triangles also need to be adjusted.
            p_Imarker(iel)=ibset(0,0)
            p_Imarker(jel)=ibset(0,0)
            p_Imarker(kel)=ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(1,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(4,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(4,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),4)
           
          CASE(STATE_TRIA_OUTERINNER)
            ! Element IEL and JEL are the result of a Quad4Tria refinement, whereby
            ! element JEL is the inner triangles  and IEL is its right neighbor.
            kel=rhadapt%p_IneighboursAtElement(2,jel)
            lel=rhadapt%p_IneighboursAtElement(3,jel)

            ! Mark the edge of the element adjacent to IEL for subdivision
            iel1=rhadapt%p_IneighboursAtElement(3,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
            
            ! Mark the edge of the element adjacent to LEL for subdivision
            iel1=rhadapt%p_IneighboursAtElement(1,lel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,lel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

            ! Now, we can physically convert the four triangles into four quadrilaterals
            CALL convert_Quad4Tria(rhadapt,iel,kel,lel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! All four elements have to be converted from triangles to quadrilaterals.
            p_Imarker(iel)=ibset(0,0)
            p_Imarker(jel)=ibset(0,0)
            p_Imarker(kel)=ibset(0,0)
            p_Imarker(lel)=ibset(0,0)
            
            ! The first edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(1,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(4,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(4,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),4)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(1,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(4,lel)
            iel2=rhadapt%p_ImidneighboursAtElement(4,lel)
            IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel)=ibset(p_Imarker(lel),4)

          CASE DEFAULT
            CALL output_line('Invalid element state!',&
                OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_refinement2D')
            CALL sys_halt()
          END SELECT

        CASE(STATE_TRIA_GREENOUTER_RIGHT)
          ! We are processing a green triangle. Here, we have to consider several cases.
          ! First, let us find out the state of the neighboring element JEL.
          jel   =rhadapt%p_IneighboursAtElement(2,iel)
          jstate=redgreen_getState(rhadapt,jel)
          
          ! What state is element JEL
          SELECT CASE(jstate)
          CASE(STATE_TRIA_GREENOUTER_LEFT)
            ! Element IEL and JEL are the result of a Tria2Tria refinement, whereby
            ! triangle IEL is located right to element JEL. We can safely convert both
            ! elements into one and perform regular refinement afterwards.
            
            ! To begin with, we need to mark the edges of the elements adjacent
            ! to IEL and JEL for subdivision.
            iel1=rhadapt%p_IneighboursAtElement(1,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)
            
            ! The same procedure must be applied to the neighbor of element JEL
            iel1=rhadapt%p_IneighboursAtElement(3,jel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)
            
            ! Now, we can physically convert the two elements IEL and JEL into four similar triangles
            CALL convert_Tria2Tria(rhadapt,jel,iel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new elements NEL0+1 and NEL0+2 have zero markers by construction.
            ! The markers for the modified elements IEL and JEL need to be adjusted
            p_Imarker(iel)=0
            p_Imarker(jel)=0
            
            iel1=rhadapt%p_IneighboursAtElement(1,jel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,jel)
            IF (ismarked_edge(iel1,iel2,jel)) p_Imarker(jel)=ibset(p_Imarker(jel),1)
            
            iel1=rhadapt%p_IneighboursAtElement(3,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),3)


          CASE(STATE_TRIA_GREENINNER)
            ! Element IEL and JEL are the result of a Quad3Tria refinement, whereby
            ! element JEL is the inner triangle and IEL is its left neighbor.

            ! To begin with, we need to mark the edges of the elements adjacent to IEL,
            ! JEL and KEL for subdivision. Let us start with the firth neighbor of IEL.
            iel1=rhadapt%p_IneighboursAtElement(1,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

            ! The same procedure must be applied to the second neighbor of element JEL
            iel1=rhadapt%p_IneighboursAtElement(2,jel)
            iel2=rhadapt%p_ImidneighboursAtElement(2,jel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(jel,iel1)

            ! And again, the same procedure must be applied to the third neighbor of 
            ! element KEL which is the third neighbor of JEL, whereby KEL needs to 
            ! be found in the dynamic data structure first
            kel=rhadapt%p_IneighboursAtElement(3,jel)

            ! Ok, now we can proceed to its first neighbor
            iel1=rhadapt%p_IneighboursAtElement(3,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,kel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(kel,iel1)

            ! Now, we can physically convert the three elements IEL,JEL and KEL
            ! into four similar quadrilaterals
            CALL convert_Quad3Tria(rhadapt,kel,iel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! The new element NEL0+1 has zero marker by construction but that of JEL must
            ! be nullified. The markers for the modified triangles also need to be adjusted.
            p_Imarker(iel)=ibset(0,0)
            p_Imarker(jel)=ibset(0,0)
            p_Imarker(kel)=ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(1,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(4,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(4,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),4)
            
          CASE(STATE_TRIA_OUTERINNER)
            ! Element IEL and JEL are the result of a Quad4Tria refinement, whereby
            ! element JEL is the inner triangles  and IEL is its left neighbor.
            kel=rhadapt%p_IneighboursAtElement(2,jel)
            lel=rhadapt%p_IneighboursAtElement(1,jel)

            ! Mark the edge of the element adjacent to IEL for subdivision
            iel1=rhadapt%p_IneighboursAtElement(1,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,iel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(iel,iel1)

            ! Mark the edge of the element adjacent to LEL for subdivision
            iel1=rhadapt%p_IneighboursAtElement(3,lel)
            iel2=rhadapt%p_ImidneighboursAtElement(3,lel)
            IF (iel1*iel2 .NE. 0 .AND. iel1 .EQ. iel2) CALL mark_edge(lel,iel1)

            ! Now, we can physically convert the four triangles into four quadrilaterals
            CALL convert_Quad4Tria(rhadapt,lel,kel,iel,jel,rcollection,fcb_hadaptCallback)
            isConform=.FALSE.
            
            ! All four elements have to be converted from triangles to quadrilaterals.
            p_Imarker(iel)=ibset(0,0)
            p_Imarker(jel)=ibset(0,0)
            p_Imarker(kel)=ibset(0,0)
            p_Imarker(lel)=ibset(0,0)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(1,lel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,lel)
            IF (ismarked_edge(iel1,iel2,lel)) p_Imarker(lel)=ibset(p_Imarker(lel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(4,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(4,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),4)

            ! The first edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(1,kel)
            iel2=rhadapt%p_ImidneighboursAtElement(1,kel)
            IF (ismarked_edge(iel1,iel2,kel)) p_Imarker(kel)=ibset(p_Imarker(kel),1)

            ! The fourth edge is only marked if it is also marked from the adjacent element
            iel1=rhadapt%p_IneighboursAtElement(4,iel)
            iel2=rhadapt%p_ImidneighboursAtElement(4,iel)
            IF (ismarked_edge(iel1,iel2,iel)) p_Imarker(iel)=ibset(p_Imarker(iel),4)

          CASE DEFAULT
            CALL output_line('Invalid element state!',&
                OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_refinement2D')
            CALL sys_halt()
          END SELECT
          
        CASE DEFAULT
          CALL output_line('Invalid element state!',&
              OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_refinement2D')
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
            jel=rhadapt%p_IneighboursAtElement(ive,iel)
            IF (jel .EQ. 0) CYCLE
            
            ! Get number of vertices per element
            mve=get_NVE(rhadapt,jel)
            
            ! Are we triangular or quadrilateral element?
            SELECT CASE(mve)
            CASE(TRIA_NVETRI2D)
              p_Imarker(jel)=ibclr(p_Imarker(jel),0)

            CASE(TRIA_NVEQUAD2D)
              p_Imarker(jel)=ibset(p_Imarker(jel),0)

            CASE DEFAULT
              CALL output_line('Invalid number of vertices per element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_refinement2D')
              CALL sys_halt()
            END SELECT

            ! Now, we need to find the local position of element IEL in the
            ! adjacency list of the nieghboring element JEL
            DO jve=1,mve
              IF (rhadapt%p_IneighboursAtElement(jve,jel) .EQ. iel) EXIT
            END DO

            IF (jve > mve) THEN
              CALL output_line('Unable to find element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'redgreen_mark_refinement2D')
              CALL sys_halt()
            END IF

            ! If the edge is already marked for refinement then we can
            ! guarantee conformity for the edge shared by IEL and JEL.
            ! Otherwise, this edge needs to be marked "looking" from JEL
            IF (.NOT.btest(p_Imarker(jel),jve)) THEN
              p_Imarker(jel)=ibset(p_Imarker(jel),jve)
              isConform     =.FALSE.
              IF (p_Imodified(jel) .NE. imodifier) p_Imodified(jel)= -imodifier
            END IF

            ! Finally, "lock" all vertices of element JEL
            p_IverticesAtElement(1:mve)=rhadapt%p_IverticesAtElement(1:mve,jel)
            rhadapt%p_IvertexAge(p_IverticesAtElement(1:mve))=&
                -ABS(rhadapt%p_IvertexAge(p_IverticesAtElement(1:mve)))
          END IF
        END DO
      END DO
      
      ! Note that the number of elements (and vertices) may have changed
      ! due to conversion of green into red elements. Hence, adjust dimensions.
      rhadapt%NEL0=rhadapt%NEL
      rhadapt%NVT0=rhadapt%NVT

      !--------------------------------------------------------------
      ! We don't want to have elements with two and/or three divided edges, 
      ! aka, blue refinement for triangles and quadrilaterals, respectively.
      ! As a remedy, these elements are filtered and marked for red refinement.
      !--------------------------------------------------------------
      DO iel=1,rhadapt%NEL0
        IF (p_Imodified(iel)*imodifier .EQ. 0) CYCLE
        
        ! What type of element are we?
        SELECT CASE(p_Imarker(iel))
        CASE(MARK_REF_TRIA3TRIA_12,MARK_REF_TRIA3TRIA_13,&
             MARK_REF_TRIA3TRIA_23)
          ! Blue refinement for triangles is not allowed.
          ! Hence, mark triangle for red refinement
          p_Imarker(iel)        = MARK_REF_TRIA4TRIA
          p_Imodified(iel)      =-imodifier
          isConform             =.FALSE.
          rhadapt%increaseNVT   = rhadapt%increaseNVT+1
          
        CASE(MARK_REF_QUADBLUE_123,MARK_REF_QUADBLUE_412,&
             MARK_REF_QUADBLUE_341,MARK_REF_QUADBLUE_234)
          ! Blue refinement for quadrilaterals is not allowed. 
          ! Hence, mark quadrilateral for red refinement
          p_Imarker(iel)        = MARK_REF_QUAD4QUAD
          p_Imodified(iel)      =-imodifier
          isConform             =.FALSE.
          rhadapt%increaseNVT   = rhadapt%increaseNVT+2

        CASE DEFAULT
          p_Imodified(iel)=(p_Imodified(iel)-imodifier)/2
        END SELECT
      END DO
      
      ! Reverse modifier
      imodifier=-imodifier
    END DO conformity

    ! As a last step, lock all vertices of those elements which are marked 
    ! for refinement and increase the number of new vertices by the number
    ! of quadrilateral elements which are marked for red refinement
    DO iel=1,SIZE(p_Imarker)
      
      ! What type of element are we?
      SELECT CASE(p_Imarker(iel))
      CASE(MARK_REF_TRIA2TRIA_1)
        i=rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        
      CASE(MARK_REF_TRIA2TRIA_2)
        i=rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_TRIA2TRIA_3)
        i=rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_QUAD3TRIA_1)
        i=rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        i=rhadapt%p_IverticesAtElement(4,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        
      CASE(MARK_REF_QUAD3TRIA_2)
        i=rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        i=rhadapt%p_IverticesAtElement(4,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_QUAD3TRIA_3)
        i=rhadapt%p_IverticesAtElement(1,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        i=rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_QUAD3TRIA_4)
        i=rhadapt%p_IverticesAtElement(2,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))
        i=rhadapt%p_IverticesAtElement(3,iel)
        rhadapt%p_IvertexAge(i)=-ABS(rhadapt%p_IvertexAge(i))

      CASE(MARK_REF_QUAD4QUAD)
        rhadapt%increaseNVT=rhadapt%increaseNVT+1
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
      INTEGER :: ive,nve

      ! Get number of vertices per element
      nve=get_NVE(rhadapt,jel)
      
      ! Find local position of element IEL in adjacency list of JEL
      SELECT CASE(nve)
      CASE(TRIA_NVETRI2D)
        ! Triangular elements
        DO ive=1,nve
          IF (rhadapt%p_IneighboursAtElement(ive,jel) .EQ. iel) THEN
            p_Imarker(jel)=ibclr(ibset(p_Imarker(jel),ive),0)
            isConform      =.FALSE.
            IF (p_Imodified(jel) .NE. imodifier) p_Imodified(jel)=-imodifier
            EXIT
          END IF
        END DO

      CASE(TRIA_NVEQUAD2D)
        ! Quadrilateral element
        DO ive=1,nve
          IF (rhadapt%p_IneighboursAtElement(ive,jel) .EQ. iel) THEN
            p_Imarker(jel)=ibset(ibset(p_Imarker(jel),ive),0)
            isConform      =.FALSE.
            IF (p_Imodified(jel) .NE. imodifier) p_Imodified(jel)=-imodifier
            EXIT
          END IF
        END DO

      CASE DEFAULT
        CALL output_line('Invalid number of vertices per element!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mark_edge')
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
      IF (iel*ielmid .EQ. 0) THEN
        bismarked=.FALSE.
        RETURN
      ELSEIF(iel .NE. ielmid) THEN
        bismarked=.TRUE.
        RETURN
      END IF
      
      ! Loop over all edges of element IEL and try to find neighbor JEL
      DO ive=1,get_NVE(rhadapt,iel)
        IF (rhadapt%p_IneighboursAtElement(ive,iel)    .EQ. jel .OR.&
            rhadapt%p_ImidneighboursAtElement(ive,iel) .EQ. jel) THEN
          bismarked=btest(p_Imarker(iel),ive)
          RETURN
        END IF
      END DO
      
      CALL output_line('Unable to find common egde!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ismarked_edge')
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
      CALL output_line('Dynamic data structures are not generated &
          &or no marker for refinement is available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'redgreen_refine')
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
    DO iel=SIZE(p_Imarker),1,-1

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
          DO ive=1,get_NVE(rhadapt,jel)
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

    INTEGER, DIMENSION(TRIA_MAXNVE) :: IvertexAge

    ! Are we in 2D or 3D?
    IF (rhadapt%ndim .EQ. NDIM2D) THEN
      
      ! Do we have quadrilaterals in the triangulation?
      IF (rhadapt%InelOfType(TRIA_NVEQUAD2D).EQ.0) THEN
        
        ! There are no quadrilaterals in the current triangulation.
        IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
            rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,iel))
        istate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
        
      ELSE

        ! There are quadrilaterals and possible also triangles in
        ! the current triangulation. If the last entry of the vertices
        ! at element list is nonzero then the current element is a 
        ! quadrilateral. Otherwise, we are dealing with  a triangle
        IF (rhadapt%p_IverticesAtElement(TRIA_NVEQUAD2D,iel).EQ.0) THEN
          IvertexAge(1:TRIA_NVETRI2D) = rhadapt%p_IvertexAge(&
            rhadapt%p_IverticesAtElement(1:TRIA_NVETRI2D,iel)) 
          istate=redgreen_getStateTria(IvertexAge(1:TRIA_NVETRI2D))
        ELSE
          IvertexAge(1:TRIA_NVEQUAD2D) = rhadapt%p_IvertexAge(&
            rhadapt%p_IverticesAtElement(1:TRIA_NVEQUAD2D,iel)) 
          istate=redgreen_getStateQuad(IvertexAge(1:TRIA_NVEQUAD2D))
        END IF

      END IF

    ELSE
      
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
    INTEGER, DIMENSION(TRIA_NVETRI2D), INTENT(IN) :: IvertexAge
!</input>

!<result>
    ! State of the triangle
    INTEGER :: istate
!</result>
!</function>

    ! local variables
    INTEGER :: ive,ipos

    ! Reset state
    istate=0
    
    ! Check if triangle belongs to the initial triangulation. Then nothing
    ! needs to be done. Otherwise, perform additional checks
    IF (SUM(ABS(IvertexAge)).NE.0) THEN
      
      ! Check if any two vertices have the same age and mark that edge.
      ! In addition, determine the largest age and its position
      ipos=1
      DO ive=1,TRIA_NVETRI2D
        IF (ABS(IvertexAge(ive)) .EQ. ABS(IvertexAge(MOD(ive,TRIA_NVETRI2D)+1))) THEN
          ! Edge connects two nodes with the same age
          istate=ibset(istate,ive)
        ELSEIF (ABS(IvertexAge(ive)) .GT. ABS(IvertexAge(ipos))) THEN
          ! The age of the "next" node is larger, so than it might be the largest
          ipos=ive
        END IF
      END DO
      
      ! If ISTATE=0, then there must be one youngest vertex which is returned
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
    istate=ibclr(istate,0)
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
    INTEGER, DIMENSION(TRIA_NVEQUAD2D), INTENT(IN) :: IvertexAge
!</input>

!<result>
    ! State of the quadrilateral
    INTEGER :: istate
!</result>
!</function>

    ! local variables
    INTEGER :: ive

    ! Reset state
    istate=0

    ! Check if quadrilateral belongs to the initial triangulation. Then nothing
    ! needs to be done. Otherwise, perform additional checks.
    IF (SUM(ABS(IvertexAge)).NE.0) THEN
      
      ! Check if any two vertices have the same age and mark that edge.
      ! After this procedure, ISTATE must be different from 0 and a unique
      ! state for the quadrilateral has been determined.
      DO ive=1,TRIA_NVEQUAD2D
        IF (ABS(IvertexAge(ive)) .EQ. ABS(IvertexAge(MOD(ive,TRIA_NVEQUAD2D)+1))) THEN
          ! Edge connects two nodes with the same age
          istate=ibset(istate,ive)
        END IF
      END DO
    END IF

    ! Mark state for quadrilateral
    istate=ibset(istate,0)
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
      inewstate=ishft(istate,-1)
      inewstate=ishftc(inewstate,irotate,4)
      inewstate=ishft(inewstate,1)
      inewstate=ibset(inewstate,0)
    ELSE
      inewstate=ishft(istate,-1)
      inewstate=ishftc(inewstate,irotate,3)
      inewstate=ishft(inewstate,1)
      inewstate=ibclr(inewstate,0)
    END IF
  END FUNCTION redgreen_rotateState
END MODULE hadaptivity

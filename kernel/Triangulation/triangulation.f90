!##############################################################################
!# ****************************************************************************
!# <name> triangulation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the geometric definition of a triangulation.
!# A triangulation consists of three elemental geometric objects:
!# - Points
!# - Lines
!# - Polygons
!# Each polygon (or just called "element" or "element primitive") consists
!# of a number of lines and a number of corner points. The polygons must not
!# overlap, but usually share common edges and/or points.
!#
!# The following routines can be found here:
!#
!#  1.) tria_duplicate
!#      -> Creates a duplicate / backup of a triangulation.
!#         Some information may be shared between two triangulation structures.
!#
!#  2.) tria_restore
!#      -> Restores a triangulation previously backed up with tria_backup.
!#
!#  3.) tria_done
!#      -> Cleans up a triangulation structure, releases memory from the heap.
!#
!#  4.) tria_createRawTria1D
!#      -> Creates a 'raw' 1D triangulation $[a,b]$ with $n$ sub-intervals 
!#         of the same length
!#
!#  5.) tria_readTriFile1D
!#      -> Reads a .TRI file and creates a 'raw' 1D mesh with only basic 
!#         information.
!#
!#  6.) tria_readTriFile2D
!#      -> Reads a .TRI file and creates a 'raw' 2D mesh with only basic 
!#         information.
!#
!#  7.) tria_readTriFile3D
!#      -> Reads a .TRI file and creates a 'raw' 3D mesh with only basic 
!#         information.
!#
!#  8.) tria_resetToRaw
!#      -> Resets a mesh to a raw mesh.
!#
!#  9.) tria_initExtendedRawMesh
!#      -> Create an extended raw mesh from a pure raw mesh.
!#
!# 10.) tria_initStandardMeshFromRaw
!#      -> Generates all standard arrays for a mesh, i.e. converts a 'raw' mesh
!#         (as set up by tria_readTriFile2D e.g.) to a standard mesh.
!#
!# 11.) tria_initMacroNodalProperty
!#      -> Attaches a macro nodal property array to an extended raw mesh.
!#
!# 12.) tria_refine2LevelOrdering
!#      -> Refines a mesh according to the 2-level ordering algorithm.
!#         Creates a 'raw' fine mesh from a 'standard' coarse mesh.
!#
!# 13.) tria_compress2LevelOrdHierarchy
!#      -> Can be used to compress a mesh hierarchy created with the 2-level
!#         ordering. Saves memory. Shares vertex coordinates between fine
!#         and coarse mesh.
!#
!# 14.) tria_quickRefine2LevelOrdering
!#      -> Refines a mesh multiple times according to the 2-level ordering 
!#         algorithm. Creates a 'raw' fine mesh from a 'raw' or 'standard' 
!#         coarse mesh.
!#
!# 15.) tria_rawGridToTri
!#      -> Converts a raw mesh into a triangular mesh.
!#
!# 16.) tria_infoStatistics
!#      -> Prints out statistics about a mesh.
!#
!# 17.) tria_exportTriFile
!#      -> Exports a triangulation structure to a .TRI file.
!#
!# 18.) tria_exportPostScript
!#      -> Exports a 2D triangulation to an (encapsulated) PostScript file
!#         (EPS) which e.g. can be imported into a LaTeX document using
!#         the \includegraphics{} command.
!#
!# 19.) tria_searchBoundaryNode
!#      -> Search for the position of a boundary vertex / edge on the boundary.
!#
!# 20.) tria_generateSubdomain
!#      -> Extract cells from a mesh and generate a subdomain from them
!#
!# 21.) tria_attachCells
!#      -> Attaches a set of cells to an existing triangulation
!#
!# 22.) tria_cellGroupGreedy
!#      -> Combines the cells of a mesh into simply connected sets of 
!#         similar size by a greedy algorithm
!#
!# 23.) tria_getNAE = tria_getNAE_direct /
!#                    tria_getNAE_indirect
!#      -> Get the number of faces on an element
!#
!# 24.) tria_getNVE = tria_getNVE_direct / 
!#                    tria_getNVE_indirect
!#      -> Get the number of vertices/edges on an element
!#
!# 25.) tria_getPointsOnEdge
!#      -> For all edges in a triangulation, calculate the coordinates 
!#         of a number of points on each edge
!#
!# 26.) tria_getNeighbourVertex
!#      -> Calculates the vertex number of the neighbour vertex of a 
!#         vertex on an edge.
!#
!# 27.) tria_getSubmeshNeighbourhood
!#      -> Calculate a list of cells in the neighbourhood of a cell set
!#
!# 28.) tria_getElementsAtMacroEdge
!#      -> Calculate a list of all elements adjacent to a macro edge
!#         inside of a macro cell.
!#
!# 29.) tria_getElementsAtMacroVertex
!#      -> Calculate a list of all elements adjacent to a vertex
!#         inside of a macro cell.
!#
!# 30.) tria_getVerticesAtFaceDirect
!#      -> Calculates the vertices on a face in 3D
!#
!# Auxiliary routines:
!#
!#  1.) tria_readRawTriangulation1D /
!#      tria_readRawTriangulation2D /
!#      tria_readRawTriangulation3D
!#      -> Reads basic data from a TRI file.
!#
!#  2.) tria_genRawBoundary1D /
!#      tria_genRawBoundary2D /
!#      tria_genRawBoundary3D
!#     -> Generates basic information about boundary vertices.
!#
!#  3.) tria_genElementsAtVertex1D2D /
!#      tria_genElementsAtVertex3D
!#      -> Generates the IelementsAtVertex array.
!#
!#  4.) tria_genNeighboursAtElement1D /
!#      tria_genNeighboursAtElement2D /
!#      tria_genNeighboursAtElement3D
!#      -> Generates the IneighboursAtElement array.
!#
!#  5.) tria_genEdgesAtElement2D /
!#      tria_genEdgesAtElement3D
!#      -> Generates the IedgesAtElement array.
!#
!#  6.) tria_genFacesAtElement3D
!#      -> Generate the arrays with the faces at an element
!#
!#  7.) tria_genElementVolume1D /
!#      tria_genElementVolume2D
!#      -> Generates the DelementVolume array.
!#
!#  8.) tria_sortBoundaryVertices1D2D
!#      -> Sorts the boundary vertices for increasing parameter values.
!#
!#  9.) tria_genElementsAtBoundary1D2D
!#      -> Generates the IelementsAtBoundary array.
!#
!# 10.) tria_genBoundaryVertexPos1D2D
!#      -> Generates the IboundaryVertexPos2D array.
!#
!# 11.) tria_genElementsAtEdge2D /
!#      tria_genElementsAtEdge3D
!#      -> Generates the IelementsAtEdge array.
!#
!# 12.) tria_genVerticesAtEdge2D /
!#      tria_genVerticesAtEdge3D
!#      -> Generates the IverticesAtEdge array.
!#
!# 13.) tria_genEdgeNodalProperty2D /
!#      tria_genEdgeNodalProperty3D
!#      -> Generates the InodalProperty array part for all edges
!#
!# 14.) tria_genFaceNodalProperty3D
!#      -> Generate the face nodal property in 3D
!#
!# 15.) tria_genEdgesAtBoundary2D /
!#      tria_genEdgesAtBoundary3D
!#      -> Generates the IedgesAtBoundary array for a 2D triangulation
!#
!# 16.) tria_genEdgeParameterValue2D
!#      -> Generates the DedgeParameterValue array for a 2D triangulatioon
!#
!# 17.) tria_genBoundaryEdgePos2D
!#      -> Generates the IboundaryEdgePos2D array for a 2D triangulatioon
!#
!# 18.) tria_genEdgesAtVertex2D
!#      -> Generates the array with the edge numbers at a vertex
!#
!# 19.) tria_genFacesAtBoundary3D
!#      -> Generate the array with face numbers on the boundary
!#
!# 20.) tria_genFacesAtVertex3D
!#      -> Generate the arrays with the faces at a vertex
!#
!# 21.) tria_genFacesAtEdge3D
!#      -> Generate the arrays with the faces at an edge
!#
!# 22.) tria_genEdgesAtFace3D
!#      -> Generate the arrays with the edges at a face
!#
!# 23.) tria_genElementsAtFace3D
!#      -> Generate the arrays with the elements at a face
!#
!# 24.) tria_genVerticesAtFace3D
!#      -> Generate the arrays with the vertices at a face
!#
!# 25.) tria_genTwistIndex
!#      -> Generates the twist index array
!#
!# 26.) tria_propMacroNodalProperty2lv
!#      -> Propagates a macro nodal property array from a coarse mesh
!#         to a fine mesh
!#
!# 27.) tria_genRefTags2lv
!#      -> Calculates refinement tags according to the uniform 2-level refinement
!#
!# 28.) tria_hangingNodeRefinement
!#      -> Performs a local hanging-node refinement of quadrilateral elements
!#
!# 29.) tria_calcBoundingBox
!#      -> Calculates a bounding box around a mesh
!#
!# Auxiliary routines for connector lists
!#
!#  1.) tria_buildConnectorList
!#      -> Builds the connector list used in the neighbours at elements routine
!#
!#  2.) tria_sortElements3D
!#      -> Calculates the lexicographical ordering for the connector list
!#
!#  3.) tria_sortElements3DInt
!#      -> Calculates the sorted numbering for the connector list
!#  
!#  4.) tria_mergesort
!#      -> Sorts a connector list
!#
!#  5.) tria_merge
!#      -> Auxiliary routine for merge sort
!#
!#
!#  FAQ - Some explainations  \\
!# -------------------------- \\
!# 1.) When reading the refinement routine, there is written something about
!#     'raw' meshes and 'standard' meshes. What does that mean?
!#
!#     A 'raw' mesh is a very basic form of a mesh with most information
!#     missing. E.g. there is no information about adjacencies.
!#     Such meshes come e.g. from .TRI files when reading them with
!#     tria_readTriFile2D or similar routines. These basic meshes can normally
!#     not be used for any computations; all the missing informationhas first
!#     to be 'extracted' or 'generated' based on them. They can only be used
!#     for very low-level modifications; e.g. the routine tria_rawGridToTri
!#     allows to convert a quad mesh in a triangular mesh, which would
!#     be much harder if all adjacency information is already computed.
!#     Another possibile thing what can be done with such a 'raw' mesh is to
!#     do quicker pre-refinement with routines like 
!#     tria_quickRefine2LevelOrdering. This routine refines the mesh without
!#     computing everything and is therefore a little bit faster than a
!#     subsequent application of tria_refine2LevelOrdering onto a 'standard'
!#     mesh.
!#
!#     If you have a 'raw' mesh, simply use tria_initStandardMeshFromRaw
!#     to convert it to a 'standard' mesh. A 'standard' mesh contains all
!#     neighbouring information and is usually the mesh you want to use
!#     for computations.
!#
!# 2.) But there is something written in the tria_readTriFileXd routines
!#     about 'extended' raw meshes? What is with that?
!#
!#     Well, the tria_readTriFileXd routines actually produce a so called
!#     'extended' raw mesh which is a mesh that contains a numbering for
!#     vertices, edges and (in 3D) faces. The generation of this extended
!#     information can be prevented by parameter, in which case the
!#     raw mesh provides a numbering only for vertices. This is usually
!#     not advisable as some refinement routines may rely on edge- and
!#     face numbering.
!#     Nevertheless, tria_initStandardMeshFromRaw will initialise even
!#     a not-extended raw mesh to a standard mesh with all information.
!#
!# 3.) And what does that mean if I want to read e.g. a 2D mesh and to refine it?
!#
!#     Well, to read a mesh and prepare it to be used as single mesh, use:
!#
!# <code>
!#       call tria_readTriFile2D (rtriangulation, 'somemesh.tri', rboundary)
!#       call tria_initStandardMeshFromRaw (rtria,rboundary)
!# </code>
!#
!#     If you want to refine a mesh let us say 4 times after reading it, use:
!#
!# <code>
!#       call tria_readTriFile2D (rtriangulation, 'somemesh.tri', rboundary)
!#       call tria_quickRefine2LevelOrdering(4,rtriangulation,rboundary)
!#       call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
!# </code>
!# 
!#     If you want to generate level 1..4 for a multiple-grid structure, use
!#
!# <code>
!#       type(t_triangulation), dimension(4) :: rtriangulation
!#
!#       call tria_readTriFile2D (rtriangulation(1), 'somemesh.tri', rboundary)
!#       call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
!#       do ilev=2,4
!#         call tria_refine2LevelOrdering (rtriangulation(ilev-1),&
!#             rtriangulation(ilev), rboundary)
!#         call tria_initStandardMeshFromRaw (rtriangulation(ilev),rboundary)
!#       end do
!# </code>
!#
!# 4.) What is the tria_compress2LevelOrdHierarchy for?
!#
!#     This routine can be used after the refinement process to save memory.
!#     As the vertex coordinates of all points in the coarse mesh are the
!#     same in the fine mesh, there is redundant data! The
!#     tria_compress2LevelOrdHierarchy now removes this redundant data by
!#     releasing the vertex coordinate array on the coarse mesh and replacing
!#     it by that of the fine mesh. The routine can be used as follows:
!#
!#     Generate the triangulation:
!#
!# <code>
!#       type(t_triangulation), dimension(4) :: rtriangulation
!#
!#       call tria_readTriFile2D (rtriangulation(1), 'somemesh.tri', rboundary)
!#       call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
!#       do ilev=2,4
!#         call tria_refine2LevelOrdering (rtriangulation(ilev-1),&
!#             rtriangulation(ilev), rboundary)
!#         call tria_initStandardMeshFromRaw (rtriangulation(ilev),rboundary)
!#       end do
!# </code>
!#
!#    Remove redundant data
!#
!# <code>
!#       do ilev=4-1,1
!#         call tria_compress2LevelOrdHierarchy (rtriangulation(ilev+1),&
!#              rtriangulation(ilev))
!#       end do
!# </code>
!#
!#    One should note that afterwards, each change of vertex coordinates on the
!#    fine grid directly affects the vertex coordinates on all coarse grid
!#    and vice versa. But that is not a bug, it is a feature :-)
!# 
!# 5.) What are twist indices?
!#
!#     !!! DEPRECATED !!!
!#     This documentation will be re-written later as it does not match
!#     the 'new' implementation of the twist indices!
!#
!#   Twist indices define a global orientation of edges and faces in a mesh.
!#   There are two different fields for this in the triangulation: ItwistIndexEdges
!#   and ItwistIndexFaces.
!#
!#   2D: Think about an edge adjacent to two elements:
!#
!# <verb>
!#     ------+------
!#           |
!#           |
!#      10   E  20
!#           |
!#           |
!#     ------+------
!# </verb>
!#
!#   We define the edge E to belong to the element with the smaller number.
!#   When 'looking' from element 10 to the edge, the edge becomes twist
!#   index 1, while when looking from element 20 to it, the edge receives
!#   twist index 0. That way, the twist index defines whether the 'global' edge
!#   is oriented counterclockwise (twist index = 1) or clockwise 
!#   (twist index = 0) relative to the current element.
!#
!#   Each cell in 2D has 4 edges. The entry ItwistIndexEdges(iel) defines a
!#   bitfield where the corresponding bit is set if the twist index is 1. So:
!# <verb>
!#    Bit 0 = twist index of edge 1
!#    Bit 1 = twist index of edge 2
!#    Bit 2 = twist index of edge 3
!#    ...
!# </verb>
!#
!#   3D: In 3D the situation is slightly more complicated. Here we have twist indices
!#   for edges as well as for faces.
!#
!#   The twist indices for edges are defined relative to the face.
!#   ItwistIndexEdges(iel) is again a bitfield that defines for each edge on
!#   each face of the element its orientation. For simplicity, we always assume
!#   that there are up to 4 edges on each face and up to 6 faces on an element
!#   (which is the case for hexahedrals). Then the bitfield looks like this:
!#
!# <verb>
!#     Bit         | 31..24 | 23..20 | 19..16 | 15..12 | 11..8  | 7..4   | 3..0
!#     ------------+--------+--------+--------+--------+--------+--------+--------
!#     Twist index | undef. | Face 6 | Face 5 | Face 4 | Face 3 | Face 2 | Face 1
!# </verb>
!#
!#   Furthermore, the array ItwistIndexFaces defines for each element the orientation
!#   of the faces. This is a 2D array of the form ItwistIndexFaces(iface,ielement).
!#   To explain the meaning, we show here an example of two adjacent elements
!#   with their local vertex numbering. (To get a better overview, we depict the two
!#   elements separately from each other; the face should be thought of being
!#   shared by the elements:)
!#
!# <verb>
!#                         
!#         ----------3              2-------
!#                  /|             /|
!#                 / |            / |
!#                /  |           /  |
!#               /   |          /   |     20
!#      --------4    |         3--------
!#              | I  |         |  J |
!#        10    |    2 <-------|--  1------
!#              |   /     =-2  |   / 
!#              |  /           |  /  
!#              | /            | /   
!#              |/             |/    
!#      --------1  --------->  4-------
!#                   =4
!#
!# </verb>
!#
!#   Every face has it is own local numbering, which is independent of the
!#   local numbering of the vertices on the cell: While the cell has local
!#   vertex numbers 1..8, the face has local vertex numbers 1..4 (and
!#   there is also a mapping between the local face numbers and the local
!#   vertex numbers, but this is not used). Again we define the face to
!#   belong to the element with the smaller element number. This element
!#   defines the 'global orientation' of the face and is by definition
!#   oriented anticlockwise.
!#   The value of the twist index now defines, which local vertex number
!#   on the face of the neighbour element belongs to local vertex number 1
!#   on the current face. In the above example, ItwistIndexFaces(I,10)=4,
!#   because local vertex number 1 on face I of element 10 maps to local
!#   vertex number 4 of face J on element 20.
!#   Similarly, ItwistIndexFaces(J,20)=-2, because local vertex number 1 of face
!#   J maps to local vertex number 2 of face I. the sign "-" signales that
!#   face J is oriented clockwise in the global orientation, as it belongs
!#   to element 10.
!#
!# 6.) What is a subdomain and how to work with it?
!#
!#   A subdomain is just a normal triangulation with the difference, that it
!#   has up to NBCT+1 boundary components. It can be generated by extracting
!#   cells from an existing triangulation using tria_generateSubdomain.
!# 
!#   The first NBCT bondary components refer to the original BC`s in the 
!#   source triangulation. The 'last' BC receives information about which 
!#   parts of the domain in rtriangulation became the boundary in the 
!#   subdomain. This last boundary component is a 'blind' boundary 
!#   component: The arrays (like IboundaryCpIdx) contain information 
!#   about one boundary component which is not counted in NBCT. This last 
!#   boundary component can even be empty! The variable NblindBCT is <> 0
!#   if there may be elements in a blind boundary component in a mesh.
!#
!#   In 2D, when edges or faces are generated for a subdomain, these may have
!#   a nodal property indicating them to be located on this blind boundary
!#   component although their position in the IedgesAtBoundary array
!#   would indicate them to not be in the blind boundary component.
!#   That is because in 2D, an edge follows every vertex. Edges in 
!#   IedgesAtBoundary on the blind boundary component can be identified
!#   either by their nodal property InodalProperty or by there parameter
!#   value, as DedgeParameterValue(iedge)=-1 for such an edge.
!#
!# 7.) What is this macro nodal property?
!#
!#   The macro nodal property array is a possibility to have information
!#   about a coarse mesh available on a fine mesh. It allows for vertices,
!#   edges, faces and elements to determine the coarse grid vertex,
!#   edge, face and element where it comes from. This array is not
!#   generated by default and must explicitly be calculated by the
!#   user calling tria_initMacroNodalProperty for a given (extended raw) mesh.
!#   Roughly said, tria_initMacroNodalProperty marks a mesh as coarse mesh.
!#   As soon as a mesh is marked as coarse mesh by tria_initMacroNodalProperty,
!#   all derived fine meshes (standard as well as extended raw meshes)
!#   will automatically have their own macro nodal property array attached
!#   when they are created by refinement. 
!#
!#   Example: The following command sequence reads a tri file,
!#   attaches a macro nodal property array to the triangulation and
!#   refines 4 times. The mesh at level 5 will also have a macro nodal
!#   property array attached that allows to determine the origin of
!#   a fine grid vertex, edge, face or element.
!#
!# <code>
!#       call tria_readTriFile2D (rtriangulation, 'somemesh.tri', rboundary)
!#       call tria_initMacroNodalProperty (rtriangulation)
!#       call tria_quickRefine2LevelOrdering(4,rtriangulation,rboundary)
!# </code>
!#
!#   Note that the macro nodal property array is part of the extended raw
!#   mesh as soon as tria_initMacroNodalProperty is called. It is not necessary
!#   to call tria_initStandardMeshFromRaw in order to have this array 
!#   available!
!#
!# </purpose>
!##############################################################################

module triangulation

  use fsystem
  use genoutput
  use storage
  use basicgeometry
  use linearalgebra
  use boundary
  use geometryaux
  use sort
  use io

  implicit none
  
  private
  
!<constants>

!<constantblock description="Triangulation constants">
  
  ! Maximum number of vertices per element supported by the triangulation
  ! module. Currently, this is 8 which is the case for 3D hexahedral meshes.
  integer, parameter, public :: TRIA_MAXNVE   = 8
  
  ! Maximum number of corner-vertices in each 1D element.
  integer, parameter, public :: TRIA_MAXNVE1D = 2

  ! Maximum number of corner-vertices in each 2D element.
  ! We set this to 4 to allow triangle and quadrilateral element shapes.
  ! May be set to higher vales in the future for the support of
  ! isoparametric elements! This is the old NNVE.
  integer, parameter, public :: TRIA_MAXNVE2D = 4

  ! Maximum number of corner-vertices in each 2D element.
  integer, parameter, public :: TRIA_MAXNVE3D = 8

  
  ! Maximum number of edges in each 2D element.
  ! We set this to 4 to allow triangle and quadrilateral element shapes.
  ! This is the old NNVE, too.
  integer, parameter, public :: TRIA_MAXNME2D = TRIA_MAXNVE2D

  ! Maximum number of faces per element supported by the triangulation
  ! module. Currently, this is 6 which is the case for 3D hexahedral meshes.
  integer, parameter, public :: TRIA_MAXNAE = 6

  ! Number of elements of a 3D connector
  integer, parameter, public :: TRIA_NCONNECT3D = 5


  ! Number of vertices per line element in 1D.
  integer, parameter, public :: TRIA_NVELINE1D = 2
  
  ! Number of vertices per element for triangular element shapes in 2D.
  integer, parameter, public :: TRIA_NVETRI2D  = 3

  ! Number of vertices per element for quadrilateral element shapes in 2D.
  integer, parameter, public :: TRIA_NVEQUAD2D = 4

  ! Number of vertices per element for tetrahedral element shapes in 3D.
  integer, parameter, public :: TRIA_NVETET3D  = 4

  ! Number of vertices per element for pyramidal element shapes in 3D.
  integer, parameter, public :: TRIA_NVEPYR3D  = 5

  ! Number of vertices per element for prismatic element shapes in 3D.
  integer, parameter, public :: TRIA_NVEPRIS3D = 6

  ! Number of vertices per element for hexahedral element shapes in 3D.
  integer, parameter, public :: TRIA_NVEHEXA3D = 8

  
  ! Number of edges per element for triangular element shapes in 2D.
  integer, parameter, public :: TRIA_NNETRI2D = 3

  ! Number of edges per element for quadrilateral element shapes in 2D.
  integer, parameter, public :: TRIA_NNEQUAD2D = 4

  ! Number of edges per element for tetrahedral element shapes in 3D.
  integer, parameter, public :: TRIA_NNETET3D = 6

  ! Number of edges per element for pyramidal element shapes in 3D.
  integer, parameter, public :: TRIA_NNEPYR3D = 8

  ! Number of edges per element for prismatic element shapes in 3D.
  integer, parameter, public :: TRIA_NNEPRIS3D = 9
  
  ! Number of edges per element for hexahedral element shapes in 3D.
  integer, parameter, public :: TRIA_NNEHEXA3D = 12


  ! Number of faces per element for tetrahedral element shapes in 3D.
  integer, parameter, public :: TRIA_NAETET3D = 4

  ! Number of faces per element for pyramidal element shapes in 3D.
  integer, parameter, public :: TRIA_NAEPYR3D = 5

  ! Number of faces per element for prismatic element shapes in 3D.
  integer, parameter, public :: TRIA_NAEPRIS3D = 5

  ! Number of faces per element for hexahedral element shapes in 3D.
  integer, parameter, public :: TRIA_NAEHEXA3D = 6
!</constantblock>

!<constantblock description="KIND values for triangulation data">

  ! Remark:
  ! All PREC_* constants are deprecated. Do not use them anymore as they
  ! will be removed in future!
  
  ! DEPRECATED: kind value for indexing the vertices in a triangulation
  integer, parameter, public :: PREC_VERTEXIDX  = I32

  ! DEPRECATED: Alternative name for PREC_VERTEXIDX
  integer, parameter, public :: PREC_POINTIDX   = PREC_VERTEXIDX

  ! DEPRECATED: kind value for indexing the edges in a triangulation
  integer, parameter, public :: PREC_EDGEIDX    = I32
  
  ! DEPRECATED: kind value for indexing the faces in a triangulation
  integer, parameter, public :: PREC_FACEIDX    = I32

  ! DEPRECATED: kind value for indexing the elements in a triangulation
  integer, parameter, public :: PREC_ELEMENTIDX = I32

!</constantblock>
  
!<constantblock description="Flags to be specified as cflags in the refinement routines.">
  
  ! After refinement, those points of quad elements on the fine mesh which 
  ! were formally element midpoints on the coarse mesh were recalculated
  ! by taking the mean of the corners.
  ! This helps avoiding tangled elements when there is a 'hole' in the domain.
  integer(I32), parameter, public :: TRIA_R2LV_AVERAGEMIDPOINTS  = 2**0

  ! After refinement, the coordinates of points on the boundary are
  ! recalculated using their parameter values (only 2D)
  integer(I32), parameter, public :: TRIA_R2LV_RECALCCOORDSONBD  = 2**1

  ! After refinement, do not produce an extended raw mesh, just produce
  ! a simple raw mesh.
  integer(I32), parameter, public :: TRIA_R2LV_NOEXTENDEDRAW     = 2**2

  ! Standard parameter settings for 2-level refinement
  integer(I32), parameter, public :: TRIA_R2LV_STANDARD = TRIA_R2LV_AVERAGEMIDPOINTS + &
                                                  TRIA_R2LV_RECALCCOORDSONBD

!</constantblock>
  
!<constantblock description="Duplication flags. Specifies which information is shared \
!                            between triangulation structures">

  integer(I32), parameter, public :: TR_SHARE_DVERTEXCOORDS          = 2** 0  ! DCORVG 
  integer(I32), parameter, public :: TR_SHARE_DFREEVERTEXCOORDINATES = 2** 1  ! DCORMG
  integer(I32), parameter, public :: TR_SHARE_IVERTICESATELEMENT     = 2** 2  ! KVERT 
  integer(I32), parameter, public :: TR_SHARE_IEDGESATELEMENT        = 2** 3  ! KMID  
  integer(I32), parameter, public :: TR_SHARE_INEIGHBOURSATELEMENT   = 2** 4  ! KADJ  
  integer(I32), parameter, public :: TR_SHARE_IELEMENTSATVERTEX      = 2** 5  ! KVEL  
  integer(I32), parameter, public :: TR_SHARE_IELEMENTSATEDGE        = 2** 6  ! KMEL  
  integer(I32), parameter, public :: TR_SHARE_INODALPROPERTY         = 2** 7  ! KNPR  
  integer(I32), parameter, public :: TR_SHARE_IVERTICESATBOUNDARY    = 2** 8  ! KVBD  
  integer(I32), parameter, public :: TR_SHARE_IELEMENTSATBOUNDARY    = 2** 9  ! KEBD  
  integer(I32), parameter, public :: TR_SHARE_IBOUNDARYCPIDX         = 2**10  ! KBCT  
  integer(I32), parameter, public :: TR_SHARE_DVERTEXPARAMETERVALUE  = 2**11  ! DVBDP 
  integer(I32), parameter, public :: TR_SHARE_DEDGEPARAMETERVALUE    = 2**12  ! DMBDP 
  integer(I32), parameter, public :: TR_SHARE_IEDGESATBOUNDARY       = 2**13  ! KMBD  
  integer(I32), parameter, public :: TR_SHARE_IVERTICESATEDGE        = 2**14  ! KEAN  
  integer(I32), parameter, public :: TR_SHARE_IBOUNDARYVERTEXPOS     = 2**15  ! KVBDI 
  integer(I32), parameter, public :: TR_SHARE_IBOUNDARYEDGEPOS       = 2**16  ! KMBDI 
  integer(I32), parameter, public :: TR_SHARE_DELEMENTAREA           = 2**17  ! DAREA 
  integer(I32), parameter, public :: TR_SHARE_IREFINEMENTPATCH       = 2**18   
  integer(I32), parameter, public :: TR_SHARE_ICOARSEGRIDELEMENT     = 2**19  
  
  integer(I32), parameter, public :: TR_SHARE_IVERTICESATFACE        = 2**20  ! KVAR 
  integer(I32), parameter, public :: TR_SHARE_IFACESATELEMENT        = 2**21  ! KAREA
  integer(I32), parameter, public :: TR_SHARE_IELEMENTSATFACE        = 2**22
  integer(I32), parameter, public :: TR_SHARE_IEDGESATFACE           = 2**23 
  integer(I32), parameter, public :: TR_SHARE_IFACESATEDGE           = 2**24
  integer(I32), parameter, public :: TR_SHARE_IFACESATVERTEX         = 2**25
  integer(I32), parameter, public :: TR_SHARE_IFACESATBOUNDARY       = 2**26
  
  integer(I32), parameter, public :: TR_SHARE_IEDGESATVERTEX         = 2**27
  integer(I32), parameter, public :: TR_SHARE_ITWISTINDEX            = 2**28
  
  integer(I32), parameter, public :: TR_SHARE_IMACRONODALPROPERTY    = 2**29  
  
  ! Share information for an extended raw mesh.
  integer(I32), parameter, public :: TR_SHARE_EXTENDEDRAW            = &
            TR_SHARE_IELEMENTSATVERTEX + TR_SHARE_INEIGHBOURSATELEMENT + &
            TR_SHARE_IEDGESATELEMENT + TR_SHARE_IFACESATELEMENT
  
  ! Share everything
  integer(I32), parameter, public :: TR_SHARE_ALL = not(0_I32)

  ! Share nothing
  integer(I32), parameter, public :: TR_SHARE_NONE = 0_I32

!</constantblock>

!<constantblock description="Generation flags. Specifies which information should \
!                            be generated in the standard mesh">

  integer(I32), parameter, public :: TR_GEN_DVERTEXCOORDS          = TR_SHARE_DVERTEXCOORDS
  integer(I32), parameter, public :: TR_GEN_DFREEVERTEXCOORDINATES = TR_SHARE_DFREEVERTEXCOORDINATES
  integer(I32), parameter, public :: TR_GEN_IVERTICESATELEMENT     = TR_SHARE_IVERTICESATELEMENT
  integer(I32), parameter, public :: TR_GEN_IEDGESATELEMENT        = TR_SHARE_IEDGESATELEMENT
  integer(I32), parameter, public :: TR_GEN_INEIGHBOURSATELEMENT   = TR_SHARE_INEIGHBOURSATELEMENT
  integer(I32), parameter, public :: TR_GEN_IELEMENTSATVERTEX      = TR_SHARE_IELEMENTSATVERTEX
  integer(I32), parameter, public :: TR_GEN_IELEMENTSATEDGE        = TR_SHARE_IELEMENTSATEDGE
  integer(I32), parameter, public :: TR_GEN_INODALPROPERTY         = TR_SHARE_INODALPROPERTY
  integer(I32), parameter, public :: TR_GEN_IVERTICESATBOUNDARY    = TR_SHARE_IVERTICESATBOUNDARY
  integer(I32), parameter, public :: TR_GEN_IELEMENTSATBOUNDARY    = TR_SHARE_IELEMENTSATBOUNDARY
  integer(I32), parameter, public :: TR_GEN_IBOUNDARYCPIDX         = TR_SHARE_IBOUNDARYCPIDX
  integer(I32), parameter, public :: TR_GEN_DVERTEXPARAMETERVALUE  = TR_SHARE_DVERTEXPARAMETERVALUE
  integer(I32), parameter, public :: TR_GEN_DEDGEPARAMETERVALUE    = TR_SHARE_DEDGEPARAMETERVALUE
  integer(I32), parameter, public :: TR_GEN_IEDGESATBOUNDARY       = TR_SHARE_IEDGESATBOUNDARY
  integer(I32), parameter, public :: TR_GEN_IVERTICESATEDGE        = TR_SHARE_IVERTICESATEDGE
  integer(I32), parameter, public :: TR_GEN_IBOUNDARYVERTEXPOS     = TR_SHARE_IBOUNDARYVERTEXPOS
  integer(I32), parameter, public :: TR_GEN_IBOUNDARYEDGEPOS       = TR_SHARE_IBOUNDARYEDGEPOS
  integer(I32), parameter, public :: TR_GEN_DELEMENTAREA           = TR_SHARE_DELEMENTAREA
  integer(I32), parameter, public :: TR_GEN_IREFINEMENTPATCH       = TR_SHARE_IREFINEMENTPATCH
  integer(I32), parameter, public :: TR_GEN_ICOARSEGRIDELEMENT     = TR_SHARE_ICOARSEGRIDELEMENT
  
  integer(I32), parameter, public :: TR_GEN_IVERTICESATFACE        = TR_SHARE_IVERTICESATFACE
  integer(I32), parameter, public :: TR_GEN_IFACESATELEMENT        = TR_SHARE_IFACESATELEMENT
  integer(I32), parameter, public :: TR_GEN_IELEMENTSATFACE        = TR_SHARE_IELEMENTSATFACE
  integer(I32), parameter, public :: TR_GEN_IEDGESATFACE           = TR_SHARE_IEDGESATFACE
  integer(I32), parameter, public :: TR_GEN_IFACESATEDGE           = TR_SHARE_IFACESATEDGE
  integer(I32), parameter, public :: TR_GEN_IFACESATVERTEX         = TR_SHARE_IFACESATVERTEX
  integer(I32), parameter, public :: TR_GEN_IFACESATBOUNDARY       = TR_SHARE_IFACESATBOUNDARY
  
  integer(I32), parameter, public :: TR_GEN_IEDGESATVERTEX         = TR_SHARE_IEDGESATVERTEX
  integer(I32), parameter, public :: TR_GEN_ITWISTINDEX            = TR_SHARE_ITWISTINDEX
  
  ! Generate information for an extended raw mesh.
  integer(I32), parameter, public :: TR_GEN_EXTENDEDRAW            = TR_SHARE_EXTENDEDRAW
  
  ! Generate everything
  integer(I32), parameter, public :: TR_GEN_ALL = not(0_I32)

!</constantblock>
  
!<constantblock description="Format tags for TRI file formats.">

  ! Standard TRI file format, compatible to FEAT1.
  ! For 2D triangulations, Vertex coordinates of boundary vertices are
  ! replaced by parameter values.
  integer(I32), parameter, public :: TRI_FMT_STANDARD          = 0
  
  ! Standard TRI file format, but the vertex coordinates are 
  ! exported 'as they are', not as parameter values.
  integer(I32), parameter, public :: TRI_FMT_NOPARAMETRISATION = 2**0
  
!</constantblock>
  
!<constantblock description="Neighbourhood specifiers for tria_getSubmeshNeighbourhood.">
  
!</constantblock>
  
  ! All elements adjacent by vertices.
  integer(i32), parameter, public :: TRI_NEIGH_VERTEXNEIGHBOURS         = 2**0
  
  ! All elements adjacent by edges (only 2D and 3D).
  integer(i32), parameter, public :: TRI_NEIGH_EDGENEIGHBOURS           = 2**1
  
  ! All elements adjacent by faces (only 3D).
  integer(i32), parameter, public :: TRI_NEIGH_FACENEIGHBOURS           = 2**2
  
  ! All adjacent elements
  integer(i32), parameter, public :: TRI_NEIGH_ALL = TRI_NEIGH_VERTEXNEIGHBOURS + &
                                             TRI_NEIGH_EDGENEIGHBOURS + &
                                             TRI_NEIGH_FACENEIGHBOURS
  
!</constants>


!<types>

!<typeblock>

  ! Each 1D element consitst of at most TRIA_MAXNVE1D points.
  ! Each point has a number, which is usually an integer value.
  type t_elementCorners1D
    integer, dimension(TRIA_MAXNVE1D) :: Icorners
  end type t_elementCorners1D
  
  public :: t_elementCorners1D
  
  ! Each 2D element consists of at most TRIA_MAXNVE2D points.
  ! Each point has a number, which is usually an integer value.
  type t_elementCorners2D
    integer, dimension(TRIA_MAXNVE2D) :: Icorners
  end type t_elementCorners2D

  public :: t_elementCorners2D

  ! Each 2D element contains at most TRIA_MAXNME2D edges.
  ! Each edge has a number, which is usually an integer value.
  type t_elementEdges2D
    integer, dimension(TRIA_MAXNME2D) :: Iedges
  end type t_elementEdges2D
  
  public :: t_elementEdges2D
  
  ! Each 1D element contains at most TRIA_MAXNVE1D neighbour elements,
  ! each meeting the element in a vertice.
  ! Each neighbour element has a number, which is usually an integer value.
  type t_elementNeighbours1D
    integer, dimension(TRIA_MAXNVE1D) :: Ineighbours
  end type t_elementNeighbours1D
  
  public :: t_elementNeighbours1D
  
  ! Each 2D element contains at most NMAXEDGES neighbour elements,
  ! each meeting the element in an edge.
  ! Each neighbour element has a number, which is usually an integer value.
  type t_elementNeighbours2D
    integer, dimension(TRIA_MAXNME2D) :: Ineighbours
  end type t_elementNeighbours2D
  
  public :: t_elementNeighbours2D
  
!<typeblock>
  
!</typeblock>
  
  ! The basic triangulation structure for a triangulation.
  type t_triangulation
  
    ! Duplication flag. Bitfield. Used by TRIDUP/TRIRST/TRIDEL to 
    ! mark which information of a triangulation structure 
    ! coincides with those of another triangulation structure and
    ! must not be released from memory when deleting a copy of a
    ! triangulation structure.
    ! When a bit is set to 1, the corresponding array is
    ! maintained by another triangulation structure and must not
    ! be deleted by TRIDEL. When the bit is 0, the array is a real
    ! copy of another array and must be deleted in TRIDEL.
    ! <verb>
    ! Bit  0: DvertexCoords          is a copy of another structure (DCORVG)
    ! Bit  1: DfreeVertexCoordinates is a copy of another structure (DCORMG)
    ! Bit  2: IverticesAtElement     is a copy of another structure (KVERT)
    ! Bit  3: IedgesAtElement        is a copy of another structure (KMID)
    ! Bit  4: IneighboursAtElement   is a copy of another structure (KADJ)
    ! Bit  5: IelementsAtVertex +
    !         IelementsAtVertexIdx   is a copy of another structure (KVEL)
    ! Bit  6: IelementsAtEdge +
    !         IelementsAtEdge3D +
    !         IelementsAtEdgeIdx3D   is a copy of another structure (KMEL)
    ! Bit  7: InodalProperty         is a copy of another structure (KNPR)
    ! Bit  8: IverticesAtBoundary    is a copy of another structure (KVBD)
    ! Bit  9: IelementsAtBoundary    is a copy of another structure (KEBD)
    ! Bit 10: IboundaryCpIdx         is a copy of another structure (KBCT)
    ! Bit 11: DvertexParameterValue  is a copy of another structure (DVBDP)
    ! Bit 12: DedgeParameterValue    is a copy of another structure (DMBDP)
    ! Bit 13: IedgesAtBoundary       is a copy of another structure (KMBD)
    ! Bit 14: IverticesAtEdge        is a copy of another structure (KEAN)
    ! Bit 15: IboundaryVertexPos     is a copy of another structure (KVBDI)
    ! Bit 16: IboundaryEdgePos       is a copy of another structure (KMBDI)
    ! Bit 17: DelementArea           is a copy of another structure (DAREA)
    ! Bit 18: IrefinementPatch +
    !         IrefinementPatchIdx    is a copy of another structure
    ! Bit 19: IcoarseGridElement     is a copy of another structure
    ! Bit 20: IverticesAtFace        is a copy of another structure (KVAR)
    ! Bit 21: IfacesAtElement        is a copy of another structure (KAREA)
    ! Bit 22: IelementsAtFace        is a copy of another structure
    ! Bit 23: IedgesAtFace           is a copy of another structure
    ! Bit 24: IfacesAtEdge +
    !         IfacesAtEdgeIdx        is a copy of another structure
    ! Bit 25: IfacesAtVertex +
    !         IfacesAtVertexIdx      is a copy of another structure
    ! Bit 26: IfacesAtBoundary +
    !         IboundaryCpFacesIdx    is a copy of another structure
    ! Bit 27: IedgesAtVertex +
    !          IedgesAtVertexIdx     is a copy of another structure
    ! Bit 28: ItwistIndexEdges +
    !         ItwistIndexFaces       is a copy of another structure
    ! Bit 29: ImacroNodalProperty    is a copy of another structure
    ! </verb>
    integer(I32)             :: iduplicationFlag = 0
  
    ! Dimension of the triangulation.
    ! NDIM1D=1D triangulation, NDIM2D=2D triangulation
    ! NDIM2D=3D triangulation, 0=not initialised
    integer                  :: ndim = 0
  
    ! Number of points in the domain corresponding to corners of elements;
    ! coincides with SIZE(RcornerCoordinates)
    integer                  :: NVT = 0
    
    ! Number of edges in the domain belonging to elements
    integer                  :: NMT = 0
    
    ! total number of faces in the domain
    integer                  :: NAT = 0
    
    
    ! Number of elements in the domain; 
    ! corresponding to SIZE(RverticesOnElement)
    integer                  :: NEL = 0
    
    ! Number of boundary components
    integer                  :: NBCT = 0
    
    ! Number of 'blind' boundary components. 'Blind' boundary components
    ! belong to the boundary of subdomains but do not count to the real boundary.
    ! In arrays like IverticesAtBoundary, IedgesAtBoundary etc., blind 
    ! boundary components are always attached to the information about the
    ! real boundary.
    integer                  :: NblindBCT = 0
    
    ! Number of vertices on the boundary.
    ! For 2D domains, this coincides with the number of edges on 
    ! the boundary: Every vertex on the boundary has an edge following
    ! the vertex in mathematical positive sense.
    integer                  :: NVBD = 0
    
    ! The number of faces on the boundary
    integer                  :: NABD = 0
    
    ! Number of edges on the boundary; coincides with NVBD
    ! for 2D domains.
    integer                  :: NMBD = 0

    ! Maximum number of vertices per element.
    ! 3 for 2D triangles, 4 for 2D quads or 3D tetrahedrons, 12 for 3D hexas
    integer                  :: NNVE = 0
    
    ! Maximum number of edges per element.
    ! 3 for 2D triangles, 4 for 2D quads or 3D tetrahedrons, 12 for 3D hexas.
    integer                  :: NNEE = 0
    
    ! Maximum number of areas per element. One hexa has e.g. 6 areas.
    integer                  :: NNAE = 0
    
    ! Maximum number of vertices per face. 3 for 3D tetraheral meshes,
    ! 4 for 3D hexahedral meshes. Unused in 2D.
    integer                  :: NNVA = 0
    
    ! Maximum number of elements adjacent to a vertex.
    integer                  :: NNelAtVertex = 0
    
    ! Maximum number of elements adjacent to an edge. 
    ! =0 in 1D, =2 in 2D, arbitrary in 3D.
    integer                  :: NNelAtEdge = 0
    
    ! Number of elements with a defined number of vertices per element.
    ! <verb>
    ! InelOfType(TRIA_NVELINE1D) = number of lines in the mesh (1D).
    ! InelOfType(TRIA_NVETRI2D)  = number of triangles in the mesh (2D).
    ! InelOfType(TRIA_NVEQUAD2D) = number of quads in the mesh (2D).
    ! InelOfType(TRIA_NVEPYR3D)  = number of pyramides in the mesh (3D).
    ! InelOfType(TRIA_NVEPRIS3D) = number of prisms in the mesh (3D).
    ! InelOfType(TRIA_NVETET3D)  = number of tetrahedra in the mesh (3D).
    ! InelOfType(TRIA_NVEHEXA3D) = number of hexahedra in the mesh (3D).
    ! </verb>
    integer, dimension(TRIA_MAXNVE) :: InelOfType = 0
  
    ! Number of vertices per edge; normally = 0.
    ! If a regular distribution of vertices on edges is given,
    ! NVPED saves the number of vertices on each edge;
    ! e.g. 1 if midpoints on edges exist in p_DfreeVertexCoordinates.
    integer                  :: nverticesPerEdge = 0
    
    ! Total number of vertices on edges; normally = 0.
    ! Total number of vertices on all edges, realized in p_DfreeVertexCoordinates. 
    ! E.g. if midpoints on edges exist, there is NVEDT=NMT.
    integer                  :: nVerticesOnAllEdges = 0
    
    ! Number of inner-element vertices; normally = 0.
    ! If a regular distribution of vertices in the inner of 
    ! each element is given, NIELV saves the number of vertices 
    ! in the inner of each element; e.g. 1 if element-midpoints 
    ! exist in p_DfreeVertexCoordinates.
    integer                  :: nverticesInEachElement = 0

    ! Total number of vertices in elements; normally = 0.
    ! Total number of vertices on all elements, realized in p_DfreeVertexCoordinates. 
    ! E.g. if element-midpoints exist in DCORMG, there is NIEVT=NEL.
    integer                  :: nverticesInAllElements = 0
    
    ! Number of additional vertices; normally = 0.
    ! Can be set <> 0 if there are any additional vertices 
    ! realized in p_DfreeVertexCoordinates, that do not belong to a regular 
    ! distribution of vertices in corners, on edges or on elements.
    integer                  :: nadditionalVertices = 0
  
    ! Minimum X/Y/Z coordinate of a bounding box around the mesh.
    ! This information is only present in a 'standard' mesh.
    real(DP), dimension(NDIM3D) :: DboundingBoxMin = (/0.0_DP,0.0_DP,0.0_DP/)

    ! Maximum X/Y/Z coordinate of a bounding box around the mesh.
    ! This information is only present in a 'standard' mesh.
    real(DP), dimension(NDIM3D) :: DboundingBoxMax = (/0.0_DP,0.0_DP,0.0_DP/)
    
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    ! for 1D meshes,
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! for 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Z-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    integer        :: h_DvertexCoords = ST_NOHANDLE
    
    ! Vertices Adjacent to an Element.
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    integer        :: h_IverticesAtElement = ST_NOHANDLE

    ! Edges Adjacent to an Element.
    ! Handle to 
    !       p_IedgesAtElement = array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the edges following the
    ! corner vertices in mathematically positive sense.
    ! This is similar to the old KMID array (but it is not exactly
    ! the old KMID of FEAT1 as KMID started the edge numbers at NVT+1
    ! while IedgesAtElement starts the edge numbering with 1).
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IedgesAtElement(4,.)=0 for a triangle in a quad mesh.
    integer        :: h_IedgesAtElement = ST_NOHANDLE
    
    ! Neighbour Elements Adjacent to an Element.
    ! Handle to 
    !       p_IneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element, the numbers of adjacent elements
    ! in mathematically positive sense, meeting the element in an edge.
    ! p_RneighbourElement(IEL)\%Ineighbours(.) describes the elements adjacent 
    ! to IEL along the edges (p_RedgesOnElement(IEL)\%Iedges(.)-NVT).
    ! This is the old KADJ array.
    !
    ! Note:  For meshes with hanging vertices, this array is slightly
    ! modified. For 'big' elements this array contains the element
    ! numbers of the 'first' adjacent element via an edge/face.
    ! Note: To access all elements adjacent to an element via a
    ! hanging vertex, calculate the vertex number of the hanging
    ! vertex via InodalProperty and access all adjacent elements.
    ! For small 'hanging' elements, this array contains as usual
    ! the number of the 'big' adjacent element(s).
    integer        :: h_IneighboursAtElement = ST_NOHANDLE
    
    ! Elements Adjacent to an Edge. Only 2D.
    ! Handle to 
    !       p_IelementsAtEdge = array [1..2,1..NMT] of integer.
    ! The numbers of the two elements adjacent to an edge IMT in 2D. 
    ! For boundary edges, p_IelementsOnEdge(2,IMT) is set to 0.
    ! This is the old KMEL array.
    integer        :: h_IelementsAtEdge = ST_NOHANDLE

    ! Vertices Adjacent to an Edge. 
    ! Handle to 
    !       p_IverticesAtEdge = array [1..2,1..NMT]
    ! The numbers of the two vertices adjacent to an edge IMT. 
    integer        :: h_IverticesAtEdge = ST_NOHANDLE
    
    ! Nodal property array. 
    ! Handle to 
    !       p_InodalProperty=array [1..NVT+NMT+NAT] of integer.
    ! p_InodalProperty(i) defines for each vertex i=(1..NVT),
    ! each edge i=(NVT+1..NVT+NMT) and face i=NVT+NMT+1..NVT+NMT+NAT
    ! its function inside of the geometry.
    ! Generally said, the range of the p_InodalProperty-array 
    ! characterizes the type of the node (=vertex/edge):
    ! = 0    : The vertex/edge is an inner vertex/edge
    ! > 0    : The vertex/edge is a boundary vertex/edge on the real
    !           boundary. KNPR(.) defines the number of the boundary
    !           component.
    ! This is the old KNPR-array, slightly modified for edges and
    ! hanging nodes.
    !
    ! In case there are hanging nodes in the mesh, this array
    ! has a special meaning for all hanging vertices and all edges
    ! containing hanging vertices.  Values < 0 indicate hanging
    ! vertices at an edge. 
    ! Let iedgeC (NVT+1..NVT+NMT) be the number of
    ! a a 'full' edge containing the hanging vertex jvertex. 
    ! Let iedge be one of the sub-edges inside of edge iedgeC.
    ! Then there is:
    !   p_InodalProperty(jvertex) = -iedgeC
    !   p_InodalProperty(iedge)   = -iedgeC
    !   p_InodalProperty(iedgeC)  = -jvertex
    ! Let kfaceC (NVT+NMT+1..NVT+NMT+NAT) be the number of a 'full' face
    ! containing the hanging vertex jvertex. 
    ! Let kface be the number of a one of the subfaces inside
    ! the face kfaceC. Let iedge be the number of one of the sub-edges
    ! inside face kfaceC.
    ! Then there is:
    !   p_InodalProperty(jvertex) = -kfaceC
    !   p_InodalProperty(kface)   = -kfaceC
    !   p_InodalProperty(iedge)   = -kfaceC
    !   p_InodalProperty(kfaceC)  = -jvertex
    ! A hanging vertex is either the midpoint of a face or of an edge,
    ! therefore this assignment is unique due to the range of the number.
    ! 'Hanging edges' (only appear in 3D) without a hanging vertex
    ! in the center of an edge/face are not supported.
    integer         :: h_InodalProperty = ST_NOHANDLE
    
    ! Macro nodal property array.
    ! Handle to 
    !       p_ImacroNodalProperty=array [1..NVT+NMT+NAT+NEL] of integer.
    ! p_InodalProperty(i) defines for each vertex i=(1..NVT),
    ! each edge i=(NVT+1..NVT+NMT), face i=(NVT+NMT+1..NVT+NMT+NAT)
    ! and element i=(NVT+NMT+NAT+1..NVT+NMT+NAT+NEL)
    ! its origin in relation to a coarse mesh, a refined mesh stems from.
    ! Let nvtC, nmtC, natC and nelC define the number of vertices, edges and
    ! faces, resp., on the coarse mesh. Then, the value of ImacroNodalProperty
    ! is as follows:
    ! 1..nvtC: The number of a coarse grid vertex, a fine grid vertex stems from.
    ! nvtC+1..nvtC+nmtC: The number of a coarse grid edge, a fine grid 
    !                    vertex/edge stems from.
    ! nvtC+nmtC+1..nvtC+nmtC+natC : The number of a coarse grid face, a fine grid
    !                    vertex/edge/face stems from.
    ! nvtC+nmtC+natC+1..nvtC+nmtC+natC+nelC : The number of a coarse grid element, a 
    !                    fine grid element stems from.
    ! The macro nodal property array can be created based on for an extended
    ! raw coarse mesh by tria_initMacroNodalProperty. Once created, this information
    ! belongs to the set of information of a raw mesh and is automatically 
    ! propagated to every finer mesh upon refinement.
    integer         :: h_ImacroNodalProperty = ST_NOHANDLE
    
    ! 2D triangulation: Array with area of each element.
    ! 3D triangulation: Array with volume of each element.
    ! Handle to 
    !       p_DelementArea = array [1..NEL+1] of double.
    ! p_DelementArea [NEL+1] gives the total area/voloume of the domain.
    integer         :: h_DelementVolume = ST_NOHANDLE
    
    ! Handle to 
    !       p_IelementsAtVertexIdx=array [1..NVT+1] of integer.
    ! Index array for p_IelementsAtVertex of length NVT+1 for describing the
    ! elements adjacent to a corner vertex. for vertex IVT, the array
    ! p_IelementsAtVertex contains the numbers of the elements around this
    ! vertex at indices 
    !     p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1.
    ! By subtracting
    !     p_IelementsAtVertexIdx(IVT+1)-p_IelementsAtVertexIdx(IVT)
    ! One can get the number of elements adjacent to a vertex IVT.
    integer        :: h_IelementsAtVertexIdx = ST_NOHANDLE
    
    ! Array containing the Elements Adjacent to a Vertex.
    ! Handle to 
    !       p_IelementsAtVertex = array(1..*) of integer
    ! p_IelementsAtVertex ( p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1 )
    ! contains the number of the adjacent element in a vertex.
    ! This replaces the old KVEL array.
    !
    ! Note: For hanging vertices, this array contains only those
    ! elements which are 'corner adjacent' to a vertex (i.e. the 'smaller' elements).
    ! The 'big' elements adjacent to the edge which the hanging vertex
    ! is a midpoint of are not part of the vertex neighbourhood
    ! in this array.
    integer        :: h_IelementsAtVertex = ST_NOHANDLE
    
    ! This array defines a patch index and is created during the
    ! refinement. For a mesh that does not stem from a refinement,
    ! this array is undefined. For a mesh that stems from a refinement,
    ! this array defines a pointer into p_IrefinementPatch.
    ! Using p_IrefinementPatchIdx in combination with p_IrefinementPatch
    ! allows to get the element numbers of the fine grid elements
    ! that were created from a coarse grid element.
    ! Handle to 
    !       p_IrefinementPatchIdx=array [1..NELcoarse+1] of integer.
    ! So if IEL is the number of a coarse grid element, the fine grid 
    ! elements are to be found in
    ! p_IrefinementPatch ( p_IrefinementPatchIdx(IEL)..p_IrefinementPatchIdx(IEL+1)-1 )
    ! By subtracting
    !     p_IrefinementPatchIdx(IVT+1)-p_IrefinementPatchIdx(IVT)
    ! one can get the number of elements that a coarse grid element
    ! was refined to.
    integer        :: h_IrefinementPatchIdx = ST_NOHANDLE
    
    ! This array defines patches that were created during the
    ! refinement. For a mesh that does not step from a refinement,
    ! this array is undefined. For a mesh that stems from a refinement,
    ! this array defines for every coarse grid element the number
    ! of the fine grid elements that were created from that coarse
    ! grid element.
    ! Handle to
    !       p_IrefinementPatches = array(1..*) of integer
    ! p_IrefinementPatch ( p_IrefinementPatchIdx(IEL)..p_IrefinementPatchIdx(IEL+1)-1 )
    ! contains the numbers of the elements, a coarse grid element
    ! was divided into.
    integer        :: h_IrefinementPatch = ST_NOHANDLE
    
    ! If a mesh stems from a refinement of a coarse mesh, this array
    ! defines for every element on the fine mesh the element
    ! number on the coarse mesh where the fine grid element comes from.
    ! Handle to
    !       p_IcoarseGridElement = array(1..NEL) of integer
    ! For a mesh that does not come from a refinement, this handle is 
    ! undefined.
    integer        :: h_IcoarseGridElement = ST_NOHANDLE
    
    ! Boundary component index vector of length NBCT+NblindBCT+1 for 
    ! p_IverticesAtBoundary / p_IedgesAtBoundary / ... arrays.
    ! For standard meshes, this is a handle to 
    !      p_IboundaryCpIdx = array [1..NBCT+1] of integer.
    ! For subdomains, this is a handle to
    !      p_IboundaryCpIdx = array [1..NBCT+2] of integer.
    ! where p_IboundaryCpIdx(NBCT+1) points to the beginning of the vertices
    ! of the 'blind' boundary component, which does not belong to the
    ! physical boundary of the domain.
    !
    ! For a (real) boundary component i all corner nodes 
    ! for that boundary component are saved in 
    !   p_IverticesAtBoundary ( p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1 ).
    ! All boundary edges of this boundary component are saved in 
    !   p_IedgesAtBoundary ( p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1 ).
    ! p_IboundaryCpIdx(NBCT+1 / +2) points to NVBD+1 for easier access to the last
    ! boundary component.
    ! This is the old KBCT array.
    integer          :: h_IboundaryCpIdx = ST_NOHANDLE
    
    ! Vertices on boundary. 
    ! Handle to 
    !       p_IverticesAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all vertices on the (real) boundary
    ! in mathematically positive sense.
    ! The boundary vertices of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KVBD array.
    integer          :: h_IverticesAtBoundary = ST_NOHANDLE

    ! In 2d this is a reference to the IboundaryCpIdx array.
    ! (Note: By duplication flag, both arrays are identical and have
    ! the same handles).
    ! In 3D, this array is in index into the p_IedgesAtBoundary array that
    ! works like the p_IboundaryCpIdx for the vertices... see above.
    integer          :: h_IboundaryCpEdgesIdx = ST_NOHANDLE
    
    ! Edges adjacent to the boundary. 
    ! Handle to 
    !       p_IedgesAtBoundary = array [1..NMBD] of integer.
    ! This array contains a list of all edges on the (real) boundary.
    ! 2D: in mathematically positive sense. 
    ! 3D: with increasing number.
    ! The boundary edges of boundary component i are saved at
    !        p_IboundaryCpEdgesIdx(i)..p_IboundaryCpEdgesIdx(i+1)-1.
    ! This is the old KMBD array.
    ! (Note: In 2D, the above index pointer coincides with
    !        p_IboundaryCpEdgesIdx(i)..p_IboundaryCpEdgesIdx(i+1)-1 ).
    integer          :: h_IedgesAtBoundary = ST_NOHANDLE
    
    ! Ihis array is an index into the p_IfacesAtBoundary array it
    ! works like the p_IboundaryCpIdx for the vertices... see above.
    integer          :: h_IboundaryCpFacesIdx = ST_NOHANDLE

    ! Faces adjacent to the boundary. Only 3D, undefined in 2D.
    ! Handle to 
    !       p_IfacesAtBoundary = array [1..NMBD] of integer.
    ! This array contains a list of all edges on the (real) boundary
    ! with increasing number.
    ! The boundary edges of boundary component i are saved at
    !        p_IboundaryCpFacesIdx(i)..p_IboundaryCpFacesIdx(i+1)-1.
    integer          :: h_IfacesAtBoundary = ST_NOHANDLE

    ! Elements Adjacent to the boundary. 
    ! Handle to 
    !       p_IelementsAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all elements on the (real) boundary
    ! in mathematically positive sense.
    ! p_IelementsAtBoundary(i) is the element adjacent to edge
    ! h_IedgesAtBoundary - therefore one element number might appear
    ! more than once in this array!
    ! The boundary elements of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KEBD array.
    integer          :: h_IelementsAtBoundary = ST_NOHANDLE
    
    ! Parameter values of vertices on the boundary.
    ! Handle to 
    !       p_DvertexParameterValue = array [1..NVBD] of real
    ! p_DvertexParameterValue(i) contains the parameter value of 
    ! boundary vertex i, which corresponds to the corner vertex 
    ! p_IverticesAtBoundary(I).
    ! There may be vertices on the boundary of the triangulation which
    ! are not part of the physical boundary of the domain (-> subdomains).
    ! Such vertices receive parameter value -1.0_DP! All such vertices
    ! are collected in the 'blind' boundary component NBCT+1.
    integer          :: h_DvertexParameterValue = ST_NOHANDLE

    ! Parameter values of edge midpoints on the boundary.
    ! Handle to 
    !       p_DedgeParameterValue = array [1..NMBD] of real
    ! p_DedgeParameterValue(i) contains the parameter value of
    ! the midpoint of boundary edge i, which corresponds to the 
    ! edge p_IedgesAtBoundary(I).
    ! There may be edges on the boundary of the triangulation which
    ! are not part of the physical boundary of the domain (-> subdomains).
    ! Such edges receive parameter value -1.0_DP!
    ! Such edges are not necessarily collected in the 'blind' boundary
    ! component NBCT+1, but can be anywhere in the IedgesAtBoundary array.
    integer          :: h_DedgeParameterValue = ST_NOHANDLE
    
    ! Inverse index array to p_IverticesAtBoundary. 
    ! Handle to 
    !       p_IboundaryVertexPos = array [1..2,1..NVBD] of integer.
    ! p_IboundaryVertexPos(1,.) contains a vertex number on the 
    ! boundary and p_IboundaryVertexPos(2,.) the appropriate index
    ! of this vertex inside of the p_IverticesAtBoundary-array. 
    ! The array is first sorted for the boundary component and therefore
    ! has a layout as defined by IboundaryCpIdx.
    ! Inside of each boundary component, the entries are sorted for
    ! the vertex number, thus allowing quick access to the index of a 
    ! vertex in p_IverticesAtBoundary.
    integer          :: h_IboundaryVertexPos = ST_NOHANDLE

    ! Inverse index array to p_IedgesAtBoundary. 
    ! Handle to 
    !       p_IboundaryEdgePos = array [1..2,1..NMBD] of integer.
    ! p_IboundaryEdgePos(1,.) contains a vertex number on the 
    ! boundary and p_IboundaryEdgePos(2,.) the appropriate index
    ! of this vertex inside of the p_IedgesAtBoundary-array. 
    ! The array is first sorted for the boundary component and therefore
    ! has a layout as defined by IboundaryCpIdx.
    ! Inside of each boundary component, the entries are sorted for
    ! the edge number, thus allowing quick access to the index of an
    ! edge in p_IedgesAtBoundary.
    integer          :: h_IboundaryEdgePos = ST_NOHANDLE
    
    ! Handle to 
    !       p_DfreeVertexCoordinates = array [1..NDIM2D,1..NVT] of double
    ! Array containing the coordinates of all vertices on edges,
    ! inner element vertices and additional nodes in the geometry.
    !
    ! p_DfreeVertexCoordinates(1..nVerticesOnAllEdges) contains the coordinates 
    ! of the regular distributed vertices on edges. 
    !
    ! p_DfreeVertexCoordinates(nVerticesOnAllEdges + 1 ..
    !                          nVerticesOnAllEdges + nVerticesOnAllElements) 
    ! contains the coordinates of the regular distributed 
    ! inner-element vertices. 
    !
    ! p_DfreeVertexCoordinates(nVerticesOnAllEdges + nVerticesOnAllElements + 1 ..
    !                          nVerticesOnAllEdges + nVerticesOnAllElements + nadditionalVertices)
    ! contains the coordinates of any additional vertices that do not 
    ! belong to regular distributed vertices on edges or on elements.
    !
    ! This is an extended version of the old DCORMG array.
    ! The DCORMG-array is originally designed to collect the midpoints 
    ! that are associated to the edges. In the new style, this behaviour 
    ! is only a special case that happens if NVEDT=NMT is set.
    ! The new style allowS p_DfreeVertexCoordinates to save all vertex 
    ! coordinates that are not corner vertices. This includes regular distributed
    ! vertices on edges (either midpoints or points at 1/3 and 2/3 of the edge, 
    ! or...), on elements (either midpoints or the 4 Gauss-points, or...)
    ! as well as additional vertices the user wants to be realized, that
    ! do not belong to any of the previous two groups.
    !
    ! The numbers of regularly distributed vertices on edges/elements can
    ! be calculated directly with an appropriate formula. E.g. let us assume,
    ! we have n regularly distributed vertices on each edge (excluding the
    ! starting/ending point). Then the corresponding vertices on 
    ! egde E (=1..NMT) have the numbers:
    !       (NVT-1)+(NMT-1) + (E-1)*n + 1
    !    .. (NVT-1)+(NMT-1) + (E-1)*n + n
    ! the formula for regular distributed vertices on elements is similar.
    integer           :: h_DfreeVertexCoordinates = ST_NOHANDLE
    
    ! Handle to h_IelementsAtEdgeIdx3D=array [1..NMT+1] of integer.
    ! Index array for h_IelementsAtEdge3D of length NMT+1 for describing the
    ! elements attached to an edge. for edge IVE, the array
    ! p_IelementsAtEdge3D contains the indices of the elements around this
    ! edge at the array positions
    !     p_IelementsAtEdgeIdx3D(IVE)..p_IelementsAtEdgeIdx3D(IVE+1)-1.
    ! By subtracting
    !     p_IelementsAtEdgeIdx3D(IVE+1)-p_IelementsAtEdgeIdx3D(IVE)
    ! One can get the number of elements attached to edge IVE. Only 3D.
    integer        :: h_IelementsAtEdgeIdx3D = ST_NOHANDLE
    
    ! Elements Adjacent to an Edge. Only 3D.
    ! Array containing the Elements Adjacent to an edge.
    ! Handle to 
    !       p_IelementsAtEdge3D = array(1..*) of integer
    ! p_IelementsAtEdge3D ( p_IelementsAtEdgeIdx3D(IVT)..p_IelementsAtEdgeIdx3D(IVT+1)-1 )
    ! contains the number of the adjacent element in an edge.
    integer        :: h_IelementsAtEdge3D = ST_NOHANDLE
    
    ! Handle to 
    !       p_IfacesAtEdgeIdx=array [1..NMT+1] of integer.
    ! Index array for p_IfacesAtEdge of length NMT+1 for describing the
    ! faces adjacent to an edge. For edge IMT, the array
    ! p_IfacesAtEdge contains the numbers of the faces around this
    ! edge at indices 
    !     p_IfacesAtEdgeIdx(IMT)..p_IfacesAtEdgeIdx(IMT+1)-1.
    ! By subtracting
    !     p_IfacesAtEdgeIdx(IMT+1)-p_IfacesAtEdgeIdx(IMT)
    ! One can get the number of faces adjacent to an edge. Only 3D.
    integer        :: h_IfacesAtEdgeIdx = ST_NOHANDLE

    ! Array containing the Faces Adjacent to an Edge.
    ! Handle to 
    !       p_IfacesAtEdge = array(1..*) of integer
    ! p_IfacesAtEdge ( p_IfacesAtEdgeIdx(IVT)..p_IfacesAtEdgeIdx(IVT+1)-1 )
    ! contains the number of the adjacent faces in an edge.
    integer        :: h_IfacesAtEdge    = ST_NOHANDLE
    
    ! Handle to 
    !       p_IfacesAtVertexIdx=array [1..NVT+1] of integer.
    ! Index array for p_IfacesAtVertex of length NVT+1 for describing the
    ! faces adjacent to an edge. For vertex IVT, the array
    ! p_IfacesAtVertex contains the numbers of the faces around this
    ! vertex at indices 
    !     p_IfacesAtVertexIdx(IVT)..p_IfacesAtVertexIdx(IVT+1)-1.
    ! By subtracting
    !     p_IfacesAtVertexIdx(IVT+1)-p_IfacesAtVertexIdx(IVT)
    ! One can get the number of faces adjacent to an edge. Only 3D.
    integer        :: h_IfacesAtVertexIdx = ST_NOHANDLE

    ! Array containing the Faces Adjacent to a Vertex. Only 3D.
    ! Handle to 
    !       p_IfacesAtVertex = array(1..*) of integer
    ! p_IfacesAtVertex ( p_IfacesAtVertexIdx(IVT)..p_IfacesAtVertexIdx(IVT+1)-1 )
    ! contains the number of the adjacent faces in a vertex.
    integer        :: h_IfacesAtVertex    = ST_NOHANDLE
     
    ! Faces adjacent to an element. Only 3D.
    ! Handle to 
    !       p_IfacesAtElement = array [1..NNAE,1..NEL] of integer
    ! For each element the node numbers of the edges following the
    ! corner vertices in mathematically positive sense.
    ! This is the old KMID array.
    ! On pure tethrahedral meshes, there is NNAE=4. 
    ! On mixed or pure hexahedral meshes, there is NNAE=6. In this case,
    ! there is IedgesAtElement(5:6,.)=0 for a tethrahedral in a hexahedral mesh.
    ! To be able to distinguish a number of an edge from a vertex number, 
    ! edges are numbered in the range NVT+NMT+1..NVT+NMT+NAT.
    integer        :: h_IfacesAtElement = ST_NOHANDLE
    
    ! Vertices Adjacent to a face.
    ! Handle to 
    !       p_IverticesAtFace = array [1..NVA,1..NAT] of integer
    ! For each face, the numbers of the vertices adjacent to that
    ! face. Ordered in mathematically positive sense when looking
    ! from the midpoint of the adjacent element with the smaller
    ! element number.
    ! On pure tetrahedral meshes, there is NVA=3. On mixed or pure 
    ! hexahedral meshes, there is NVA=4. In this case, there is 
    ! IverticesAtFace(4,.)=0 for a tetrahedral in a hexahedral mesh.
    integer        :: h_IverticesAtFace = ST_NOHANDLE 
    
    ! Elements Adjacent to a Face. Only 3D.
    ! Handle to 
    !       p_IelementsAtEdge = array [1..2,1..NAT] of integer.
    ! The numbers of the two elements adjacent to an edge IMT in 2D. 
    ! For boundary edges, p_IelementsOnEdge(2,IMT) is set to 0.
    ! This is the old KMEL array.
    integer        :: h_IelementsAtFace = ST_NOHANDLE

    ! Edges Adjacent to a Face. Only 3D.
    ! Handle to 
    !       p_IedgesAtFace = array [1..NNVA,1..NAT] of integer.
    ! The numbers of the edges adjacent to a face in 3D. 
    integer        :: h_IedgesAtFace = ST_NOHANDLE
    
    ! handle to an index array that is used to
    ! acces the IedgesAtVertices array
    ! this way we can get the edges attached
    ! to a vertex
    integer        :: h_IedgesAtVertexIdx = ST_NOHANDLE
    
    ! here we can store the edges adjacent
    ! to a vertex, to access this array
    ! use the IedgesAtVertexIdx array.
    ! Edge numbers in this array are in the range 1..NMT!
    integer        :: h_IedgesAtVertex = ST_NOHANDLE
    
    ! Handle to the twist index array:
    !
    !       p_ItwistIndexEdges = array [1..NEL] of integer(I32).
    !
    ! Each entry is a bitfield that specifies the orientation of the
    ! faces/edges of the element.
    ! Warning:
    ! This is a handle to an integer(I32) array and not an integer array!
    integer        :: h_ItwistIndex = ST_NOHANDLE

  end type t_triangulation
  
  public :: t_triangulation

!</typeblock>

!<typeblock>

  ! Defines a set of cells (coordinates, connectivity) that can be
  ! attached to a mesh.
  type t_cellSet
    
    ! Number of vertices in the set
    integer :: NVT = 0
    
    ! Number of elements in the set
    integer :: NEL = 0
    
    ! Array with the coordinates of the vertices defining the cells.
    ! dimension(#dimensions,#vertices)
    real(DP), dimension(:,:), pointer :: p_DvertexCoords => null()
    
    ! Array defining the connectivity.
    ! dimension(max. #vertices per element, #elements)
    integer, dimension(:,:), pointer :: p_IverticesAtElement => null()
    
    ! Array defining the nodal property of all vertices in the set.
    integer, dimension(:), pointer :: p_InodalProperty => null()
    
    ! Array with parameter values for all vertices in the
    ! cell set that are located on the physical boundary.
    ! dimension(#vertices).
    ! Vertices not on the boundary are identified by DvertexPar(.) = -1.
    real(DP), dimension(:), pointer :: p_DallVerticesParameterValue => null()
    
  end type t_cellSet
  
  public :: t_cellSet

!</typeblock>

!<typeblock>

  ! a connector connects to adjacent cells (i.e. a face in 3D)
  ! structure used to generate 3D connectivity
  type t_connector3D
    ! the array stores at 1-4 the vertices of a face
    ! at 5 the element the face belongs to
    ! at 6 it stores the local face number
    integer, dimension(6) :: I_conData
  end type t_connector3D
    
  public :: t_connector3D

!</typeblock>
  
!</types>

  interface tria_getNVE
    module procedure tria_getNVE_direct
    module procedure tria_getNVE_indirect
  end interface
  
  public :: tria_getNVE

  interface tria_getNAE
    module procedure tria_getNAE_direct
    module procedure tria_getNAE_indirect
  end interface
  
  public :: tria_getNAE

  public :: tria_duplicate
  public :: tria_restore
  public :: tria_done
  public :: tria_createRawTria1D
  public :: tria_readTriFile1D
  public :: tria_readTriFile2D
  public :: tria_readTriFile3D
  public :: tria_resetToRaw
  public :: tria_initExtendedRawMesh
  public :: tria_initStandardMeshFromRaw
  public :: tria_initMacroNodalProperty
  public :: tria_refine2LevelOrdering
  public :: tria_compress2LevelOrdHierarchy
  public :: tria_quickRefine2LevelOrdering
  public :: tria_rawGridToTri
  public :: tria_infoStatistics
  public :: tria_exportTriFile
  public :: tria_exportPostScript
  public :: tria_searchBoundaryNode
  public :: tria_generateSubdomain
  public :: tria_attachCells
  public :: tria_cellGroupGreedy
  public :: tria_getPointsOnEdge
  public :: tria_getNeighbourVertex
  public :: tria_getSubmeshNeighbourhood
  public :: tria_getElementsAtMacroEdge
  public :: tria_getElementsAtMacroVertex
  public :: tria_getVerticesAtFaceDirect
  public :: tria_readRawTriangulation1D 
  public :: tria_readRawTriangulation2D 
  public :: tria_readRawTriangulation3D
  public :: tria_genRawBoundary1D 
  public :: tria_genRawBoundary2D
  public :: tria_genRawBoundary3D
  public :: tria_genElementsAtVertex1D2D 
  public :: tria_genElementsAtVertex3D
  public :: tria_genNeighboursAtElement1D 
  public :: tria_genNeighboursAtElement2D 
  public :: tria_genNeighboursAtElement3D
  public :: tria_genEdgesAtElement2D 
  public :: tria_genEdgesAtElement3D
  public :: tria_genFacesAtElement3D
  public :: tria_genElementVolume1D 
  public :: tria_genElementVolume2D
  public :: tria_genElementVolume3D
  public :: tria_sortBoundaryVertices1D2D
  public :: tria_genElementsAtBoundary1D2D
  public :: tria_genBoundaryVertexPos1D2D
  public :: tria_genElementsAtEdge2D 
  public :: tria_genElementsAtEdge3D
  public :: tria_genVerticesAtEdge2D 
  public :: tria_genVerticesAtEdge3D
  public :: tria_genEdgeNodalProperty2D 
  public :: tria_genEdgeNodalProperty3D
  public :: tria_genFaceNodalProperty3D
  public :: tria_genEdgesAtBoundary2D 
  public :: tria_genEdgesAtBoundary3D
  public :: tria_genEdgeParameterValue2D
  public :: tria_genBoundaryEdgePos2D
  public :: tria_genEdgesAtVertex2D
  public :: tria_genFacesAtBoundary3D
  public :: tria_genFacesAtVertex3D
  public :: tria_genFacesAtEdge3D
  public :: tria_genEdgesAtFace3D
  public :: tria_genElementsAtFace3D
  public :: tria_genVerticesAtFace3D
  public :: tria_genTwistIndex
  public :: tria_propMacroNodalProperty2lv
  public :: tria_genRefTags2lv
  public :: tria_hangingNodeRefinement
  public :: tria_buildConnectorList
  public :: tria_sortElements3D
  public :: tria_sortElements3DInt
  public :: tria_mergesort
  public :: tria_merge

contains

  ! ***************************************************************************

!<subroutine>

  subroutine tria_duplicate (rtriangulation, rbackupTriangulation,&
                             iduplicationFlag, bupdate)
  
!<description>
  ! This routine makes a copy of a triangulation structure in memory.
  ! The variable iduplicationFlag decides on which arrays are copied in memory
  ! and which not.
  ! 
  ! By setting the corresponding bit in iduplicationFlag to 0, the array is
  ! duplicated in memory, and any change to the new array will not harm
  ! the original one.
  ! By setting a flag TR_SHARE_xxxx in iduplicationFlag, the corresponding array is
  ! not duplicated. The handle of the original structure is simply put into
  ! the new structure to make the information accessable. Then this array
  ! is shared between two triangulation structures!
!</description>

!<input>
  ! The "source" discretisation structure that provides the information.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Bitfield that decides which handles are a copy of another
  ! structure, thus which arrays are shared between the new
  ! and the old structure. 
  ! The bitfield contains a combination of TR_SHARE_xxxx canstants.
  ! Every triangulation information whose flag is set in iduplicationFlag
  ! is shared between rtriangulation and rbackupTriangulation.
  ! Therefore e.g., iduplicationFlag=0 copies all arrays from
  ! rtriangulation in memory, while TR_SHARE_ALL will copy
  ! nothing, but will hare everything between rtriangulation
  ! and rbackupTriangulation.
  integer(I32), intent(in)          :: iduplicationFlag
  
  ! OPTIONAL. Defines how to create the backup.
  ! = .FALSE.: Treat rbackupTriangulation as empty destination structure.
  !    If necessary, information in rbackupTriangulation is released.
  !    rbackupTriangulation is rebuild according to rtriangulation 
  !    and iduplicationFlag. 
  !    This is the standard setting if bupdate is not specified.
  ! = .TRUE. : Treat rbackupTriangulation as existing copy of rtriangulation
  !    which has to be updated. Recover all vectors of rbackupTriangulation 
  !    by those of rtriangulation.
  !    (I.e. those arrays which were duplicated by a previous
  !    call to iduplicationFlag with IUPD=0.)
  !    iduplicationFlag can used to specify which data do copy
  !    from rtriangulation to rbackupTriangulation. It is OR`ed
  !    with rbackupTriangulation%iduplicationFlag to get the actual
  !    duplication flag. This leads to the following interpretation
  !    of iduplicationFlag:
  !     =0:  Copy all data that was copied previously. This is the usual
  !          setting.
  !    <>0:  Copy all arrays where the corresponding flag
  !          in iduplicationFlag is not set and which exist as
  !          duplicates. Arrays corresponding to flags which
  !          are set or where the handles in rtriangulation and
  !          rbackupTriangulation coincide are not touched.
  logical, intent(in), optional     :: bupdate
!</input>

!<inputoutput>
  ! The "destination" discretisation structure receives the information.
  ! Depending on iduplicationFlag, the arrays are duplicated or shared
  ! between rtriangulation and rbackupTriangulation
  type(t_triangulation), intent(inout), target :: rbackupTriangulation
!</inputoutput>
  
!</subroutine>
  
    ! local variables
    integer(I32) :: idupFlag
    
    logical :: bupd
    
    bupd = .false.
    if (present(bupdate)) bupd = bupdate

    if (.not. bupd) then
      ! Release any old data.
      call tria_done (rbackupTriangulation)
      
      rbackupTriangulation%ndim                   = rtriangulation%ndim
      rbackupTriangulation%NVT                    = rtriangulation%NVT
      rbackupTriangulation%NMT                    = rtriangulation%NMT
      rbackupTriangulation%NAT                    = rtriangulation%NAT
      rbackupTriangulation%NEL                    = rtriangulation%NEL
      rbackupTriangulation%NBCT                   = rtriangulation%NBCT
      rbackupTriangulation%NblindBCT              = rtriangulation%NblindBCT
      rbackupTriangulation%NVBD                   = rtriangulation%NVBD
      rbackupTriangulation%NABD                   = rtriangulation%NABD
      rbackupTriangulation%NMBD                   = rtriangulation%NMBD
      rbackupTriangulation%NNVE                   = rtriangulation%NNVE
      rbackupTriangulation%NNEE                   = rtriangulation%NNEE
      rbackupTriangulation%NNAE                   = rtriangulation%NNAE
      rbackupTriangulation%NNVA                   = rtriangulation%NNVA
      rbackupTriangulation%NNelAtVertex           = rtriangulation%NNelAtVertex
      rbackupTriangulation%NNelAtEdge             = rtriangulation%NNelAtEdge
      rbackupTriangulation%InelOfType(:)          = rtriangulation%InelOfType(:)
      rbackupTriangulation%nverticesPerEdge       = rtriangulation%nverticesPerEdge
      rbackupTriangulation%nVerticesOnAllEdges    = rtriangulation%nVerticesOnAllEdges
      rbackupTriangulation%nverticesInEachElement = rtriangulation%nverticesInEachElement
      rbackupTriangulation%nverticesInAllElements = rtriangulation%nverticesInAllElements
      rbackupTriangulation%nadditionalVertices    = rtriangulation%nadditionalVertices   
      rbackupTriangulation%DboundingBoxMin        = rtriangulation%DboundingBoxMin   
      rbackupTriangulation%DboundingBoxMax        = rtriangulation%DboundingBoxMax
      
      ! Decide on IDPFLG which arrays to copy
      rbackupTriangulation%iduplicationFlag = iduplicationFlag
      idupFlag = iduplicationFlag
      
    else

      ! Create a bitfield what to copy by ORing iduplicationFlag with what
      ! we have in rbackupTriangulation. That way, only arrays that exist as
      ! real duplicates are copied from rtriangulation to rbackupTriangulation.
      
      idupFlag = ior(iduplicationFlag,rbackupTriangulation%iduplicationFlag)
    
    end if
    
    ! Call checkAndCopy for all the arrays. this will either copy the handle
    ! or allocate new memory and copy the content of the array.
    
    ! Bit  0: DCORVG
    call checkAndCopy(idupflag, TR_SHARE_DVERTEXCOORDS,&
          rtriangulation%h_DvertexCoords, &
          rbackupTriangulation%h_DvertexCoords)

    ! Bit  1: DCORMG 
    call checkAndCopy(idupflag, TR_SHARE_DFREEVERTEXCOORDINATES,&
          rtriangulation%h_DfreeVertexCoordinates, &
          rbackupTriangulation%h_DfreeVertexCoordinates)

    ! Bit  2: KVERT  
    call checkAndCopy(idupflag, TR_SHARE_IVERTICESATELEMENT,&
          rtriangulation%h_IverticesAtElement, &
          rbackupTriangulation%h_IverticesAtElement)

    ! Bit  3: KMID   
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATELEMENT,&
          rtriangulation%h_IedgesAtElement, &
          rbackupTriangulation%h_IedgesAtElement)
    
    ! Bit  4: KADJ   
    call checkAndCopy(idupflag, TR_SHARE_INEIGHBOURSATELEMENT,&
          rtriangulation%h_IneighboursAtElement, &
          rbackupTriangulation%h_IneighboursAtElement)

    ! Bit  5: KVEL   
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertexIdx, &
          rbackupTriangulation%h_IelementsAtVertexIdx)
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertex, &
          rbackupTriangulation%h_IelementsAtVertex)

    ! Bit  6: KMEL   
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge, &
          rbackupTriangulation%h_IelementsAtEdge)
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge3D, &
          rbackupTriangulation%h_IelementsAtEdge3D)
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdgeIdx3D, &
          rbackupTriangulation%h_IelementsAtEdgeIdx3D)

    ! Bit  7: KNPR   
    call checkAndCopy(idupflag, TR_SHARE_INODALPROPERTY,&
          rtriangulation%h_InodalProperty, &
          rbackupTriangulation%h_InodalProperty)

    ! Bit  8: KVBD   
    call checkAndCopy(idupflag, TR_SHARE_IVERTICESATBOUNDARY,&
          rtriangulation%h_IverticesAtBoundary, &
          rbackupTriangulation%h_IverticesAtBoundary)

    ! Bit  9: KEBD   
    call checkAndCopy(idupflag,TR_SHARE_IELEMENTSATBOUNDARY,&
          rtriangulation%h_IelementsAtBoundary, &
          rbackupTriangulation%h_IelementsAtBoundary)

    ! Bit 10: KBCT   
    call checkAndCopy(idupflag,TR_SHARE_IBOUNDARYCPIDX,&
          rtriangulation%h_IboundaryCpIdx, &
          rbackupTriangulation%h_IboundaryCpIdx)

    ! Bit 11: DVBDP  
    call checkAndCopy(idupflag,TR_SHARE_DVERTEXPARAMETERVALUE,&
          rtriangulation%h_DvertexParameterValue, &
          rbackupTriangulation%h_DvertexParameterValue)

    ! Bit 12: DMBDP  
    call checkAndCopy(idupflag,TR_SHARE_DEDGEPARAMETERVALUE,&
          rtriangulation%h_DedgeParameterValue, &
          rbackupTriangulation%h_DedgeParameterValue)

    ! Bit 13: KMBD   
    call checkAndCopy(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IedgesAtBoundary, &
          rbackupTriangulation%h_IedgesAtBoundary)
    call checkAndCopy(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IboundaryCpEdgesIdx, &
          rbackupTriangulation%h_IboundaryCpEdgesIdx)

    ! Bit 14: KEAN   
    call checkAndCopy(idupflag,TR_SHARE_IVERTICESATEDGE,&
          rtriangulation%h_IverticesAtEdge, &
          rbackupTriangulation%h_IverticesAtEdge)

    ! Bit 15: KVBDI  
    call checkAndCopy(idupflag,TR_SHARE_IBOUNDARYVERTEXPOS,&
          rtriangulation%h_IboundaryVertexPos, &
          rbackupTriangulation%h_IboundaryVertexPos)

    ! Bit 16: KMBDI  
    call checkAndCopy(idupflag,TR_SHARE_IBOUNDARYEDGEPOS,&
          rtriangulation%h_IboundaryEdgePos, &
          rbackupTriangulation%h_IboundaryEdgePos)

    ! Bit 17: DAREA  
    call checkAndCopy(idupflag,TR_SHARE_DELEMENTAREA,&
          rtriangulation%h_DelementVolume, &
          rbackupTriangulation%h_DelementVolume)

    ! Bit 18: IrefinementPatch
    call checkAndCopy(idupflag, TR_SHARE_IREFINEMENTPATCH,&
          rtriangulation%h_IrefinementPatch, &
          rbackupTriangulation%h_IrefinementPatch)
    call checkAndCopy(idupflag, TR_SHARE_IREFINEMENTPATCH,&
          rtriangulation%h_IrefinementPatchIdx, &
          rbackupTriangulation%h_IrefinementPatchIdx)

    ! Bit 19: IcoarseGridElement
    call checkAndCopy(idupflag, TR_SHARE_ICOARSEGRIDELEMENT,&
          rtriangulation%h_IcoarseGridElement, &
          rbackupTriangulation%h_IcoarseGridElement)

    ! Bit 20: KVAR
    call checkAndCopy(idupflag, TR_SHARE_IVERTICESATFACE,&
          rtriangulation%h_IverticesAtFace, &
          rbackupTriangulation%h_IverticesAtFace)

    ! Bit 21: KAREA
    call checkAndCopy(idupflag, TR_SHARE_IFACESATELEMENT,&
          rtriangulation%h_IfacesAtElement, &
          rbackupTriangulation%h_IfacesAtElement)

    ! Bit 22: IelementsAtFace
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATFACE,&
          rtriangulation%h_IelementsAtFace, &
          rbackupTriangulation%h_IelementsAtFace)

    ! Bit 23: IedgesAtFace
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATFACE,&
          rtriangulation%h_IedgesAtFace, &
          rbackupTriangulation%h_IedgesAtFace)

    ! Bit 24: IfacesAtEdge
    call checkAndCopy(idupflag, TR_SHARE_IFACESATEDGE,&
          rtriangulation%h_IfacesAtEdge, &
          rbackupTriangulation%h_IfacesAtEdge)
    call checkAndCopy(idupflag, TR_SHARE_IFACESATEDGE,&
          rtriangulation%h_IfacesAtEdgeIdx, &
          rbackupTriangulation%h_IfacesAtEdgeIdx)

    ! Bit 25: IfacesAtVertex
    call checkAndCopy(idupflag, TR_SHARE_IFACESATVERTEX,&
          rtriangulation%h_IfacesAtVertex, &
          rbackupTriangulation%h_IfacesAtVertex)
    call checkAndCopy(idupflag, TR_SHARE_IFACESATVERTEX,&
          rtriangulation%h_IfacesAtVertexIdx, &
          rbackupTriangulation%h_IfacesAtVertexIdx)

    ! Bit 26: IfacesAtBoundary
    call checkAndCopy(idupflag, TR_SHARE_IFACESATBOUNDARY,&
          rtriangulation%h_IfacesAtBoundary, &
          rbackupTriangulation%h_IfacesAtBoundary)
    call checkAndCopy(idupflag, TR_SHARE_IFACESATBOUNDARY,&
          rtriangulation%h_IboundaryCpFacesIdx, &
          rbackupTriangulation%h_IboundaryCpFacesIdx)

    ! Bit 27: IedgesAtVertex
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertex, &
          rbackupTriangulation%h_IedgesAtVertex)
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertexIdx, &
          rbackupTriangulation%h_IedgesAtVertexIdx)

    ! Bit 28: ItwistIndex
    call checkAndCopy(idupflag, TR_SHARE_ITWISTINDEX,&
          rtriangulation%h_ItwistIndex, &
          rbackupTriangulation%h_ItwistIndex)
    
    ! Bit 29: ImacroNodalProperty
    call checkAndCopy(idupflag, TR_SHARE_IMACRONODALPROPERTY,&
          rtriangulation%h_ImacroNodalProperty, &
          rbackupTriangulation%h_ImacroNodalProperty)
    
  contains
  
    subroutine checkAndCopy (idupFlag,ibitfield,isourcehandle,idesthandle)
    
    ! Checks if idupFlag has all bits ibitfield set.
    ! If yes, idesthandle is set to isourcehandle.
    ! Otherwise, the memory behind isourcehandle is duplicated in memory
    ! and idesthandle receives the handle to the new memory block.
    
    integer(I32), intent(in) :: ibitfield
    integer(I32), intent(in) :: idupFlag
    integer, intent(in) :: isourcehandle
    integer, intent(inout) :: idesthandle
    
      if (iand(idupFlag,ibitfield) .ne. ibitfield) then
        if (isourcehandle .ne. ST_NOHANDLE) then
          call storage_copy(isourcehandle,idesthandle)
        end if
      else
        idesthandle = isourcehandle
      end if
      
    end subroutine checkAndCopy

  end subroutine tria_duplicate

  ! ***************************************************************************

!<subroutine>

  subroutine tria_restore (rbackupTriangulation,rtriangulation)
  
!<description>
  ! This routine restores data of a triangulation structure. All
  ! information arrays not shared between rbackupTriangulation and another 
  ! triangulation structure is copied (back) into the rtriangulation.
!</description>

!<input>
  ! Backup of a triangulation structure.
  type(t_triangulation), intent(in) :: rbackupTriangulation
!</input>

!<inputoutput>
  ! Destination triangulation.
  ! All arrays where a duplicates exist in rtriangulation are copied
  ! to rtriangulation, overwriting the old information arrays.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer(I32) :: idupFlag
  
    idupFlag = rtriangulation%iduplicationFlag
      
    ! Call checkAndCopy for all the arrays. this will either copy the handle
    ! or copy the content of the array.
    
    ! Bit  0: DCORVG 
    call checkAndCopy(idupflag, TR_SHARE_DVERTEXCOORDS,&
          rtriangulation%h_DvertexCoords, &
          rbackupTriangulation%h_DvertexCoords)

    ! Bit  1: DCORMG 
    call checkAndCopy(idupflag, TR_SHARE_DFREEVERTEXCOORDINATES,&
          rtriangulation%h_DfreeVertexCoordinates, &
          rbackupTriangulation%h_DfreeVertexCoordinates)

    ! Bit  2: KVERT  
    call checkAndCopy(idupflag, TR_SHARE_IVERTICESATELEMENT,&
          rtriangulation%h_IverticesAtElement, &
          rbackupTriangulation%h_IverticesAtElement)

    ! Bit  3: KMID   
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATELEMENT,&
          rtriangulation%h_IedgesAtElement, &
          rbackupTriangulation%h_IedgesAtElement)
    
    ! Bit  4: KADJ   
    call checkAndCopy(idupflag, TR_SHARE_INEIGHBOURSATELEMENT,&
          rtriangulation%h_IneighboursAtElement, &
          rbackupTriangulation%h_IneighboursAtElement)

    ! Bit  5: KVEL   
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertexIdx, &
          rbackupTriangulation%h_IelementsAtVertexIdx)
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertex, &
          rbackupTriangulation%h_IelementsAtVertex)

    ! Bit  6: KMEL   
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge, &
          rbackupTriangulation%h_IelementsAtEdge)
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge3D, &
          rbackupTriangulation%h_IelementsAtEdge3D)
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdgeIdx3D, &
          rbackupTriangulation%h_IelementsAtEdgeIdx3D)

    ! Bit  7: KNPR   
    call checkAndCopy(idupflag, TR_SHARE_INODALPROPERTY,&
          rtriangulation%h_InodalProperty, &
          rbackupTriangulation%h_InodalProperty)

    ! Bit  8: KVBD   
    call checkAndCopy(idupflag, TR_SHARE_IVERTICESATBOUNDARY,&
          rtriangulation%h_IverticesAtBoundary, &
          rbackupTriangulation%h_IverticesAtBoundary)

    ! Bit  9: KEBD   
    call checkAndCopy(idupflag,TR_SHARE_IELEMENTSATBOUNDARY,&
          rtriangulation%h_IelementsAtBoundary, &
          rbackupTriangulation%h_IelementsAtBoundary)

    ! Bit 10: KBCT
    call checkAndCopy(idupflag,TR_SHARE_IBOUNDARYCPIDX,&
          rtriangulation%h_IboundaryCpIdx, &
          rbackupTriangulation%h_IboundaryCpIdx)

    ! Bit 11: DVBDP  
    call checkAndCopy(idupflag,TR_SHARE_DVERTEXPARAMETERVALUE,&
          rtriangulation%h_DvertexParameterValue, &
          rbackupTriangulation%h_DvertexParameterValue)

    ! Bit 12: DMBDP  
    call checkAndCopy(idupflag,TR_SHARE_DEDGEPARAMETERVALUE,&
          rtriangulation%h_DedgeParameterValue, &
          rbackupTriangulation%h_DedgeParameterValue)

    ! Bit 13: KMBD   
    call checkAndCopy(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IedgesAtBoundary, &
          rbackupTriangulation%h_IedgesAtBoundary)
    call checkAndCopy(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IboundaryCpEdgesIdx, &
          rbackupTriangulation%h_IboundaryCpEdgesIdx)
    
    ! Bit 14: KEAN   
    call checkAndCopy(idupflag,TR_SHARE_IVERTICESATEDGE,&
          rtriangulation%h_IverticesAtEdge, &
          rbackupTriangulation%h_IverticesAtEdge)

    ! Bit 15: KVBDI  
    call checkAndCopy(idupflag,TR_SHARE_IBOUNDARYVERTEXPOS,&
          rtriangulation%h_IboundaryVertexPos, &
          rbackupTriangulation%h_IboundaryVertexPos)

    ! Bit 16: KMBDI  
    call checkAndCopy(idupflag,TR_SHARE_IBOUNDARYEDGEPOS,&
          rtriangulation%h_IboundaryEdgePos, &
          rbackupTriangulation%h_IboundaryEdgePos)

    ! Bit 17: DAREA  
    call checkAndCopy(idupflag,TR_SHARE_DELEMENTAREA,&
          rtriangulation%h_DelementVolume, &
          rbackupTriangulation%h_DelementVolume)

    ! Bit 18: IrefinementPatch
    call checkAndCopy(idupflag, TR_SHARE_IREFINEMENTPATCH,&
          rtriangulation%h_IrefinementPatch, &
          rbackupTriangulation%h_IrefinementPatch)
    call checkAndCopy(idupflag, TR_SHARE_IREFINEMENTPATCH,&
          rtriangulation%h_IrefinementPatchIdx, &
          rbackupTriangulation%h_IrefinementPatchIdx)

    ! Bit 19: IcoarseGridElement
    call checkAndCopy(idupflag, TR_SHARE_ICOARSEGRIDELEMENT,&
          rtriangulation%h_IcoarseGridElement, &
          rbackupTriangulation%h_IcoarseGridElement)

    ! Bit 20: KVAR
    call checkAndCopy(idupflag, TR_SHARE_IVERTICESATFACE,&
          rtriangulation%h_IverticesAtFace, &
          rbackupTriangulation%h_IverticesAtFace)

    ! Bit 21: KAREA
    call checkAndCopy(idupflag, TR_SHARE_IFACESATELEMENT,&
          rtriangulation%h_IfacesAtElement, &
          rbackupTriangulation%h_IfacesAtElement)

    ! Bit 22: IelementsAtFace
    call checkAndCopy(idupflag, TR_SHARE_IELEMENTSATFACE,&
          rtriangulation%h_IelementsAtFace, &
          rbackupTriangulation%h_IelementsAtFace)
    
    ! Bit 23: IedgesAtFace
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATFACE,&
          rtriangulation%h_IedgesAtFace, &
          rbackupTriangulation%h_IedgesAtFace)

    ! Bit 24: IfacesAtEdge
    call checkAndCopy(idupflag, TR_SHARE_IFACESATEDGE,&
          rtriangulation%h_IfacesAtEdge, &
          rbackupTriangulation%h_IfacesAtEdge)
    call checkAndCopy(idupflag, TR_SHARE_IFACESATEDGE,&
          rtriangulation%h_IfacesAtEdgeIdx, &
          rbackupTriangulation%h_IfacesAtEdgeIdx)

    ! Bit 25: IfacesAtVertex
    call checkAndCopy(idupflag, TR_SHARE_IFACESATVERTEX,&
          rtriangulation%h_IfacesAtVertex, &
          rbackupTriangulation%h_IfacesAtVertex)
    call checkAndCopy(idupflag, TR_SHARE_IFACESATVERTEX,&
          rtriangulation%h_IfacesAtVertexIdx, &
          rbackupTriangulation%h_IfacesAtVertexIdx)

    ! Bit 26: IfacesAtBoundary
    call checkAndCopy(idupflag, TR_SHARE_IFACESATBOUNDARY,&
          rtriangulation%h_IfacesAtBoundary, &
          rbackupTriangulation%h_IfacesAtBoundary)
    call checkAndCopy(idupflag, TR_SHARE_IFACESATBOUNDARY,&
          rtriangulation%h_IboundaryCpFacesIdx, &
          rbackupTriangulation%h_IboundaryCpFacesIdx)

    ! Bit 27: IedgesAtVertex
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertex, &
          rbackupTriangulation%h_IedgesAtVertex)
    call checkAndCopy(idupflag, TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertexIdx, &
          rbackupTriangulation%h_IedgesAtVertexIdx)

    ! Bit 28: ItwistIndex
    call checkAndCopy(idupflag, TR_SHARE_ITWISTINDEX,&
          rtriangulation%h_ItwistIndex, &
          rbackupTriangulation%h_ItwistIndex)

    ! Bit 29: ImacroNodalProperty
    call checkAndCopy(idupflag, TR_SHARE_IMACRONODALPROPERTY,&
          rtriangulation%h_ImacroNodalProperty, &
          rbackupTriangulation%h_ImacroNodalProperty)
    
  contains
  
    subroutine checkAndCopy (idupFlag,ibitfield,idesthandle,isourcehandle)
    
    ! Checks if idupFlag has all bits ibitfield set.
    ! If not, the memory behind isourcehandle is copied to idesthandle
    ! overwriting all previous information.
    
    integer(I32), intent(in) :: ibitfield
    integer(I32), intent(in) :: idupFlag
    integer, intent(in) :: isourcehandle
    integer, intent(inout) :: idesthandle
    
      if (iand(idupFlag,ibitfield) .ne. ibitfield) then
        if (isourcehandle .ne. ST_NOHANDLE) then
          call storage_copy(isourcehandle,idesthandle)
        end if
      end if
      
    end subroutine checkAndCopy

  end subroutine tria_restore

  ! ***************************************************************************

!<subroutine>

  subroutine tria_done (rtriangulation)
  
!<description>
  ! This routine cleans up a triangulation structure.
  ! All memory allocated by handles in the structure is released from the heap.
  ! The routine does not clean up old FEAT 1.x legacy arrays in the
  ! rtriangulation%Itria substructure!
!</description>

!<inputoutput>
  ! The triangulation structure to be cleaned up.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    integer(I32) :: idupflag
    
    if (rtriangulation%ndim .eq. 0) return
    
    idupflag = rtriangulation%iduplicationFlag
    
    ! Just release all allocated handles.
    ! Take care of which handles are duplicates from other structures - 
    ! these must not be released, as we are not the owner of them!
    
    ! Bit  0: DCORVG is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_DVERTEXCOORDS,&
          rtriangulation%h_DvertexCoords)

    ! Bit  1: DCORMG is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_DFREEVERTEXCOORDINATES,&
          rtriangulation%h_DfreeVertexCoordinates)

    ! Bit  2: KVERT  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IVERTICESATELEMENT,&
          rtriangulation%h_IverticesAtElement)

    ! Bit  3: KMID   is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IEDGESATELEMENT,&
          rtriangulation%h_IedgesAtElement)
    
    ! Bit  4: KADJ   is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_INEIGHBOURSATELEMENT,&
          rtriangulation%h_IneighboursAtElement)

    ! Bit  5: KVEL   is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertexIdx)
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertex)

    ! Bit  6: KMEL   is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge)
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdgeIdx3D)
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge3D)

    ! Bit  7: KNPR   is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_INODALPROPERTY,&
          rtriangulation%h_InodalProperty)

    ! Bit  8: KVBD   is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IVERTICESATBOUNDARY,&
          rtriangulation%h_IverticesAtBoundary)

    ! Bit  9: KEBD   is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IELEMENTSATBOUNDARY,&
          rtriangulation%h_IelementsAtBoundary)

    ! Bit 10: KBCT   is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IBOUNDARYCPIDX,&
          rtriangulation%h_IboundaryCpIdx)

    ! Bit 11: DVBDP  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_DVERTEXPARAMETERVALUE,&
          rtriangulation%h_DvertexParameterValue)

    ! Bit 12: DMBDP  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_DEDGEPARAMETERVALUE,&
          rtriangulation%h_DedgeParameterValue)

    ! Bit 13: KMBD   is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IedgesAtBoundary)
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IboundaryCpEdgesIdx)

    ! Bit 14: KEAN   is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IVERTICESATEDGE,&
          rtriangulation%h_IverticesAtEdge)

    ! Bit 15: KVBDI  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IBOUNDARYVERTEXPOS,&
          rtriangulation%h_IboundaryVertexPos)

    ! Bit 16: KMBDI  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IBOUNDARYEDGEPOS,&
          rtriangulation%h_IboundaryEdgePos)

    ! Bit 17: DAREA  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_DELEMENTAREA,&
          rtriangulation%h_DelementVolume)

    ! Bit 18: IrefinementPatch  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IREFINEMENTPATCH, &
           rtriangulation%h_IrefinementPatch)
    call checkAndRelease(idupflag, TR_SHARE_IREFINEMENTPATCH, &
           rtriangulation%h_IrefinementPatchIdx)

    ! Bit 19: IcoarseGridElement  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_ICOARSEGRIDELEMENT, &
           rtriangulation%h_IcoarseGridElement)

    ! Bit 20: KVAR  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IVERTICESATFACE, &
           rtriangulation%h_IverticesAtFace)
           
    ! Bit 21: KAREA  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IFACESATELEMENT, &
           rtriangulation%h_IFacesAtElement)
    
    ! Bit 22: IelementsAtFace  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATFACE, &
           rtriangulation%h_IelementsAtFace)
    
    ! Bit 23: IedgesAtEdge  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IEDGESATFACE, &
           rtriangulation%h_IedgesAtFace)
           
    ! Bit 24: IfacesAtEdge  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IFACESATEDGE, &
           rtriangulation%h_IfacesAtEdgeIdx)
    call checkAndRelease(idupflag, TR_SHARE_IFACESATEDGE, &
           rtriangulation%h_IfacesAtEdge)

    ! Bit 25: IfacesAtVertex  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IFACESATVERTEX, &
           rtriangulation%h_IfacesAtVertexIdx)
    call checkAndRelease(idupflag, TR_SHARE_IFACESATVERTEX, &
           rtriangulation%h_IfacesAtVertex)
    
    ! Bit 26: IfacesAtBoundary  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IFACESATBOUNDARY, &
           rtriangulation%h_IfacesAtBoundary)
    call checkAndRelease(idupflag,TR_SHARE_IFACESATBOUNDARY,&
          rtriangulation%h_IboundaryCpFacesIdx)
     
    ! Bit 27: IedgesAtVertex  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertexIdx)
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertex)

    ! Bit 28: ItwistIndex  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_ITWISTINDEX,&
          rtriangulation%h_ItwistIndex)
    
    ! Bit 29: ImacroNodalProperty  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IMACRONODALPROPERTY,&
          rtriangulation%h_InodalProperty)

    ! Clean up the rest of the structure

    rtriangulation%iduplicationFlag = 0
    rtriangulation%ndim = 0
    rtriangulation%NVT = 0
    rtriangulation%NNVE = 0
    rtriangulation%NNAE = 0
    rtriangulation%NNVA = 0
    rtriangulation%NNEE = 0
    rtriangulation%NMT = 0
    rtriangulation%NAT = 0
    rtriangulation%NEL = 0
    rtriangulation%NBCT = 0
    rtriangulation%NblindBCT = 0
    rtriangulation%NVBD = 0
    rtriangulation%NMBD = 0
    rtriangulation%nverticesPerEdge = 0
    rtriangulation%nVerticesOnAllEdges = 0
    rtriangulation%nverticesInEachElement = 0
    rtriangulation%nverticesInAllElements = 0
    rtriangulation%nadditionalVertices = 0
    rtriangulation%InelOfType(:) = 0
    rtriangulation%NNelAtVertex = 0
    rtriangulation%NNelAtEdge = 0
    rtriangulation%DboundingBoxMin(:) = 0.0_DP
    rtriangulation%DboundingBoxMax(:) = 0.0_DP

    ! That is it...

  contains
  
    ! **********************************************************
    ! Release handle ihandle if bitfield ibitfield in idubFlag is not set.
    ! Otherwise, ihandle is set to ST_NOHANDLE.    
    subroutine checkAndRelease (idupFlag,ibitfield,ihandle)
    
    integer(I32), intent(in) :: ibitfield
    integer(I32), intent(in) :: idupFlag
    integer, intent(inout) :: ihandle
    
      if (iand(idupFlag,ibitfield) .ne. ibitfield) then
        if (ihandle .ne. ST_NOHANDLE) call storage_free(ihandle)
      else
        ihandle = ST_NOHANDLE
      end if
      
    end subroutine checkAndRelease

  end subroutine tria_done

  ! ***************************************************************************

!<subroutine>

  subroutine tria_resetToRaw (rtriangulation,bkeepExtendedRaw)
  
!<description>
  ! This routine resets a mesh to a raw mesh by deleting all not necessary
  ! information.
!</description>

!<inputoutput>
  ! The triangulation structure to be resetted.
  type(t_triangulation), intent(inout) :: rtriangulation
  
  ! OPTIONAL: If not specified or set to TRUE, arrays of the extended
  ! raw mesh are not removed. When set to FALSE, the mesh is reset to pure
  ! RAW state.
  logical, intent(in), optional :: bkeepExtendedRaw
!</inputoutput>
  
!</subroutine>

    integer(I32) :: idupflag
    logical :: bext
    
    bext = .true.
    if (present(bkeepExtendedRaw)) bext = bkeepExtendedRaw
    
    if (rtriangulation%ndim .eq. 0) return
    
    idupflag = rtriangulation%iduplicationFlag
    
    ! Just release all allocated handles.
    ! Take care of which handles are duplicates from other structures - 
    ! these must not be released, as we are not the owner of them!
    
    ! Bit  3: KMID   is a copy of another structure
    if (.not. bext) then
      call checkAndRelease(idupflag, TR_SHARE_IEDGESATELEMENT,&
            rtriangulation%h_IedgesAtElement)
    end if
    
    ! Bit  4: KADJ   is a copy of another structure
    if (.not. bext) then
      call checkAndRelease(idupflag, TR_SHARE_INEIGHBOURSATELEMENT,&
            rtriangulation%h_IneighboursAtElement)
    end if
    
    if (.not. bext) then
      ! Bit  5: KVEL   is a copy of another structure
      call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
            rtriangulation%h_IelementsAtVertexIdx)
      call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
            rtriangulation%h_IelementsAtVertex)
    end if

    ! Bit  6: KMEL   is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge)
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdgeIdx3D)
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge3D)

    ! Bit  9: KEBD   is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IELEMENTSATBOUNDARY,&
          rtriangulation%h_IelementsAtBoundary)

    ! Bit 12: DMBDP  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_DEDGEPARAMETERVALUE,&
          rtriangulation%h_DedgeParameterValue)

    ! Bit 13: KMBD   is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IedgesAtBoundary)
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IboundaryCpEdgesIdx)

    ! Bit 14: KEAN   is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IVERTICESATEDGE,&
          rtriangulation%h_IverticesAtEdge)

    ! Bit 15: KVBDI  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IBOUNDARYVERTEXPOS,&
          rtriangulation%h_IboundaryVertexPos)

    ! Bit 16: KMBDI  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IBOUNDARYEDGEPOS,&
          rtriangulation%h_IboundaryEdgePos)

    ! Bit 17: DAREA  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_DELEMENTAREA,&
          rtriangulation%h_DelementVolume)

    ! Bit 18: IrefinementPatch  is a copy of another structure      
    call checkAndRelease(idupflag, TR_SHARE_IREFINEMENTPATCH, &
           rtriangulation%h_IrefinementPatch)
    call checkAndRelease(idupflag, TR_SHARE_IREFINEMENTPATCH, &
           rtriangulation%h_IrefinementPatchIdx)

    ! Bit 19: IcoarseGridElement  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_ICOARSEGRIDELEMENT, &
           rtriangulation%h_IcoarseGridElement)

    ! Bit 20: KVAR  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IVERTICESATFACE, &
           rtriangulation%h_IverticesAtFace)

    if (.not. bext) then           
      ! Bit 21: KAREA  is a copy of another structure
      call checkAndRelease(idupflag, TR_SHARE_IFACESATELEMENT, &
            rtriangulation%h_IFacesAtElement)
    end if
    
    ! Bit 22: IelementsAtFace  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATFACE, &
           rtriangulation%h_IelementsAtFace)
           
    ! Bit 23: IedgesAtFace  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IEDGESATFACE, &
           rtriangulation%h_IedgesAtFace)

    ! Bit 24: IfacesAtEdge  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IFACESATEDGE, &
           rtriangulation%h_IfacesAtEdge)
    call checkAndRelease(idupflag, TR_SHARE_IFACESATEDGE, &
           rtriangulation%h_IfacesAtEdgeIdx)
   
    ! Bit 25: IfacesAtVertex  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IFACESATVERTEX, &
           rtriangulation%h_IfacesAtVertex)
    call checkAndRelease(idupflag, TR_SHARE_IFACESATVERTEX, &
           rtriangulation%h_IfacesAtVertexIdx)

    ! Bit 26: IfacesAtBoundary  is a copy of another structure
    call checkAndRelease(idupflag, TR_SHARE_IFACESATBOUNDARY, &
           rtriangulation%h_IfacesAtBoundary)
    call checkAndRelease(idupflag,TR_SHARE_IFACESATBOUNDARY,&
          rtriangulation%h_IboundaryCpFacesIdx)

    ! Bit 27: IedgesAtVertex  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertexIdx)
    call checkAndRelease(idupflag,TR_SHARE_IEDGESATVERTEX,&
          rtriangulation%h_IedgesAtVertex)

    ! Bit 28: ItwistIndex  is a copy of another structure
    call checkAndRelease(idupflag,TR_SHARE_ITWISTINDEX,&
          rtriangulation%h_ItwistIndex)
           
    if (.not. bext) then           
      ! Bit 29: ImacroNodalProperty  is a copy of another structure
      call checkAndRelease(idupflag, TR_SHARE_IMACRONODALPROPERTY,&
          rtriangulation%h_ImacroNodalProperty)
    end if

    ! Clean up the rest of the structure

    rtriangulation%DboundingBoxMin(:) = 0.0_DP
    rtriangulation%DboundingBoxMax(:) = 0.0_DP
    if (.not. bext) then  
      rtriangulation%NMT = 0
      rtriangulation%NAT = 0
      rtriangulation%NNelAtVertex = 0
      rtriangulation%NNelAtEdge = 0
      rtriangulation%iduplicationFlag = iand(rtriangulation%iduplicationFlag,&
        not(TR_SHARE_DVERTEXCOORDS+TR_SHARE_IVERTICESATELEMENT+&
        TR_SHARE_INODALPROPERTY+TR_SHARE_IVERTICESATBOUNDARY))
    else
      rtriangulation%iduplicationFlag = iand(rtriangulation%iduplicationFlag,&
        not(TR_SHARE_EXTENDEDRAW+TR_SHARE_DVERTEXCOORDS+TR_SHARE_IVERTICESATELEMENT+&
        TR_SHARE_INODALPROPERTY+TR_SHARE_IVERTICESATBOUNDARY))
    end if

    ! That is it...

  contains
  
    ! **********************************************************
    ! Release handle ihandle if bitfield ibitfield in idubFlag is not set.
    ! Otherwise, ihandle is set to ST_NOHANDLE.    
    subroutine checkAndRelease (idupFlag,ibitfield,ihandle)
    
    integer(I32), intent(in) :: ibitfield
    integer(I32), intent(in) :: idupFlag
    integer, intent(inout) :: ihandle
    
      if (iand(idupFlag,ibitfield) .ne. ibitfield) then
        if (ihandle .ne. ST_NOHANDLE) call storage_free(ihandle)
      else
        ihandle = ST_NOHANDLE
      end if
      
    end subroutine checkAndRelease

  end subroutine tria_resetToRaw

  ! ***************************************************************************

!<subroutine>

  subroutine tria_initExtendedRawMesh (rtriangulation)

!<description>
  ! This routine initialises an 'extended raw mesh'. An extended raw mesh
  ! provides a numbering for edges and faces (in 3D) but does not have
  ! any adjacency structures attached.
  ! A standard mesh can be formed from an extended raw mesh by the
  ! usual call to tria_initStandardMeshFromRaw.
!</description>
  
!<inputoutput>
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>

    select case (rtriangulation%ndim)
    case (NDIM1D)
      ! Generate basic arrays for 1D meshes.
      call tria_genElementsAtVertex1D2D  (rtriangulation)
      call tria_genNeighboursAtElement1D (rtriangulation)

    case (NDIM2D)
      ! Generate all standard arrays for 2D meshes.
      call tria_genElementsAtVertex1D2D  (rtriangulation)
      call tria_genNeighboursAtElement2D (rtriangulation)
      call tria_genEdgesAtElement2D      (rtriangulation)

    case (NDIM3D)
      ! Generate all standard arrays for 3D meshes.
      call tria_genElementsAtVertex3D    (rtriangulation)
      call tria_genNeighboursAtElement3D (rtriangulation)
      call tria_genEdgesAtElement3D      (rtriangulation)
      call tria_genFacesAtElement3D      (rtriangulation)

    case DEFAULT
      call output_line ('Triangulation structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_initExtendedRawMesh')
      call sys_halt()
    end select

  end subroutine tria_initExtendedRawMesh

  !************************************************************************

!<subroutine>

  subroutine tria_initStandardMeshFromRaw(rtriangulation, rboundary, igenflag)

!<description>
  ! This routine creates all 'standard' arrays based on a 'raw'
  ! triangulation.
  ! rtriangulation is a 'raw' triangulation as provided by preprocessing
  ! routines like tria_readTriFile2D. Based on this raw triangulation,
  ! all standard arrays will be created (adjacencies, edge information,...)
  ! such that the application can work with the triangulation as usual.
  !
  ! If rtriangulation is already a 'standard' triangulation, all
  ! information in the mesh is regenerated. Missing information is added.
!</description>

!<input>
  ! OPTIONAL: Boundary structure that defines the domain.
  ! If not specified, information about boundary vertices (e.g. 
  ! parameter values of edge midpoints in 2D) are not initialised.
  type(t_boundary), intent(in), optional :: rboundary

  ! OPTIONAL: Flags which can be used to prevent the creation of some structure.
  ! Note that this routine does not check whether arrays already exist,
  ! but regenerates all information indicated by this flags.
  integer(I32), intent(in), optional     :: igenflag
!</input>

!<inputoutput>
  ! Triangulation structure to be initialised.
  type(t_triangulation), intent(inout)   :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer(I32) :: iflag

    iflag = iand(TR_GEN_ALL,not(TR_GEN_EXTENDEDRAW))
    if (present(igenflag)) iflag = igenflag
    
    ! If the mesh is not an extended raw mesh, we have to initialise it
    ! at first.
    if ((rtriangulation%h_IelementsAtVertexIdx .eq. ST_NOHANDLE) .or. &
        (iand(iflag,TR_GEN_EXTENDEDRAW) .eq. TR_GEN_EXTENDEDRAW)) then
      
      ! Generate an extended raw mesh
      call tria_initExtendedRawMesh (rtriangulation)
      
      ! Do not regenerate just generated arrays in the following
      ! select-case statement.
      iflag = iand(iflag,not(TR_GEN_EXTENDEDRAW))
    end if
    
    ! Get the bounding box of the mesh.
    call tria_calcBoundingBox(rtriangulation,&
        rtriangulation%DboundingBoxMin,rtriangulation%DboundingBoxMax)
    
    select case (rtriangulation%ndim)
    case (NDIM1D)
      ! Generate all standard arrays for 1D meshes.
      if (checkGen(iflag, TR_GEN_IELEMENTSATVERTEX))&
          call tria_genElementsAtVertex1D2D (rtriangulation)
      
      if (checkGen(iflag, TR_GEN_INEIGHBOURSATELEMENT))&
          call tria_genNeighboursAtElement1D (rtriangulation)

      if (checkGen(iflag, TR_GEN_DELEMENTAREA))&
          call tria_genElementVolume1D (rtriangulation)
      
      if (checkGen(iflag, TR_GEN_IELEMENTSATBOUNDARY) .or.&
          checkGen(iflag, TR_GEN_IBOUNDARYVERTEXPOS)) then
        call tria_sortBoundaryVertices1D2D  (rtriangulation)
        call tria_genElementsAtBoundary1D2D (rtriangulation)
        call tria_genBoundaryVertexPos1D2D  (rtriangulation)
      end if
       
    case (NDIM2D)
      ! Generate all standard arrays for 2D meshes.
      if (checkGen(iflag, TR_GEN_IELEMENTSATVERTEX))&
          call tria_genElementsAtVertex1D2D (rtriangulation)

      if (checkGen(iflag, TR_GEN_INEIGHBOURSATELEMENT))&
          call tria_genNeighboursAtElement2D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IEDGESATELEMENT))&
          call tria_genEdgesAtElement2D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IELEMENTSATEDGE))&
          call tria_genElementsAtEdge2D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IVERTICESATEDGE))&
          call tria_genVerticesAtEdge2D (rtriangulation)

      if (checkGen(iflag, TR_GEN_INODALPROPERTY))then
            call tria_genEdgeNodalProperty2D (rtriangulation)
      end if

      if (checkGen(iflag, TR_GEN_DELEMENTAREA))&
          call tria_genElementVolume2D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IELEMENTSATBOUNDARY) .or.&
          checkGen(iflag, TR_GEN_IEDGESATBOUNDARY)    .or.&
          checkGen(iflag, TR_GEN_DEDGEPARAMETERVALUE) .or.&
          checkGen(iflag, TR_GEN_IBOUNDARYVERTEXPOS)  .or.&
          checkGen(iflag, TR_GEN_IBOUNDARYEDGEPOS)) then
        call tria_sortBoundaryVertices1D2D (rtriangulation)

        if (checkGen(iflag, TR_GEN_IELEMENTSATBOUNDARY))&
            call tria_genElementsAtBoundary1D2D (rtriangulation)

        if (checkGen(iflag, TR_GEN_IEDGESATBOUNDARY))&
            call tria_genEdgesAtBoundary2D (rtriangulation)
        
        if (present(rboundary) .and. &
            checkGen(iflag, TR_GEN_DEDGEPARAMETERVALUE))&
            call tria_genEdgeParameterValue2D (rtriangulation, rboundary)
        
        if (checkGen(iflag, TR_GEN_IBOUNDARYVERTEXPOS))&
            call tria_genBoundaryVertexPos1D2D (rtriangulation)
        
        if (checkGen(iflag, TR_GEN_IBOUNDARYEDGEPOS))&
            call tria_genBoundaryEdgePos2D (rtriangulation)
      end if
      
      if (checkGen(iflag, TR_GEN_IEDGESATVERTEX))&
          call tria_genEdgesAtVertex2D (rtriangulation)

      if (checkGen(iflag, TR_GEN_ITWISTINDEX))&
          call tria_genTwistIndex (rtriangulation)
      
    case (NDIM3D)
      ! vertices at element info provided by tri-File
      if (checkGen(iflag, TR_GEN_IELEMENTSATVERTEX))&
          call tria_genElementsAtVertex3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_INEIGHBOURSATELEMENT))&
          call tria_genNeighboursAtElement3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IEDGESATELEMENT))&
          call tria_genEdgesAtElement3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IELEMENTSATEDGE))&
          call tria_genElementsAtEdge3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IVERTICESATEDGE))&
          call tria_genVerticesAtEdge3D (rtriangulation)
      
      ! faces have global numbers
      ! nvt+nmt+1 = first face
      if (checkGen(iflag, TR_GEN_IFACESATELEMENT))&
          call tria_genFacesAtElement3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IVERTICESATFACE))&
          call tria_genVerticesAtFace3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IELEMENTSATFACE))&
          call tria_genElementsAtFace3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IEDGESATFACE))&
          call tria_genEdgesAtFace3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IFACESATEDGE))&
          call tria_genFacesAtEdge3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IFACESATVERTEX))&
          call tria_genFacesAtVertex3D (rtriangulation)
      
      !----BOUNDARY------
      if (checkGen(iflag, TR_GEN_IFACESATBOUNDARY))&
          call tria_genFacesAtBoundary3D (rtriangulation)

      if (checkGen(iflag, TR_GEN_IEDGESATBOUNDARY))&
          call tria_genEdgesAtBoundary3D (rtriangulation)
      
      !----Properties----!
      if (checkGen(iflag, TR_GEN_INODALPROPERTY)) then
        call tria_genEdgeNodalProperty3D (rtriangulation)
        call tria_genFaceNodalProperty3D (rtriangulation)
      end if

      if (checkGen(iflag, TR_GEN_ITWISTINDEX))&
          call tria_genTwistIndex (rtriangulation)
      
      ! call tria_genElementVolume3D (rtriangulation)
    
    case DEFAULT
      call output_line ('Triangulation structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_generateStandardMeshFromRaw')
      call sys_halt()
    end select

  contains
    
    function checkGen (igenFlag,ibitfield)
      
      ! Checks if igenFlag has all bits ibitfield set.
      integer(I32), intent(in) :: igenFlag
      integer(I32), intent(in) :: ibitfield
      
      logical                  :: checkGen
      
      checkGen = (iand(igenFlag,ibitfield) .eq. ibitfield)

    end function checkGen
  end subroutine tria_initStandardMeshFromRaw

  !************************************************************************

!<subroutine>
  
  subroutine tria_initMacroNodalProperty (rtriangulation)
  
!<description>
  ! Calculates and attaches a macro nodal property array to the
  ! triangulation rtriangulation. This information is propagated
  ! upon refinement as part of an extended raw mesh to all finer 
  ! meshes. It allows to determine the origin of vertices, edges and
  ! faces in relation to a 'coarse' mesh.
!</description>

!<inputoutput>
  ! The triangulation structure where a macro nodal property array
  ! should be attached to.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: isize

    ! Allocate memory for that array if it does not exist.
    if (rtriangulation%h_ImacroNodalProperty .eq. ST_NOHANDLE) then
      call storage_new ('tria_initMacroNodalProperty', 'KMCPR', &
          rtriangulation%NVT+rtriangulation%NMT+&
          rtriangulation%NAT+rtriangulation%NEL, &
          ST_INT, rtriangulation%h_ImacroNodalProperty, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_ImacroNodalProperty, isize)
      if (isize .lt. rtriangulation%NVT+rtriangulation%NMT+&
                     rtriangulation%NAT+rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        ! Copy the old content as we must not destroy the old nodal property
        ! tags of the vertices.
        call storage_realloc ('tria_initMacroNodalProperty', &
            rtriangulation%NVT+rtriangulation%NMT+&
            rtriangulation%NAT+rtriangulation%NEL, &
            rtriangulation%h_ImacroNodalProperty, &
            ST_NEWBLOCK_NOINIT, .true.)
      end if
    end if
    
    ! Initialise the array by increasíng numbers.
    ! Vertices get numbers 1..NVT, edges NVT+1..NMT, faces NVT+NMT+1,NVT+NMT+NAT,
    ! and elements NVT+NMT+NAT+1..NVT+NMT+NAT+NEL
    ! so simply initialise by a sequence 1..NVT+NMT+NAT+NEL.
    call storage_initialiseBlock (&
        rtriangulation%h_ImacroNodalProperty, ST_NEWBLOCK_ORDERED)

  end subroutine tria_initMacroNodalProperty

  !************************************************************************

!<subroutine>

  subroutine tria_refine2LevelOrdering(rsourceTriangulation,&
                                       rdestTriangulation, rboundary, cflags)

!<description>
  ! This routine refines the given mesh rsourceTriangulation according to
  ! the 2-level ordering algorithm. The refined mesh is saved in 
  ! rdestTriangulation.
  !
  ! rsourceTriangulation must be a 'standard' mesh. The resulting refined
  ! mesh in rdestTriangulation will be an 'extended raw' mesh (like as being read 
  ! from a .TRI file, except if this generation is prevented by the
  ! appropriate cflags setting), i.e. the caller must add further information 
  ! to it with routines like tria_initStandardMeshFromRaw!
!</description>

!<input>
  ! OPTIONAL: Boundary structure that defines the parametrisation of the boundary.
  ! If specified, the coordinates of the new boundary vertices are
  ! recomputed according to the analytic boundary.
  type(t_boundary), intent(in), optional :: rboundary
  
  ! OPTIONAL: Bitfield of TRIA_R2LV_xxxx constants that allow to specify
  ! options for the refinement. If not specified, TRIA_R2LV_STANDARD
  ! is used as default.
  integer(I32), intent(in), optional :: cflags
!</input>

!<inputoutput>
  ! The source triangulation to be refined; must be a 'standard' mesh.
  type(t_triangulation), intent(inout) :: rsourceTriangulation
!</inputoutput>

!<output>
  ! OPTIONAL: Destination triangulation structure that receives the refined mesh.
  ! The refined mesh will be a 'raw' mesh as being read from a .TRI file e.g..
  ! If not specified, the source triangulation rsourceTriangulation is replaced
  ! by the refined mesh.
  type(t_triangulation), intent(out), optional :: rdestTriangulation
!</output>
  
!</subroutine>
 
    type(t_triangulation) :: rdestTria
    integer(I32) :: cflagsAct
    
    cflagsAct = TRIA_R2LV_STANDARD
    if (present(cflags)) cflagsAct = cflags
    
    ! Call the correct submethod depending on the dimension.
    select case (rsourceTriangulation%ndim)
    case (NDIM1D)
      ! Refine the basic mesh
      call tria_refineMesh2lv1D(rsourceTriangulation, rdestTria)

    case (NDIM2D)
      ! Refine the basic mesh
      call tria_refineMesh2lv2D(rsourceTriangulation, rdestTria)
      
      ! Refine the boundary
      call tria_refineBdry2lv2D(rsourceTriangulation, rdestTria,&
          iand(cflagsAct,TRIA_R2LV_RECALCCOORDSONBD) .ne. 0,rboundary)
      
      if (iand(cflagsAct,TRIA_R2LV_AVERAGEMIDPOINTS) .ne. 0) then
        ! Recalculate corner points of quads that were element midpoints
        ! on the coarse mesh.
        call tria_averageMidpoints2D(rsourceTriangulation, rdestTria)
      end if
      
    case (NDIM3D)
      ! refine the basic mesh
      call tria_refineMesh2lv3D(rsourceTriangulation, rdestTria)

      ! Refine the boundary
      call tria_refineBdry2lv3D(rsourceTriangulation, rdestTria, rboundary)
      
    case DEFAULT
      call output_line ('Triangulation structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_refine2LevelOrdering')
      call sys_halt()
    end select
    
    ! Generate an extended raw mesh if not prohibited.
    if (iand(cflagsAct,TRIA_R2LV_NOEXTENDEDRAW) .eq. 0) then
      call tria_initExtendedRawMesh (rdestTria)
      
      ! If there is a macro nodal property array attached to the coarse
      ! mesh, create one for the fine mesh. This array is an additional
      ! and optional array and will only be created on the fine mesh
      ! if there is a 'parent' one available to be propagated.
      if (rsourceTriangulation%h_ImacroNodalProperty .ne. ST_NOHANDLE) then
        call tria_propMacroNodalProperty2lv (rsourceTriangulation, rdestTria)
      end if
    end if
    
    ! Either copy rdestTria to rdestTriangulation or overwrite the source
    ! triangulation.
    if (present(rdestTriangulation)) then
      rdestTriangulation = rdestTria
    else
      call tria_done (rsourceTriangulation)
      rsourceTriangulation = rdestTria
    end if

  contains

    ! ---------------------------------------------------------------
  
    subroutine tria_refineMesh2lv1D(rsourceTriangulation, rdestTriangulation)

    ! This routine refines the given 1D mesh rsourceTriangulation according to
    ! the 2-level ordering algorithm. The refined mesh is saved in 
    ! rdestTriangulation. There will be no correction of boundary
    ! vertices. Boundary parameter values are not handled here!

    ! The source triangulation to be refined
    type(t_triangulation), intent(in) :: rsourceTriangulation

    ! Destination triangulation structure that receives the refined mesh. 
    type(t_triangulation), intent(out) :: rdestTriangulation
    
      ! local variables
      real(DP), dimension(:,:), pointer :: p_DcoordSource
      real(DP), dimension(:,:), pointer :: p_DcoordDest
      integer, dimension(:,:), pointer :: p_IvertAtElementSource
      integer, dimension(:,:), pointer :: p_IvertAtElementDest
      integer, dimension(:), pointer :: p_IrefinementPatchIdx
      integer, dimension(:), pointer :: p_IrefinementPatch
      integer, dimension(:), pointer :: p_IcoarseGridElement
      integer, dimension(:), pointer :: p_InodalPropSource, p_InodalPropDest
      integer :: iel,iel2
      integer :: ivt1,ivt2, ivtoffset, ivt
      integer, dimension(2) :: Isize
      integer :: nnve, ive
      real(DP) :: x,y
      
      ! Get the arrays with information of the source mesh.
      call storage_getbase_double2D (rsourceTriangulation%h_DvertexCoords,&
          p_DcoordSource)
      call storage_getbase_int2D (rsourceTriangulation%h_IverticesAtElement,&
          p_IvertAtElementSource)
      
      ! The 2-level ordering has the following properties:
      !
      ! - Vertices in the coarse mesh are vertices in the fine mesh 
      !   with the same number. Coordinates of coarse mesh vertices
      !   stay unchanged in the fine mesh.
      ! - Line midpoints in the coarse mesh become vertices in the
      !   fine mesh. They are appended to the vertices from the
      !   coarse mesh
      
      nnve = rsourceTriangulation%NNVE
      
      ! Initialise the basic mesh data in rdestTriangulation:
      
      ! 1D mesh
      rdestTriangulation%ndim = NDIM1D
      
      rdestTriangulation%NNVE = nnve
      
      ! Every element is divided into 2 subelements
      rdestTriangulation%NEL = 2 * rsourceTriangulation%NEL
      rdestTriangulation%InelOfType(:) = 2 * rsourceTriangulation%InelOfType(:)

      ! We expect NVT+NEL new points.
      rdestTriangulation%NVT = rsourceTriangulation%NVT + &
                               rsourceTriangulation%NEL
      
      ! Allocate memory for the new vertex coordinates and
      ! get the pointers to the coordinate array
      Isize = (/NDIM1D,rdestTriangulation%NVT/)
      call storage_new ('tria_refineMesh2lv1D', 'DCORVG',&
          Isize, ST_DOUBLE, &
          rdestTriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2D(&
          rdestTriangulation%h_DvertexCoords,p_DcoordDest)
      
      ! Allocate memory for the refinement information arrays.
      ! These arrays define for every coarse grid element the fine
      ! grid elements and for every fine grid element the coarse
      ! grid element where it comes from.
      call storage_new ('tria_refineMesh2lv1D', 'h_IrefinementPatchIdx', &
          rsourceTriangulation%NEL+1, ST_INT, &
          rdestTriangulation%h_IrefinementPatchIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatchIdx,p_IrefinementPatchIdx)

      call storage_new ('tria_refineMesh2lv1D', 'h_IrefinementPatch', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IrefinementPatch, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatch,p_IrefinementPatch)
    
      call storage_new ('tria_refineMesh2lv1D', 'h_IcoarseGridElement', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IcoarseGridElement, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IcoarseGridElement,p_IcoarseGridElement)

      ! The p_IrefinementPatchIdx array can directly be initialised
      ! as every coarse grid element gets 4 fine grid elements.
      do iel = 0,rsourceTriangulation%NEL
        p_IrefinementPatchIdx(1+iel) = 1 + iel*2
      end do
      
      ! Also the p_IcoarseGridElement array can directly be initialised.
      ! The coarse grid element number of element 1..NEL(coarse) stays
      ! the same. The number of the other coarse grid elements can
      ! be calculated by a formula.
      do iel = 1,rsourceTriangulation%NEL
        p_IcoarseGridElement(iel) = iel
      end do

      do iel = rsourceTriangulation%NEL+1,rdestTriangulation%NEL
        p_IcoarseGridElement(iel) = iel-rsourceTriangulation%NEL
      end do

      ! Ok, let us start the refinement. In the first step, we copy the
      ! corner coordinates of the coarse mesh to the fine mesh; they
      ! do not change during the refinement.
      do ivt=1,rsourceTriangulation%NVT
        p_DcoordDest(1,ivt) = p_DcoordSource(1,ivt)
      end do
      
      ! Each line produces an line midpoint which is stored as new
      ! point in the fine mesh. To calculate the coordinates, take
      ! the mean of the coordinates in the coarse mesh.
      ivtoffset = rsourceTriangulation%NVT

      do iel=1, rsourceTriangulation%NEL
        ivt1 = p_IvertAtElementSource (1,iel)
        ivt2 = p_IvertAtElementSource (2,iel)
        p_DcoordDest(1,ivtoffset+iel) = &
            0.5_DP * ( p_DcoordSource (1,ivt1) + p_DcoordSource (1,ivt2) )
      end do
      
      ! Allocate memory for IverticesAtElement and get a pointer to it.
      ! Fill the array with zero, so we will not have problems when mixing
      ! triangles into a quad mesh.
      Isize = (/nnve,rdestTriangulation%NEL/)
      call storage_new ('tria_refineMesh2lv1D', 'KVERT', Isize, ST_INT, &
          rdestTriangulation%h_IverticesAtElement, ST_NEWBLOCK_ZERO)
      call storage_getbase_int2D(&
          rdestTriangulation%h_IverticesAtElement,p_IvertAtElementDest)
    
      ! Look at the following line element with local and global vertex numbers:
      !
      !   1                   2
      !  X---------------------X
      ! 101        IEL        102
      !
      ! The two-level ordering refinement strategy says now:
      !  - The vertex numbers in the coarse grid are transferred to the fine grid
      !    without change.
      !  - The element numbers in the coarse grid define the vertex numbers of the
      !    element midpoints in the fine grid.
      !  - The line element number IEL in the coarse grid is transferred
      !    to the subelement at local vertex 1. The other element
      !    gets the element number NEL+IEL.
      !
      ! In the above picture, we therefore have:
      !
      !   1        2 1        2
      !  X----------X----------X
      ! 101  IEL1  201  IEL2  102
      !
      ! IEL1 = IEL
      ! IEL2 = NEL+IEL
      do iel = 1,rsourceTriangulation%NEL
    
        ! Determine number of subelements.
        iel2 = rsourceTriangulation%NEL + iel
        
        ! Save them
        p_IrefinementPatch(1+(iel-1)*2) = iel
        p_IrefinementPatch(1+(iel-1)*2+1) = iel2

        ! IEL1
        p_IvertAtElementDest(1,iel) = p_IvertAtElementSource(1, iel)
        p_IvertAtElementDest(2,iel) = ivtoffset + iel
        
        ! IEL2
        p_IvertAtElementDest(1,iel2) = ivtoffset + iel
        p_IvertAtElementDest(2,iel2) = p_IvertAtElementSource(2, iel)
      
      end do
      
      ! Create a new nodal property array
      call storage_new('tria_refineMesh2lv1D', 'h_InodalProperty',&
          rdestTriangulation%NVT, ST_INT,&
          rdestTriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
      
      ! Get the nodal property arrays
      call storage_getbase_int(rsourceTriangulation%h_InodalProperty, &
          p_InodalPropSource)
      call storage_getbase_int(rdestTriangulation%h_InodalProperty, &
          p_InodalPropDest)
      
      ! Copy the nodal property of the coarse mesh
      do ivt=1, rsourceTriangulation%NVT
        p_InodalPropDest(ivt) = p_InodalPropSource(ivt)
      end do
      
      ! We also need to copy the boundary vertices array to the finer level.
      rdestTriangulation%NVBD      = rsourceTriangulation%NVBD
      rdestTriangulation%NBCT      = rsourceTriangulation%NBCT
      rdestTriangulation%NblindBCT = rsourceTriangulation%NblindBCT

      call storage_copy (rsourceTriangulation%h_IverticesAtBoundary,&
          rdestTriangulation%h_IverticesAtBoundary)
      call storage_copy (rsourceTriangulation%h_IboundaryCpIdx,&
          rdestTriangulation%h_IboundaryCpIdx)

      ! Let us see if parameter values of boundary vertices are available.
      if (rsourceTriangulation%h_DvertexParameterValue .ne. ST_NOHANDLE) then
        call storage_copy (rsourceTriangulation%h_DvertexParameterValue,&
            rdestTriangulation%h_DvertexParameterValue)
      end if
      
    end subroutine tria_refineMesh2lv1D
    
    ! ---------------------------------------------------------------
  
    subroutine tria_refineMesh2lv2D(rsourceTriangulation, rdestTriangulation)

    ! This routine refines the given 2D mesh rsourceTriangulation according to
    ! the 2-level ordering algorithm. The refined mesh is saved in 
    ! rdestTriangulation. There will be no correction of boundary
    ! vertices. Boundary parameter values are not handled here!

    ! The source triangulation to be refined
    type(t_triangulation), intent(in) :: rsourceTriangulation

    ! Destination triangulation structure that receives the refined mesh
    type(t_triangulation), intent(out) :: rdestTriangulation
    
      ! local variables
      real(DP), dimension(:,:), pointer :: p_DcoordSource
      real(DP), dimension(:,:), pointer :: p_DcoordDest
      integer, dimension(:,:), pointer :: p_IvertAtElementSource
      integer, dimension(:,:), pointer :: p_IvertAtElementDest
      integer, dimension(:,:), pointer :: p_IedgesAtElementSource
      integer, dimension(:,:), pointer :: p_IvertAtEdgeSource
      integer, dimension(:), pointer :: p_IrefinementPatchIdx
      integer, dimension(:), pointer :: p_IrefinementPatch
      integer, dimension(:), pointer :: p_IcoarseGridElement
      
      integer :: nquads,iel,iel1,iel2,iel3
      integer :: ivt1,ivt2,ivtoffset,ivt,imt
      integer, dimension(2) :: Isize
      integer :: nnve,ive,NVTsrc,NMTsrc
      real(DP) :: x,y
      
      ! Get the arrays with information of the source mesh.
      call storage_getbase_double2D (rsourceTriangulation%h_DvertexCoords,&
          p_DcoordSource)
      call storage_getbase_int2D (rsourceTriangulation%h_IverticesAtElement,&
          p_IvertAtElementSource)
      call storage_getbase_int2D (rsourceTriangulation%h_IedgesAtElement,&
          p_IedgesAtElementSource)
      call storage_getbase_int2D (rsourceTriangulation%h_IverticesAtEdge,&
          p_IvertAtEdgeSource)
      
      ! Get the total number of vertices and edges
      NVTsrc = rsourceTriangulation%NVT
      NMTsrc = rsourceTriangulation%NMT
      
      ! The 2-level ordering has the following properties:
      !
      ! - Vertices in the coarse mesh are vertices in the fine mesh 
      !   with the same number. Coordinates of coarse mesh vertices
      !   stay unchanged in the fine mesh.
      ! - Edges in the coarse mesh produce vertices in the fine mesh
      !   at the midpoints of the edges in the coarse mesh.
      !   Edge numbers in the coarse mesh get vertex numbers in
      !   the fine mesh.
      ! - For quad meshes: Element midpoints in the coarse mesh
      !   get vertices in the fine mesh. They are appended to the
      !   vertices generated by the edges.
      !
      ! So at first, get the number of quads in the mesh (may be =0 which is ok).
      
      nquads = rsourceTriangulation%InelOfType (TRIA_NVEQUAD2D)
      nnve = rsourceTriangulation%NNVE
      
      if ((nnve .lt. TRIA_NVETRI2D) .or. (nnve .gt. TRIA_NVEQUAD2D)) then
      
        call output_line (&
            '2-level refinement supports only triangular and quad meshes!', &
            OU_CLASS_ERROR,OU_MODE_STD,'tria_refineMesh2lv2D')
        call sys_halt()
        
      end if
      
      ! Initialise the basic mesh data in rdestTriangulation:
      
      ! 2D mesh
      rdestTriangulation%ndim = NDIM2D
      
      rdestTriangulation%nnve = nnve
      rdestTriangulation%NNEE = rsourceTriangulation%NNEE
      
      ! Every element is divided into 4 subelements
      rdestTriangulation%NEL           = 4 * rsourceTriangulation%NEL
      rdestTriangulation%InelOfType(:) = 4 * rsourceTriangulation%InelOfType(:)

      ! We expect NVT+NMT+nquads new points.
      rdestTriangulation%NVT = rsourceTriangulation%NVT + &
                               rsourceTriangulation%NMT + nquads
      
      
      ! Allocate memory for the new vertex coordinates and
      ! get the pointers to the coordinate array
      Isize = (/NDIM2D,rdestTriangulation%NVT/)
      call storage_new ('tria_refineMesh2lv2D', 'DCORVG',&
          Isize, ST_DOUBLE, &
          rdestTriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2D(&
          rdestTriangulation%h_DvertexCoords,p_DcoordDest)
      
      ! Ok, let us start the refinement. In the first step, we copy the
      ! corner coordinates of the coarse mesh to the fine mesh; they
      ! do not change during the refinement.
      do ivt = 1, rsourceTriangulation%NVT
        p_DcoordDest(1,ivt) = p_DcoordSource(1,ivt)
        p_DcoordDest(2,ivt) = p_DcoordSource(2,ivt)
      end do
      
      ! Each edge produces an edge midpoint which is stored as new
      ! point in the fine mesh. To calculate the coordinates, take
      ! the mean of the coordinates in the coarse mesh.
      ivtoffset = rsourceTriangulation%NVT
      do imt = 1, rsourceTriangulation%NMT
        ivt1 = p_IvertAtEdgeSource (1,imt)
        ivt2 = p_IvertAtEdgeSource (2,imt)
        p_DcoordDest(1,ivtoffset+imt) = &
            0.5_DP * ( p_DcoordSource (1,ivt1) + p_DcoordSource (1,ivt2) )
        p_DcoordDest(2,ivtoffset+imt) = &
            0.5_DP * ( p_DcoordSource (2,ivt1) + p_DcoordSource (2,ivt2) )
      end do
      
      ! Allocate memory for IverticesAtElement and get a pointer to it.
      ! Fill the array with zero, so we will not have problems when mixing
      ! triangles into a quad mesh.
      Isize = (/nnve,rdestTriangulation%NEL/)
      call storage_new ('tria_refineMesh2lv2D', 'KVERT', Isize, ST_INT, &
          rdestTriangulation%h_IverticesAtElement, ST_NEWBLOCK_ZERO)
      call storage_getbase_int2D(&
          rdestTriangulation%h_IverticesAtElement,p_IvertAtElementDest)

      ! Allocate memory for the refinement information arrays.
      ! These arrays define for every coarse grid element the fine
      ! grid elements and for every fine grid element the coarse
      ! grid element where it comes from.
      call storage_new ('tria_refineMesh2lv2D', 'h_IrefinementPatchIdx', &
          rsourceTriangulation%NEL+1, ST_INT, &
          rdestTriangulation%h_IrefinementPatchIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatchIdx,p_IrefinementPatchIdx)

      call storage_new ('tria_refineMesh2lv2D', 'h_IrefinementPatch', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IrefinementPatch, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatch,p_IrefinementPatch)
    
      call storage_new ('tria_refineMesh2lv2D', 'h_IcoarseGridElement', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IcoarseGridElement, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IcoarseGridElement,p_IcoarseGridElement)
    
      ! The p_IrefinementPatchIdx array can directly be initialised
      ! as every coarse grid element gets 4 fine grid elements.
      do iel = 0, rsourceTriangulation%NEL
        p_IrefinementPatchIdx(1+iel) = 1 + iel*4
      end do
      
      ! Also the p_IcoarseGridElement array can directly be initialised.
      ! The coarse grid element number of element 1..NEL(coarse) stays
      ! the same. The number of the other coarse grid elements can
      ! be calculated by a formula.
      do iel = 1,rsourceTriangulation%NEL
        p_IcoarseGridElement(iel) = iel
      end do

      do iel = rsourceTriangulation%NEL+1, rdestTriangulation%NEL
        p_IcoarseGridElement(iel) = (iel-rsourceTriangulation%NEL-1)/3 + 1
      end do
    
      ! Are there quads in the mesh? They produce additional element midpoints
      ! which get numbers NVT+NMT+1..*
      if (nnve .gt. TRIA_NVETRI2D) then
      
        ! Loop over the elements to find the quads.
        ivtoffset = rsourceTriangulation%NVT + rsourceTriangulation%NMT
        do iel = 1, rsourceTriangulation%NEL
        
          if (p_IvertAtElementSource (TRIA_NVEQUAD2D,iel) .ne. 0) then
          
            ! New element midpoint
            ivtoffset = ivtoffset+1
            
            ! Sum up the coordinates of the corners to get the midpoint.
            x = 0.0_DP
            y = 0.0_DP
            do ive = 1,TRIA_NVEQUAD2D
              x = x + p_DcoordSource (1,p_IvertAtElementSource(ive,iel))
              y = y + p_DcoordSource (2,p_IvertAtElementSource(ive,iel))
            end do

            ! Store the midpoint
            p_DcoordDest(1,ivtoffset) = 0.25_DP*x
            p_DcoordDest(2,ivtoffset) = 0.25_DP*y
          
          end if
        
        end do
      
      end if
      
      ! Ok, that was easy up to here. Now the interesting part: 
      ! the two-level ordering.
      !
      ! Look at the following triangular and quad element with local
      ! and global vertex and edge numbers:
      !
      !    102                               104           203           103
      !       X__                               X-----------------------X
      !       |2 \__                            |4                     3|
      !       |     \__                         |                       |
      !       |        \__201                   |                       |
      !    202|           \__                204|          IEL          |202
      !       |     IEL      \__                |                       |
      !       |                 \__             |                       |
      !       | 3                 1\_           |1                     2|
      !       X----------------------+X         X-----------------------X
      !    103          203           101    101           201           102
      !
      ! The two-level ordering refinement strategy says now:
      !  - Every edge produces an edge midpoint. Connect opposite midpoints to get
      !    the fine grid.
      !  - The vertex numbers in the coarse grid are transferred to the fine grid
      !    without change.
      !  - The edge numbers in the coarse grid define the vertex numbers of the
      !    edge midpoints in the fine grid.
      !  - For triangles: The element number IEL in the coarse grid is transferred
      !    to the inner element of the three subelements in the fine grid.
      !    The midpoint of edge 1 gets local vertex 1 of the inner triangle.
      !    The other three elements get the numbers NEL+3*(iel-1)+1, +2 and +3
      !    according to which edge of the inner triangle they touch.
      !    The first vertex of each subtriangle is the vertex with
      !    local number 1, 2 or 3, respectively, in the coarse grid.
      !  - For quads: The element number IEL in the coarse grid is transferred
      !    to the subelement at local vertex 1. The other elements
      !    get element numbers NEL+3*(iel-1)+1, +2 and +3 in counterclockwise
      !    order. The first vertex on every subelement is the old (corner) vertex
      !    in the coarse mesh.
      !
      ! In the above picture, we therefore have:
      !
      !    102                               104           203           103
      !       X__                               X-----------X-----------X
      !       |1 \__                            |1          |          1|
      !       | IEL2\__                         |    IEL3   |   IEL2    |
      !       |        \_ 201                   |           |           |
      !    202X----------X___                204X-----------X-----------X202
      !       | \__ IEL 1|   \__                |           |           |
      !       |    \__   |      \__             |    IEL    |   IEL1    |
      !       |1 IEL3 \__|  IEL1  1\_           |1          |          1|
      !       X----------X-----------X          X-----------X-----------X
      !    103          203           101    101           201           102
      !
      ! IEL1 = NEL+3*(IEL-1)+1
      ! IEL2 = NEL+3*(IEL-1)+2
      ! IEL3 = NEL+3*(IEL-1)+3
      !
      ! Ok, to produce all that, we have to set up IverticesOnElement correctly!
      ! Let us loop over the vertices on the coarse grid. Everyone produces
      ! four elements on the fine grid:
      !
      ! In nquads we count the number of quads +NVT+NMT we reach.
      ! That is the number of the midpoint of that element!
      nquads = rsourceTriangulation%NVT+rsourceTriangulation%NMT
      
      if (nnve .eq. TRIA_NVETRI2D) then
      
        ! Pure triangle mesh
        do iel = 1, rsourceTriangulation%NEL
      
          ! Determine number of subelements.
          iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
          iel2 = iel1+1
          iel3 = iel1+2
          
          ! Save them
          p_IrefinementPatch(1+(iel-1)*4)   = iel
          p_IrefinementPatch(1+(iel-1)*4+1) = iel1
          p_IrefinementPatch(1+(iel-1)*4+2) = iel2
          p_IrefinementPatch(1+(iel-1)*4+3) = iel3
          
          ! Step 1: Initialise IverticesOnElement for element IEL
          p_IvertAtElementDest(1,iel) = p_IedgesAtElementSource (1,iel)+NVTsrc
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (2,iel)+NVTsrc
          p_IvertAtElementDest(3,iel) = p_IedgesAtElementSource (3,iel)+NVTsrc
          
          ! Step 2: Initialise IverticesOnElement for element IEL1
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (1,iel)
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (1,iel)+NVTsrc
          p_IvertAtElementDest(3,iel1) = p_IedgesAtElementSource (3,iel)+NVTsrc

          ! Step 3: Initialise IverticesOnElement for element IEL2
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (2,iel)
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (2,iel)+NVTsrc
          p_IvertAtElementDest(3,iel2) = p_IedgesAtElementSource (1,iel)+NVTsrc

          ! Step 4: Initialise IverticesOnElement for element IEL3
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (3,iel)
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (3,iel)+NVTsrc
          p_IvertAtElementDest(3,iel3) = p_IedgesAtElementSource (2,iel)+NVTsrc
        
        end do
      
      elseif (nquads .eq. rsourceTriangulation%NEL) then
        
        ! Pure QUAD mesh
        do iel = 1, rsourceTriangulation%NEL
      
          ! Determine number of subelements.
          iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
          iel2 = iel1+1
          iel3 = iel1+2
          
          ! Save them
          p_IrefinementPatch(1+(iel-1)*4)   = iel
          p_IrefinementPatch(1+(iel-1)*4+1) = iel1
          p_IrefinementPatch(1+(iel-1)*4+2) = iel2
          p_IrefinementPatch(1+(iel-1)*4+3) = iel3

          ! In nquads we count the number of the quad +NVT+NMT we process.
          ! That is the number of the midpoint of that element!
          ! As we reached a new quad, we increase nquads
          
          !nquads = nquads+1 !replaced because of OpenMP directive
          nquads = rsourceTriangulation%NVT + rsourceTriangulation%NMT + iel
          
          ! Step 1: Initialise IverticesOnElement for element IEL
          p_IvertAtElementDest(1,iel) = p_IvertAtElementSource (1,iel)
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (1,iel)+NVTsrc
          p_IvertAtElementDest(3,iel) = nquads
          p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource (4,iel)+NVTsrc
          
          ! Step 2: Initialise IverticesOnElement for element IEL1
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (2,iel)
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (2,iel)+NVTsrc
          p_IvertAtElementDest(3,iel1) = nquads
          p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource (1,iel)+NVTsrc
        
          ! Step 3: Initialise IverticesOnElement for element IEL2
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (3,iel)
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (3,iel)+NVTsrc
          p_IvertAtElementDest(3,iel2) = nquads
          p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource (2,iel)+NVTsrc

          ! Step 4: Initialise IverticesOnElement for element IEL3
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (4,iel)
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (4,iel)+NVTsrc
          p_IvertAtElementDest(3,iel3) = nquads
          p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource (3,iel)+NVTsrc
        
        end do ! iel
        
      else
      
        ! Triangles and quads mixed.
      
        do iel = 1, rsourceTriangulation%NEL
      
          ! Is that a triangle or a quad?  
          if (p_IvertAtElementSource(TRIA_NVEQUAD2D,iel) .eq. 0) then
          
            ! Triangular element.
            !
            ! Determine number of subelements.
            iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
            iel2 = iel1+1
            iel3 = iel1+2
            
            ! Save them
            p_IrefinementPatch(1+(iel-1)*4)   = iel
            p_IrefinementPatch(1+(iel-1)*4+1) = iel1
            p_IrefinementPatch(1+(iel-1)*4+2) = iel2
            p_IrefinementPatch(1+(iel-1)*4+3) = iel3

            ! Step 1: Initialise IverticesOnElement for element IEL
            p_IvertAtElementDest(1,iel) = p_IedgesAtElementSource (1,iel)+NVTsrc
            p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (2,iel)+NVTsrc
            p_IvertAtElementDest(3,iel) = p_IedgesAtElementSource (3,iel)+NVTsrc
            
            ! Step 2: Initialise IverticesOnElement for element IEL1
            p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (1,iel)
            p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (1,iel)+NVTsrc
            p_IvertAtElementDest(3,iel1) = p_IedgesAtElementSource (3,iel)+NVTsrc

            ! Step 3: Initialise IverticesOnElement for element IEL2
            p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (2,iel)
            p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (2,iel)+NVTsrc
            p_IvertAtElementDest(3,iel2) = p_IedgesAtElementSource (1,iel)+NVTsrc

            ! Step 4: Initialise IverticesOnElement for element IEL3
            p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (3,iel)
            p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (3,iel)+NVTsrc
            p_IvertAtElementDest(3,iel3) = p_IedgesAtElementSource (2,iel)+NVTsrc
          
          else
          
            ! Quadrilateral element
            !
            ! Determine number of subelements.
            iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
            iel2 = iel1+1
            iel3 = iel1+2
            
            ! Save them
            p_IrefinementPatch(1+(iel-1)*4)   = iel
            p_IrefinementPatch(1+(iel-1)*4+1) = iel1
            p_IrefinementPatch(1+(iel-1)*4+2) = iel2
            p_IrefinementPatch(1+(iel-1)*4+3) = iel3

            ! In nquads we count the number of the quad +NVT+NMT we process.
            ! That is the number of the midpoint of that element!
            ! As we reached a new quad, we increase nquads
            nquads = nquads+1
            
            ! Step 1: Initialise IverticesOnElement for element IEL
            p_IvertAtElementDest(1,iel) = p_IvertAtElementSource (1,iel)
            p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (1,iel)+NVTsrc
            p_IvertAtElementDest(3,iel) = nquads
            p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource (4,iel)+NVTsrc
            
            ! Step 2: Initialise IverticesOnElement for element IEL1
            p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (2,iel)
            p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (2,iel)+NVTsrc
            p_IvertAtElementDest(3,iel1) = nquads
            p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource (1,iel)+NVTsrc
          
            ! Step 3: Initialise IverticesOnElement for element IEL2
            p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (3,iel)
            p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (3,iel)+NVTsrc
            p_IvertAtElementDest(3,iel2) = nquads
            p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource (2,iel)+NVTsrc

            ! Step 4: Initialise IverticesOnElement for element IEL3
            p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (4,iel)
            p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (4,iel)+NVTsrc
            p_IvertAtElementDest(3,iel3) = nquads
            p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource (3,iel)+NVTsrc
          
          end if
        
        end do ! iel
        
      end if
      
      ! The last step of setting up the raw mesh on the finer level:
      ! Set up InodalProperty. But that is the most easiest thing: Simply
      ! copy the nodal property array from the coarse mesh to the fine mesh.
      ! The nodal information of the edges (number NVT+1..NVT+NMT)
      ! that way converts into the nodal information about the new vertices
      ! on the fine mesh!
      call storage_copy (rsourceTriangulation%h_InodalProperty,&
                         rdestTriangulation%h_InodalProperty)
    
    end subroutine tria_refineMesh2lv2D

    ! ---------------------------------------------------------------
  
    subroutine tria_refineBdry2lv2D(rsourceTriangulation, rdestTriangulation,&
                                    brecalcBoundaryCoords, rboundary)

    ! This routine refines the boundary definition of rsourceTriangulation
    ! according to the 2-level ordering algorithm to generate a new 
    ! IverticesAtBoundary. 
    ! If rboundary is specified, the parameter values of the boundary vertices are 
    ! updated and the coordinates of the boundary points are corrected according 
    ! to the analytic boundarty. the boundary vertices are not sorted for their
    ! parameter value!

    ! The source triangulation to be refined
    type(t_triangulation), intent(in) :: rsourceTriangulation

    ! Destination triangulation structure that receives the refined mesh
    type(t_triangulation), intent(inout) :: rdestTriangulation
    
    ! Recalculate the coordinates of all boundary vertices according to
    ! their parameter value.
    logical, intent(in) :: brecalcBoundaryCoords
    
    ! OPTIONAL: Defintion of analytic boundary.
    ! If specified, the coordinates of the new boundary vertices are
    ! recomputed according to the analytic boundary.
    type(t_boundary), intent(in), optional :: rboundary
    
      ! local variables
      real(DP), dimension(:), pointer :: p_DvertParamsSource
      real(DP), dimension(:), pointer :: p_DedgeParamsSource
      real(DP), dimension(:), pointer :: p_DvertParamsDest
      real(DP), dimension(:,:), pointer :: p_DcornerCoordDest
      integer, dimension(:), pointer :: p_IvertAtBoundartySource
      integer, dimension(:), pointer :: p_IedgesAtBoundartySource
      integer, dimension(:), pointer :: p_IvertAtBoundartyDest
      integer, dimension(:), pointer :: p_IboundaryCpIdxSource
      integer, dimension(:), pointer :: p_IboundaryCpIdxDest
      integer, dimension(:), allocatable :: IverticesAtBoundaryTmp
      integer :: ivbd,ibct,ivtpos,ivtsource,ivtdest,ivtstart,ibc
      
      ! Get the definition of the boundary vertices and -edges.
      call storage_getbase_int(rsourceTriangulation%h_IverticesAtBoundary,&
          p_IvertAtBoundartySource)
      call storage_getbase_int(rsourceTriangulation%h_IedgesAtBoundary,&
          p_IedgesAtBoundartySource)
      call storage_getbase_int(rsourceTriangulation%h_IboundaryCpIdx,&
          p_IboundaryCpIdxSource)

      ! The number of boundary vertices is doubled.
      rdestTriangulation%NVBD = 2*rsourceTriangulation%NVBD
      rdestTriangulation%NBCT = rsourceTriangulation%NBCT 
      rdestTriangulation%NblindBCT = rsourceTriangulation%NblindBCT 
          
      ! Create new arrays in the fine grid for the vertices and indices.
      call storage_new('tria_refineBdry2lv2D', 'KVBD',&
          rdestTriangulation%NVBD,  ST_INT,&
          rdestTriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)

      call storage_new('tria_generateBasicBoundary', 'KBCT',&
          size(p_IboundaryCpIdxSource), ST_INT,&
          rdestTriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT)

      call storage_getbase_int(&
          rdestTriangulation%h_IverticesAtBoundary, p_IvertAtBoundartyDest)
      call storage_getbase_int(&
          rdestTriangulation%h_IboundaryCpIdx, p_IboundaryCpIdxDest)

      ! The 2-level ordering algorithm says:
      !  - Every edge produces an edge midpoint. Connect opposite midpoints to get
      !    the fine grid.
      !  - The vertex numbers in the coarse grid are transferred to the fine grid
      !    without change.
      !  - The edge numbers in the coarse grid define the vertex numbers of the
      !    edge midpoints in the fine grid.
      ! Therefore, we have to interleave the vertex and edge numbers to get the
      ! new IverticesAtBoundary. The IboundaryCpIdx can be transferred, but we have
      ! to multiply each index by 2 as we have twice as many vertices per boundary
      ! component now.
      
      do ibct = 1, size(p_IboundaryCpIdxDest)
        p_IboundaryCpIdxDest(ibct) = 2*(p_IboundaryCpIdxSource(ibct)-1) + 1
      end do
      do ivbd = 0, size(p_IvertAtBoundartySource)-1
        p_IvertAtBoundartyDest(1+2*ivbd) = p_IvertAtBoundartySource(1+ivbd)
        p_IvertAtBoundartyDest(2+2*ivbd) = p_IedgesAtBoundartySource(1+ivbd)&
                                         + rsourceTriangulation%NVT
      end do
      
      ! Let us see if parameter values of boundary vertices are available.
      if ((rsourceTriangulation%h_DvertexParameterValue .ne. ST_NOHANDLE) .and. &
          (rsourceTriangulation%h_DedgeParameterValue .ne. ST_NOHANDLE)) then

        ! Also interleave the parameter values of the vertices and edge midpoints.
        call storage_getbase_double(&
            rsourceTriangulation%h_DvertexParameterValue, p_DvertParamsSource)
        call storage_getbase_double(&
            rsourceTriangulation%h_DedgeParameterValue, p_DedgeParamsSource)
        
        ! Create a new array for the boundary vertex parameter values
        ! and fill it with the data from the coarse grid.
        call storage_new ('tria_refineBdry2lv2D', 'DVBDP',&
            rdestTriangulation%NVBD, ST_DOUBLE,&
            rdestTriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double(&
            rdestTriangulation%h_DvertexParameterValue, p_DvertParamsDest)
            
        do ivbd = 0,size(p_IvertAtBoundartySource)-1
          p_DvertParamsDest(1+2*ivbd) = p_DvertParamsSource(1+ivbd)
          p_DvertParamsDest(2+2*ivbd) = p_DedgeParamsSource(1+ivbd)
        end do
        
        ! Now we have an array like
        !  1.0 2.0 -1.0 3.0    1.0 -1.0 2.0  ...
        ! With some -1 inbetween for all vertices that do not belong to the physical
        ! boundary. We now 'shift' all these '-1'-nodes to the end of the array
        ! and that way assign them to the 'blind' boundary component.
        
        allocate(IverticesAtBoundaryTmp(size(p_IvertAtBoundartyDest)))
        ivtdest = 1
        ivtpos = 1
        ivtstart = p_IboundaryCpIdxDest(1)
        
        ! Loop over all physical-boundary-BC`s; ignore any existing 'blind'
        ! boundary component.
        do ibc = 1, rdestTriangulation%NBCT
        
          ! Loop over all boundary vertices in that BC. Compress the BC.
          ! All boundary vertices not on the physical boundary are
          ! extracted to IverticesAtBoundaryTmp.
          
          do ivtsource = ivtstart, p_IboundaryCpIdxDest(ibc+1)-1
            if (p_DvertParamsDest(ivtsource) .ne. -1.0_DP) then
              ! This is a vertex on the boundary. Copy it to the destination position.
              p_DvertParamsDest(ivtdest) = p_DvertParamsDest(ivtsource)
              p_IvertAtBoundartyDest(ivtdest) = p_IvertAtBoundartyDest(ivtsource)
              ivtdest = ivtdest + 1
            else
              ! Extract that DOF, so we can overwrite with the forthcoming DOF`s.
              IverticesAtBoundaryTmp(ivtpos) = p_IvertAtBoundartyDest(ivtsource)
              ivtpos = ivtpos + 1
            end if
          end do
          
          ! Remember the new start address where the DOF`s of the next
          ! boundary component are found; it is immediately overwritten
          ! by ivtdest!
          ivtstart = p_IboundaryCpIdxDest(ibc+1)
          
          ! ivtdest points now behind the DOF`s of boundary component ibc.
          ! Save the pointer to the p_IboundaryCpIdx array.
          p_IboundaryCpIdxDest(ibc+1) = ivtdest
          
        end do
        
        ! Ok, the array is compressed now. The only thing that remains is
        ! to copy the 'blind' vertices to the 'blind' boundary component.
        do ivtsource = 1, ivtpos-1
          p_IvertAtBoundartyDest(ivtdest) = IverticesAtBoundaryTmp(ivtsource)
          p_DvertParamsDest(ivtdest) = -1.0_DP
          ivtdest = ivtdest + 1
        end do
        
        deallocate(IverticesAtBoundaryTmp)
        
        ! The last entry in p_IboundaryCpIdx does not have to be changed since
        ! it counts the total number of vertices on the boundary -- which
        ! has not changed.        
        
        ! If the analytic boundary is given, compute the coordinates of the
        ! boundary vertices from that.
        if (present(rboundary) .and. brecalcBoundaryCoords) then
          
          ! Get the array with the vertex coordinates.
          ! We want to correct the coordinates of the boundary points
          ! according to the analytic boundary.
          call storage_getbase_double2d(&
              rdestTriangulation%h_DvertexCoords, p_DcornerCoordDest)
            
          ! Loop through the boundary points and calculate the correct
          ! coordinates. Points that are not on the physical bondary but
          ! only on the boundary of a subdomain (identified by parameter value -1)
          ! must stay where they are. 
          do ibct = 1, rdestTriangulation%NBCT
            do ivbd = p_IboundaryCpIdxDest(ibct), p_IboundaryCpIdxDest(ibct+1)-1
              if (p_DvertParamsDest(ivbd) .ge. 0.0_DP) then
                call boundary_getCoords(rboundary, ibct, p_DvertParamsDest(ivbd),&
                    p_DcornerCoordDest(1,p_IvertAtBoundartyDest(ivbd)),&
                    p_DcornerCoordDest(2,p_IvertAtBoundartyDest(ivbd)))
              end if
            end do
          end do
          
        end if
        
      end if
    
    end subroutine tria_refineBdry2lv2D

    ! ---------------------------------------------------------------
  
    subroutine tria_averageMidpoints2D(rsourceTriangulation,rdestTriangulation)

    ! The follwing function corrects the grid on domains where non-linear/
    ! curved boundary segments are used. 
    ! When a boundary segment of the domain is a line, new boundary nodes
    ! are automatically positioned on the boundary. But when the boundary
    ! is a curve, a new boundary vertex is not automatically on the
    ! boundary, but it first has tobe moved to there. This is already
    ! performed in the refinement routine XSB0X. Unfortunately this
    ! procedure can lead to very anisotropic elements near the boundary,
    ! depending on how sharp the curved boundary is.
    !
    ! tria_averageMidpoints2D now tries to reduce these effects of anisotropy. 
    ! The element midpoint of the coarser element (which is at the same time 
    ! the vertex where all the four finer elements meet) is taken as the 
    ! average of the four edge midpoints that arise from natural refinement.

    ! The source triangulation specifying the coarse mesh
    type(t_triangulation), intent(in) :: rsourceTriangulation

    ! Destination triangulation that specifies the fine mesh after
    ! 2-level refinement
    type(t_triangulation), intent(inout) :: rdestTriangulation

      ! local variables
      integer :: i
      real(DP) :: dx,dy
      integer :: iel,ivtoffset !,ielbd
      !INTEGER, dimension(:), POINTER :: p_IelementsAtBoundaryCoarse
      integer, dimension(:,:), pointer :: p_IverticesAtElementCoarse
      integer, dimension(:,:), pointer :: p_IverticesAtElementFine
      integer, dimension(:,:), pointer :: p_IedgesAtElementCoarse
      real(DP), dimension(:,:), pointer :: p_DvertexCoordsFine
      integer :: ipt

      ! The approach is simple:
      !    
      ! o----------o----------o
      ! \          \          |
      ! |\         \          | 
      ! | \      -> \         |
      ! |  o---------o--------o   Averaging of the element midpoint by 
      ! | /      -> /         |   interpolation
      ! |/         /          |
      ! /          /          |
      ! o----------o----------o
      !
      ! Actually, we only have to tackle the elements on the boundary of the coarse
      ! mesh. Loop over them. Some midpoints may be calculated
      ! twice that way, but we do not care...
      !
      ! Unfortunately, due to the incompleteness of the triangulation structure,
      ! this approach does not work for mixed mesh. We have to loop through
      ! all elements on the coarse mesh to find the numbers of the quads :-(
      
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IverticesAtElement, p_IverticesAtElementCoarse)
          
      ! If this is a pure triangle mesh, there is nothing to do.
      if (ubound(p_IverticesAtElementCoarse,1) .le. TRIA_NVETRI2D) return
      
      ! call storage_getbase_int (rsourceTriangulation%h_IelementsAtBoundary,&
      !     p_IelementsAtBoundaryCoarse)

      call storage_getbase_int2d(&
          rsourceTriangulation%h_IedgesAtElement, p_IedgesAtElementCoarse)
      call storage_getbase_int2d(&
          rdestTriangulation%h_IverticesAtElement, p_IverticesAtElementFine)
      call storage_getbase_double2d(&
          rdestTriangulation%h_DvertexCoords, p_DvertexCoordsFine)

      ivtoffset = rsourceTriangulation%NVT + rsourceTriangulation%NMT
      
      do iel = 1, rsourceTriangulation%NEL

        ! Get the element number of the coarse mesh element
        ! iel = p_IelementsAtBoundaryCoarse(ielbd)
      
        ! Triangle or quad?
        if (p_IverticesAtElementCoarse(TRIA_NVEQUAD2D,iel) .ne. 0) then
        
          ! New quad. Increase ivtoffset, this is then the number of the
          ! point that was the midpoint in the coarse mesh.
          ivtoffset = ivtoffset + 1
        
          ! Ok, that is a quad. The midpoint must be recalculated.
          ! We do this based on the 'edge midpoints' which are now
          ! vertices in the fine mesh. The coordinates of these points
          ! may differ from the coordinates calculated by linear 
          ! interpolation due to boundary adjustment!
          dx = 0.0_DP
          dy = 0.0_DP
          do i = 1, TRIA_NVEQUAD2D
            ipt = p_IedgesAtElementCoarse(i,iel)+rsourceTriangulation%NVT
            dx = dx + p_DvertexCoordsFine(1,ipt)
            dy = dy + p_DvertexCoordsFine(2,ipt)
          end do
          
          ! Save the vertex coordinates. 
          p_DvertexCoordsFine(1,ivtoffset) = 0.25_DP*dx
          p_DvertexCoordsFine(2,ivtoffset) = 0.25_DP*dy
        
        end if
      
      end do
    
    end subroutine tria_averageMidpoints2D

    ! ---------------------------------------------------------------
    
    subroutine tria_refineMesh2lv3D(rsourceTriangulation,rdestTriangulation)

      ! This routine refines the given 3D mesh rsourceTriangulation
      ! according to the 2-level ordering algorithm.
      ! The refined mesh is saved in rdestTriangulation.
      ! There will be no correction of boundary vertices.
      ! Boundary parameter values are not handled here!
      
      ! The source triangulation to be refined
      type(t_triangulation), intent(in) :: rsourceTriangulation

      ! Destination triangulation structure that receives the refined mesh
      type(t_triangulation), intent(out) :: rdestTriangulation
  
      ! local variables
      
      real(DP), dimension(:,:), pointer :: p_DcoordSource
      real(DP), dimension(:,:), pointer :: p_DcoordDest
      integer, dimension(:,:), pointer :: p_IvertAtElementSource
      integer, dimension(:,:), pointer :: p_IvertAtElementDest
      integer, dimension(:,:), pointer :: p_IedgesAtElementSource
      integer, dimension(:,:), pointer :: p_IvertAtEdgeSource
      integer, dimension(:,:), pointer :: p_IverticesAtFace
      integer, dimension(:,:), pointer :: p_IfacesAtElement
      integer, dimension(:), pointer :: p_IfacesAtBoundary
      integer, dimension(:), pointer :: p_InodalPropertySource
      integer, dimension(:), pointer :: p_InodalPropertyDest
      integer, dimension(:), pointer :: p_IrefinementPatch
      integer, dimension(:), pointer :: p_IrefinementPatchIdx
      integer, dimension(:), pointer :: p_IcoarseGridElement
      integer, dimension(:), pointer :: p_ImidpointAtFace
      integer, dimension(2) :: Isize
      integer :: iel,iel1,iel2,iel3,iel4,iel5,iel6,iel7,iel8,iel9
      integer :: ihex,nhexs,ntets,npris,npyrs,nquads
      integer :: ivt1,ivt2,ivtoffset,ivt,iae,imt,idx,iee,nee
      integer :: ieloffset,idxoffset,midPointOfIel,iabd
      integer :: midpointFaceA,midpointFaceB,midpointFaceC
      integer :: midpointFaceD,midpointFaceE,midpointFaceF
      integer :: ive,nve,nnve,NVTsrc,NMTsrc,NELsrc
      integer :: h_ImidpointAtFace
      real(DP) :: x,y,z
      
      ! Get the arrays with information of the source mesh.
      call storage_getbase_double2d(&
          rsourceTriangulation%h_DvertexCoords, p_DcoordSource)
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IverticesAtElement, p_IvertAtElementSource)
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IedgesAtElement, p_IedgesAtElementSource)
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IverticesAtEdge, p_IvertAtEdgeSource)
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IverticesAtFace, p_IverticesAtFace)
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IfacesAtElement, p_IfacesAtElement)
      call storage_getbase_int(&
          rsourceTriangulation%h_IfacesAtBoundary, p_IfacesAtBoundary)
      
      ! Get the total number of vertices, edges, and elements
      NVTsrc = rsourceTriangulation%NVT
      NMTsrc = rsourceTriangulation%NMT
      NELsrc = rsourceTriangulation%NEL

      ! The 2-level ordering has the following properties:
      !
      ! - Vertices in the coarse mesh are vertices in the fine mesh 
      !   with the same number. Coordinates of coarse mesh vertices
      !   stay unchanged in the fine mesh.
      ! - Edges in the coarse mesh produce vertices in the fine mesh
      !   at the midpoints of the edges in the coarse mesh.
      !   Edge numbers in the coarse mesh get vertex numbers in
      !   the fine mesh.
      ! - For hexahedral meshes: Element midpoints in the coarse mesh
      !   get vertices in the fine mesh. They are appended to the
      !   vertices generated by the edges.
      
      nhexs = rsourceTriangulation%InelOfType(TRIA_NVEHEXA3D)
      ntets = rsourceTriangulation%InelOfType(TRIA_NVETET3D)
      npris = rsourceTriangulation%InelOfType(TRIA_NVEPRIS3D)
      npyrs = rsourceTriangulation%InelOfType(TRIA_NVEPYR3D)
      nnve  = rsourceTriangulation%NNVE

      ! Initialise the basic mesh data in rdestTriangulation:
      
      ! 3D mesh
      rdestTriangulation%ndim = NDIM3D
      
      rdestTriangulation%nnve = nnve
      rdestTriangulation%NNEE = rsourceTriangulation%NNEE
     
      ! Calculate the number of elements after regular subdivision:
      ! - hexahedra, tetrahedra, and prisms
      !   are subdivided into 8 elements
      ! - pyramides are subdivided into 10 elements
      !   (6 pyramids + 4 tetrahedra)
      rdestTriangulation%NEL = 8 * (nhexs+ntets+npris) + 10 * npyrs

      rdestTriangulation%InelOfType(TRIA_NVEHEXA3D) = &
          8 * rsourceTriangulation%InelOfType(TRIA_NVEHEXA3D)
      rdestTriangulation%InelOfType(TRIA_NVETET3D) = &
          8 * rsourceTriangulation%InelOfType(TRIA_NVETET3D) +&
          4 * rsourceTriangulation%InelOfType(TRIA_NVEPYR3D)
      rdestTriangulation%InelOfType(TRIA_NVEPRIS3D) = &
          8 * rsourceTriangulation%InelOfType(TRIA_NVEPRIS3D)
      rdestTriangulation%InelOfType(TRIA_NVEPYR3D) = &
          6 * rsourceTriangulation%InelOfType(TRIA_NVEPYR3D)
      
      ! We have now found rdestTriangulation%NEL, so we can
      ! allocate memory for the refinement information arrays.
      ! These arrays define for every coarse grid element the fine
      ! grid elements and for every fine grid element the coarse
      ! grid element where it comes from.
      call storage_new('tria_refineMesh2lv3D', 'h_IrefinementPatchIdx', &
          NELsrc+1, ST_INT, &
          rdestTriangulation%h_IrefinementPatchIdx, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatchIdx,p_IrefinementPatchIdx)
      
      call storage_new('tria_refineMesh2lv3D', 'h_IrefinementPatch', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IrefinementPatch, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatch,p_IrefinementPatch)
      
      call storage_new('tria_refineMesh2lv3D', 'h_IcoarseGridElement', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IcoarseGridElement, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(&
          rdestTriangulation%h_IcoarseGridElement,p_IcoarseGridElement)
      
      ! Tetrahedra, hexahedra and prisms are refined into 8 elements of
      ! the same type so that no distinction is required in this phase
      ! of the refinement process. Pyramids are refined into 6 pyramids
      ! and 4 tetrahedral elements so that we have to consider this 
      ! special case separately and compute all data step-by-step.
      if (npyrs .eq. 0) then
        
        ! The p_IrefinementPatchIdx array can be directly initialised
        ! as every coarse grid element gets 8 fine grid elements.
        do iel = 0, NELsrc
          p_IrefinementPatchIdx(1+iel) = 1 + iel*8
        end do
        
        ! Also the p_IcoarseGridElement array can directly be initialised.
        ! The coarse grid element number of element 1..NEL(coarse) stays
        ! the same. The number of the other coarse grid elements can
        ! be calculated by a formula.
        do iel = 1, NELsrc
          p_IcoarseGridElement(iel) = iel
        end do
        
        ! assign the other element numbers
        ! the number for an element is iel1 = rsourceTriangulation%NEL+7*(iel-1)+1  
        do iel = NELsrc+1, rdestTriangulation%NEL
          p_IcoarseGridElement(iel) = (iel-NELsrc-1)/7 + 1
        end do
        
      else

        ! The p_IrefinementPatchIdx array has to be initialised
        p_IrefinementPatchIdx(1) = 1

        ! The index offset has to be initialised
        idxoffset = 1

        ! The element offset has to be initialised
        ieloffset = NELsrc
        
        ! and updated step-by-step depending on the type of element
        do iel = 1, NELsrc
          
          ! Get number of vertices in element and
          ! increase the element offset accordingly
          nve = tria_getNVE(p_IvertAtElementSource, iel)
          nee = merge(10, 8, nve .eq. TRIA_NVEPYR3D)
          idxoffset = idxoffset + nee
          
          p_IrefinementPatchIdx(1+iel) = idxoffset

          ! The coarse grid element number of element 
          ! 1..NEL(coarse) stays the same.
          p_IcoarseGridElement(iel) = iel

          ! The number of the other coarse grid elements can
          ! be calculated depending on the type of element
          do iee = 1, nee-1
            p_IcoarseGridElement(ieloffset+iee) = iel
          end do

          ! Increase element offset accordingly
          ieloffset = ieloffset + nee-1

        end do

      end if

      ! Calculate the number of vertices after regular subdivision.
      ! That is a little bit tricky since only quadrilateral faces
      ! give rise to new vertices. Here is the topology of elements:
      ! - hexahedra:  8 vertices + 12 edges + 6 faces + 1 = 27 vertices
      ! - tetrahedra: 4 vertices +  6 edges               = 10 vertices
      ! - prisms:     6 vertices +  9 edges               = 15 vertices
      ! - pyramids:   5 vertices +  8 edges + 1 face      = 14 vertices
      !
      ! Note that faces are addressed from both adjacent elements if they
      ! are located in the interior of the domain and only once if they
      ! are located at the boundary. Thus the total number of faces equals:
      !   NAT = (6*NHEXs + 5*NPRIs + 5*NPYRs + 4*NTETs + NABD) / 2
      ! The same formula holds for triangular and/or quadrilateral faces
      ! if the coefficients are modified and NABD is replaced by the number
      ! of triangular and/or quadrilateral faces at the boundary. Hence,
      ! we have to count the number of quadrilateral faces at the boundary.
      nquads = 0
      do iabd = 1, rsourceTriangulation%NABD
        ! Get the global face number
        iae = p_IfacesAtBoundary(iabd)

        ! Check if 4th vertex is not empty
        if (p_IverticesAtFace(4,iae) .ne. 0) nquads = nquads+1
      end do

      ! Then we can compute the number of vertices
      rdestTriangulation%NVT = NVTsrc + NMTsrc + nhexs +&
                               (6*nhexs + 3*npris + npyrs + nquads)/2

      ! Allocate memory for the new vertex coordinates and
      ! get the pointers to the coordinate array
      Isize = (/NDIM3D,rdestTriangulation%NVT/)
      call storage_new('tria_refineMesh2lv3D', 'DCORVG',&
          Isize, ST_DOUBLE, &
          rdestTriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2D(&
          rdestTriangulation%h_DvertexCoords,p_DcoordDest)

      ! Allocate auxiliary memory for the vertex numbers at face midpoints
      h_ImidpointAtFace = ST_NOHANDLE
      call storage_new('tria_refineMesh2lv3D', 'ImidpointAtFace',&
          rsourceTriangulation%NAT, ST_INT,&
          h_ImidpointAtFace, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(h_ImidpointAtFace, p_ImidpointAtFace)
      
      ! Ok, let us start the refinement. In the first step, we copy the
      ! corner coordinates of the coarse mesh to the fine mesh; they
      ! do not change during the refinement.
      do ivt = 1, NVTsrc
        p_DcoordDest(1,ivt) = p_DcoordSource(1,ivt)
        p_DcoordDest(2,ivt) = p_DcoordSource(2,ivt)
        p_DcoordDest(3,ivt) = p_DcoordSource(3,ivt)
      end do
      
      ! Each edge produces an edge midpoint which is stored as new
      ! point in the fine mesh. To calculate the coordinates, take
      ! the mean of the coordinates in the coarse mesh.
      
      ! This applies to ALL types of elements.
      do imt = 1, NMTsrc
        ivt1 = p_IvertAtEdgeSource (1,imt)
        ivt2 = p_IvertAtEdgeSource (2,imt)
        p_DcoordDest(1,NVTsrc+imt) = &
            0.5_DP * ( p_DcoordSource (1,ivt1) + p_DcoordSource (1,ivt2) )
        p_DcoordDest(2,NVTsrc+imt) = &
            0.5_DP * ( p_DcoordSource (2,ivt1) + p_DcoordSource (2,ivt2) )
        p_DcoordDest(3,NVTsrc+imt) = &
            0.5_DP * ( p_DcoordSource (3,ivt1) + p_DcoordSource (3,ivt2) )     
      end do
      
      ! Initialise the offset for vertex numbers
      ivtoffset = NVTsrc + NMTsrc
      
      ! Each midpoint of a quadrilateral face produces a vertex on
      ! the next level. This applies to hexahedral elements and to
      ! the single quadrilateral bottom face of pyramidal elements
      do iae = 1, rsourceTriangulation%NAT

        ! Check if the current face is quadrilateral, 
        ! that is, it features four corner vertices
        if (p_IverticesAtFace(4,iae) .eq. 0) cycle

        ! Increase vertex offset by one
        ivtoffset = ivtoffset+1

        ! Store the vertex number at midpoint of the face
        p_ImidpointAtFace(iae) = ivtoffset

        ! Initialize vertex coordinates
        p_DcoordDest(:,ivtoffset) = 0.0_DP
        
        ! Loop over all four corner vertices and 
        ! sum up the x-, y-, and z-coordinates
        do ive = 1, 4
          ! Get global vertex number
          ivt1 = p_IverticesAtFace(ive,iae)
          
          p_DcoordDest(1,ivtoffset) = p_DcoordDest(1,ivtoffset) + &
                                      p_DcoordSource(1,ivt1)
          
          p_DcoordDest(2,ivtoffset) = p_DcoordDest(2,ivtoffset) + &
                                      p_DcoordSource(2,ivt1)
          
          p_DcoordDest(3,ivtoffset) = p_DcoordDest(3,ivtoffset) + &
                                      p_DcoordSource(3,ivt1)
        end do ! end ive
        
        ! Scale by factor 1/4
        p_DcoordDest(1,ivtoffset) = p_DcoordDest(1,ivtoffset) * 0.25_DP
        p_DcoordDest(2,ivtoffset) = p_DcoordDest(2,ivtoffset) * 0.25_DP
        p_DcoordDest(3,ivtoffset) = p_DcoordDest(3,ivtoffset) * 0.25_DP
      end do ! end iae
       
      ! Each hexahedral element gives rise to one additional
      ! vertex located in the interior/center of the element.
      do iel = 1, NELsrc
        
        ! Get number of vertices per element and do not create an
        ! vertex vertex if the element is not a hexahedral element.
        nve = tria_getNVE(p_IvertAtElementSource, iel)

        if (nve .ne. TRIA_NVEHEXA3D) cycle
        
        ! Increase vertex offset by one
        ivtoffset = ivtoffset+1
        
        ! Sum up the coordinates of the corners to get the midpoint.
        x = 0.0_DP
        y = 0.0_DP
        z = 0.0_DP
        do ive = 1, nve
          x = x + p_DcoordSource (1,p_IvertAtElementSource(ive,iel))
          y = y + p_DcoordSource (2,p_IvertAtElementSource(ive,iel))
          z = z + p_DcoordSource (3,p_IvertAtElementSource(ive,iel))
        end do
        
        ! Store the midpoint scaled by factor 1/8
        p_DcoordDest(1,ivtoffset) = x * 0.125_DP
        p_DcoordDest(2,ivtoffset) = y * 0.125_DP
        p_DcoordDest(3,ivtoffset) = z * 0.125_DP
        
      end do
      
      ! Allocate memory for IverticesAtElement and get a pointer to it.
      ! Fill the array with zero, so we will not have problems when mixing
      ! different elements featuring less than NNVE corner vertices.
      Isize = (/nnve,rdestTriangulation%NEL/)
      call storage_new('tria_refineMesh2lv3D', 'KVERT', Isize, ST_INT, &
          rdestTriangulation%h_IverticesAtElement, ST_NEWBLOCK_ZERO)
      call storage_getbase_int2D(&
          rdestTriangulation%h_IverticesAtElement, p_IvertAtElementDest)
      
      ! Initialise the element offset
      ieloffset = NELsrc

      ! Initialise the number of hexahedral elements
      ihex = 0

      ! Now the interesting part: the two-level ordering.
      do iel = 1, NELsrc
        
        ! Get number of vertices per element
        nve = tria_getNVE(p_IvertAtElementSource, iel)

        ! What type of element are we?
        select case(nve)

        case (TRIA_NVEHEXA3D)
          !========================================================= 
          ! Increase number of hexahedral elements
          ihex = ihex+1
          
          ! Determine elements numbers of subelements ...
          iel1 = ieloffset+1
          iel2 = ieloffset+2
          iel3 = ieloffset+3
          iel4 = ieloffset+4
          iel5 = ieloffset+5
          iel6 = ieloffset+6
          iel7 = ieloffset+7
          
          ! ... and store them in the refinement patch
          idx = p_IrefinementPatchIdx(iel)
          
          p_IrefinementPatch(idx)   = iel
          p_IrefinementPatch(idx+1) = iel1
          p_IrefinementPatch(idx+2) = iel2
          p_IrefinementPatch(idx+3) = iel3
          p_IrefinementPatch(idx+4) = iel4
          p_IrefinementPatch(idx+5) = iel5
          p_IrefinementPatch(idx+6) = iel6
          p_IrefinementPatch(idx+7) = iel7
          
          ! Determine the vertex number of face midpoints
          midpointFaceA = p_ImidpointAtFace(p_IfacesAtElement(1,iel))
          midpointFaceB = p_ImidpointAtFace(p_IfacesAtElement(2,iel))
          midpointFaceC = p_ImidpointAtFace(p_IfacesAtElement(3,iel))
          midpointFaceD = p_ImidpointAtFace(p_IfacesAtElement(4,iel))
          midpointFaceE = p_ImidpointAtFace(p_IfacesAtElement(5,iel))
          midpointFaceF = p_ImidpointAtFace(p_IfacesAtElement(6,iel))
          
          ! The vertex number of the midpoint in the interior
          ! of the element is rdestTriangulation%NVT - nhexs + ihex
          midPointOfIel = rdestTriangulation%NVT - nhexs + ihex

          ! This regular refinement gives rise to 7 new elements
          ieloffset=ieloffset+7
          
          ! Step 1: Initialise IverticesOnElement for element IEL
          !
          ! the 1st vertex of the old hexahedron
          p_IvertAtElementDest(1,iel) = p_IvertAtElementSource(1,iel)
          ! the index of the 1st edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex 
          p_IvertAtElementDest(3,iel) = midpointFaceA
          ! the index of the 4th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource(4,iel)+NVTsrc
          
          ! the index of the 5th edge is the index of the 5th vertex     
          p_IvertAtElementDest(5,iel) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the midpoint of face B is the index of the 6th vertex
          p_IvertAtElementDest(6,iel) = midpointFaceB
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel) = midPointOfIel
          ! the index of the midpoint of face E is the index of the 8th vertex
          p_IvertAtElementDest(8,iel) = midpointFaceE
          
          
          ! Step 2: Initialise IverticesOnElement for element IEL1
          !
          ! the 2nd vertex of the old hexahedron
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource(2,iel)
          ! the index of the 2nd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel1) = midPointFaceA
          ! the index of the 1st edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource(1,iel)+NVTsrc
          
          ! the index of the 6th edge is the index of the 5th vertex     
          p_IvertAtElementDest(5,iel1) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the index of the midpoint of face C is the index of the 6th vertex
          p_IvertAtElementDest(6,iel1) = midpointFaceC
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel1) = midPointOfIel
          ! the index of the midpoint of face B is the index of the 8th vertex
          p_IvertAtElementDest(8,iel1) = midpointFaceB
          

          ! Step 3: Initialise IverticesOnElement for element IEL2
          !
          ! the 3rd vertex of the old hexahedron
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource(3,iel)
          ! the index of the 3rd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel2) = midpointFaceA
          ! the index of the 2nd edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource(2,iel)+NVTsrc
          
          ! the index of the 7th edge is the index of the 5th vertex     
          p_IvertAtElementDest(5,iel2) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the index of the midpoint of face D is the index of the 6th vertex
          p_IvertAtElementDest(6,iel2) = midpointFaceD
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel2) = midPointOfIel
          ! the index of the midpoint of face C is the index of the 8th vertex
          p_IvertAtElementDest(8,iel2) = midpointFaceC
          

          ! Step 4: Initialise IverticesOnElement for element IEL3
          !
          ! the 4th vertex of the old hexahedron
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource(4,iel)
          ! the index of the 4th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel3) = midpointFaceA
          ! the index of the 3rd edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource(3,iel)+NVTsrc
          
          ! the index of the 8th edge is the index of the 5th vertex  
          p_IvertAtElementDest(5,iel3) = p_IedgesAtElementSource(8,iel)+NVTsrc
          ! the index of the midpoint of face E is the index of the 6th vertex
          p_IvertAtElementDest(6,iel3) = midpointFaceE
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel3) = midPointOfIel
          ! the index of the midpoint of face D is the index of the 8th vertex
          p_IvertAtElementDest(8,iel3) = midpointFaceD
          

          ! Step 5: Initialise IverticesOnElement for element IEL4
          !
          ! the 5th vertex of the old hexahedron
          p_IvertAtElementDest(1,iel4) = p_IvertAtElementSource(5,iel)
          ! the index of the 9th edge is the index of the 2nd vertex
          !p_IvertAtElementDest(2,iel4) = p_IedgesAtElementSource(9,iel)+NVTsrc
          p_IvertAtElementDest(2,iel4) = p_IedgesAtElementSource(12,iel)+NVTsrc
          ! the midpoint of face F is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel4) = midpointFaceF
          ! the index of the 12th edge is the index of the 4th vertex
          !p_IvertAtElementDest(4,iel4) = p_IedgesAtElementSource(12,iel)+NVTsrc
          p_IvertAtElementDest(4,iel4) = p_IedgesAtElementSource(9,iel)+NVTsrc
          
          ! the index of the 5th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel4) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the midpoint of face B is the index of the 6th vertex
          !p_IvertAtElementDest(6,iel4) = midpointFaceB
          p_IvertAtElementDest(6,iel4) = midpointFaceE
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel4) = midPointOfIel
          ! the index of the midpoint of face E is the index of the 8th vertex
          !p_IvertAtElementDest(8,iel4) = midpointFaceE
          p_IvertAtElementDest(8,iel4) = midpointFaceB
          

          ! Step 6: Initialise IverticesOnElement for element IEL5
          !
          ! the 6th vertex of the old hexahedron
          p_IvertAtElementDest(1,iel5) = p_IvertAtElementSource(6,iel)
          ! the index of the 10th edge is the index of the 2nd vertex
          !p_IvertAtElementDest(2,iel5) = p_IedgesAtElementSource(10,iel)+NVTsrc
          p_IvertAtElementDest(2,iel5) = p_IedgesAtElementSource(9,iel)+NVTsrc
          ! the midpoint of face F is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel5) = midpointFaceF
          ! the index of the 9th edge is the index of the 4th vertex
          !p_IvertAtElementDest(4,iel5) = p_IedgesAtElementSource(9,iel)+NVTsrc
          p_IvertAtElementDest(4,iel5) = p_IedgesAtElementSource(10,iel)+NVTsrc
          
          ! the index of the 6th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel5) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the index of the midpoint of face C is the index of the 6th vertex
          !p_IvertAtElementDest(6,iel5) = midpointFaceC
          p_IvertAtElementDest(6,iel5) = midpointFaceB
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel5) = midPointOfIel
          ! the index of the midpoint of face B is the index of the 8th vertex
          !p_IvertAtElementDest(8,iel5) = midpointFaceB
          p_IvertAtElementDest(8,iel5) = midpointFaceC
          
          
          ! Step 7: Initialise IverticesOnElement for element IEL6
          !
          ! the 7th vertex of the old hexahedron
          p_IvertAtElementDest(1,iel6) = p_IvertAtElementSource(7,iel)
          ! the index of the 11th edge is the index of the 2nd vertex
          !p_IvertAtElementDest(2,iel6) = p_IedgesAtElementSource(11,iel)+NVTsrc
          p_IvertAtElementDest(2,iel6) = p_IedgesAtElementSource(10,iel)+NVTsrc
          ! the midpoint of face F is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel6) = midpointFaceF
          ! the index of the 10th edge is the index of the 4th vertex
          !p_IvertAtElementDest(4,iel6) = p_IedgesAtElementSource(10,iel)+NVTsrc
          p_IvertAtElementDest(4,iel6) = p_IedgesAtElementSource(11,iel)+NVTsrc
          
          ! the index of the 7th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel6) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the index of the midpoint of face D is the index of the 6th vertex
          !p_IvertAtElementDest(6,iel6) = midpointFaceD
          p_IvertAtElementDest(6,iel6) = midpointFaceC
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel6) = midPointOfIel
          ! the index of the midpoint of face C is the index of the 8th vertex
          !p_IvertAtElementDest(8,iel6) = midpointFaceC
          p_IvertAtElementDest(8,iel6) = midpointFaceD
          

          ! Step 8: Initialise IverticesOnElement for element IEL7
          !
          ! the 8th vertex of the old hexahedron
          p_IvertAtElementDest(1,iel7) = p_IvertAtElementSource(8,iel)
          ! the index of the 12th edge is the index of the 2nd vertex
          !p_IvertAtElementDest(2,iel7) = p_IedgesAtElementSource(12,iel)+NVTsrc
          p_IvertAtElementDest(2,iel7) = p_IedgesAtElementSource(11,iel)+NVTsrc
          ! the midpoint of face F is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel7) = midpointFaceF
          ! the index of the 11th edge is the index of the 4th vertex
          !p_IvertAtElementDest(4,iel7) = p_IedgesAtElementSource(11,iel)+NVTsrc
          p_IvertAtElementDest(4,iel7) = p_IedgesAtElementSource(12,iel)+NVTsrc
          
          ! the index of the 8th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel7) = p_IedgesAtElementSource(8,iel)+NVTsrc
          ! the index of the midpoint of face E is the index of the 6th vertex
          !p_IvertAtElementDest(6,iel7) = midpointFaceE
          p_IvertAtElementDest(6,iel7) = midpointFaceD
          ! the index of the midpoint of the old hexahedron is the index of the 7th vertex
          p_IvertAtElementDest(7,iel7) = midPointOfIel
          ! the index of the midpoint of face D is the index of the 8th vertex
          !p_IvertAtElementDest(8,iel7) = midpointFaceD
          p_IvertAtElementDest(8,iel7) = midpointFaceE
          

        case (TRIA_NVEPRIS3D)
          !========================================================= 
          ! Determine elements numbers of subelements ...
          iel1 = ieloffset+1
          iel2 = ieloffset+2
          iel3 = ieloffset+3
          iel4 = ieloffset+4
          iel5 = ieloffset+5
          iel6 = ieloffset+6
          iel7 = ieloffset+7
          
          ! ... and store them in the refinement patch
          idx = p_IrefinementPatchIdx(iel)
          
          p_IrefinementPatch(idx)   = iel
          p_IrefinementPatch(idx+1) = iel1
          p_IrefinementPatch(idx+2) = iel2
          p_IrefinementPatch(idx+3) = iel3
          p_IrefinementPatch(idx+4) = iel4
          p_IrefinementPatch(idx+5) = iel5
          p_IrefinementPatch(idx+6) = iel6
          p_IrefinementPatch(idx+7) = iel7
          
          ! Determine the vertex number of face midpoints
          midpointFaceB = p_ImidpointAtFace(p_IfacesAtElement(2,iel))
          midpointFaceC = p_ImidpointAtFace(p_IfacesAtElement(3,iel))
          midpointFaceD = p_ImidpointAtFace(p_IfacesAtElement(4,iel))

          ! This regular refinement gives rise to 7 new elements
          ieloffset=ieloffset+7

          ! Step 1: Initialise IverticesOnElement for element IEL
          !
          ! the 1st vertex of the old prism
          p_IvertAtElementDest(1,iel) = p_IvertAtElementSource(1,iel)
          ! the index of the 1st edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 3rd edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel) = p_IedgesAtElementSource(3,iel)+NVTsrc

          ! the index of the 4th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the index of the midpoint of face B is the index of the 5th vertex
          p_IvertAtElementDest(5,iel) = midpointFaceB
          ! the index of the midpoint of face D is the index of the 6th vertex
          p_IvertAtElementDest(6,iel) = midpointFaceD


          ! Step 2: Initialise IverticesOnElement for element IEL1
          !
          ! the 2nd vertex of the old prism
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource(2,iel)
          ! the index of the 2nd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 1st edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel1) = p_IedgesAtElementSource(1,iel)+NVTsrc

          ! the index of the 5th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the midpoint of face C is the index of the 5th vertex
          p_IvertAtElementDest(5,iel1) = midpointFaceC
          ! the index of the midpoint of face B is the index of the 6th vertex
          p_IvertAtElementDest(6,iel1) = midpointFaceB

          
          ! Step 3: Initialise IverticesOnElement for element IEL2
          !
          ! the 3rd vertex of the old prism
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource(3,iel)
          ! the index of the 3rd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the index of the 2nd edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel2) = p_IedgesAtElementSource(2,iel)+NVTsrc

          ! the index of the 6th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the index of the midpoint of face D is the index of the 5th vertex
          p_IvertAtElementDest(5,iel2) = midpointFaceD
          ! the index of the midpoint of face C is the index of the 6th vertex
          p_IvertAtElementDest(6,iel2) = midpointFaceC


          ! Step 4: Initialise IverticesOnElement for element IEL3
          !
          ! the 4th vertex of the old prism
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource(4,iel)
          ! the index of the 9th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource(9,iel)+NVTsrc
          ! the index of the 7th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel3) = p_IedgesAtElementSource(7,iel)+NVTsrc

          ! the index of the 4th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the midpoint of face D is the index of the 5th vertex
          p_IvertAtElementDest(5,iel3) = midpointFaceD
          ! the midpoint of face B is the index of the 6th vertex
          p_IvertAtElementDest(6,iel3) = midpointFaceB
          

          ! Step 5: Initialise IverticesOnElement for element IEL4
          !
          ! the 5th vertex of the old prism
          p_IvertAtElementDest(1,iel4) = p_IvertAtElementSource(5,iel)
          ! the index of the 7th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel4) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the index of the 8th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel4) = p_IedgesAtElementSource(8,iel)+NVTsrc

          ! the index of the 5th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel4) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the midpoint of face B is the index of the 5th vertex
          p_IvertAtElementDest(5,iel4) = midpointFaceB
          ! the index of the midpoint of face C is the index of the 6th vertex
          p_IvertAtElementDest(6,iel4) = midpointFaceC

          
          ! Step 6: Initialise IverticesOnElement for element IEL5
          !
          ! the 6th vertex of the old prism
          p_IvertAtElementDest(1,iel5) = p_IvertAtElementSource(6,iel)
          ! the index of the 8th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel5) = p_IedgesAtElementSource(8,iel)+NVTsrc
          ! the index of the 9th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel5) = p_IedgesAtElementSource(9,iel)+NVTsrc

          ! the index of the 6th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel5) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the index of the midpoint of face C is the index of the 5th vertex
          p_IvertAtElementDest(5,iel5) = midpointFaceC
          ! the index of the midpoint of face D is the index of the 6th vertex
          p_IvertAtElementDest(6,iel5) = midpointFaceD

          
          ! Step 7: Initialise IverticesOnElement for element IEL6
          !
          ! the index of the 1st edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel6) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 2nd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel6) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 3rd edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel6) = p_IedgesAtElementSource(3,iel)+NVTsrc

          ! the index of the midpoint of face B is the index of the 4th vertex
          p_IvertAtElementDest(4,iel6) = midpointFaceB
          ! the index of the midpoint of face C is the index of the 5th vertex
          p_IvertAtElementDest(5,iel6) = midpointFaceC
          ! the index of the midpoint of face D is the index of the 6th vertex
          p_IvertAtElementDest(6,iel6) = midpointFaceD

          
          ! Step 8: Initialise IverticesOnElement for element IEL7
          !
          ! the index of the 7th edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel7) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the index of the 9th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel7) = p_IedgesAtElementSource(9,iel)+NVTsrc
          ! the index of the 8th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel7) = p_IedgesAtElementSource(8,iel)+NVTsrc

          ! the index of the midpoint of face B is the index of the 4th vertex
          p_IvertAtElementDest(4,iel7) = midpointFaceB
          ! the index of the midpoint of face D is the index of the 5th vertex
          p_IvertAtElementDest(5,iel7) = midpointFaceD
          ! the index of the midpoint of face C is the index of the 6th vertex
          p_IvertAtElementDest(6,iel7) = midpointFaceC


        case (TRIA_NVEPYR3D)
          !========================================================= 
          ! Determine elements numbers of subelements ...
          iel1 = ieloffset+1
          iel2 = ieloffset+2
          iel3 = ieloffset+3
          iel4 = ieloffset+4
          iel5 = ieloffset+5
          iel6 = ieloffset+6
          iel7 = ieloffset+7
          iel8 = ieloffset+8
          iel9 = ieloffset+9

          ! ... and store them in the refinement patch
          idx = p_IrefinementPatchIdx(iel)
          
          p_IrefinementPatch(idx)   = iel
          p_IrefinementPatch(idx+1) = iel1
          p_IrefinementPatch(idx+2) = iel2
          p_IrefinementPatch(idx+3) = iel3
          p_IrefinementPatch(idx+4) = iel4
          p_IrefinementPatch(idx+5) = iel5
          p_IrefinementPatch(idx+6) = iel6
          p_IrefinementPatch(idx+7) = iel7
          p_IrefinementPatch(idx+8) = iel8
          p_IrefinementPatch(idx+9) = iel9

          ! The vertex number of midpoints located on the 
          ! element faces is the old face number +NVT+NMT
          midpointFaceA = p_ImidpointAtFace(p_IfacesAtElement(1,iel))

          ! This regular refinement gives rise to 9 new elements
          ieloffset=ieloffset+9

          ! Step 1: Initialise IverticesOnElement for element IEL
          !
          ! the 1st vertex of the old pyramid
          p_IvertAtElementDest(1,iel) = p_IvertAtElementSource(1,iel)
          ! the index of the 1st edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex 
          p_IvertAtElementDest(3,iel) = midpointFaceA
          ! the index of the 4th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel) = p_IedgesAtElementSource(5,iel)+NVTsrc


          ! Step 2: Initialise IverticesOnElement for element IEL1
          !
          ! the 2nd vertex of the old pyramid
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource(2,iel)
          ! the index of the 2nd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel1) = midPointFaceA
          ! the index of the 1st edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel1) = p_IedgesAtElementSource(6,iel)+NVTsrc


          ! Step 3: Initialise IverticesOnElement for element IEL2
          !
          ! the 3rd vertex of the old pyramid
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource(3,iel)
          ! the index of the 3rd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel2) = midpointFaceA
          ! the index of the 2nd edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 7th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel2) = p_IedgesAtElementSource(7,iel)+NVTsrc


          ! Step 4: Initialise IverticesOnElement for element IEL3
          !
          ! the 4th vertex of the old pyramid
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource(4,iel)
          ! the index of the 4th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the midpoint of face A is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel3) = midpointFaceA
          ! the index of the 3rd edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the index of the 8th edge is the index of the 5th vertex
          p_IvertAtElementDest(5,iel3) = p_IedgesAtElementSource(8,iel)+NVTsrc


          ! Step 5: Initialise IverticesOnElement for element IEL4
          !
          ! the index of the 5th edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel4) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel4) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the index of the 7th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel4) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the index of the 8th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel4) = p_IedgesAtElementSource(8,iel)+NVTsrc
          ! the 5th vertex of the old pyramid
          p_IvertAtElementDest(5,iel4) = p_IvertAtElementSource(5,iel)


          ! Step 6: Initialise IverticesOnElement for element IEL5
          !
          ! the index of the 8th edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel5) = p_IedgesAtElementSource(8,iel)+NVTsrc
          ! the index of the 7th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel5) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel5) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel5) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the midpoint of face A is the index of the 5th vertex
          p_IvertAtElementDest(5,iel5) = midpointFaceA

          
          ! Step 7: Initialise IverticesOnElement for element IEL6
          !
          ! the index of the 1st edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel6) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel6) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel6) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the midpoint of face A is the index of the 4th vertex
          p_IvertAtElementDest(4,iel6) = midpointFaceA


          ! Step 8: Initialise IverticesOnElement for element IEL7
          !
          ! the index of the 2nd edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel7) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel7) = p_IedgesAtElementSource(6,iel)+NVTsrc
          ! the index of the 7th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel7) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the midpoint of face A is the index of the 4th vertex
          p_IvertAtElementDest(4,iel7) = midpointFaceA


          ! Step 9: Initialise IverticesOnElement for element IEL8
          !
          ! the index of the 3rd edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel8) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the index of the 7th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel8) = p_IedgesAtElementSource(7,iel)+NVTsrc
          ! the index of the 8th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel8) = p_IedgesAtElementSource(8,iel)+NVTsrc
          ! the midpoint of face A is the index of the 4th vertex
          p_IvertAtElementDest(4,iel8) = midpointFaceA
          

          ! Step 10: Initialise IverticesOnElement for element IEL9
          !
          ! the index of the 4th edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel9) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the index of the 8th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel9) = p_IedgesAtElementSource(8,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel9) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the midpoint of face A is the index of the 4th vertex
          p_IvertAtElementDest(4,iel9) = midpointFaceA
          
          
        case (TRIA_NVETET3D)
          !========================================================= 
          ! Determine elements numbers of subelements ...
          iel1 = ieloffset+1
          iel2 = ieloffset+2
          iel3 = ieloffset+3
          iel4 = ieloffset+4
          iel5 = ieloffset+5
          iel6 = ieloffset+6
          iel7 = ieloffset+7
          
          ! ... and store them in the refinement patch
          idx = p_IrefinementPatchIdx(iel)
          
          p_IrefinementPatch(idx)   = iel
          p_IrefinementPatch(idx+1) = iel1
          p_IrefinementPatch(idx+2) = iel2
          p_IrefinementPatch(idx+3) = iel3
          p_IrefinementPatch(idx+4) = iel4
          p_IrefinementPatch(idx+5) = iel5
          p_IrefinementPatch(idx+6) = iel6
          p_IrefinementPatch(idx+7) = iel7

          ! This regular refinement gives rise to 7 new elements
          ieloffset=ieloffset+7

          ! Step 1: Initialise IverticesOnElement for element IEL
          !
          ! the 1st vertex of the old tetrahedron
          p_IvertAtElementDest(1,iel) = p_IvertAtElementSource(1,iel)
          ! the index of the 1st edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 1st edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the index of the 1st edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource(4,iel)+NVTsrc


          ! Step 2: Initialise IverticesOnElement for element IEL1
          !
          ! the 2nd vertex of the old tetrahedron
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource(2,iel)
          ! the index of the 2nd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 1st edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel1) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource(5,iel)+NVTsrc


          ! Step 3: Initialise IverticesOnElement for element IEL2
          !
          ! the 3rd vertex of the old tetrahedron
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource(3,iel)
          ! the index of the 3rd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the index of the 2nd edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel2) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource(6,iel)+NVTsrc


          ! Step 4: Initialise IverticesOnElement for element IEL3
          !
          ! the 4th vertex of the old tetrahedron
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource(4,iel)
          ! the index of the 4th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel3) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource(6,iel)+NVTsrc
          

          ! Step 5: Initialise IverticesOnElement for element IEL4
          !
          ! the index of the 1st edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel4) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 2nd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel4) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 3rd edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel4) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel4) = p_IedgesAtElementSource(6,iel)+NVTsrc


          ! Step 6: Initialise IverticesOnElement for element IEL5
          !
          ! the index of the 1st edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel5) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 3rd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel5) = p_IedgesAtElementSource(3,iel)+NVTsrc
          ! the index of the 4th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel5) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel5) = p_IedgesAtElementSource(6,iel)+NVTsrc

          
          ! Step 7: Initialise IverticesOnElement for element IEL6
          !
          ! the index of the 1st edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel6) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 2nd edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel6) = p_IedgesAtElementSource(2,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel6) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel6) = p_IedgesAtElementSource(6,iel)+NVTsrc


          ! Step 8: Initialise IverticesOnElement for element IEL7
          !
          ! the index of the 1st edge is the index of the 1st vertex
          p_IvertAtElementDest(1,iel7) = p_IedgesAtElementSource(1,iel)+NVTsrc
          ! the index of the 4th edge is the index of the 2nd vertex
          p_IvertAtElementDest(2,iel7) = p_IedgesAtElementSource(4,iel)+NVTsrc
          ! the index of the 5th edge is the index of the 3rd vertex
          p_IvertAtElementDest(3,iel7) = p_IedgesAtElementSource(5,iel)+NVTsrc
          ! the index of the 6th edge is the index of the 4th vertex
          p_IvertAtElementDest(4,iel7) = p_IedgesAtElementSource(6,iel)+NVTsrc


        case DEFAULT
          call output_line('Unsupported type of element shape',&
                           OU_CLASS_ERROR,OU_MODE_STD,'tria_refineMesh2lv3D')
          call sys_halt()
        end select
      end do ! iel
      
      ! The last step of setting up the raw mesh on the finer level:
      call storage_new ('tria_refineMesh2lv3D', 'KNPR',&
          rdestTriangulation%NVT, ST_INT,&
          rdestTriangulation%h_InodalProperty, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(&
          rdestTriangulation%h_InodalProperty, p_InodalPropertyDest)
      call storage_getbase_int(&
          rsourceTriangulation%h_InodalProperty, p_InodalPropertySource)

p_InodalPropertyDest = -4711

      ! Copy the nodal information of vertices and edges from the
      ! coarse mesh to the nodal information about the new vertices
      ! on the fine mesh!
      call lalg_copyVectorInt(p_InodalPropertySource,&
                              p_InodalPropertyDest, NVTsrc+NMTsrc)

      ! The treatment of faces is slightly more difficult. For pure
      ! hexahedral meshes the nodal information of the faces can be
      ! directly converted into the nodal information about the new
      ! vertices on the fine mesh! Note that for mixed meshes, only
      ! quadrilateral faces create new vertices in the fine mesh.
      do iae = 1, rsourceTriangulation%NAT

        ! Get global vertex number of face midpoint
        ivt = p_ImidpointAtFace(iae)

        ! Check if vertex is not empty
        if (ivt > 0) p_InodalPropertyDest(ivt) =&
            p_InodalPropertySource(NVTsrc+NMTsrc+iae)
      end do
      
      ! There are still NHEXS vertices created in the interior of hexahedra
      do ivt = rdestTriangulation%NVT-nhexs+1, rdestTriangulation%NVT
        p_InodalPropertyDest(ivt) = 0
      end do

      rdestTriangulation%NNEE = rsourceTriangulation%NNEE
      rdestTriangulation%NNAE = rsourceTriangulation%NNAE
      rdestTriangulation%NNVA = rsourceTriangulation%NNVA
      
      ! In principle, this subroutine does not refine the boundary. 
      ! However, the number of quadrilaterals at the boundary has 
      ! already been compute so that is it expedient to set NVBD here.
      rdestTriangulation%NVBD =  rsourceTriangulation%NVBD + &
                                 rsourceTriangulation%NMBD + nquads

      ! Release auxiliary memory
      call storage_free(h_ImidpointAtFace)

    end subroutine tria_refineMesh2lv3D
    
    ! ---------------------------------------------------------------

    subroutine tria_refineBdry2lv3D(rsourceTriangulation, rdestTriangulation, rboundary)
      
      ! This routine refines the boundary definition of rsourceTriangulation
      ! according to the 2-level ordering algorithm to generate a new 
      ! IverticesAtBoundary. 
      
      ! The source triangulation to be refined
      type(t_triangulation), intent(in) :: rsourceTriangulation
      
      ! Destination triangulation structure that receives the refined mesh
      type(t_triangulation), intent(inout) :: rdestTriangulation
      
      ! OPTIONAL: Defintion of analytic boundary.
      ! If specified, the coordinates of the new boundary vertices are
      ! recomputed according to the analytic boundary.
      type(t_boundary), intent(in), optional :: rboundary


      ! local variables
      
      integer, dimension(:), pointer :: p_IvertAtBoundartySource
      integer, dimension(:), pointer :: p_IvertAtBoundartyDest
      integer, dimension(:), pointer :: p_IedgesAtBoundartySource
      integer, dimension(:), pointer :: p_IedgesAtBoundartyDest
      integer, dimension(:), pointer :: p_InodalPropertyDest
      integer, dimension(:), pointer :: p_IboundaryCpIdxSource
      integer, dimension(:), pointer :: p_IboundaryCpIdxDest
      integer, dimension(:), pointer :: p_IfacesAtBoundary
      integer, dimension(:), pointer :: p_IedgesAtBoundary
      integer, dimension(:,:), pointer :: p_IverticesAtEdge
      integer, dimension(:,:), pointer :: p_IverticesAtFace
      
      integer :: ivt,ivbd,ibct,isize
      
      ! Set pointers
      call storage_getbase_int(&
          rsourceTriangulation%h_IverticesAtBoundary, p_IvertAtBoundartySource)
      call storage_getbase_int(&
          rsourceTriangulation%h_IedgesAtBoundary, p_IedgesAtBoundartySource)
      call storage_getbase_int(&
          rsourceTriangulation%h_IboundaryCpIdx, p_IboundaryCpIdxSource)
      call storage_getbase_int(&
          rsourceTriangulation%h_IfacesAtBoundary, p_IfacesAtBoundary)
      call storage_getbase_int(&
          rsourceTriangulation%h_IedgesAtBoundary, p_IedgesAtBoundary)
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IverticesAtEdge, p_IverticesAtEdge)
      call storage_getbase_int2d(&
          rsourceTriangulation%h_IverticesAtFace, p_IverticesAtFace)
      
      ! Set number of boundary components
      rdestTriangulation%NBCT      = rsourceTriangulation%NBCT 
      rdestTriangulation%NblindBCT = rsourceTriangulation%NblindBCT 
          
      ! Create new arrays in the fine grid for the vertices and indices.
      call storage_new ('tria_refineBdry2lv3D', 'KVBD', rdestTriangulation%NVBD, &
          ST_INT, rdestTriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
      
      call storage_getbase_int(&
          rdestTriangulation%h_IverticesAtBoundary, p_IvertAtBoundartyDest) 
      
      call storage_getbase_int(&
          rdestTriangulation%h_InodalProperty, p_InodalPropertyDest)
      
      if(rdestTriangulation%h_IboundaryCpIdx .eq. ST_NOHANDLE) then
        call storage_new('tria_refineBdry2lv3D', 'IboundaryCpIdx', &
            rdestTriangulation%NBCT+rdestTriangulation%NblindBCT+1, ST_INT, &
            rdestTriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT)
      end if
      
      call storage_getbase_int(&
          rdestTriangulation%h_IboundaryCpIdx, p_IboundaryCpIdxDest)
      
      ! Regenerate the array of vertices at the boundary
      ivbd = 1
      do ivt = 1, rdestTriangulation%NVT
        if(p_InodalPropertyDest(ivt) > 0) then
          p_IvertAtBoundartyDest(ivbd) = ivt
          ivbd = ivbd + 1
        end if
      end do
      
      ! Initialise p_IboundaryCpIdx
      call lalg_clearVector(p_IboundaryCpIdxDest)
      p_IboundaryCpIdxDest(1) = 1
      
      ! Assign the indices of the boundary vertices:
      ! Step 1: save the number of vertices in each boundary 
      !         component in p_IboundaryCpIdx(2:NBCT+1)
      do ivt = 1, rdestTriangulation%NVT
        if(p_InodalPropertyDest(ivt) .ne. 0) then
          ibct = p_InodalPropertyDest(ivt)
          p_IboundaryCpIdxDest(ibct+1) = p_IboundaryCpIdxDest(ibct+1) + 1
        end if
      end do
      
      ! Step 2: create the actual index array
      do ibct = 2, rdestTriangulation%NBCT+rdestTriangulation%NblindBCT+1
        p_IboundaryCpIdxDest(ibct) = p_IboundaryCpIdxDest(ibct)+ &
                                     p_IboundaryCpIdxDest(ibct-1)
      end do
      
    end subroutine tria_refineBdry2lv3D

  end subroutine tria_refine2LevelOrdering

  ! ***************************************************************************

!<subroutine>

  subroutine tria_compress2LevelOrdHierarchy (rtriangulationFine, rtriangulationCoarse)
  
!<description>
  ! This routine can be used to save memory after a refinement with the
  ! 2-level-ordering algorithm was applied. It is typically applied 'backwards'
  ! from the maximum level to the minimum one, releases redundant data on the
  ! coarse mesh and replaces it by the corresponding data on the fine mesh.
  !
  ! In detail, the routine does the following:
  ! 1.) Release the DvertexCoords array in rtriangulationCoarse.
  ! 2.) Shares the DvertexCoords array in rtriangulationCoarse with that one
  !     in rtriangulationFine
  ! This technique is possible since the vertex coordinates of the 
  ! (rtriangulationCoarse%) NVT vertices are the same in the coarse as well 
  ! as in the fine mesh by definition of the 2-level ordering!
  !
  ! The effect is the following: 1.) It saves memory (as the
  ! vertex coordinates exist only once) and 2.) every change of the coordinates
  ! in the fine mesh b the application will also affect the coordinates 
  ! on the coarse mesh and vice versa!
  !
  ! Example: 
  !
  !   do I=NLMAX-1,NLMIN,-1             
  !     call tria_compress2LevelOrdHierarchy (rtria(i+1),rtria(i))
  !   end do
  !
  ! Afterwards, the finest mesh contains all coordinates and the coarse grid
  ! coordinates are part of the fine grid coordinates.
  !
  ! WARNING: Use this routine with care. It does not check whether
  ! rtriangulationCoarse,rtriangulationFine are 'compatible' to each other.
!</description>

!<input>
  ! Fine grid triangulation.
  type(t_triangulation), intent(in) :: rtriangulationFine
!</input>

!<inputoutput>
  ! Coarse grid triangulation where redundant data should be removed from.
  type(t_triangulation), intent(inout) :: rtriangulationCoarse
!</inputoutput>
  
!</subroutine>
  
    ! Release the vertex coordinates array in the coarse triangulation -- as long
    ! as it is not a copy of someone else...
    if (iand(rtriangulationCoarse%iduplicationFlag,TR_SHARE_DVERTEXCOORDS) .eq. 0) then
      call storage_free (rtriangulationCoarse%h_DvertexCoords)
    end if 
    
    ! Share the same handle.
    rtriangulationCoarse%h_DvertexCoords = rtriangulationFine%h_DvertexCoords
    
    ! Mark the DvertexCoords array in the source mesh as `being a copy
    ! of someone else`, so the array is not released when the coarse
    ! mesh is released!
    rtriangulationCoarse%iduplicationFlag = &
      ior(rtriangulationCoarse%iduplicationFlag,TR_SHARE_DVERTEXCOORDS)

  end subroutine tria_compress2LevelOrdHierarchy

  !************************************************************************

!<subroutine>

  subroutine tria_quickRefine2LevelOrdering(nfine, rtriangulation, rboundary, cflags)

!<description>
  ! This routine refines the given mesh rsourceTriangulation according to
  ! the 2-level ordering algorithm nfine times. The refined mesh is saved in 
  ! rdestTriangulation.
  !
  ! rtriangulation can be either a 'raw' mesh (as read from a .TRI
  ! file) or a 'standard' mesh. The mesh is overwritten by the refined
  ! one. The resulting mesh will be a 'raw' mesh, i.e. the caller must 
  ! add further information to it with routines like 
  ! tria_initStandardMeshFromRaw!
!</description>

!<input>
  ! Number of refinements that should be applied to rtriangulation. > 0.
  ! For a value .le. 0, nothing will happen.
  integer, intent(in) :: nfine

  ! OPTIONAL: Boundary structure that defines the parametrisation of the boundary.
  ! If specified, the coordinates of the new boundary vertices are
  ! recomputed according to the analytic boundary.
  type(t_boundary), intent(in), optional :: rboundary

  ! OPTIONAL: Bitfield of TRIA_R2LV_xxxx constants that allow to specify
  ! options for the refinement. If not specified, TRIA_R2LV_STANDARD
  ! is used as default.
  integer(I32), intent(in), optional :: cflags
!</input>

!<inputoutput>
  ! The triangulation to be refined; can be a 'raw' or a 'standard' mesh.
  ! Is overwritten by the refined mesh.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>
 
    integer :: ifine

    select case(rtriangulation%ndim)
    case (NDIM1D)
      ! Refine nfine times:
      do ifine = 1, nfine

        ! Create missing arrays in the source mesh.
        ! Create only those arrays we need to make the refinement process
        ! as fast as possible.
        if (rtriangulation%h_IelementsAtVertex .eq. ST_NOHANDLE) &
            call tria_genElementsAtVertex1D2D (rtriangulation)
        
        if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) &
            call tria_genNeighboursAtElement2D (rtriangulation)
        
        call tria_sortBoundaryVertices1D2D (rtriangulation)
        
        if (rtriangulation%h_IelementsAtBoundary .eq. ST_NOHANDLE) &
            call tria_genElementsAtBoundary1D2D (rtriangulation)
        
        ! Refine the mesh, replace the source mesh.
        call tria_refine2LevelOrdering(rtriangulation,rboundary=rboundary,&
                                       cflags=cflags)
           
      end do
      
    case (NDIM2D)
      ! Refine nfine times:
      do ifine = 1, nfine
      
        ! Create missing arrays in the source mesh.
        ! Create only those arrays we need to make the refinement process
        ! as fast as possible.
        if (rtriangulation%h_IelementsAtVertex .eq. ST_NOHANDLE) &
            call tria_genElementsAtVertex1D2D (rtriangulation)
        
        if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) &
            call tria_genNeighboursAtElement2D (rtriangulation)
        
        if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) &
            call tria_genEdgesAtElement2D (rtriangulation)
        
        if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) &
            call tria_genVerticesAtEdge2D (rtriangulation)
        
        ! Reallocate the nodal property array if necessary.
        ! IMPORTANT NOTE:
        !   Reallocation introduces a very serious trick!
        !   The reallocated memory contains the old nodal property array
        !   in the beginning and is filled with zero. Thus all
        !   new entries (for the new edge/element midpoints) are by
        !   default set to be inner nodes!
        !   The successive call to tria_genEdgeNodalProperty2D will
        !   generate the edge nodal property information for the
        !   coarse mesh which then defines the vertex nodal property
        !   for the fine mesh. The nodal property for the vertices
        !   stemming from element midpoints are not touched, so they
        !   stay at 0, identifying an inner vertex!
        if (reallocRefinedInodalProperty2D (rtriangulation)) &
            call tria_genEdgeNodalProperty2D (rtriangulation)
        
        call tria_sortBoundaryVertices1D2D (rtriangulation)
        
        if (rtriangulation%h_IelementsAtBoundary .eq. ST_NOHANDLE) &
            call tria_genElementsAtBoundary1D2D (rtriangulation)
        
        if (rtriangulation%h_IedgesAtBoundary .eq. ST_NOHANDLE) &
            call tria_genEdgesAtBoundary2D (rtriangulation)
        
        if (present(rboundary) .and. &
            (rtriangulation%h_DedgeParameterValue .eq. ST_NOHANDLE)) then
          call tria_genEdgeParameterValue2D (rtriangulation,rboundary)
        end if
       
        ! Refine the mesh, replace the source mesh.
        call tria_refine2LevelOrdering(rtriangulation,rboundary=rboundary,&
                                       cflags=cflags)

      end do
      
    case (NDIM3D)
      do ifine = 1, nfine     
        
        if (rtriangulation%h_IelementsAtVertex .eq. ST_NOHANDLE) &
            call tria_genElementsAtVertex3D (rtriangulation)
        
        if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) &
            call tria_genNeighboursAtElement3D (rtriangulation)
        
        if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) &
            call tria_genEdgesAtElement3D (rtriangulation)
        
        if (rtriangulation%h_IelementsAtEdge3D .eq. ST_NOHANDLE) &
            call tria_genElementsAtEdge3D (rtriangulation)
        
        if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) &
            call tria_genVerticesAtEdge3D (rtriangulation)
        
        if (rtriangulation%h_IfacesAtElement .eq. ST_NOHANDLE) &
            call tria_genFacesAtElement3D (rtriangulation)
        
        if (rtriangulation%h_IverticesAtFace .eq. ST_NOHANDLE) &
            call tria_genVerticesAtFace3D (rtriangulation)          
        
        if (rtriangulation%h_IelementsAtFace .eq. ST_NOHANDLE) &
            call tria_genElementsAtFace3D (rtriangulation)             
        
        if (rtriangulation%h_IfacesAtBoundary .eq. ST_NOHANDLE) &
            call tria_genFacesAtBoundary3D (rtriangulation)             
        
        if (rtriangulation%h_IedgesAtFace .eq. ST_NOHANDLE) &
            call tria_genEdgesAtFace3D (rtriangulation)             
        
        if (rtriangulation%h_IfacesAtEdge .eq. ST_NOHANDLE) &
            call tria_genFacesAtEdge3D (rtriangulation)             
        
        if (rtriangulation%h_IedgesAtBoundary .eq. ST_NOHANDLE) &
            call tria_genEdgesAtBoundary3D (rtriangulation)
        
        call tria_genEdgeNodalProperty3D (rtriangulation)

        call tria_genFaceNodalProperty3D (rtriangulation)
        
        ! Refine the mesh, replace the source mesh.
        call tria_refine2LevelOrdering(rtriangulation,rboundary=rboundary,&
                                       cflags=cflags)
      end do
    end select
    
  contains
  
    logical function reallocRefinedInodalProperty2D (rtriangulation)
    
    ! Reallocates the InodalProperty-array in rtriangulation such that
    ! is provides enough space for the refined mesh.
    
    type(t_triangulation), intent(inout) :: rtriangulation
    
    ! Return value: whether the memory was reallocated or not.
    
      ! local variables
      integer :: nnodes,isize
      integer, dimension(:), pointer :: p_InodalProperty
      
      ! Calculate the number of nodes in the mesh.
      ! This is: #vertices 
      !         +#edges (as every edge generates a new vertex)
      !         +#quads (as every quad generates a midpoint)
      nnodes = rtriangulation%NVT + rtriangulation%NMT + &
          rtriangulation%InelOfType(TRIA_NVEQUAD2D)
          
      ! Reallocate the memory if necessary.
      ! Copy the old content as we must not destroy the old nodal 
      ! property tags of the vertices.
      ! New elements are filled with zero = specify inner vertices.
      call storage_getsize (rtriangulation%h_InodalProperty, isize)
      if (isize .lt. nnodes) then
        call storage_realloc ('tria_genEdgeNodalProperty2D', &
            nnodes, rtriangulation%h_InodalProperty, &
            ST_NEWBLOCK_NOINIT, .true.)
        
        if (rtriangulation%InelOfType(TRIA_NVEQUAD2D) .ne. 0) then
          ! Fill the last InelOfType(TRIA_NVEQUAD2D) entries of
          ! the array by 0. These elements will generate the
          ! nodal property for the element midpoints.
          ! Edge midpoints are recomputed anyway, so we do not
          ! have to initialise that part of the array!
          call storage_getbase_int (&
              rtriangulation%h_InodalProperty,p_InodalProperty)
          call lalg_clearVectorInt (&
              p_InodalProperty(nnodes-rtriangulation%InelOfType(TRIA_NVEQUAD2D)+1:))
        end if
        
        reallocRefinedInodalProperty2D = .true.
      else
        reallocRefinedInodalProperty2D = .false.
      end if
      
    end function reallocRefinedInodalProperty2D
    
  end subroutine tria_quickRefine2LevelOrdering

  ! ***************************************************************************

!<subroutine>

  subroutine tria_infoStatistics (rtriangulation,bheadline,ilevel)
  
!<description>
  ! Prints out statistical information of the given triangulation to the
  ! terminal. The output is formatted as a table with an optional headline
  ! and an optional level identifier.
!</description>

!<input>
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! OPTIONAL: Print out a headline above the statistical data.
  ! =FALSE: do not print = standard.
  ! =TRUE: print a headline.
  logical, intent(in), optional :: bheadline
  
  ! OPTIONAL: Level identifier.
  ! If specified, an additional column 'Level' is added to the front of
  ! the statistics table. ilevel is printed to this column.
  integer, intent(in), optional :: ilevel
!</input>

!</subroutine>

    select case (rtriangulation%NDIM)
    case (NDIM1D)
      if (present(bheadline)) then
        if (bheadline) then
          ! Print a headline
          if (present(ilevel)) then
            call output_line(&
              'Lv. dim.        NVT        NMT        NEL    NBCT' &
            //'  NblindBCT    NVBD     #lines')
            call output_line(&
              '--------------------------------------------' &
            //'------------------------------------')
          else
            call output_line(&
              'dim.        NVT        NMT        NEL    NBCT' &
            //'  NblindBCT    NVBD     #lines')
            call output_line(&
              '--------------------------------------------' &
            //'-------------------------------')
          end if
        end if
      end if

      ! Print out the statistics
      if (present(ilevel)) call output_line (trim(sys_si(ilevel,3))//' ',&
          bnolinebreak=.true.)
      
      call output_line (&
          trim(sys_si(rtriangulation%NDIM,4)) &
        //trim(sys_si(rtriangulation%NVT,11)) &
        //trim(sys_si(rtriangulation%NMT,11)) &
        //trim(sys_si(rtriangulation%NEL,11)) &
        //trim(sys_si(rtriangulation%NBCT,8)) &
        //trim(sys_si(rtriangulation%NblindBCT,11)) &
        //trim(sys_si(rtriangulation%NVBD,8)) &
        //trim(sys_si(rtriangulation%InelOfType(TRIA_NVELINE1D),11)),cdateTimeLogPolicy=OU_DTP_NONE )

    case (NDIM2D)
      if (present(bheadline)) then
        if (bheadline) then
          ! Print a headline
          if (present(ilevel)) then
            call output_line(&
              'Lv. dim.        NVT        NMT        NEL    NBCT' &
            //'  NblindBCT    NVBD     #trias     #quads')
            call output_line(&
              '--------------------------------------------' &
            //'----------------------------------------------')
          else
            call output_line(&
              'dim.        NVT        NMT        NEL    NBCT' &
            //'  NblindBCT    NVBD     #trias     #quads')
            call output_line(&
              '--------------------------------------------' &
            //'------------------------------------------')
          end if
        end if
      end if

      ! Print out the statistics
      if (present(ilevel)) call output_line (trim(sys_si(ilevel,3))//' ',&
          bnolinebreak=.true.)
      
      call output_line (&
          trim(sys_si(rtriangulation%NDIM,4)) &
        //trim(sys_si(rtriangulation%NVT,11)) &
        //trim(sys_si(rtriangulation%NMT,11)) &
        //trim(sys_si(rtriangulation%NEL,11)) &
        //trim(sys_si(rtriangulation%NBCT,8)) &
        //trim(sys_si(rtriangulation%NblindBCT,11)) &
        //trim(sys_si(rtriangulation%NVBD,8)) &
        //trim(sys_si(rtriangulation%InelOfType(TRIA_NVETRI2D),11)) &
        //trim(sys_si(rtriangulation%InelOfType(TRIA_NVEQUAD2D),11)),cdateTimeLogPolicy=OU_DTP_NONE )
    end select

  end subroutine tria_infoStatistics

  !************************************************************************

!<subroutine>

  subroutine tria_exportTriFile(rtriangulation, sfilename, ctriFormat)

!<description>
    ! This routine exports a triangulation into a .TRI file.
!</description>

!<input>
    ! Triangulation structure, to be exported
    type(t_triangulation), intent(inout) :: rtriangulation

    ! The name of the .tri file to write.
    character(LEN=*), intent(in) :: sfilename
    
    ! OPTIONAL: Format tag of the TRI file to export.
    ! TRI_FMT_STANDARD: Standard TRI file format, compatible to FEAT1.
    !    Vertex coordinates of boundary vertices are in 2D replaced
    !    by parameter values.
    ! TRI_FMT_NOPARAMETRISATION: Standard TRI file format, but the 
    !    vertex coordinates are exported 'as they are', not as parameter
    !    values.
    ! If not specified, TRI_FMT_STANDARD is assumed.
    integer(I32), intent(in), optional :: ctriFormat
!</input>
!</subroutine>

    ! Depending on the dimension of the triangulation, call the corresponding
    ! export routine!
    select case(rtriangulation%ndim)
    case (NDIM1D)
      call output_line ('1D TRI file export not implemented!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_exportTriFile')
      call sys_halt()
      
    case (NDIM2D)
      call tria_exportTriFile2D(rtriangulation, sfilename, ctriFormat)
      
    case (NDIM3D)
      call output_line ('3D TRI file export not implemented!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_exportTriFile')
      call sys_halt()
      
    case DEFAULT
      call output_line ('Triangulation structure not properly initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_exportTriFile')
      call sys_halt()
    end select

  contains

    ! ---------------------------------------------------------------

    subroutine tria_exportTriFile2D(rtriangulation, sfilename, ctriFormat)

      ! Auxiliary routine. This routine exports a 2D triangulation into a .TRI file.
      
      ! Triangulation structure, to be exported
      type(t_triangulation), intent(inout) :: rtriangulation
      
      ! The name of the .tri file to write.
      character(LEN=*), intent(in) :: sfilename
      
      ! OPTIONAL: Format tag of the TRI file to export.
      ! TRI_FMT_STANDARD: Standard TRI file format, compatible to FEAT1.
      !    Vertex coordinates of boundary vertices are in 2D replaced
      !    by parameter values.
      ! TRI_FMT_NOPARAMETRISATION: Standard TRI file format, but the 
      !    vertex coordinates are exported 'as they are', not as parameter
      !    values.
      ! If not specified, TRI_FMT_STANDARD is assumed.
      integer(I32), intent(in), optional :: ctriFormat

      ! Local variables
      real(DP), dimension(:,:), pointer     :: p_Ddata2D
      real(DP), dimension(:), pointer     :: p_Ddata
      integer, dimension(:,:), pointer :: p_Idata2D
      integer, dimension(:), pointer   :: p_Idata
      integer, dimension(:), pointer   :: p_IverticesAtBoundary
      integer, dimension(:), pointer   :: p_InodalProperty
      integer, dimension(:), pointer   :: p_IboundaryCpIdx
      character(SYS_STRLEN) :: ckmmstr
      integer :: ivt, iel, ivbd
      integer      :: idim, ive, ibct
      integer      :: iunit
      logical      :: bnoParameters
      
      bnoParameters = .false.
      if (present(ctriFormat)) &
          bnoParameters = ctriFormat .eq. TRI_FMT_NOPARAMETRISATION
      
      ! Open the file
      call io_openFileForWriting(sfilename, iunit, SYS_REPLACE)
      
      ! Comment: Header
      write (iunit,*) 'Coarse mesh exported by FeatFlow2 exporter'
      if (.not. bnoParameters) then
        write (iunit,*) 'Parametrisation PARXC, PARYC, TMAXC'
        !WRITE (iunit,*) 'Parametrisierung PARXC, PARYC, TMAXC'
      else
        write (iunit,*) 'No parametrisation'
      end if
      
      ! Write NEL,NVT,NMT,NVE,NBCT to the file
      write (iunit,*) rtriangulation%NEL,rtriangulation%NVT,rtriangulation%NMT,&
          rtriangulation%NNVE,rtriangulation%NBCT,'NEL NVT NMT NVE NBCT'
      
      ! Write: 'DCORVG'
      write (iunit,*) 'DCORVG'
      
      ! Get the pointers to the coordinate array
      call storage_getbase_double2D(&
          rtriangulation%h_DvertexCoords,p_Ddata2D)
      
      ! Write the data to the file
      !
      if (.not. bnoParameters) then
        ! Standard FEAT1 format. For boundary vertices, the parameter value is
        ! exported instead of the vertex coordinate!
        call storage_getbase_int(&
            rtriangulation%h_InodalProperty,p_InodalProperty)
        call storage_getbase_double(&
            rtriangulation%h_DvertexParameterValue,p_Ddata)
        
        do ivt = 1, rtriangulation%NVT
          if (p_InodalProperty(ivt) .gt. 0) then
            ! Get the index of the vertex in the IverticesAtBoundary array.
            call tria_searchBoundaryNode(ivt,rtriangulation,ivbd)
            
            ! Write the parameter value. 2nd entry is =0.
            write (iunit,*) p_Ddata(ivbd),0.0_DP
          else
            write (iunit,*) (p_Ddata2D(idim,ivt),idim=1,NDIM2D)
          end if
        end do
      else
        ! Parameterless format.
        ! The coordinates of the vertices are exportes 'as they are' independent
        ! of whether they are on the boundary or not.
        do ivt = 1, rtriangulation%NVT
          write (iunit,*) (p_Ddata2D(idim,ivt),idim=1,NDIM2D)
        end do
      end if
      
      ! Write: 'KVERT'
      write (iunit,*) 'KVERT'
      
      ! Get the pointer to the IverticesAtElement array and read the array
      call storage_getbase_int2D(&
          rtriangulation%h_IverticesAtElement,p_Idata2D)
      
      ! Write the data to the file
      do iel = 1, rtriangulation%NEL
        write (iunit,*) (p_Idata2D(ive,iel),ive=1,size(p_Idata2D,1))
      end do
      
      ! Write: 'KNPR'
      write (iunit,*) 'KNPR'
      
      ! Get the pointer to the InodalProperty array
      call storage_getbase_int(&
          rtriangulation%h_InodalProperty,p_Idata)
      
      ! Write the data
      do ivt = 1, rtriangulation%NVT
        write (iunit,*) p_Idata(ivt)
      end do
      
      ! Write: 'KMM'
      write (iunit,*) 'KMM'
      
      ! Get the pointer to the IboundaryCpIdx and IverticesAtBoundary arrays
      call storage_getbase_int(&
          rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
      call storage_getbase_int(&
          rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
      
      ! Write the data
      ckmmstr = ''
      do ibct = 1, rtriangulation%NBCT
        ckmmstr = adjustl(trim(ckmmstr))//' '//&
            trim(sys_siL(p_IverticesAtBoundary(p_IboundaryCpIdx(ibct)),10))//' '//&
            trim(sys_siL(p_IverticesAtBoundary(p_IboundaryCpIdx(ibct+1)-1),10))
      end do
      write (iunit,*) trim(ckmmstr)
      
      ! Close the file, finish
      close(iunit)
      
    end subroutine tria_exportTriFile2D
    
  end subroutine tria_exportTriFile

  !****************************************************************************

!<subroutine>

  subroutine tria_exportPostScript(rtria, sfilename, Dbox, Dtrafo, &
                                   dlineWidth, bkeepAR)

!<description>
  ! Exports a 2D triangulation into a PostScript file.
  ! The created encapsulated PostScript file can be included into a LaTeX
  ! document by using the \includegraphics{} command.
  ! The mesh will be aligned at the bottom left corner of the drawing box.
!</description>

!<input>

  ! The triangulation that is to be exported. Must be a 2D mesh.
  type(t_triangulation), intent(in) :: rtria
  
  ! The filename of the PostScript file that is to be written.
  character(len=*), intent(in) :: sfilename
  
  ! OPTIONAL:
  ! The dimensions of the drawing box into which the mesh is to be exported.
  ! Dbox(1) = width of drawing box in millimeters (default: 100 mm)
  ! Dbox(2) = height of drawing box in millimeters (default: 100 mm)
  real(DP), dimension(2), optional, intent(in) :: Dbox
  
  ! OPTIONAL:
  ! A transformation matrix that the vertice coordinates are to be multiplied
  ! with. This matrix can be used to e.g. rotate a mesh by 90 degrees by
  ! setting Dtrafo to the corresponding rotation matrix.
  ! If not given, the transformation matrix is the identity matrix.
  real(DP), dimension(2,2), optional, intent(in) :: Dtrafo
  
  ! OPTIONAL:
  ! The line width for the PostScript file in millimeters.
  ! If not given, 0.1 mm is used.
  real(DP), optional, intent(in) :: dlineWidth
  
  ! OPTIONAL:
  ! If set to .true. (default), the aspect ratio of the mesh is kept.
  ! If set to .false., the mesh is stretched such that the bounding box of
  ! the mesh is equal to the drawing box.
  logical, optional, intent(in) :: bkeepAR
  
!</input>

!</subroutine>

  ! local variables
  real(DP) :: ddet,dwidth
  real(DP), dimension(2) :: Dv, Db, DbboxPS
  real(DP), dimension(2,2) :: Dt
  real(DP), dimension(:,:), pointer :: p_Dcoords
  integer, dimension(:,:), pointer :: p_Iedges
  integer :: iunit, i, NVT, NMT
  logical :: bkeepAspectRatio
  
  ! Bounding box of mesh
  real(DP) :: dbboxMinX, dbboxMinY, dbboxMaxX, dbboxMaxY,&
              dbboxWidth, dbboxHeight
  
  ! Scaling factors
  real(DP) :: dscaleX, dscaleY

  ! PostScript units are 'points'
  ! 1 inch = 25.4 millimeters
  ! 1 inch = 72 points
  ! => 1 millimeter = 72 / 25.4 points
  real(DP), parameter :: MM2PTS = 72.0_DP / 25.4_DP
  
    ! Intitialise default values
    Db(1) = 100.0_DP            ! drawing box dimensions in millimeters
    Db(2) = 100.0_DP
    Dt(1,1) = 1.0_DP            ! transformation matrix
    Dt(1,2) = 0.0_DP
    Dt(2,1) = 0.0_DP
    Dt(2,2) = 1.0_DP
    dwidth = 0.1_DP             ! line width in millimeters
    bkeepAspectRatio = .true.   ! self explaining
  
    ! First of all, let us make sure the drawing box is not empty.
    if(present(Dbox)) then
    
      if((Dbox(1) .le. SYS_EPSREAL) .or. (Dbox(2) .le. SYS_EPSREAL)) then
        call output_line('Drawing box is invalid!', OU_CLASS_ERROR,&
                         OU_MODE_STD, 'tria_exportPostScript')
        call sys_halt()
      end if
      
      Db = Dbox
      
    end if
      
    ! And let us make sure the transformation matrix is regular - otherwise
    ! we would divide by zero later!
    if(present(Dtrafo)) then
    
      ! Calculate determinant of trafo matrix
      ddet = Dtrafo(1,1)*Dtrafo(2,2) - Dtrafo(1,2)*Dtrafo(2,1)
      
      if(abs(ddet) .le. SYS_EPSREAL) then
        call output_line('Transformation matrix is singular!', OU_CLASS_ERROR,&
                         OU_MODE_STD, 'tria_exportPostScript')
        call sys_halt()
      end if
      
      Dt = Dtrafo
      
    end if
    
    ! And make sure the line width is positive.
    if(present(dlineWidth)) then
    
      if(dlineWidth .le. SYS_EPSREAL) then
        call output_line('Line width must be positive!', OU_CLASS_ERROR,&
                         OU_MODE_STD, 'tria_exportPostScript')
        call sys_halt()
      end if
      
      dwidth = dlineWidth
      
    end if
    
    ! There is no way the caller can mess up with this parameter...
    if(present(bkeepAR)) bkeepAspectRatio = bkeepAR
    
    ! Now make sure the triangulation is a 2D mesh.
    if(rtria%ndim .ne. NDIM2D) then
      call output_line('Only 2D triangulations supported!', OU_CLASS_ERROR,&
                       OU_MODE_STD, 'tria_exportPostScript')
      call sys_halt()
    end if
    
    ! Okay, get the necessary information from the mesh.
    NVT = rtria%NVT
    NMT = rtria%NMT
    call storage_getbase_double2D(rtria%h_DvertexCoords, p_Dcoords)
    call storage_getbase_int2D(rtria%h_IverticesAtEdge, p_Iedges)
    
    ! We now need to calculate the bounding box of the mesh.
    ! So get the first vertice and initialise the bounding box to it.
    Dv = matmul(Dt, p_Dcoords(1:2,1))
    dbboxMinX = Dv(1)
    dbboxMaxX = Dv(1)
    dbboxMinY = Dv(2)
    dbboxMaxY = Dv(2)
    
    ! And loop through the rest of the vertices.
    do i = 2, NVT
      
      ! Calculate transformed vertice
      Dv = matmul(Dt, p_Dcoords(1:2,i))
      
      ! And update the bounding box.
      dbboxMinX = min(dbboxMinX, Dv(1))
      dbboxMaxX = max(dbboxMaxX, Dv(1))
      dbboxMinY = min(dbboxMinY, Dv(2))
      dbboxMaxY = max(dbboxMaxY, Dv(2))
      
    end do
    
    ! Calculate the dimensions of the bounding box
    dbboxWidth  = dbboxMaxX - dbboxMinX
    dbboxHeight = dbboxMaxY - dbboxMinY
    
    ! Make sure the bounding box is fully-dimensional
    if((dbboxWidth .le. SYS_EPSREAL) .or. (dbboxHeight .le. SYS_EPSREAL)) then
      call output_line('Triangulation is not a 2D domain!', OU_CLASS_ERROR,&
                       OU_MODE_STD, 'tria_exportPostScript')
      call sys_halt()
    end if
    
    ! Calculate the scaling parameters:
    dscaleX = (Db(1) / dbboxWidth) * MM2PTS
    dscaleY = (Db(2) / dbboxHeight) * MM2PTS
    DbboxPS(1:2) = Db(1:2)
    
    ! Do we have to keep the aspect ratio?
    if(bkeepAspectRatio) then
      
      ! Yes, so choose the minimum scaling factor.
      dscaleX = min(dscaleX, dscaleY)
      dscaleY = dscaleX
    
      ! And calculate the bounding box for EPS
      DbboxPS(1) = dscaleX * dbboxWidth
      DbboxPS(2) = dscaleY * dbboxHeight
      
    else
    
      ! The bounding box is equal to the drawing box in this case
      DbboxPS(1) = Db(1)
      DbboxPS(2) = Db(2)
    
    end if
    
    ! Okay, open a file for writing
    call io_openFileForWriting(sfilename,iunit,SYS_REPLACE,bformatted=.true.)
    
    ! Fail?
    if(iunit .le. 0) then
      call output_line('Failed to open file for writing!', OU_CLASS_ERROR,&
                       OU_MODE_STD, 'tria_exportPostScript')
      call sys_halt()
    end if
    
    ! Okay, write PostScript header
    write(iunit,'(A)') '%!!PS-Adobe-3.0 EPSF-3.0'
    
    ! Write bounding box
    write(iunit,'(A,F12.6,F12.6,F12.6,F12.6)') '%%BoundingBox: ', &
      0.0_DP, 0.0_DP, DbboxPS(1), DbboxPS(2)
    
    ! Write the line width
    write(iunit,'(F12.6,A)') real(dwidth*MM2PTS,dp), ' setlinewidth'
    
    ! Begin a new path
    write(iunit,'(A)') 'newpath'
    
    ! Now go through all edges
    do i = 1, NMT
    
      ! Get the first vertice
      Dv = matmul(Dt, p_Dcoords(1:2, p_Iedges(1,i)))
      
      ! Transform the coordinates
      Dv(1) = (Dv(1) - dbboxMinX) * dscaleX
      Dv(2) = (Dv(2) - dbboxMinY) * dscaleY
      
      ! Write first vertice
      write(iunit,'(F12.6,F12.6,A)') Dv(1), Dv(2), ' moveto'
      
      ! Get the second vertice
      Dv = matmul(Dt, p_Dcoords(1:2, p_Iedges(2,i)))
      
      ! Transform the coordinates
      Dv(1) = (Dv(1) - dbboxMinX) * dscaleX
      Dv(2) = (Dv(2) - dbboxMinY) * dscaleY
      
      ! Write second vertice
      write(iunit,'(F12.6,F12.6,A)') Dv(1), Dv(2), ' lineto'
    
    end do
    
    ! Draw the path
    write(iunit,'(A)') 'stroke'
    
    ! And show the page
    write(iunit,'(A)') 'showpage'
    
    ! Close the file
    close(iunit)
    
    ! That is it

  end subroutine
  
  !************************************************************************

!<subroutine>

  subroutine tria_searchBoundaryNode(inode, rtriangulation, iindex)

!<description>
  ! This routine accepts as inode a vertex or an edge number of a boundary 
  ! vertex/edge (1..NVT=vertex, NVT+1..NVT+NMT=edge) and determines the 
  ! appropriate index of that vertex in the IverticesAtBoundary/
  ! IedgesAtBoundary-array.
!</description>

!<input>
  ! The boundary node to search for. A node number 1..NVT will search for
  ! a boundary vertex. A boundary node number NVT+1..NVT+NMT will search
  ! for boundary edges.
  integer, intent(in) :: inode

  ! The triangulation structure where to search the boundary node.
  type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<output>
  ! If inode is a boundary vertex: The index of the inode in IboundaryVertexPos.
  ! If inode is a boundary edge: The index of the inode in IboundaryEdgePos.
  ! =0 if inode was not found (e.g. because inode is not on the boundary e.g.).
  integer, intent(out) :: iindex
!</output>

!</subroutine>

    ! local variables
    integer :: ibct
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:,:), pointer :: p_InodePos
    integer :: ipos,ileft,iright
    
    ! Get the boundary component of the vertex
    call storage_getbase_int (rtriangulation%h_InodalProperty,p_InodalProperty)
    ibct = p_InodalProperty(inode)

    ! We can quit if ibct=0: The node is not on the boundary.
    if (ibct .eq. 0) then
      iindex = 0
      return
    end if
    
    ! Do we have a vertex or an edge number?
    if (inode .le. rtriangulation%NVT) then
      ! Search in the IboundaryVertexPos array
      call storage_getbase_int2d (rtriangulation%h_IboundaryVertexPos,p_InodePos)
    else
      ! Search in the IboundaryEdgePos array
      call storage_getbase_int2d (rtriangulation%h_IboundaryEdgePos,p_InodePos)
    end if

    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Use bisection search to find the node in the array.
    ileft = p_IboundaryCpIdx(ibct)
    iright = p_IboundaryCpIdx(ibct+1)-1
    
    if (p_InodePos(1,ileft) .eq. inode) then
      ! Return the index in the node array.
      iindex = p_InodePos(2,ileft)
    elseif (p_InodePos(1,iright) .eq. inode) then
      ! Return the index in the node array.
      iindex = p_InodePos(2,iright)
    else
      do while (ileft .lt. iright)
        ipos = (ileft+iright)/2
        if (p_InodePos(1,ipos) .gt. inode) then
          iright = ipos
        elseif (p_InodePos(1,ipos) .lt. inode) then
          ileft = ipos
        else
          ! We found the node. Return the index in the node array.
          iindex = p_InodePos(2,ipos)
          exit
        end if
      end do
    end if
    
  end subroutine tria_searchBoundaryNode

  ! ***************************************************************************

!<subroutine>

  subroutine tria_generateSubdomain(rtriangulation, Ielements, rtriaDest, rboundary)
  
!<description>
  ! Generates a triangulation structure for a subdomain of a larger domain.
  !
  ! The routine accepts a set of elements Ielements on the mesh rtriangulation,
  ! extracts all these cells and forms a new triangulation structure from it.
!</description>

!<input>
  ! Source triangulation; provides the 'parent' domain. A subdomain will
  ! be extracted from this.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A list of elements in rtriangulation that form the subdomain.
  integer, dimension(:), intent(in) :: Ielements
  
  ! OPTIONAL: A boundary structure that defines the parametrisation of the boumdary.
  type(t_boundary), intent(in), optional :: rboundary
!</input>

!<output>
  ! Destination triangulation structure. Receives the subdomain that consists
  ! only of the elements in Ielements. This will be a 'raw' mesh!
  type(t_triangulation), intent(out) :: rtriaDest
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoordsSrc
    real(DP), dimension(:,:), pointer :: p_DvertexCoordsDest
    real(DP), dimension(:), pointer :: p_DvertexParSrc
    real(DP), dimension(:), pointer :: p_DvertexParDest
    integer, dimension(:,:), pointer :: p_IverticesAtElementSrc
    integer, dimension(:,:), pointer :: p_IverticesAtElementDest
    integer, dimension(:,:), pointer :: p_IedgesAtElementSrc
    integer, dimension(:,:), pointer :: p_IedgesAtElementDest
    integer, dimension(:,:), pointer :: p_IfacesAtElementSrc
    integer, dimension(:,:), pointer :: p_IfacesAtElementDest
    integer, dimension(:), pointer :: p_IedgesAtVertex,p_IedgesAtVertexIdx
    integer, dimension(:,:), pointer :: p_IelementsAtEdge2D
    integer, dimension(:), pointer :: p_InodalPropertySrc
    integer, dimension(:), pointer :: p_InodalPropertyDest
    integer, dimension(:), pointer :: p_Idata,p_IdataInverse
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), allocatable :: IverticesAtBoundaryTmp
    integer, dimension(:), pointer :: p_ImacroNodalPropertySrc
    integer, dimension(:), pointer :: p_ImacroNodalPropertyDest
    integer :: ivt, ivt2, iel,ielsrc,ivtpos,ivtsource,ivtdest, ivtstart
    integer :: iedge, imt
    integer :: idim, ive, htemp, htempinverse, ibc,NVTsrc,NMTsrc
    integer, dimension(2) :: Isize
    integer :: ivertexAtElement

    NVTsrc = rtriangulation%NVT
    NMTsrc = rtriangulation%NMT

    ! Set up basic information
    rtriaDest%ndim = rtriangulation%ndim
    rtriaDest%NEL = size(Ielements)
    rtriaDest%NMT = 0
    rtriaDest%NNVE = rtriangulation%NNVE
    rtriaDest%NNEE = rtriangulation%NNEE
    
    ! We basically have NBCT boundary components like in the original triangulation.
    ! But there is still one 'blind' boundary component more that collects vertices/
    ! edges from inside of the domain rtriangulation that got boundary in rtriaDest.
    ! This is not counted in NBCT!
    rtriaDest%NBCT = rtriangulation%NBCT
    
    ! And we assume that this gives now a 'blind' boundary component.
    rtriaDest%NblindBCT = 1

    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriangulation%NNVE,rtriaDest%NEL/)
    call storage_new ('tria_generateSubdomain', 'KVERT', Isize, ST_INT, &
        rtriaDest%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)

    ! Get the pointer to the IverticesAtElement array
    call storage_getbase_int2D(&
        rtriaDest%h_IverticesAtElement,p_IverticesAtElementDest)
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_IverticesAtElementSrc)

    ! Temp array for regeneration of vertex numbers.
    call storage_new ('tria_generateSubdomain', 'itemp', &
        rtriangulation%NVT, ST_INT, htemp, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(htemp,p_Idata)
        
    ! Find all the vertices in the element set Ielements.
    ! For every vertex that was found, but bit 0 in the corresponding
    ! tag in Idata to 1, so the value gets 1. Later, we will sum up
    ! all these ones to get the number of vertices in our subdomain.
    do iel=1,size(Ielements)
      ielsrc = Ielements(iel)
      do ive=1,rtriaDest%NNVE
        ivertexAtElement = p_IverticesAtElementSrc(ive,ielsrc)
        p_Idata(ivertexAtElement) = 1
      end do
    end do
    
    ! Sum up all the ones to get NVT and the new vertex numbering.
    ! Afterwards, the array may have this form:
    !   0 0 0 0 0 1 1 1 1 1 2 3 3 3 3 4 ...
    ! The number changes in that moment, where a vertex is 'new'.
    ! This is then a mapping from the 'global' vertices to the 
    ! subdomain vertices. We will later access only a very small
    ! part of that array, namely the part ' 1 2 3 4 ...', forgetting
    ! the intermediate 'repetition' of the numbers. Nevertheless,
    ! this is a very easy method of calculating NVT and the new numbering
    ! simultaneously!
    do ivt=2,size(p_Idata)
      p_Idata(ivt) = p_Idata(ivt-1)+p_Idata(ivt)
    end do
    rtriaDest%NVT = p_Idata(size(p_Idata))
    
    ! Set up the inverse of the temp array -- to calculate from a 'local'
    ! vertex the corresponding 'global' vertex number.
    call storage_new ('tria_generateSubdomain', 'itempinverse', &
        rtriaDest%NVT, ST_INT, htempinverse, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(htempinverse,p_IdataInverse)
    
    ivt2 = 0
    do ivt = 1,size(p_Idata)
      ! When the number changes, we found a new vertex
      if (p_Idata(ivt) .ne. ivt2) then
        ivt2 = p_Idata(ivt)
        p_IdataInverse (p_Idata(ivt)) = ivt
      end if
    end do
    
    ! Allocate memory for the basic arrays on the heap
    ! array of size(dimension, NVT)
    Isize = (/rtriangulation%ndim,rtriaDest%NVT/)
    call storage_new ('tria_generateSubdomain', 'DCORVG',&
        Isize, ST_DOUBLE, &
        rtriaDest%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2Ddest is the pointer to the coordinate array
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_DvertexCoordsSrc)
    call storage_getbase_double2D(&
        rtriaDest%h_DvertexCoords,p_DvertexCoordsDest)

    ! Allocate memory for InodalProperty 
    call storage_new ('tria_generateSubdomain', 'KNPR', &
        rtriaDest%NVT, ST_INT, &
        rtriaDest%h_InodalProperty, ST_NEWBLOCK_ZERO)

    ! Get pointers to the nodal property array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_InodalPropertySrc)
    call storage_getbase_int(&
        rtriaDest%h_InodalProperty,p_InodalPropertyDest)
    
    ! Transfer the content of the IverticesAtElement-array for the
    ! elements in the list. Use Itemp to renumber the vertices.
    ! Get the coordinates of the points as well.
    do iel=1,size(Ielements)
      ielsrc = Ielements(iel)
      do ive=1,rtriaDest%NNVE
      
        ivt = p_IverticesAtElementSrc(ive,ielsrc)
        ivt2 = p_Idata(ivt)
        p_IverticesAtElementDest(ive,iel) = ivt2
        
        do idim = 1,Isize(1)
          p_DvertexCoordsDest(idim,ivt2) = p_DvertexCoordsSrc(idim,ivt)
        end do
        
        ! Transfer the nodal property from the vertex.
        ! The information is probably overwritten multiple times,
        ! but we do not care here...
        p_InodalPropertyDest(ivt2) = p_InodalPropertySrc(ivt)
        
      end do
    end do

    ! Initialise InelOfType.
    !    
    ! Loop through the elements and determine how many elements
    ! of each element type we have.
    rtriaDest%InelOfType(:) = 0
    do iel=1,rtriaDest%NEL
      ! start at the last index of element iel down to the first
      do ive=rtriaDest%NNVE,1,-1
        if (p_IverticesAtElementDest(ive,iel) .ne. 0) then
          rtriaDest%InelOfType(ive) = rtriaDest%InelOfType(ive)+1
          exit
        end if
      end do
    end do
    
    ! Dimension-dependent part: What is this for a mesh?
    select case(rtriangulation%ndim)
    case (NDIM1D)
      ! Nothing implemented here; probably noting to do.
      ! Has to be checked...
      
    case (NDIM2D)
    
      ! We need some additional arrays to proceed with the following tasks.
      call tria_genElementsAtVertex1D2D (rtriaDest)
      call tria_genNeighboursAtElement2D (rtriaDest)
      call tria_genEdgesAtElement2D (rtriaDest)
      call tria_genElementsAtEdge2D (rtriaDest)
      call tria_genVerticesAtEdge2D (rtriaDest)
      call tria_genEdgesAtVertex2D (rtriaDest)
    
      ! Loop through all edges adjacent to every vertex of our
      ! new triangulation. If there is one edge with only one
      ! adjacent element, the edge is a boundary edge -- and
      ! so it is the vertex!
      call storage_getbase_int2D (rtriaDest%h_IelementsAtEdge,&
          p_IelementsAtEdge2D)
      call storage_getbase_int(rtriaDest%h_IedgesAtVertexIdx, &
          p_IedgesAtVertexIdx)
      call storage_getbase_int(rtriaDest%h_IedgesAtVertex, &
          p_IedgesAtVertex)  
      
      vertexloop: do ivt=1,rtriaDest%NVT
      
        do imt = p_IedgesAtVertexIdx(ivt),p_IedgesAtVertexIdx(ivt+1)-1
        
          iedge = p_IedgesAtVertex(imt)
          
          ! Find the inner edges that are going to be boundary edges
          if ((p_IelementsAtEdge2D(2,iedge) .eq. 0) .and. &
              (p_InodalPropertyDest(ivt) .eq. 0)) then
          
            ! That is a boundary edge. So the vertex is also one.
            ! Remember #NBCT+1 as nodal property for such nodes.
            p_InodalPropertyDest(ivt) = rtriaDest%NBCT+1
            
            ! Next vertex
            cycle vertexloop
          
          end if
        
        end do
      
      end do vertexloop
      
      ! Generate boundary information
      call genRawBoundary2D (rtriaDest)

      ! If we have boundary information, we can extract the parameter values
      ! of the vertices on the physical boundary.    
      if (present(rboundary)) then

        ! Allocate memory for DvertexParameterValue
        call storage_new ('tria_generateSubdomain', &
            'DVBDP', rtriaDest%NVBD, &
            ST_DOUBLE, rtriaDest%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
        
        call storage_getbase_double (&
            rtriangulation%h_DvertexParameterValue,p_DvertexParSrc)
            
        call storage_getbase_double (&
            rtriaDest%h_DvertexParameterValue,p_DvertexParDest)

        call storage_getbase_int (&
            rtriaDest%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
        if (rtriangulation%h_IboundaryVertexPos .eq. ST_NOHANDLE) then
          call output_line ('Boundary search arrays not initialised!.', &
                            OU_CLASS_ERROR,OU_MODE_STD,'tria_generateSubdomain')
          call sys_halt()
        end if
        
        ! Loop through all vertices on the boundary. Find out their position
        ! in the original boundary-vertex array and get their parameter values.
        do ivt2=1,size(p_IverticesAtBoundary)
        
          ivt = p_IverticesAtBoundary(ivt2)
        
          ! If this is a vertex on the real boundary...
          if (p_InodalPropertyDest(ivt) .le. rtriaDest%NBCT) then
        
            ! Search the vertex position
            call tria_searchBoundaryNode(p_IdataInverse(ivt),rtriangulation,ivtpos)
            
            ! Get the parameter value.
            p_DvertexParDest(ivt2) = p_DvertexParSrc(ivtpos)
            
          else
          
            ! Othewise, save -1.
            p_DvertexParDest(ivt2) = -1.0_DP
          
          end if
          
        end do
        
        ! Now we have an array like
        !  1.0 2.0 -1.0 3.0    1.0 -1.0 2.0  ...
        ! With some -1 inbetween for all vertices that do not belong to the physical
        ! boundary. We now 'shift' all these '-1'-nodes to the end of the array
        ! and that way assign them to the 'blind' boundary component.
        
        call storage_getbase_int(rtriaDest%h_IboundaryCpIdx, &
            p_IboundaryCpIdx)  
        
        allocate(IverticesAtBoundaryTmp(size(p_IverticesAtBoundary)))
        ivtdest = 1
        ivtpos = 1
        ivtstart = p_IboundaryCpIdx(1)
        
        ! Loop over all physical-boundary-BC`s; ignore any existing 'blind'
        ! boundary component.
        do ibc = 1,rtriangulation%NBCT
        
          ! Loop over all boundary vertices in that BC. Compress the BC.
          ! All boundary vertices not on the physical boundary are
          ! extracted to IverticesAtBoundaryTmp.
          
          do ivtsource = ivtstart,p_IboundaryCpIdx(ibc+1)-1
            if (p_DvertexParDest(ivtsource) .ne. -1.0_DP) then
              ! This is a vertex on the boundary. Copy it to the destination position.
              p_DvertexParSrc(ivtdest) = p_DvertexParSrc(ivtsource)
              p_IverticesAtBoundary(ivtdest) = p_IverticesAtBoundary(ivtsource)
              ivtdest = ivtdest + 1
            else
              ! Extract that DOF, so we can overwrite with the forthcoming DOF`s.
              IverticesAtBoundaryTmp(ivtdest) = p_IverticesAtBoundary(ivtsource)
              ivtpos = ivtpos + 1
            end if
          end do
          
          ! Remember the new start address where the DOF`s of the next
          ! boundary component are found; it is immediately overwritten
          ! by ivtdest!
          ivtstart = p_IboundaryCpIdx(ibc+1)
          
          ! ivtdest points now behind the DOF`s of boundary component ibc.
          ! Save the pointer to the p_IboundaryCpIdx array.
          p_IboundaryCpIdx(ibc+1) = ivtdest
          
        end do
        
        ! Ok, the array is compressed now. The only thing that remains is
        ! to copy the 'blind' vertices to the 'blind' boundary component.
        do ivtsource = 1,ivtpos-1
          p_IverticesAtBoundary(ivtdest) = IverticesAtBoundaryTmp(ivtsource)
          ivtdest = ivtdest + 1
        end do
        
        deallocate(IverticesAtBoundaryTmp)
        
        ! The last entry in p_IboundaryCpIdx does not have to be changed since
        ! it counts the total number of vertices on the boundary -- which
        ! has not changed.
        
      end if
    
    case (NDIM3D)
      ! Nothing implemented here; probably noting to do.
      ! Has to be checked...
  
    end select

    ! Now we can release the temp arrays.
    call storage_free (htemp)
    call storage_free (htempinverse)

    ! Propagate the macro nodal property array to the submesh if available.
    if (rtriangulation%h_ImacroNodalProperty .ne. ST_NOHANDLE) then
    
      call storage_getbase_int (rtriangulation%h_ImacroNodalProperty,&
          p_ImacroNodalPropertySrc)

      ! Allocate memory for that array
      call storage_new ('tria_generateSubdomain', &
          'KMCPR', rtriaDest%NVT+rtriaDest%NMT+rtriaDest%NAT+rtriaDest%NEL, &
          ST_INT, rtriaDest%h_ImacroNodalProperty, ST_NEWBLOCK_NOINIT)

      call storage_getbase_int (rtriaDest%h_ImacroNodalProperty,&
          p_ImacroNodalPropertyDest)
    
      ! Dimension-dependent part: What is this for a mesh?
      select case(rtriangulation%ndim)
      case (NDIM1D)
        ! Loop through the elements. For every element, copy the macro nodal property
        ! entry to the new mesh.
        do iel=1,rtriaDest%NEL
          
          ! Vertices -- in every mesh
          do ive = 1,ubound(p_IverticesAtElementSrc,1)
            if (p_IverticesAtElementSrc(ive,iel) .ne. 0) then
              p_ImacroNodalPropertyDest(p_IverticesAtElementSrc(ive,iel)) = &
                  p_ImacroNodalPropertySrc(p_IverticesAtElementDest(ive,Ielements(iel))) 
            end if
          end do
          
        end do

      case (NDIM2D)

        call storage_getbase_int2d (rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElementSrc)
        call storage_getbase_int2d (rtriaDest%h_IedgesAtElement,&
            p_IedgesAtElementDest)

        ! Loop through the elements. For every element, copy the macro nodal property
        ! entry to the new mesh.
        do iel=1,rtriaDest%NEL
          
          do ive = 1,ubound(p_IverticesAtElementSrc,1)
            if (p_IverticesAtElementSrc(ive,iel) .ne. 0) then
              ! Vertices
              p_ImacroNodalPropertyDest(p_IverticesAtElementDest(ive,iel)) = &
                  p_ImacroNodalPropertySrc(p_IverticesAtElementSrc(ive,Ielements(iel))) 
                  
              ! Edges
              p_ImacroNodalPropertyDest(p_IedgesAtElementDest(ive,iel)+rtriaDest%NVT) = &
                  p_ImacroNodalPropertySrc(p_IedgesAtElementSrc(ive,Ielements(iel))+rtriangulation%NVT) 
            end if
          end do
          
          ! Elements
          p_ImacroNodalPropertyDest(rtriaDest%NVT+rtriaDest%NMT+iel) = &
              p_ImacroNodalPropertySrc(rtriangulation%NVT+rtriangulation%NMT+Ielements(iel)) 
              
        end do
          
      case (NDIM3D)

        call storage_getbase_int2d (rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElementSrc)
        call storage_getbase_int2d (rtriaDest%h_IedgesAtElement,&
            p_IedgesAtElementDest)

        call storage_getbase_int2d (rtriangulation%h_IfacesAtElement,&
            p_IfacesAtElementSrc)
        call storage_getbase_int2d (rtriaDest%h_IfacesAtElement,&
            p_IfacesAtElementDest)

        ! Loop through the elements. For every element, copy the macro nodal property
        ! entry to the new mesh.
        do iel=1,rtriaDest%NEL

          ! Vertices          
          do ive = 1,ubound(p_IverticesAtElementSrc,1)
            if (p_IverticesAtElementSrc(ive,iel) .ne. 0) then
              p_ImacroNodalPropertyDest(p_IverticesAtElementDest(ive,iel)) = &
                  p_ImacroNodalPropertySrc(p_IverticesAtElementSrc(ive,Ielements(iel))) 
            end if
          end do

          ! Edges
          do ive = 1,ubound(p_IedgesAtElementSrc,1)
            if (p_IedgesAtElementSrc(ive,iel) .ne. 0) then
              p_ImacroNodalPropertyDest(p_IedgesAtElementDest(ive,iel)+rtriaDest%NVT) = &
                  p_ImacroNodalPropertySrc(p_IedgesAtElementSrc(ive,Ielements(iel))+rtriangulation%NVT) 
            end if
          end do

          ! Faces
          do ive = 1,ubound(p_IfacesAtElementSrc,1)
            if (p_IfacesAtElementSrc(ive,iel) .ne. 0) then
              p_ImacroNodalPropertyDest(p_IfacesAtElementDest(ive,iel)+rtriaDest%NVT+rtriaDest%NMT) = &
                  p_ImacroNodalPropertySrc(p_IfacesAtElementSrc(ive,Ielements(iel))&
                    +rtriangulation%NVT+rtriangulation%NMT) 
            end if
          end do
          
          ! Elements
          p_ImacroNodalPropertyDest(rtriaDest%NVT+rtriaDest%NMT+rtriaDest%NAT+iel) = &
              p_ImacroNodalPropertySrc(&
                  rtriangulation%NVT+rtriangulation%NMT+rtriangulation%NAT+Ielements(iel)) 
              
        end do
        
      end select
          
    end if

  contains

    ! ---------------------------------------------------------------

    subroutine genRawBoundary2D (rtriangulation)

    ! Auxiliary routine.
    ! This routine initialises basic boundary arrays and cleans up 
    ! a basic triangulation. That means:
    ! -> NVBD is calculated
    ! -> IboundaryCpIdx is created and generated
    ! -> IverticesAtBoundary is created and generated.
    !    The vertices are ordered for the boundary component according
    !    to IboundaryCpIdx but not ordered for their parameter value.
    
    ! Triangulation to be initialised with basic data.
    type(t_triangulation), intent(inout) :: rtriangulation
    
      ! local variables
      real(DP), dimension(:,:), pointer :: p_DvertexCoords
      real(DP), dimension(:), pointer :: p_DvertexParameterValue
      integer, dimension(:), pointer :: p_IboundaryCpIdx
      integer, dimension(:), pointer :: p_IverticesAtBoundary
      integer :: ivbd,ivt
      integer :: ibct
      integer, dimension(:), pointer :: p_InodalProperty

      ! Get the pointer to the InodalProperty array
      call storage_getbase_int(&
          rtriangulation%h_InodalProperty,p_InodalProperty)

      ! Calculate NVBD by simply counting how many elements
      ! in p_InodalProperty are <> 0.
      rtriangulation%NVBD = 0
      ivbd = 0
      do ivt=1,rtriangulation%NVT
        !IF (p_InodalProperty(ivt) .NE. 0) rtriangulation%NVBD = rtriangulation%NVBD+1
        if (p_InodalProperty(ivt) .gt. 0) ivbd = ivbd+1
      end do
      rtriangulation%NVBD = ivbd

      ! Allocate memory for IverticesAtBoundary.
      call storage_new ('genRawBoundary2D', &
          'KVBD', rtriangulation%NVBD, &
          ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
          
      ! Allocate memory for the boundary component index vector.
      ! We reserve NBCT+1(+1) elements here, where the NBCT+1`th element
      ! corresponds to the 'blind' boundary which came from inside of
      ! the domain.
      ! Initialise everything with zero!
      call storage_new ('genRawBoundary2D', &
          'KBCT', rtriangulation%NBCT+2, &
          ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
      
      ! Get pointers to the arrays
      call storage_getbase_int (&
          rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
          
      call storage_getbase_double2D (&
          rtriangulation%h_DvertexCoords,p_DvertexCoords)
          
      call storage_getbase_int (&
          rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
      
      ! The first element in p_IboundaryCpIdx is (as the head) always =1.
      p_IboundaryCpIdx(1) = 1

      ! Perform a first loop over all vertices to check which are on the
      ! boundary. Count them. This way, we at first set up the boundary
      ! component index vector p_IboundaryCpIdx.
      ! In this first step, we save the number of vertices on each 
      ! boundary component in p_IboundaryCpIdx(2:NBCT+1)!
      do ivt=1,rtriangulation%NVT
        if (p_InodalProperty(ivt) .gt. 0) then
          ibct = p_InodalProperty(ivt)
          
          ! Increase the number of vertices in that boundary component by 1.
          ! The number of vertices on boundary component i is saved here
          ! at p_IboundaryCpIdx(i+1) for later.
          ! Note that the array was initialised with zero during the creation
          ! process!
          p_IboundaryCpIdx(ibct+1) = p_IboundaryCpIdx(ibct+1)+1

        end if
      end do
      
      ! Sum up the number of vertices on each boundary component to get the
      ! actual index vector.
      do ibct = 2,rtriangulation%NBCT+2
        p_IboundaryCpIdx(ibct) = p_IboundaryCpIdx(ibct) + &
                                 p_IboundaryCpIdx(ibct-1)
      end do
      
      ! Shift the p_IboundaryCpIdx array by one position. That is a little trick in
      ! the use of p_IboundaryCpIdx!
      ! Imagine, we have 3 boundary components with 8,6 and 4 edges. 
      ! Before summing the entries up, p_IboundaryCpIdx may have looked like this:
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   8   6   4
      !
      ! Summing that up gave the actual indices:
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   9  15  19
      !
      ! which is already the correct index array.
      ! Shifting p_IboundaryCpIdx by 1 will give:
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   1   9  15
      
      p_IboundaryCpIdx(2:rtriangulation%NBCT+2) = p_IboundaryCpIdx(1:rtriangulation%NBCT+1)
      
      ! Then, we again loop through all vertices and collect those on the
      ! boundary. In that loop, we use p_IboundaryCpIdx(2:NBCT+2) as pointer and
      ! increase them for every point we find. The loop will behave like
      ! the first one, i.e. we will again find 8,6 and 4 vertices on the boundary components
      ! 1,2 and 3. Increasing the p_IboundaryCpIdx(2:NBCT+1) for boundary components
      ! 1..NBCT the same way as done in the first loop will therefore again lead to
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   9  15  19
      !
      ! Ok, let us catch the actual vertices.
      !      
      ! Check all vertices to find out, which vertices are on the boundary.
      !%OMP PARALLEL do PRIVATE(ibct,ivbd)
      do ivt = 1, rtriangulation%NVT
        if (p_InodalProperty(ivt) .gt. 0) then
          ! id of the boundary component
          ibct = p_InodalProperty(ivt)
          
          ! set ivbd to the number of vertices on that boundary component
          ! thus ivbd holds the current number of vertices found for
          ! boundary component ibct and ivbd represents the current
          ! position in the p_IverticesAtBoundary array
          ivbd = p_IboundaryCpIdx(ibct+1)
          ! we have found a new point on that boundary component
          ! so increase the number of points by one
          p_IboundaryCpIdx(ibct+1) = ivbd+1
          
          ! Store the vertex as boundary vertex
          p_IverticesAtBoundary (ivbd) = ivt
        end if
      end do
      !%OMP END PARALLEL do
        
    end subroutine genRawBoundary2D

  end subroutine tria_generateSubdomain

  ! ***************************************************************************

!<subroutine>

  subroutine tria_attachCells (rtriangulation,rcellSet,IvertexMapping,rtriaDest)
  
!<description>
  ! Extends a mesh by a set of cells. rcellSet defines a set of cells.
  ! These cells are attached to the mesh rtriangulation. The new mesh is
  ! created in rtriaDest and will be a 'raw' mesh.
!</description>

!<input>
  ! A mesh which is to be extended.
  type(t_triangulation), intent(inout) :: rtriangulation
  
  ! A set of cells which is to be attached to rtriangulation.
  ! At least, IverticesAtElement in this structure must be defined to create
  ! a valid mesh. If DvertexCoords is undefined, all new vertices are set to
  ! the origin. If InodalProperty is undefined, new vertices are treated
  ! as 'inner' vertices.
  type(t_cellSet), intent(in) :: rcellSet
  
  ! A vertex mapping. This mapping defines for every vertex in rcellSet
  ! the number of that vertex in rtriangulation which coincides with it.
  ! If a vertex in rcellSet is not related to a vertex in rtriangulation,
  ! i.e. the vertex is a new vertex, the corresponding entry in IvertexMapping
  ! must be =0.
  integer, dimension(:), intent(in) :: IvertexMapping
!</input>

!<output>
  ! Receives the extended mesh.
  type(t_triangulation), intent(out) :: rtriaDest
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoordsNew
    integer, dimension(:,:), pointer :: p_IverticesAtElementNew
    integer, dimension(:), pointer :: p_InodalPropertyNew
    real(DP), dimension(:), pointer :: p_DvertexParNew
    integer :: h_IvertexMappingAux
    integer, dimension(:), pointer :: p_IvertexMappingAux
    integer, dimension(2) :: Isize
    integer :: i,j,idestpos,ivt,ivt2
    integer :: iel
    integer :: NVT,ipoint,ivtpos
    integer :: h_DvertParamTmp
    real(DP), dimension(:), pointer :: p_DvertParamTmp
    
    integer, dimension(:,:), pointer :: p_IverticesAtElementSrc
    integer, dimension(:,:), pointer :: p_IverticesAtElementDest
    real(DP), dimension(:,:), pointer :: p_DvertexCoordsSrc,p_DvertexCoordsDest
    integer, dimension(:), pointer :: p_InodalPropertySrc,p_InodalPropertyDest
    real(DP), dimension(:), pointer :: p_DvertexParSrc,p_DvertexParDest
    integer, dimension(:), pointer :: p_IverticesAtBoundary

    ! Get the arrays from the cell set
    p_DvertexCoordsNew => rcellSet%p_DvertexCoords
    p_IverticesAtElementNew => rcellSet%p_IverticesAtElement
    p_InodalPropertyNew => rcellSet%p_InodalProperty
    p_DvertexParNew => rcellSet%p_DallVerticesParameterValue
    if (.not. associated(p_IverticesAtElementNew)) then
      call output_line('IverticesAtElementNew in the cell set undefined!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'tria_attachCells')
      call sys_halt()
    end if
    
    ! Set up basic information
    rtriaDest%ndim = rtriangulation%ndim
    rtriaDest%NMT = 0
    rtriaDest%NNVE = max(rtriangulation%NNVE,ubound(p_IverticesAtElementNew,1))
    rtriaDest%NNEE = rtriangulation%NNEE
    rtriaDest%NBCT = rtriangulation%NBCT
    rtriaDest%NblindBCT = rtriangulation%NblindBCT
    rtriaDest%NEL = rtriangulation%NEL + rcellSet%NEL

    ! Allocate memory for IverticesAtElement.
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriaDest%NNVE,rtriaDest%NEL/)
    call storage_new ('tria_attachCells', 'KVERT', Isize, ST_INT, &
        rtriaDest%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)

    ! Get the pointer to the IverticesAtElement array
    call storage_getbase_int2D(&
        rtriaDest%h_IverticesAtElement,p_IverticesAtElementDest)
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_IverticesAtElementSrc)

    ! Copy the first part of IverticesAtElement to the destination array; 
    ! old cells are new cells.
    do j=1,rtriangulation%NEL
      do i=1,ubound(p_IverticesAtElementSrc,1)
        p_IverticesAtElementDest(i,j) = p_IverticesAtElementSrc(i,j)
      end do
    end do
    
    ! Now we have to correct the vertex numbers. The new vertices are
    ! expected in the numbering 1,2,3,4,... and must be modified to
    ! fulfil the numbering NVT+1,NVT+2,... to be new vertices.
    !
    ! For that purpose, we create a new vertex mapping array. This array basically
    ! receives the original vertex mapping array for all vertices that coincide
    ! with existing vertices. For new vertices, we create here their
    ! 'destination' number.
    call storage_new ('tria_attachCells', 'IvertexMapping', size(IvertexMapping), &
        ST_INT, h_IvertexMappingAux, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (h_IvertexMappingAux,p_IvertexMappingAux)
    
    ipoint = rtriangulation%NVT
    do i=1,size(p_IvertexMappingAux)
      if (IvertexMapping(i) .eq. 0) then
        ! New number, following NVT.
        ipoint = ipoint+1
        p_IvertexMappingAux(i) = ipoint
      else
        ! Existing point
        p_IvertexMappingAux(i) = IvertexMapping(i)
      end if
    end do
    
    ! ipoint is now the total number of vertices in the new mesh.
    rtriadest%NVT = ipoint
    
    ! Position of the new cells
    idestPos = rtriangulation%NEL
    
    ! Now, attach the cell connectivity and simultaneously renumber the vertices
    do j=1,rcellSet%NEL
      do i=1,ubound(p_IverticesAtElementNew,1)
        p_IverticesAtElementDest(i,j+idestPos) = &
            p_IvertexMappingAux(p_IverticesAtElementNew(i,j))
      end do
    end do
    
    ! Allocate memory for the basic arrays on the heap.
    ! Initialise with zero, so new points are originally at the origin.
    ! array of size(dimension, NVT)
    Isize = (/rtriangulation%ndim,rtriaDest%NVT/)
    call storage_new ('tria_generateSubdomain', 'DCORVG',&
        Isize, ST_DOUBLE, &
        rtriaDest%h_DvertexCoords, ST_NEWBLOCK_ZERO)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2Ddest is the pointer to the coordinate array
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_DvertexCoordsSrc)
    call storage_getbase_double2D(&
        rtriaDest%h_DvertexCoords,p_DvertexCoordsDest)
    
    ! Copy the first part of DvertexCoords to the destination array; 
    ! old vertices are new vertices.
    do j=1,rtriangulation%NVT
      do i=1,ubound(p_DvertexCoordsSrc,1)
        p_DvertexCoordsDest(i,j) = p_DvertexCoordsSrc(i,j)
      end do
    end do
    
    ! If new vertex coordinates are present, initialise the new coordinates.
    if (associated(p_DvertexCoordsNew)) then
    
      do j=1,rcellSet%NVT
        ivt = p_IvertexMappingAux(j)
        do i=1,ubound(p_DvertexCoordsNew,1)
          p_DvertexCoordsDest(i,ivt) = p_DvertexCoordsNew(i,j)
        end do
      end do
      
    end if
    
    ! Allocate memory for InodalProperty 
    call storage_new ('tria_generateSubdomain', 'KNPR', &
        rtriaDest%NVT, ST_INT, &
        rtriaDest%h_InodalProperty, ST_NEWBLOCK_ZERO)

    ! Get pointers to the nodal property array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_InodalPropertySrc)
    call storage_getbase_int(&
        rtriaDest%h_InodalProperty,p_InodalPropertyDest)

    ! Copy the first part of InodalProperty
    call lalg_copyVectorInt(p_InodalPropertySrc(1:rtriangulation%NVT),&
        p_InodalPropertyDest(1:rtriangulation%NVT))
    
    ! If new nodal property information tags are present, initialise the nodal property
    ! array for the new vertices.
    if (associated(p_InodalPropertyNew)) then

      do j=1,rcellSet%NVT
        p_InodalPropertyDest(p_IvertexMappingAux(j)) = p_InodalPropertyNew(j)
      end do
    
    end if
    
    ! Initialise InelOfType.
    !    
    ! Loop through the elements and determine how many elements
    ! of each element type we have.
    rtriaDest%InelOfType(:) = 0
    do iel=1,rtriaDest%NEL
      ! start at the last index of element iel down to the first
      do i=rtriaDest%NNVE,1,-1
        if (p_IverticesAtElementDest(i,iel) .ne. 0) then
          rtriaDest%InelOfType(i) = rtriaDest%InelOfType(i)+1
          exit
        end if
      end do
    end do
    
    ! Generate basic boundary information
    call genRawBoundary2D (rtriaDest)
  
    ! If we have boundary information, we can extract the parameter values
    ! of the vertices on the physical boundary.    
    if (rtriangulation%h_DvertexParameterValue .ne. ST_NOHANDLE) then

      ! Allocate memory for DvertexParameterValue
      call storage_new ('tria_generateSubdomain', &
          'DVBDP', rtriaDest%NVBD, &
          ST_DOUBLE, rtriaDest%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
      
      call storage_getbase_double (&
          rtriangulation%h_DvertexParameterValue,p_DvertexParSrc)
          
      call storage_getbase_double (&
          rtriaDest%h_DvertexParameterValue,p_DvertexParDest)

      call storage_getbase_int (&
          rtriaDest%h_IverticesAtBoundary,p_IverticesAtBoundary)
      
      if (rtriangulation%h_IboundaryVertexPos .eq. ST_NOHANDLE) then
        call output_line ('Boundary search arrays not initialised!.', &
                          OU_CLASS_ERROR,OU_MODE_STD,'tria_generateSubdomain')
        call sys_halt()
      end if
      
      NVT = rtriangulation%NVT

      ! Create temp array that holds the parameter values for all new
      ! vertices on the boundary and is undefined for all other new vertices.      
      call storage_new ('tria_attachCells', 'DvertParamTmp', ipoint-NVT, &
          ST_DOUBLE, h_DvertParamTmp, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double (h_DvertParamTmp,p_DvertParamTmp)
      
      if (associated(p_DvertexParNew)) then
        ! Copy the parameter values, reordered according to the vertex mapping.
        ! Copy only the parameter values of the 'new' vertices.
        do i=1,rcellSet%NVT
          if (p_IvertexMappingAux(i) .gt. NVT) &
            p_DvertParamTmp(p_IvertexMappingAux(i)-NVT) = p_DvertexParNew(i)
        end do
      else
        ! We do not know anything about the vertices; set their parameter value
        ! to -1.
        do i=1,ipoint-NVT
          p_DvertParamTmp(p_IvertexMappingAux(i)-NVT) = -1.0_DP
        end do
      end if
      
      ! Loop through all vertices on the boundary. Find out their position
      ! in the original boundary-vertex array and get their parameter values.
      do ivt2=1,size(p_IverticesAtBoundary)
      
        ivt = p_IverticesAtBoundary(ivt2)
        
        ! If this is a vertex on the real boundary...
        if (p_InodalPropertyDest(ivt) .le. rtriaDest%NBCT) then
      
          ! Old or new vertex. ivt=vertex number in new mesh, NVT=#vertices in old mesh.
          if (ivt .le. NVT) then
      
            ! Search the vertex position
            call tria_searchBoundaryNode(ivt,rtriangulation,ivtpos)
            
            ! Get the parameter value.
            p_DvertexParDest(ivt2) = p_DvertexParSrc(ivtpos)
            
          else
           
            ! New vertex; take the parameter value from the cell set.
            ! Take the parameter value from p_DvertParamTmp which collects the new
            ! parameter values in reordered order, not including the old vertices.
            p_DvertexParDest(ivt2) = p_DvertParamTmp(ivt-NVT)
          
          end if
          
        else
        
          ! Othewise, save -1.
          p_DvertexParDest(ivt2) = -1.0_DP
        
        end if
        
      end do
      
      ! Release memory.
      call storage_free (h_DvertParamTmp)
      
      ! Something`s still missing here!
      ! We have to delete all vertices from the boundary arrays that returned
      ! from 'blind vertex' state to 'inner vertex' state, e.g. that do not belong
      ! to the 'blind' boundary component anymore...
      
    end if
  
    ! Release temp memory
    call storage_free (h_IvertexMappingAux)
    
  contains

    ! ---------------------------------------------------------------

    subroutine genRawBoundary2D (rtriangulation)

    ! Auxiliary routine.
    ! This routine initialises basic boundary arrays and cleans up 
    ! a basic triangulation. That means:
    ! -> NVBD is calculated
    ! -> IboundaryCpIdx is created and generated
    ! -> IverticesAtBoundary is created and generated.
    !    The vertices are ordered for the boundary component according
    !    to IboundaryCpIdx but not ordered for their parameter value.
    
    ! Triangulation to be initialised with basic data.
    type(t_triangulation), intent(inout) :: rtriangulation
    
      ! local variables
      real(DP), dimension(:,:), pointer :: p_DvertexCoords
      real(DP), dimension(:), pointer :: p_DvertexParameterValue
      integer, dimension(:), pointer :: p_IboundaryCpIdx
      integer, dimension(:), pointer :: p_IverticesAtBoundary
      integer :: ivbd,ivt
      integer :: ibct
      integer, dimension(:), pointer :: p_InodalProperty

      ! Get the pointer to the InodalProperty array
      call storage_getbase_int(&
          rtriangulation%h_InodalProperty,p_InodalProperty)

      ! Calculate NVBD by simply counting how many elements
      ! in p_InodalProperty are <> 0.
      rtriangulation%NVBD = 0
      ivbd = 0
      do ivt=1,rtriangulation%NVT
        !IF (p_InodalProperty(ivt) .NE. 0) rtriangulation%NVBD = rtriangulation%NVBD+1
        if (p_InodalProperty(ivt) .gt. 0) ivbd = ivbd+1
      end do
      rtriangulation%NVBD = ivbd

      ! Allocate memory for IverticesAtBoundary.
      call storage_new ('genRawBoundary2D', &
          'KVBD', rtriangulation%NVBD, &
          ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
          
      ! Allocate memory for the boundary component index vector.
      ! We reserve NBCT+1(+1) elements here, where the NBCT+1`th element
      ! corresponds to the 'blind' boundary which came from inside of
      ! the domain.
      ! Initialise everything with zero!
      call storage_new ('genRawBoundary2D', &
          'KBCT', rtriangulation%NBCT+2, &
          ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
      
      ! Get pointers to the arrays
      call storage_getbase_int (&
          rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
          
      call storage_getbase_double2D (&
          rtriangulation%h_DvertexCoords,p_DvertexCoords)
          
      call storage_getbase_int (&
          rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
      
      ! The first element in p_IboundaryCpIdx is (as the head) always =1.
      p_IboundaryCpIdx(1) = 1

      ! Perform a first loop over all vertices to check which are on the
      ! boundary. Count them. This way, we at first set up the boundary
      ! component index vector p_IboundaryCpIdx.
      ! In this first step, we save the number of vertices on each 
      ! boundary component in p_IboundaryCpIdx(2:NBCT+1)!
      do ivt=1,rtriangulation%NVT
        if (p_InodalProperty(ivt) .gt. 0) then
          ibct = p_InodalProperty(ivt)
          
          ! Increase the number of vertices in that boundary component by 1.
          ! The number of vertices on boundary component i is saved here
          ! at p_IboundaryCpIdx(i+1) for later.
          ! Note that the array was initialised with zero during the creation
          ! process!
          p_IboundaryCpIdx(ibct+1) = p_IboundaryCpIdx(ibct+1)+1

        end if
      end do
      
      ! Sum up the number of vertices on each boundary component to get the
      ! actual index vector.
      do ibct = 2,rtriangulation%NBCT+2
        p_IboundaryCpIdx(ibct) = p_IboundaryCpIdx(ibct)+p_IboundaryCpIdx(ibct-1)
      end do
      
      ! Shift the p_IboundaryCpIdx array by one position. That is a little trick in
      ! the use of p_IboundaryCpIdx!
      ! Imagine, we have 3 boundary components with 8,6 and 4 edges. 
      ! Before summing the entries up, p_IboundaryCpIdx may have looked like this:
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   8   6   4
      !
      ! Summing that up gave the actual indices:
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   9  15  19
      !
      ! which is already the correct index array.
      ! Shifting p_IboundaryCpIdx by 1 will give:
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   1   9  15
      
      p_IboundaryCpIdx(2:rtriangulation%NBCT+2) = p_IboundaryCpIdx(1:rtriangulation%NBCT+1)
      
      ! Then, we again loop through all vertices and collect those on the
      ! boundary. In that loop, we use p_IboundaryCpIdx(2:NBCT+2) as pointer and
      ! increase them for every point we find. The loop will behave like
      ! the first one, i.e. we will again find 8,6 and 4 vertices on the boundary components
      ! 1,2 and 3. Increasing the p_IboundaryCpIdx(2:NBCT+1) for boundary components
      ! 1..NBCT the same way as done in the first loop will therefore again lead to
      !
      !         i            1   2   3   4
      ! p_IboundaryCpIdx(i)  1   9  15  19
      !
      ! Ok, let us catch the actual vertices.
      !      
      ! Check all vertices to find out, which vertices are on the boundary.
      !%OMP PARALLEL do PRIVATE(ibct,ivbd)
      do ivt=1,rtriangulation%NVT
        if (p_InodalProperty(ivt) .gt. 0) then
          ! id of the boundary component
          ibct = p_InodalProperty(ivt)
          
          ! set ivbd to the number of vertices on that boundary component
          ! thus ivbd holds the current number of vertices found for
          ! boundary component ibct and ivbd represents the current
          ! position in the p_IverticesAtBoundary array
          ivbd = p_IboundaryCpIdx(ibct+1)
          ! we have found a new point on that boundary component
          ! so increase the number of points by one
          p_IboundaryCpIdx(ibct+1) = ivbd+1
          
          ! Store the vertex as boundary vertex
          p_IverticesAtBoundary (ivbd) = ivt
        end if
      end do
      !%OMP END PARALLEL do
     
    end subroutine genRawBoundary2D

  end subroutine tria_attachCells

  ! ***************************************************************************************

!<subroutine>

  subroutine tria_cellGroupGreedy (rtriangulation, ngroups, Icells, IcellIndex, IelementGroup)
  
!<description>
  ! Divides all cells from a given triangulation structure into ngroups cell groups
  ! with all cells in a group being simply connected.
  !
  ! Uses a simple greedy algorithm to partition the cells into groups of similar
  ! size.
!</description>
  
!<input>
  ! Triangulation structure defining the cells to be decomposed into groups
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Number of groups, the cells should be divided into.
  integer, intent(in) :: ngroups
!</input>

!<output>
  ! Array with cell numbers, sorted for the groups. dimension(NEL)
  integer, dimension(:), intent(out)  :: Icells
  
  ! Array with indices for the groups. IcellIndex(i) defines the starting position
  ! of group i in Icells. dimension(1:ngroups+1)
  integer, dimension(:), intent(out) :: IcellIndex
  
  ! Array of dimension(NEL). For every element, the corresponding entry receives
  ! the ID of the group, that element belongs to.
  integer, dimension(:), intent(out) :: IelementGroup
!</output>

!</subroutine>

    ! local variables
    integer :: igroup, ibdelement, ifreeelement, iel, ineigh, nfreeelements, ive
    integer :: ngroupsact,ngroupsmax,i
    integer :: ielneigh
    
    !  Triangulation arrays
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    
    ! Element queue that saves non-processed element
    integer, dimension(:), allocatable :: IelementQueue
    integer :: iqptrRead, iqptrWrite
    
    ! Number of elements in a group
    integer, dimension(:), allocatable :: InelPerGroup,IgroupInc,IgroupMap
    
    ! Allocate memory and get some pointers.
    call storage_getbase_int2d(&
        rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
    call storage_getbase_int2d(&
        rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtBoundary,p_IelementsAtBoundary)
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtVertexIdx,p_IelementsAtVertexIdx)
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
    
    allocate(IelementQueue(max(rtriangulation%NEL,size(IcellIndex))))
    ! At the beginning, no element is assigned to a group.
    IelementGroup(:) = 0
    IcellIndex(:) = 0
    
    ibdElement = 1
    ifreeelement = 1
    nfreeelements = rtriangulation%NEL
    
    ! Loop through the queues we have to form
    !do igroup = 1,ngroups
    
    ! Loop as long as we have free elements
    igroup = 0
    do while (nfreeelements .gt. 0)
    
      ! New group
      igroup = igroup + 1
    
      ! Initialise the queue of non-processed elements
      iqptrRead = 1
      iqptrWrite = 1
    
      ! Put the first non-processed element into the queue.
      ! We try to find this on the boundary as long as we have elements
      ! there. If all boundary elements are processed and there are still
      ! groups left, we try to find new start elements in the inner.
      do while (ibdElement .le. size(p_IelementsAtBoundary))
        if (IelementGroup(p_IelementsAtBoundary(ibdElement)) .eq. 0) exit
        ibdElement = ibdElement + 1
      end do
      
      if (ibdElement .le. size(p_IelementsAtBoundary)) then
        ! Ok, a boundary element will be our next start element
        IelementQueue(iqptrWrite) = p_IelementsAtBoundary(ibdElement)
        iqptrWrite = iqptrWrite + 1
      else
        ! No element found. Try to find another free element in the domain
        ! From this point on, we use the variable ifreeelement to save
        ! which elements we already processed.
        do while (ifreeelement .le. rtriangulation%NEL)
          if (IelementGroup(ifreeelement) .le. 0) exit
          ifreeelement = ifreeelement + 1
        end do
        
        if (ifreeelement .le. rtriangulation%NEL) then
          IelementQueue(iqptrWrite) = ifreeelement
          iqptrWrite = iqptrWrite + 1
        else
          ! No free elements anymore; may happen if there are more
          ! groups than elements. Ok, then we can stop assigning groups here.
          exit
        end if
      end if
  
      ! Ok, we should now have at least one element in the queue.
      ! Mark that element as belonging to us and save its neighbours to the
      ! queue, so we can continue with marking these.
      !
      ! In IcellIndex(2:), we count how many elements we assigned to each group.
      ! The index is shifted by one to allow easier assignment later in the 
      ! collection phase.
      do
      
        ! Cancel if we have enough elements.
        ! The groups 1..ngroups can have at most nel/ngroups elements.
        if (igroup .le. ngroups) then
          if (IcellIndex(1+igroup) .ge. (rtriangulation%NEL+ngroups-1)/ngroups) exit
  
          ! We will get a new element in this group now...
          IcellIndex(1+igroup) = IcellIndex(1+igroup) + 1
          nfreeelements = nfreeelements - 1
        end if
      
        ! Get the element and put it to our current group
        iel = IelementQueue(iqptrRead)
        iqptrRead = iqptrRead + 1
        
        ! Assign a new group 
        IelementGroup (iel) = igroup
        
        ! Remember this group as the last one we found elements in.
        ngroupsmax = igroup
        
        ! Get all neighbours of the element (which are not assigned to
        ! a group) and put them to our queue.
        
        ! The following piece of code would produce some kind of 'diagonal'
        ! partitioning, where the partitions have diagonal form.
        
!        do ineigh = 1,UBOUND(p_IneighboursAtElement,1)
!          if (p_IneighboursAtElement(ineigh,iel) .ne. 0) then
!            ! Check if the element is already assigned.
!            if ((IelementGroup(p_IneighboursAtElement(ineigh,iel)) .gt. -igroup) .and. &
!                (IelementGroup(p_IneighboursAtElement(ineigh,iel)) .le. 0)) then
!              IelementQueue(iqptrWrite) = p_IneighboursAtElement(ineigh,iel)
!              iqptrWrite = iqptrWrite + 1
!              
!              ! Mark the element as to be processed in the current group.
!              ! We assign the negative group ID to the element.
!              ! The above IF statement ensures that the element is only once
!              ! processed in this group and then firstly processed again in
!              ! the next group.
!              IelementGroup(p_IneighboursAtElement(ineigh,iel)) = -igroup
!            end if
!          end if
!        end do

        ! The following piece of code would produce some kind of 'quadratic'
        ! partitioning, i.e. where the partitions have quadratic shape.

        ! Loop through all vertices on the element
        do ive = 1,ubound(p_IverticesAtElement,1)
          ! Cancel if this is a triangle in a quad mesh (or similar in 3D
          if (p_IverticesAtElement(ive,iel) .eq. 0) exit
          
          ! Loop through all neighbour elements adjacent to that vertex
          do ineigh = p_IelementsAtVertexIdx(p_IverticesAtElement(ive,iel)), &
                      p_IelementsAtVertexIdx(p_IverticesAtElement(ive,iel)+1)-1
                      
            ! Get the neighbour elements adjacent to that vertex.
            !
            ! Check if the element is already assigned.
            if ((IelementGroup(p_IelementsAtVertex(ineigh)) .gt. -igroup) .and. &
                (IelementGroup(p_IelementsAtVertex(ineigh)) .le. 0)) then
              IelementQueue(iqptrWrite) = p_IelementsAtVertex(ineigh)
              iqptrWrite = iqptrWrite + 1
              
              ! Mark the element as to be processed in the current group.
              ! We assign the negative group ID to the element.
              ! The above IF statement ensures that the element is only once
              ! processed in this group and then firstly processed again in
              ! the next group.
              IelementGroup(p_IelementsAtVertex(ineigh)) = -igroup
            end if
                      
          end do
          
        end do
        
        ! Proceed with the next element in the queue -- if we have any.
        if (iqptrRead .ge. iqptrWrite) exit
      end do
      
      ! The loop is left 
      ! a) if there are too many elements in the group or
      ! b) if no more elements are found.
      ! Continue to form a new group
      
    end do
    
    ! At this point, all elements are assigned to groups, although we may have more
    ! groups than we are allowed to have. More precisely, the number of groups
    ! we have is...
    ngroupsact = ngroupsmax
    
    if (ngroupsact .gt. ngroups) then
    
      ! Now we have a small adventure: more groups than we are allowed to have!
      !
      ! Note that all of these groups are not connected! So we have to reduce the
      ! number of groups by combining them. For that purpose, we at first search
      ! for 'small' groups and combine them with larger groups.
      !
      ! Loop through the groups and figure out how many elements each group has
      ! and which group is connected to which one.
      
      allocate(InelPerGroup(ngroupsact), IgroupInc(ngroupsact), IgroupMap(ngroupsact))
      IgroupInc(:) = 0
      InelPerGroup(:) = 0
      
      ! IgroupMap defines a group mapping for the later correction of group ID`s.
      do i=1,ngroupsact
        IgroupMap(i) = i
      end do
      
      do iel=1,rtriangulation%NEL
      
        InelPerGroup(IelementGroup(iel)) = InelPerGroup(IelementGroup(iel)) + 1
        
        ! Try to find neigbouring element groups with as few elements
        ! as possible.
        !
        ! Analyse the (edge/face-) neighbours of the element to find
        ! a neighbour group.
        do ineigh = 1,ubound(p_IneighboursAtElement,1)
        
          if (p_IneighboursAtElement(ineigh,iel) .ne. 0) then
          
            ielneigh = p_IneighboursAtElement(ineigh,iel)
            if (IelementGroup(ielneigh) .lt. IelementGroup(iel)) then
            
              ! Neighbour group with a smaller number. (Probably no group associated
              ! up to now.) Does that group have less elements than the 
              ! previous neighbour?
              if (IgroupInc(IelementGroup(iel)) .eq. 0) then
                IgroupInc(IelementGroup(iel)) = IelementGroup(ielneigh)
              else
                  if (InelPerGroup(IelementGroup(ielneigh)) .lt. &
                    InelPerGroup(IgroupInc(IelementGroup(iel)))) then
                  IgroupInc(IelementGroup(iel)) = IelementGroup(ielneigh)
                end if
              end if
              
            end if
          end if
          
        end do
        
      end do
      
      ! Now we know how many elements are in each group and (roughly)
      ! which group is connected to which one. Search for the element
      ! group with the smallest number of elements and combine it
      ! with the incident group until we have only ngroups groups.
      do while (ngroupsact .gt. ngroups)
      
        igroup = 1
        do i=2,ngroupsact
          if (InelPerGroup(i) .lt. InelPerGroup(igroup)) igroup = i
        end do
        
        ! Combine group igroup with group IgroupInc(igroup)
        InelPerGroup(IgroupInc(igroup)) = InelPerGroup(IgroupInc(igroup)) + InelPerGroup(igroup)
        
        ! Remember how this group is mapped...
        IgroupMap(igroup) = IgroupInc(igroup)
        
        ! Do not forget to correct incidental groups!
        do i=1,ngroupsact
          if (IgroupInc(i) .eq. igroup) &
            IgroupInc(i) = IgroupInc(igroup)
        end do
        
        ! Remove the group from InelPerGroup by setting it to a high value
        InelPerGroup(igroup) = rtriangulation%NEL+1
        
        ! Number of groups reduced.
        ngroupsact = ngroupsact - 1
      
      end do
      
      ! Ok, now the number of groups are ok, but the group ID`s in IelementGroup not!
      ! We have to compress them...
      ! Build in IcellIndex (2:) the number of elements per group.
      igroup = 0
      do i = 1,ngroupsmax
        if (InelPerGroup(i) .ne. rtriangulation%NEL+1) then
          igroup = igroup + 1
          IcellIndex(1+igroup) = InelPerGroup(i)
          
          ! Update the group map array for those groups that still have elements.
          IgroupMap(i) = igroup
        end if
      end do
      
      ! Now a final loop through the elements to correct the group ID`s.
      do iel=1,rtriangulation%NEL
        IelementGroup(iel) = IgroupMap(IelementGroup(iel))
      end do
      
      deallocate(InelPerGroup,IgroupInc,IgroupMap)
      
    end if

    ! Now calculate the start positions of all groups.
    IcellIndex(1) = 1
    do igroup = 2,ngroups+1
      IcellIndex(igroup) = IcellIndex(igroup) + IcellIndex(igroup-1)
    end do
  
    ! Abuse the IelementQueue array as group pointer, it is large enough.
    call lalg_copyVectorInt(IcellIndex,IelementQueue)
  
    ! Loop through the elements and collect them.
    do iel = 1,rtriangulation%NEL
      Icells(IelementQueue(IelementGroup(iel))) = iel
      IelementQueue(IelementGroup(iel)) = IelementQueue(IelementGroup(iel)) + 1
    end do
    
    ! Deallocate memory, finish.
    deallocate(IelementQueue)
  
  end subroutine tria_cellGroupGreedy

  ! ***************************************************************************

!<function>

  integer function tria_getNAE_direct (rtriangulation, iel)
  
!<description>
  ! This routine calculates the number of faces on element iel.
!</description>

!<input>
  ! Triangulation structure
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Number of the element whose NAE should be calculated
  integer, intent(in) :: iel
!</input>

!<result>
  ! Number of faces on element iel.
!</result>

!</function>

    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: i,nve
    
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)

    do i = ubound(p_IverticesAtElement,1), 1, -1
      if (p_IverticesAtElement(i,iel) .ne. 0) then
        nve = i
        exit
      end if
    end do
  
    select case(nve)
    case (TRIA_NVETET3D)
      tria_getNAE_direct = TRIA_NAETET3D

    case (TRIA_NVEPYR3D)
      tria_getNAE_direct = TRIA_NAEPYR3D
      
    case (TRIA_NVEPRIS3D)
      tria_getNAE_direct = TRIA_NAEPRIS3D

    case (TRIA_NVEHEXA3D)
      tria_getNAE_direct = TRIA_NAEHEXA3D

    case DEFAULT
      tria_getNAE_direct = 0
    end select

  end function tria_getNAE_direct

  ! ***************************************************************************

!<function>

  pure integer function tria_getNAE_indirect (IverticesAtElement, iel)
  
!<description>
  ! This routine calculates the number of faces on element iel.
!</description>

!<input>
  ! This is the IverticesAtElement array of a triangulation
  integer, dimension(:,:), intent(in) :: IverticesAtElement
  
  ! Number of the element whose NAE should be calculated
  integer, intent(in) :: iel
!</input>

!<result>
  ! Number of faces on element iel.
!</result>

!</function>

    integer :: i,nve

    do i = ubound(IverticesAtElement,1), 1, -1
      if (IverticesAtElement(i,iel) .ne. 0) then
        nve = i
        exit
      end if
    end do
    
    select case(nve)
    case (TRIA_NVETET3D)
      tria_getNAE_indirect = TRIA_NAETET3D

    case (TRIA_NVEPYR3D)
      tria_getNAE_indirect = TRIA_NAEPYR3D
      
    case (TRIA_NVEPRIS3D)
      tria_getNAE_indirect = TRIA_NAEPRIS3D

    case (TRIA_NVEHEXA3D)
      tria_getNAE_indirect = TRIA_NAEHEXA3D

    case DEFAULT
      tria_getNAE_indirect = 0
    end select
  
  end function tria_getNAE_indirect

  ! ***************************************************************************

!<function>

  integer function tria_getNVE_direct (rtriangulation, iel)
  
!<description>
  ! This routine calculates the number of vertices/edges on element iel.
!</description>

!<input>
  ! Triangulation structure
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Number of the element whose NVE should be calculated
  integer, intent(in) :: iel
!</input>

!<result>
  ! Number of vertices on element iel.
!</result>

!</function>

    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: i
    
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)

    do i = ubound(p_IverticesAtElement,1), 1, -1
      if (p_IverticesAtElement(i,iel) .ne. 0) then
        tria_getNVE_direct = i
        return
      end if
    end do
  
  end function tria_getNVE_direct

  ! ***************************************************************************

!<function>

  pure integer function tria_getNVE_indirect (IvertEdgAtElement,iel)
  
!<description>
  ! This routine calculates the number of vertices/edges on element iel.
!</description>

!<input>
  ! This may be either the IverticesAtElement or the IedgesAtElement
  ! array of a triangulation
  integer, dimension(:,:), intent(in) :: IvertEdgAtElement
  
  ! Number of the element whose NVE should be calculated
  integer, intent(in) :: iel
!</input>

!<result>
  ! Number of vertices on element iel.
!</result>

!</function>

    integer :: i

    do i=ubound(IvertEdgAtElement,1),1,-1
      if (IvertEdgAtElement(i,iel) .ne. 0) then
        tria_getNVE_indirect = i
        return
      end if
    end do
  
  end function tria_getNVE_indirect

  ! ***************************************************************************

!<subroutine>

  subroutine tria_getPointsOnEdge (rtriangulation, Dcoords, npointsPerEdge, &
                                   DparValue, rboundary)
  
!<description>
  ! This routine can be used to calculate the coordinates of a fixed number
  ! of regularly distributed (inner-edge) points per edge.
  !
  ! npointsPerEdge specifies how many points per edge should be calculated.
  ! E.g. if npointsPerEdge=3, the routine will place 3 points
  ! on each edge and calculate their coordinates. By default, the points will
  ! be distributed regularly, so in this example at the 'relative position'
  ! (0.25, 0.5, 0.75) of the edge. If the caller wants to specify a different
  ! position of the points manually, the DparValue array can be specified.
  !
  ! Another example: By calling this routine with npointsPerEdge=1, the 
  ! routine will calculate the midpoints of all edges.
!</description>

!<input>
  ! Triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Number of points per edge to generate
  integer, intent(in) :: npointsPerEdge
  
  ! OPTIONAL: Array with parameter values of the points on the edge.
  ! dimension(npointsPerEdge). DparValue is a value in the range [0,1]
  ! and specifies the 'relative position' or 'parameter value' of each 
  ! point on the edge. If not specified, tria_getPointsOnEdge assumes a regular
  ! distribution of points, defined by the number of points on the edge.
  ! E.g. if npointsPerEdge=2 and DparValue is not specified, 
  ! the routine assumes 3 inner points on the edge corresponding to 
  ! DparValue=/(0.333333,0.666666)/.
  ! If specified, the caller can specify the exact parameter values of
  ! the three points on the edge, e.g. DparValue=/(0.25,0.75)/.
  real(DP), dimension(:), intent(in), optional :: DparValue
  
  ! OPTIONAL: Definition of the domain.
  ! If specified, the routine will calculate the points on the boundary
  ! edges using the definition on the boundary. If not specified, the routine
  ! calculates the points on the edges based only on the triangulation.
  type(t_boundary), intent(in), optional :: rboundary
!</input>

!<output>
  ! Coordinates of the regularly distributed points on all the edges
  ! of the triangulation.
  !   dimension(space-dimension, npointsPerEdge * #edges)
  ! Here, Dcoords(:,1..npointsPerEdge) receives the coordinates of the 
  ! points on the first edge, Dcoords(:,npointsPerEdge+1:2*npointsPerEdge)
  ! the coordinates on edge number 2 etc.
  real(DP), dimension(:,:), intent(out) :: Dcoords
!</output>

!</subroutine>

    ! local variables      
    real(DP), dimension(npointsPerEdge) :: Dparameters
    real(DP), dimension(:,:), pointer :: p_Dcoords
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer :: iedge
    integer :: ipointpos,ipoint1,ipoint2
    integer :: idim,ipoint
    integer :: ibdc,ibdedge
    real(DP) :: dpar1,dpar2
    
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    
    ! If DparValue is specified, take that. Otherwise, create the parameter
    ! values of the (inner-edge) points manually.
    if (present(DparValue)) then
      if (size(DparValue) .lt. npointsPerEdge) then
        call output_line ('DparValue not large enough!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'tria_getPointsOnEdge')
        call sys_halt()
      end if
      Dparameters(:) = DparValue(1:npointsPerEdge)
    else
      ! Equal distributiuon of points
      do ipoint1=1,npointsPerEdge
        Dparameters(ipoint1) = real(ipoint1,DP)/real(npointsPerEdge+1,DP)
      end do
    end if
    
    ! Get the triangulation stuff
    if (rtriangulation%ndim .eq. 0) then
      call output_line ('Triangulation not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_getPointsOnEdge')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtEdge .eq. 0) then
      call output_line ('IverticesAtEdge not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_getPointsOnEdge')
      call sys_halt()
    end if
    
    call storage_getbase_double2d(rtriangulation%h_DvertexCoords,p_Dcoords)
    call storage_getbase_int2d(rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
  
    ! Loop through all edges
    do iedge = 1,rtriangulation%NMT
    
      ! Where do the set of points start in Dcoords?
      ipointpos = (iedge-1)*npointsPerEdge
      
      ! Endpoints of that edge?
      ipoint1 = p_IverticesAtEdge(1,iedge)
      ipoint2 = p_IverticesAtEdge(2,iedge)
    
      ! Calculate the points on the edge
      do ipoint = 1,npointsPerEdge
    
        ! Calculate the point coordinates
        do idim = 1,ubound(Dcoords,1)
          
          Dcoords(idim,ipoint+ipointpos) = &
            p_Dcoords(idim,ipoint1) * Dparameters(ipoint) + &
            p_Dcoords(idim,ipoint2) * (1.0_DP-Dparameters(ipoint))
        
        end do
      
      end do
      
    end do
    
    ! Is the boundary structure given?
    if (present(rboundary)) then
    
      ! 2D?
      if (rtriangulation%ndim .eq. NDIM2D) then
      
        ! Ok, we can calculate the correct position of the points on the boundary!
        ! Get the array with the boundary vertices and their parameter values.
        call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
            p_IboundaryCpIdx)
        call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
            p_IverticesAtBoundary)
        call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
            p_IedgesAtBoundary)
        call storage_getbase_double (rtriangulation%h_DvertexParameterValue,&
            p_DvertexParameterValue)
            
        ! Loop through the boundary components and the points in each component
        do ibdc = 1,rtriangulation%NBCT
        
          do ibdedge = p_IboundaryCpIdx(ibdc),p_IboundaryCpIdx(ibdc+1)-1
          
            ! Get the boundary edge
            iedge = p_IedgesAtBoundary(ibdedge)
            
            ! Get the parameter value of the points adjacent to that edge.
            dpar1 = p_DvertexParameterValue(ibdedge)
            if (ibdedge .ne. p_IboundaryCpIdx(ibdc+1)-1) then
              dpar2 = p_DvertexParameterValue(ibdedge+1)
            else
              ! Last edge ends with maximum parameter value on that boundary component
              dpar2 = boundary_dgetMaxParVal(rboundary,ibdc)
            end if
            
            ! Where do the set of points start in Dcoords?
            ipointpos = iedge*npointsPerEdge
            
            ! Calculate the points on the edge
            do ipoint = 1,npointsPerEdge
          
              ! Calculate the point coordinates
              call boundary_getCoords(rboundary, ibdc, &
                  dpar1 * Dparameters(ipoint) + dpar2 * (1.0_DP-Dparameters(ipoint)),  &
                  Dcoords(1,ipoint+ipointpos), Dcoords(2,ipoint+ipointpos))
                  
            end do
              
          end do
        
        end do
      
      end if
    
    end if
  
  end subroutine tria_getPointsOnEdge

  !************************************************************************

!<subroutine>

  elemental subroutine tria_getNeighbourVertex(ivertex,ivt1,ivt2,ineighbour)

!<description>
  ! Calculates the vertex number of the neighbour vertex of a vertex on an
  ! edge.
  !
  ! ivt1 and ivt2 are the vertex numbers of two vertices connected
  ! by an edge. ivertex is either ivt1 or ivt2.
  ! The result of this routine is the number of the neighbour vertex
  ! of ivertex, i.e. if ivertex=ivt1, ineighbour=ivt2 and if
  ! ivertex=ivt2, ineighbour=ivt1.
!</description>

!<input>
  ! Vertex number. Either ivt1 or ivt2.
  integer, intent(in) :: ivertex
  
  ! Vertex number of one vertex adjacent to an edge.
  integer, intent(in) :: ivt1

  ! Vertex number of the other vertex adjacent to that edge.
  integer, intent(in) :: ivt2
!</input>

!<output>
  ! Vertex number of the neighbour of ivertex.
  integer, intent(out) :: ineighbour
!</output>

!</subroutine>

    ! Note: Directly implementing this formula into the program
    ! code brings more speed :-)
    ! But to have a reference not to forget the formula, we have
    ! this routine...
    
    ineighbour = ivt1 + ivt2 - ivertex
    
  end subroutine tria_getNeighbourVertex

  !************************************************************************

!<subroutine>

  subroutine tria_getSubmeshNeighbourhood (rtriangulation, IsubmeshElements,&
                                           cneighbourhood, p_IsubmeshNeighbourhood,&
                                           p_IneighbourhoodType)

!<description>
  ! Creates a list of all elements which are directly adjacent to a set
  ! of elements. The routine returns all elements in one cell layer 
  ! around a given cell set.
!</description>

!<input>
  ! A given mesh.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A list of all elements from rtriangulation. The neighbourhood of
  ! this cell set is to be determined.
  integer, dimension(:), intent(in) :: IsubmeshElements
  
  ! A combination of TRIA_NEIGH_xxxx flags that specify the neighbourhood
  ! to be computed. If TRI_NEIGH_VERTEXNEIGHBOURS is set, all vertex neighbours
  ! of IsubmeshElements are found, TRI_NEIGH_EDGENEIGHBOURS will find all
  ! edge adjacent elements and TRI_NEIGH_FACENEIGHBOURS all face
  ! adjacent elements.
  integer(I32), intent(in) :: cneighbourhood
!</input>

!<output>
  ! A list of all elements in a cell layer around the given submesh.
  ! This is a pointer and will be allocated in this routine.
  integer, dimension(:), pointer :: p_IsubmeshNeighbourhood
  
  ! OPTIONAL: A list of neighbourhood specifiers. For every element
  ! in p_IsubmeshNeighbourhood, the corresponding value in this
  ! array specifies the type of neighbourhood of that element in relation
  ! to the subset. A value i in this array has the following meaning:
  ! i=1..NVT: The element is adjacent to the subset by vertex i.
  ! i=NVT+1..NVT+NMT: The element is adjacent to the subset by edge i-NVT.
  ! i=NVT+NMT+1..NVT+NMT+NAT: The element is adjacent to the subset by
  !                           face i-NVT-NMT.
  ! Elements may be adjacent by different types of neighbourhoods; the
  ! higher the connectivity of a neighbour element, the higher its
  ! priority and thus the higher the value i. That means: If an element
  ! is a face neighbour, i will be the face. If it is not a face neighbour
  ! but an edge neighbout, i will be the edge number. If it is neighter
  ! face nor edge but only vertex neighbour, i will be the vertex number.
  integer, dimension(:), pointer, optional :: p_IneighbourhoodType
!</output>

    ! local variables
    integer :: iel,nel,ielidx
    integer :: ive
    integer :: ivt
    integer, dimension(:), allocatable :: IelementFlag
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    
    ! To compute the neighbours, we simply mark them and collect them later.
    ! First allocate a flag array for all elements in the mesh.
    allocate(IelementFlag(rtriangulation%NEL))
    call lalg_clearVectorInt(IelementFlag,rtriangulation%NEL)
    
    ! Mark all elements in the submesh with a -1. This prevents these elements
    ! from being added to the element list below.
    do iel=1,size(IsubmeshElements)
      IelementFlag(IsubmeshElements(iel)) = -1
    end do
    
    ! Now loop through all elements. This is a dimension dependent loop,
    ! as the neighbourhood is dimension dependent.
    select case (rtriangulation%ndim)
    case (NDIM2D)
    
      ! Find all elements adjacent to the edges and mark them with a tag --
      ! as long as they are not already marked.
      call storage_getbase_int2d(rtriangulation%h_IneighboursAtElement,&
          p_IneighboursAtElement)
      call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,&
          p_IverticesAtElement)
      call storage_getbase_int2d(rtriangulation%h_IedgesAtElement,&
          p_IedgesAtElement)
      call storage_getbase_int(rtriangulation%h_IelementsAtVertex,&
          p_IelementsAtVertex)
      call storage_getbase_int(rtriangulation%h_IelementsAtVertexIdx,&
          p_IelementsAtVertexIdx)
      
      ! nel counts the elements in the vincinity of our element list.
      nel = 0
      
      ! Find edge-adjacent elements elements at first
      if (iand(cneighbourhood,TRI_NEIGH_VERTEXNEIGHBOURS) .ne. 0) then
      
        do iel=1,size(IsubmeshElements)
          ! Loop through all edges on this element
          do ive = 1,ubound(p_IneighboursAtElement,1)
            if (p_IneighboursAtElement(ive,IsubmeshElements(iel)) .ne. 0) then
              if (IelementFlag(p_IneighboursAtElement(ive,IsubmeshElements(iel))) .eq. 0) then
                IelementFlag(p_IneighboursAtElement(ive,IsubmeshElements(iel))) = &
                  p_IedgesAtElement(ive,IsubmeshElements(iel)) 
                nel = nel+1
              end if
            end if
          end do
        end do
        
      end if
      
      ! Additionally find vertex-adjacent elements
      if (iand(cneighbourhood,TRI_NEIGH_VERTEXNEIGHBOURS) .ne. 0) then
      
        do iel=1,size(IsubmeshElements)
          ! Loop through all all vertices on this element
          do ive = 1,ubound(p_IneighboursAtElement,1)
            ivt = p_IverticesAtElement(ive,IsubmeshElements(iel))
            do ielidx = p_IelementsAtVertexIdx(ivt),p_IelementsAtVertexIdx(ivt+1)-1
              if (IelementFlag(p_IelementsAtVertex(ielidx)) .eq. 0) then
                IelementFlag(p_IelementsAtVertex(ielidx)) = ivt
                nel = nel+1
              end if
            end do
          end do
        end do
        
      end if
      
      ! Allocate memory for p_IsubmeshNeighbourhood and collect the elements.
      allocate (p_IsubmeshNeighbourhood(nel))

      if (.not. present(p_IneighbourhoodType)) then
        nel = 0
        do iel = 1,rtriangulation%NEL
          if (IelementFlag(iel) .gt. 0) then
            nel = nel+1
            p_IsubmeshNeighbourhood(nel) = iel
          end if
        end do
      else
        ! Also collect the neighbourhood type.
        allocate (p_IneighbourhoodType(nel))
        nel = 0
        do iel = 1,rtriangulation%NEL
          if (IelementFlag(iel) .gt. 0) then
            nel = nel+1
            p_IsubmeshNeighbourhood(nel) = iel
            p_IneighbourhoodType(nel) = IelementFlag(iel)
          end if
        end do
      end if
      
      ! Deallocate memory, that is it.
      deallocate(IelementFlag)
    
    case DEFAULT
      call output_line ('tria_getSubmeshNeighbourhood: Dimension not implemented!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
      call sys_halt()
    end select

  end subroutine tria_getSubmeshNeighbourhood

  !************************************************************************

!<subroutine>

  subroutine tria_getElementsAtMacroEdge (rtriangulation, rmacroMesh,&
                                          ielcoarse, imtcoarse, p_Ineighbourhood)

!<description>
  ! Returns a list of all elements in the mesh rtriangulation that
  ! stem from element ielcoarse in the macro mesh and are adjacent
  ! to edge imtcoarse on the macro mesh.
  ! Note: The macro nodal property array must be present in rtriangulation
  ! to refer to the macro mesh!
!</description>

!<input>
  ! The mesh where to extract information from.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! The mesh defining the macro cells. rtriangulation must have been derived
  ! from this by refinement.
  type(t_triangulation), intent(in) :: rmacroMesh
  
  ! Number of the element on the coarse mesh which should contain
  ! all elements on the fine mesh adjacent to imtcoarse.
  integer, intent(in) :: ielcoarse
  
  ! Edge number (1..rmacroMesh%NMT) on the coarse mesh where to search for
  ! adjacent cells on the fine mesh.
  integer, intent(in) :: imtcoarse
!</input>

!<output>
  ! List of all elements adjacent to edge imtcoarse in element 
  ! ielcoarse on the macro mesh.
  ! This is a pointer and will be allocated in this routine.
  integer, dimension(:), pointer :: p_Ineighbourhood
!</output>

!</subroutine>

    ! local variables
    integer :: imt
    integer :: ivt,ielidx,ivtidx,irel
    integer :: nel,iel
    integer(I32), dimension(:), allocatable :: IelementTag
    integer, dimension(:), pointer :: p_ImacroNodalProperty
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:,:), pointer :: p_IelementsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    ! Primitive implementation. Mark all elements and collect them.
    !
    ! Allocate memory for element tags
    allocate(IelementTag(rtriangulation%NEL))
    call lalg_clearVectorint(IelementTag,rtriangulation%NEL)
    
    ! Get the macro nodal property array of the mesh
    call storage_getbase_int(rtriangulation%h_ImacroNodalProperty,&
        p_ImacroNodalProperty)
    call storage_getbase_int(rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)
    call storage_getbase_int(rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
    call storage_getbase_int2d(rtriangulation%h_IelementsAtEdge,&
        p_IelementsAtEdge)
    call storage_getbase_int2d(rtriangulation%h_IverticesAtEdge,&
        p_IverticesAtEdge)
    
    ! Now the dimension-dependent part...
    select case (rtriangulation%ndim)
    case (NDIM2D)
    
      ! NEL counts the number of found elements
      nel = 0
    
      ! Loop through the edges, find those that coincide with imtcoarse
      ! on the coarse mesh.
      irel = rmacroMesh%NVT+rmacroMesh%NMT
      
      do imt = 1,rtriangulation%NMT
      
        if (p_ImacroNodalProperty(rtriangulation%NVT+imt) .eq. rmacroMesh%NVT+imtcoarse) then
        
          ! Process the elements adjacent to the vertices of the edge.
          ! This will find all elements adjacent to the endpoints
          ! of the edge which includes the element adjacent to the edge itself.
          do ivtidx = 1,2
            ivt = p_IverticesAtEdge(ivtidx,imt)
            do ielidx = p_IelementsAtVertexIdx(ivt),p_IelementsAtVertexIdx(ivt+1)-1
              if (p_ImacroNodalProperty(irel+ p_IelementsAtVertex(ielidx)) .eq. &
                  irel+ielcoarse) then
                if (IelementTag(p_IelementsAtVertex(ielidx)) .eq. 0) then
                  IelementTag(p_IelementsAtVertex(ielidx)) = 1
                  nel = nel + 1
                end if
              end if            
            end do
          end do
        
        end if
      
      end do
      
      ! Now collect the elements on that edge.
      allocate(p_Ineighbourhood(nel))
      nel = 0
      do iel = 1,rtriangulation%NEL
        if (IelementTag(iel) .ne. 0) then
          nel = nel + 1
          p_Ineighbourhood(nel) = iel
        end if
      end do
      
      deallocate (IelementTag)

    case default
      call output_line ('tria_getElementsAtMacroEdge: Dimension not implemented!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
      call sys_halt()
    end select

  end subroutine tria_getElementsAtMacroEdge
  
  !************************************************************************

!<subroutine>

  subroutine tria_getElementsAtMacroVertex (rtriangulation, rmacroMesh,&
                                            ielcoarse, ivt, Ineighbourhood,&
                                            nelements)

!<description>
  ! Returns a list of all elements in the mesh rtriangulation that
  ! stem from element ielcoarse in the macro mesh and are adjacent
  ! to vertex ivt on the fine mesh.
  ! Note: The macro nodal property array must be present in rtriangulation
  ! to refer to the macro mesh!
!</description>

!<input>
  ! The mesh where to extract information from.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! The mesh defining the macro cells. rtriangulation must have been derived
  ! from this by refinement.
  type(t_triangulation), intent(in) :: rmacroMesh

  ! Number of the element on the coarse mesh which should contain
  ! all elements on the fine mesh adjacent to imtcoarse.
  integer, intent(in) :: ielcoarse
  
  ! Vertex number (1..rtriangulation%NVT) on the fine mesh where to search for
  ! adjacent cells on the fine mesh.
  integer, intent(in) :: ivt
!</input>

!<output>
  ! List of all elements adjacent to vertex ivtcoarse in element 
  ! ielcoarse on the macro mesh.
  ! The buffer must be large enough; a safe size is the maximum number
  ! of elements adjacent to a vertex (NNelAtVertex).
  integer, dimension(:), intent(out) :: Ineighbourhood
  
  ! Number of elements found and written to Ineighbourhood.
  integer, intent(out) :: nelements
!</output>

!</subroutine>

    ! local variables
    integer :: irel,irelcoarse
    integer :: ielidx
    integer, dimension(:), pointer :: p_ImacroNodalProperty
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    
    ! Get the macro nodal property array of the mesh
    call storage_getbase_int(rtriangulation%h_ImacroNodalProperty,&
        p_ImacroNodalProperty)
    call storage_getbase_int(rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)
    call storage_getbase_int(rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)

    ! NEL counts the number of found elements
    nelements = 0
    
    irel = rtriangulation%NVT+rtriangulation%NMT+rtriangulation%NAT
    irelcoarse = rmacroMesh%NVT+rmacroMesh%NMT+rmacroMesh%NAT
    
    ! Loop through all elements adjacent to ivt and collect them.
    do ielidx = p_IelementsAtVertexIdx(ivt),p_IelementsAtVertexIdx(ivt+1)-1
      if (p_ImacroNodalProperty(irel+p_IelementsAtVertex(ielidx)) .eq. &
          irelcoarse+ielcoarse) then
        nelements = nelements + 1
        Ineighbourhood(nelements) = p_IelementsAtVertex(ielidx)
      end if            
    end do
    
  end subroutine tria_getElementsAtMacroVertex

  ! ***************************************************************************

!<subroutine>  

  subroutine tria_genTwistIndex(rtriangulation)

!<description>
  ! Generates the twist indices for the elements.
!</description>
  
!<inputoutput>
    ! The triangulation structure to be updated.
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    integer :: isize

    ! 1D does not have twist indices.
    if (rtriangulation%ndim .eq. NDIM1D) return
    
    ! Allocate memory if necessary
    if (rtriangulation%h_ItwistIndex .eq. ST_NOHANDLE) then
      call storage_new ('tria_genTwistIndex', 'ItwistIndex', &
          rtriangulation%NEL, ST_INT32, &
          rtriangulation%h_ItwistIndex, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_ItwistIndex, isize)
      if (isize .ne. rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genTwistIndex', &
            rtriangulation%NEL, rtriangulation%h_ItwistIndex, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Now we can start to generate the twist index.
    ! Call the appropriate generation routine.
    select case(rtriangulation%ndim)
    case (NDIM2D)
      call genTwistIndex2D(rtriangulation)
    
    case (NDIM3D)
      call genTwistIndex3D(rtriangulation)
      
    end select


  contains
  
    ! ---------------------------------------------------------------
    
    subroutine genTwistIndex2D(rtriangulation)
    type(t_triangulation), intent(inout) :: rtriangulation
    
    ! local variables
    integer, dimension(:,:), pointer :: p_IedgesAtElement, &
                     p_IverticesAtElement, p_IverticesAtEdge
    integer(I32), dimension(:), pointer :: p_ItwistIndex
    integer :: iel, imt, iedge, ivt
    integer(I32) :: itwist
    
      ! Get all arrays from the storage
      call storage_getbase_int2D(rtriangulation%h_IedgesAtElement,&
                                 p_IedgesAtElement)
      call storage_getbase_int2D(rtriangulation%h_IverticesAtElement,&
                                 p_IverticesAtElement)
      call storage_getbase_int2D(rtriangulation%h_IverticesAtEdge,&
                                 p_IverticesAtEdge)
      call storage_getbase_int32(rtriangulation%h_ItwistIndex,&
                                 p_ItwistIndex)
                                 
      ! Okay, let us loop over all elements
      do iel = 1, rtriangulation%NEL
        
        ! Format the twist index entry
        itwist = 0
        
        ! Loop over all edges of the element
        do imt = 1, ubound(p_IedgesAtElement,1)
        
          ! Get the index of the edge
          iedge = p_IedgesAtElement(imt,iel)
          
          ! Jump out of the inner loop if iedge is zero - this can happen if
          ! the triangulation is a mixed triangle/quadrilateral mesh!
          if(iedge .eq. 0) exit
          
          ! Get the index of the first vertice on the edge
          ivt = p_IverticesAtEdge(1,iedge)
          
          ! If the index of the vertice is equal to the index of the vertice
          ! at which our local edge starts, then the corresponding twist bit
          ! is 0, otherwise it is 1.
          if(ivt .ne. p_IverticesAtelement(imt,iel)) then
            itwist = ior(itwist, int(ishft(1,imt-1),I32))
          end if
        
        end do ! imt
        
        ! Store the twist index entry
        p_ItwistIndex(iel) = itwist
        
      end do ! iel
      
      ! That is it
    
    end subroutine ! genTwistIndex2D
    
    ! ---------------------------------------------------------------
    
    subroutine genTwistIndex3D(rtriangulation)
    type(t_triangulation), intent(inout) :: rtriangulation
    
    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement, &
        p_IedgesAtElement, p_IfacesAtElement, p_IverticesAtEdge, &
        p_IverticesAtFace, p_IedgesAtFace
    integer(I32), dimension(:), pointer :: p_ItwistIndex
    integer :: iel, imt, iat, ivt, iedge, iface, nve, maxnve, maxnva
    integer(I32) :: itwist
    
    ! Tetrahedron
    integer, dimension(TRIA_NNETET3D), parameter :: Itev = (/1,2,3,1,2,3/)
    integer, dimension(TRIA_NAETET3D), parameter :: Itfv = (/1,1,2,3/)
    integer, dimension(TRIA_NAETET3D), parameter :: Itfn = (/2,4,4,4/)
    integer, dimension(TRIA_NAETET3D), parameter :: Itav = (/3,3,3,3/)

    ! Pyramid
    integer, dimension(TRIA_NNEPYR3D), parameter :: Iyev = (/1,2,3,4,1,2,3,4/)
    integer, dimension(TRIA_NAEPYR3D), parameter :: Iyfv = (/1,1,2,3,4/)
    integer, dimension(TRIA_NAEPYR3D), parameter :: Iyfn = (/2,4,5,6,6/)
    integer, dimension(TRIA_NAEPYR3D), parameter :: Iyav = (/4,3,3,3,3/)
    
    ! Prism
    integer, dimension(TRIA_NNEPRIS3D), parameter :: Irev = (/1,2,3,1,2,3,4,5,6/)
    integer, dimension(TRIA_NAEPRIS3D), parameter :: Irfv = (/1,1,2,3,4/)
    integer, dimension(TRIA_NAEPRIS3D), parameter :: Irfn = (/2,5,5,5,5/)
    integer, dimension(TRIA_NAEPRIS3D), parameter :: Irav = (/3,4,4,4,3/)
    
    ! Hexahedron
    integer, dimension(TRIA_NNEHEXA3D), parameter :: Ihev = (/1,2,3,4,1,2,3,4,5,6,7,8/)
    !integer, dimension(TRIA_NAEHEXA3D), parameter :: Ihfv = (/1,1,2,3,4,5/)
    !integer, dimension(TRIA_NAEHEXA3D), parameter :: Ihfn = (/2,5,6,7,8,8/)
    integer, dimension(TRIA_NAEHEXA3D), parameter :: Ihfv = (/1,1,7,7,1,7/)
    integer, dimension(TRIA_NAEHEXA3D), parameter :: Ihfn = (/2,5,3,8,4,6/)
    integer, dimension(TRIA_NAEHEXA3D), parameter :: Ihav = (/4,4,4,4,4,4/)
    

      ! Get all arrays from the storage
      call storage_getbase_int2D(rtriangulation%h_IverticesAtElement,&
                                 p_IverticesAtElement)
      call storage_getbase_int2D(rtriangulation%h_IedgesAtElement,&
                                 p_IedgesAtElement)
      call storage_getbase_int2D(rtriangulation%h_IfacesAtElement,&
                                 p_IfacesAtElement)
      call storage_getbase_int2D(rtriangulation%h_IverticesAtEdge,&
                                 p_IverticesAtEdge)
      call storage_getbase_int2D(rtriangulation%h_IverticesAtFace,&
                                 p_IverticesAtFace)
      call storage_getbase_int2D(rtriangulation%h_IedgesAtFace,&
                                 p_IedgesAtFace)
      call storage_getbase_int32(rtriangulation%h_ItwistIndex,&
                                 p_ItwistIndex)
      
      ! Get the maximum number of vertices adjacent to ...
      maxnve = ubound(p_IverticesAtElement,1)
      maxnva = ubound(p_IverticesAtFace,1)

      ! Okay, let us loop over all elements
      do iel = 1, rtriangulation%NEL
        
        ! Format the twist index entry
        itwist = 0
        
        ! Determine the number of vertices adjacent to this element
        do nve = 1, maxnve
          if(p_IverticesAtElement(nve,iel) .le. 0) exit
        end do
        nve = nve - 1
        
        ! Okay, what type of element do we have here?
        select case(nve)
        case(TRIA_NVETET3D)
          ! It is a tetrahedron
          
          ! Calculate edge-twist-indices
          do imt = 1, TRIA_NNETET3D
            if(p_IverticesAtEdge(1,p_IedgesAtElement(imt,iel)) .ne. &
               p_IverticesAtElement(Itev(imt),iel)) then
               itwist = ior(itwist,int(ishft(1,imt-1),I32))
             end if
          end do ! imt
          
          ! Calculate face-twist indices
          do iat = 1, TRIA_NAETET3D
          
            ! Get face index
            iface = p_IfacesAtElement(iat,iel)
            
            ! Calculate shift
            do ivt = 1, Itav(iat)
              if(p_IverticesAtFace(ivt,iface) .eq. &
                p_IverticesAtElement(Itfv(iat),iel)) exit
            end do ! ivt
            
            ! Calculate orientation
            if(p_IverticesAtFace(mod(ivt,Itav(iat))+1,iface) .ne. &
              p_IverticesAtElement(Itfn(iat),iel)) then
              itwist = ior(itwist, int(ishft(ivt+3,9+3*iat),I32))
            else
              itwist = ior(itwist, int(ishft(ivt-1,9+3*iat),I32))
            end if
          end do ! iat
        
        case(TRIA_NVEPYR3D)
          ! It is a pyramid
          
          ! Calculate edge-twist-indices
          do imt = 1, TRIA_NNEPYR3D
            if(p_IverticesAtEdge(1,p_IedgesAtElement(imt,iel)) .ne. &
               p_IverticesAtElement(Iyev(imt),iel)) then
               itwist = ior(itwist,int(ishft(1,imt-1),I32))
             end if
          end do ! imt

          ! Calculate face-twist indices
          do iat = 1, TRIA_NAEPYR3D
          
            ! Get face index
            iface = p_IfacesAtElement(iat,iel)
            
            ! Calculate shift
            do ivt = 1, Iyav(iat)
              if(p_IverticesAtFace(ivt,iface) .eq. &
                p_IverticesAtElement(Iyfv(iat),iel)) exit
            end do ! ivt
            
            ! Calculate orientation
            if(p_IverticesAtFace(mod(ivt,Iyav(iat))+1,iface) .ne. &
              p_IverticesAtElement(Iyfn(iat),iel)) then
              itwist = ior(itwist, int(ishft(ivt+3,9+3*iat),I32))
            else
              itwist = ior(itwist, int(ishft(ivt-1,9+3*iat),I32))
            end if
          end do ! iat

        case(TRIA_NVEPRIS3D)
          ! It is a prism
        
          ! Calculate edge-twist-indices
          do imt = 1, TRIA_NNEPRIS3D
            if(p_IverticesAtEdge(1,p_IedgesAtElement(imt,iel)) .ne. &
               p_IverticesAtElement(Irev(imt),iel)) then
               itwist = ior(itwist,int(ishft(1,imt-1),I32))
             end if
          end do ! imt
          
          ! Calculate face-twist indices
          do iat = 1, TRIA_NAEPRIS3D
          
            ! Get face index
            iface = p_IfacesAtElement(iat,iel)
            
            ! Calculate shift
            do ivt = 1, Irav(iat)
              if(p_IverticesAtFace(ivt,iface) .eq. &
                p_IverticesAtElement(Irfv(iat),iel)) exit
            end do ! ivt
            
            ! Calculate orientation
            if(p_IverticesAtFace(mod(ivt,Irav(iat))+1,iface) .ne. &
              p_IverticesAtElement(Irfn(iat),iel)) then
              itwist = ior(itwist, int(ishft(ivt+3,9+3*iat),I32))
            else
              itwist = ior(itwist, int(ishft(ivt-1,9+3*iat),I32))
            end if
          end do ! iat

        case(TRIA_NVEHEXA3D)
          ! It is a hexahedron
        
          ! Calculate edge-twist-indices
          do imt = 1, TRIA_NNEHEXA3D
            if(p_IverticesAtEdge(1,p_IedgesAtElement(imt,iel)) .ne. &
               p_IverticesAtElement(Ihev(imt),iel)) then
               itwist = ior(itwist,int(ishft(1,imt-1),I32))
             end if
          end do ! imt
          
          ! Calculate face-twist indices
          do iat = 1, TRIA_NAEHEXA3D
          
            ! Get face index
            iface = p_IfacesAtElement(iat,iel)
            
            ! Calculate shift
            do ivt = 1, Ihav(iat)
              if(p_IverticesAtFace(ivt,iface) .eq. &
                p_IverticesAtElement(Ihfv(iat),iel)) exit
            end do ! ivt
            
            ! Calculate orientation
            if(p_IverticesAtFace(mod(ivt,Ihav(iat))+1,iface) .ne. &
              p_IverticesAtElement(Ihfn(iat),iel)) then
              itwist = ior(itwist, int(ishft(ivt+3,9+3*iat),I32))
            else
              itwist = ior(itwist, int(ishft(ivt-1,9+3*iat),I32))
            end if
          end do ! iat

        case default
          ! Unknown cell type...
          itwist = 0
        
        end select
        
        ! Store the twist index entry
        p_ItwistIndex(iel) = itwist
      
      end do ! iel
      
      ! That is it
    
    end subroutine genTwistIndex3D

!    ! 'Old implementation'
!    subroutine genTwistIndexEdges2D(rtriangulation)
!    
!    ! Generate the twist index for 2D meshes
!    
!    ! The triangulation where the twist index should be generated.
!    type(t_triangulation), intent(inout) :: rtriangulation
!    
!      ! local variables
!      integer(I32), dimension(:,:), pointer :: p_IneighboursAtElement
!      integer :: iel,ielneighbour
!      integer(I32), dimension(:), pointer :: p_ItwistIndex
!      integer :: iedge
!      integer :: imt
!      integer(I32) :: itwistindex
!    
!      if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) then
!        call output_line ('IverticesAtEdge not available!', &
!                          OU_CLASS_ERROR,OU_MODE_STD,'genTwistIndexEdges2D')
!        call sys_halt()
!      end if
!    
!      if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
!        call output_line ('IedgesAtElement not available!', &
!                          OU_CLASS_ERROR,OU_MODE_STD,'genTwistIndexEdges2D')
!        call sys_halt()
!      end if
!      
!      call storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,&
!          p_IneighboursAtElement)
!      call storage_getbase_int (rtriangulation%h_ItwistIndexEdges,&
!          p_ItwistIndex)
!          
!      ! Chech the edge on each each element.
!      ! If the first vertex has the higher number, the twist index is
!      ! one. In that case, set the appropriate bit in the twist index field.
!      do iel=1,rtriangulation%NEL
!        itwistindex = 0
!        do imt=1,ubound(p_IneighboursAtElement,1)
!          ielneighbour = p_IneighboursAtElement(imt,iel)
!          if ((iel .lt. ielneighbour) .or. (ielneighbour .eq. 0)) then
!            itwistindex = ior(itwistindex,ishft(1,imt-1))
!          end if
!        end do
!        p_ItwistIndex(iel) = itwistindex
!      end do
!          
!    end subroutine
    
!    ! ---------------------------------------------------------------
!
!    subroutine genTwistIndexFaces3D(rtriangulation)
!    
!    ! Generate the twist index for 3D meshes
!
!    ! The triangulation where the twist index should be generated.
!    type(t_triangulation), intent(inout) :: rtriangulation
!    
!    ! Implementation must be done again.
!    ! The calculation works, but the target array is not properly installed
!    ! in the triangulation structure!
!    
!      ! local variables
!      integer(I32), dimension(:,:), pointer :: p_IverticesAtEdge
!      integer(I32), dimension(:,:), pointer :: p_IfacesAtElement
!      integer(I32), dimension(:,:), pointer :: p_IverticesAtFace
!      integer(I32), dimension(:,:), pointer :: p_IneighboursAtElement
!      integer(I32), dimension(:,:), pointer :: p_IverticesAtElement
!      integer(I32), dimension(:,:), pointer :: p_ItwistIndex
!      integer :: iedge,iface
!      integer :: ivt,iidx
!      integer :: itwistsgn,imyface,ineighbourface,nva
!      integer :: iel,ielneighbour
!      integer, dimension(TRIA_MAXNVE2D) :: IverticesAtFace1
!      integer, dimension(TRIA_MAXNVE2D) :: IverticesAtFace2
!    
!      if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) then
!        call output_line ('IverticesAtEdge not available!', &
!                          OU_CLASS_ERROR,OU_MODE_STD,'genTwistIndexFaces3D')
!        call sys_halt()
!      end if
!
!      if (rtriangulation%h_IverticesAtFace .eq. ST_NOHANDLE) then
!        call output_line ('IverticesAtFace not available!', &
!                          OU_CLASS_ERROR,OU_MODE_STD,'genTwistIndexFaces3D')
!        call sys_halt()
!      end if
!    
!      if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
!        call output_line ('IneighboursAtElement not available!', &
!                          OU_CLASS_ERROR,OU_MODE_STD,'genTwistIndexFaces3D')
!        call sys_halt()
!      end if
!    
!      if (rtriangulation%h_IfacesAtElement .eq. ST_NOHANDLE) then
!        call output_line ('IfacesAtElement not available!', &
!                          OU_CLASS_ERROR,OU_MODE_STD,'genTwistIndexFaces3D')
!        call sys_halt()
!      end if
!
!      if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
!        call output_line ('IverticesAtElement not available!', &
!                          OU_CLASS_ERROR,OU_MODE_STD,'genTwistIndexFaces3D')
!        call sys_halt()
!      end if
!      
!      call storage_getbase_int2d (rtriangulation%h_IverticesAtEdge,&
!          p_IverticesAtEdge)
!      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
!          p_IverticesAtElement)
!      call storage_getbase_int2d (rtriangulation%h_IverticesAtFace,&
!          p_IverticesAtFace)
!      call storage_getbase_int2d (rtriangulation%h_IfacesAtElement,&
!          p_IfacesAtElement)
!      call storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,&
!          p_IneighboursAtElement)
!      call storage_getbase_int2d (rtriangulation%h_ItwistIndexFaces,&
!          p_ItwistIndex)
!          
!      nva = rtriangulation%NNVA
!          
!      ! Generate the twist index for all elements.
!      do iel=1,rtriangulation%NEL
!      
!        ! Loop over the faces
!        do imyface = 1,ubound(p_IfacesAtElement,1)
!      
!          ! Get the neighbour, face, vertex with local vertex number 1
!          ielneighbour = p_IneighboursAtElement(imyface,iel)
!          
!          ! Figure out which vertex on the face of the neighbour element
!          ! corresponds to vertex 1 of the current face.
!          ! If there is no neighbour, we set the index to 1, that is it.
!          if (ielneighbour .eq. 0) then
!            
!            p_ItwistIndex(imyface,iel) = 1
!            
!          else
!
!            ! We only have to do something if there is no neighbour with
!            ! a larger element number.
!            if (iel .lt. ielneighbour) then
!
!              iface = p_IfacesAtElement (imyface,iel)
!
!              ! There is a neighbour. which face is the current one?
!              nsearch: do ineighbourface = 1,ubound(p_IfacesAtElement,1)
!                
!                if (p_IfacesAtElement(ineighbourface,ielneighbour) .eq. iface) then
!                
!                  ! Get the vertex numbers on the current face...
!                  ! From the one element as well as from the other.
!                  call tria_getVerticesAtFaceDirect(&
!                      imyface,nva,p_IverticesAtElement(:,iel),.true.,&
!                      IverticesAtFace1)
!                  call tria_getVerticesAtFaceDirect(&
!                      ineighbourface,nva,p_IverticesAtElement(:,ielneighbour),.true.,&
!                      IverticesAtFace2)
!
!                  ! First vertex on the face
!                  ivt = IverticesAtFace1 (1)
!                
!                  ! Which vertex on the neighbour corresponds to vertex 1
!                  ! of the current element?
!                  do iidx = 1,ubound(p_IverticesAtFace,1)
!                  
!                    if (IverticesAtFace2(iidx) .eq. ivt) then
!                    
!                      ! This value is the twist index. Save it for our 
!                      ! element as well as for the neighbour (with inverse sign
!                      ! since the orientation is changed on the other
!                      ! element, although the number is the same).
!                      p_ItwistIndex(imyface,iel) = iidx
!                      p_ItwistIndex(ineighbourface,ielneighbour) = -iidx
!                      
!                      ! Next face
!                      exit nsearch
!                      
!                    end if
!                  
!                  end do
!                
!                end if
!                
!              end do nsearch
!            
!            end if
!            
!          end if
!          
!        end do
!        
!      end do
!    
!    end subroutine genTwistIndexFaces3D

  end subroutine tria_genTwistIndex

  !************************************************************************

!<subroutine>
  
  subroutine tria_propMacroNodalProperty2lv (rtriaCoarse,rtriangulation)
  
!<description>
  ! Propagates a macro nodal property array from a coarse mesh to a fine mesh
  ! according to the uniform 2-level refinement strategy.
!</description>

!<input>
  ! The triangulation structure of the coarse mesh.
  type(t_triangulation), intent(in) :: rtriaCoarse
!</input>

!<inputoutput>
  ! The triangulation structure of the fine mesh which is to recevie
  ! a macro nodal property array
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: isize
    integer, dimension(:), pointer :: p_ImacroNodPropSource
    integer, dimension(:), pointer :: p_ImacroNodPropDest

    ! Allocate memory for that array if it does not exist.
    if (rtriangulation%h_ImacroNodalProperty .eq. ST_NOHANDLE) then
      call storage_new ('propMacroNodalProperty2lv', 'KMCPR', &
          rtriangulation%NVT+rtriangulation%NMT+&
          rtriangulation%NAT+rtriangulation%NEL, &
          ST_INT, rtriangulation%h_ImacroNodalProperty, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_ImacroNodalProperty, isize)
      if (isize .lt. rtriangulation%NVT+rtriangulation%NMT+rtriangulation%NAT) then
        ! If the size is wrong, reallocate memory.
        ! Copy the old content as we must not destroy the old nodal property
        ! tags of the vertices.
        call storage_realloc ('propMacroNodalProperty2lv', &
            rtriangulation%NVT+rtriangulation%NMT+&
            rtriangulation%NAT+rtriangulation%NEL, &
            rtriangulation%h_ImacroNodalProperty, &
            ST_NEWBLOCK_NOINIT, .true.)
      end if
    end if
    
    ! Get the macro nodal property array for the coarse and fine grid
    call storage_getbase_int(rtriaCoarse%h_ImacroNodalProperty,&
        p_ImacroNodPropSource)
    call storage_getbase_int(rtriangulation%h_ImacroNodalProperty,&
        p_ImacroNodPropDest)

    ! Propagate the macro nodal property as refinement tag to the
    ! fine mesh.
    call tria_genRefTags2lv(rtriaCoarse,rtriangulation,&
      p_ImacroNodPropSource,p_ImacroNodPropDest)

  end subroutine tria_propMacroNodalProperty2lv

!************************************************************************

!<subroutine>

  subroutine tria_genRefTags2lv(rtriaCoarse,rtriaFine,&
      IrefTagsCoarse,IrefTagsFine)

!<description>
  ! Calculates refinement tags for a mesh that was constructed by
  ! 2-level ordering.
  !
  ! Refinement tags represent a user-defined possibility to automaticall 
  ! classify elements, faces, edges and vertices during the refinement
  ! process. Using this technique allows to transfer coarse grid information
  ! to finegrid information.
  !
  ! Before the refinement starts, the caller may assign each a user defined 
  ! number to each vertex, edge, face and element; this number is called
  ! 'refinement tag'. Upon refinement, the refinement tag is inherited by
  ! new geometric elements in the following sense:
  ! -> Vertices, edges, faces and elements that appear in the inner of a 
  !    coarse grid cell receive the refinement tag of the coarse grid cell.
  ! -> Vertices, edges and faces that appear in the inner of a coarse grid 
  !    face receive the refinement tag of the coarse grid face.
  ! -> Vertices and edges that appear in the inner of a coarse grid edge
  !    receive the refinement tag of the coarse grid edge.
  ! -> Vertices that coincide with coarse grid vertices receive the
  !    refinement tag of the coarse grid vertices.
!</description>

!<input>
  ! Coarse mesh.
  ! Shall be a standard mesh.
  type(t_triangulation), intent(in) :: rtriaCoarse
  
  ! Fine mesh that was constructed from rtriaCoarse by 2-level refinement.
  ! Shall be a standard mesh.
  type(t_triangulation), intent(in) :: rtriaFine
  
  ! Refinement tags for the coarse mesh that should be propagated to
  ! the fine mesh. The size and shape of the array is depending 
  ! on the dimension.
  ! 1D: dimension(1:NVT+NEL).
  !     IrefTagsCoarse(1:NVT) = ref. tags for vertices.
  !     IrefTagsCoarse(NVT+1:NVT+NEL) = ref. tags for elements.
  ! 2D: dimension(1:NVT+NMT+NEL)
  !     IrefTagsCoarse(1:NVT) = ref. tags for vertices.
  !     IrefTagsCoarse(NVT+1:NVT+NMT) = ref. tags for edges.
  !     IrefTagsCoarse(NVT+NMT+1:NVT+NMT+NEL) = ref. tags for elements.
  ! 3D: dimension(1:NVT+NMT+NAT+NEL)
  !     IrefTagsCoarse(1:NVT) = ref. tags for vertices.
  !     IrefTagsCoarse(NVT+1:NVT+NMT) = ref. tags for edges.
  !     IrefTagsCoarse(NVT+NMT+1:NVT+NMT+NAT) = ref. tags for faces.
  !     IrefTagsCoarse(NVT+NMT+NAT+1:NVT+NMT+NAT+NEL) = ref. tags for elements.
  integer, dimension(:), intent(in) :: IrefTagsCoarse
!</input>

!</output>
  ! Refinement tags for the fine mesh calculated from coarse mesh 
  ! information. The size and shape of the array is depending 
  ! on the dimension.
  ! 1D: dimension(1:NVT+NEL).
  !     IrefTagsFine(1:NVT) = ref. tags for vertices.
  !     IrefTagsFine(NVT+1:NVT+NEL) = ref. tags for elements.
  ! 2D: dimension(1:NVT+NMT+NEL)
  !     IrefTagsFine(1:NVT) = ref. tags for vertices.
  !     IrefTagsFine(NVT+1:NVT+NMT) = ref. tags for edges.
  !     IrefTagsFine(NVT+NMT+1:NVT+NMT+NEL) = ref. tags for elements.
  ! 3D: dimension(1:NVT+NMT+NAT+NEL)
  !     IrefTagsFine(1:NVT) = ref. tags for vertices.
  !     IrefTagsFine(NVT+1:NVT+NMT) = ref. tags for edges.
  !     IrefTagsFine(NVT+NMT+1:NVT+NMT+NAT) = ref. tags for faces.
  !     IrefTagsFine(NVT+NMT+NAT+1:NVT+NMT+NAT+NEL) = ref. tags for elements.
  integer, dimension(:), intent(out) :: IrefTagsFine
!</output>

    ! Select the correct calculation routine depending on the dimension.
    select case (rtriaCoarse%ndim)
    case (NDIM1D)
      call output_line ('Refinement tags not implemented for 1D.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genRefTags2lv')
      call sys_halt()
    case (NDIM2D)
      call calcRefTags2D (rtriaCoarse,rtriaFine,IrefTagsCoarse,IrefTagsFine)
    case (NDIM3D)
      call output_line ('Refinement tags not implemented for 3D.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genRefTags2lv')
      call sys_halt()
    end select

  contains

    ! ---------------------------------------------------------------

    subroutine calcRefTags2D(rtriaCoarse, rtriaFine,&
                             IrefTagsCoarse, IrefTagsFine)

    type(t_triangulation), intent(in) :: rtriaCoarse
    type(t_triangulation), intent(in) :: rtriaFine
    integer, dimension(:), intent(in) :: IrefTagsCoarse
    integer, dimension(:), intent(inout) :: IrefTagsFine
    
      ! local variables
      integer :: i
      integer :: nvt,nvtfine
      integer :: nmt,nmtfine
      integer :: nve
      integer :: nel,iel,ielidx
      integer, dimension(5) :: IrefTag
      integer, dimension(4) :: iellocal
      integer, dimension(:,:), pointer :: p_IedgesAtElementC
      integer, dimension(:,:), pointer :: p_IedgesAtElementF
      integer, dimension(:,:), pointer :: p_IverticesAtElement
      
      nvt = rtriaCoarse%NVT
      nmt = rtriaCoarse%NMT
      nel = rtriaCoarse%NEL

      nvtfine = rtriaFine%NVT
      nmtfine = rtriaFine%NMT
      
      ! Fetch some arrays
      call storage_getbase_int2d (rtriaCoarse%h_IedgesAtElement,&
          p_IedgesAtElementC)
      call storage_getbase_int2d (rtriaFine%h_IedgesAtElement,&
          p_IedgesAtElementF)
      
      ! Keep the following picture in mind:
      !
      !    4---7---3      2          
      !    |   |   |      |`.        
      !    |   |   |      |  `.      
      !    8---9---6      5----4.    
      !    |   |   |      |`.  | `.  
      !    |   |   |      |  `.|   `.
      !    1---5---2      3----6-----1
      !
      ! We have vertices 1..4 (1..3) coming from old vertices, vertices 5..8 (4..6)
      ! coming from edges and vertex 9 coming from the element on quad-, resp. tri-, 
      ! meshes. This is exactly the order of the tags on the coarse grid.
      !
      ! Copy the refinement tags of the coarse grid vertices to the
      ! fine grid vertices that stem from coarse grid vertices.
      ! Copy the refinement tags of the old edges to the new vertices.
      ! Copy the refinement tags of the old quad elements to the new vertices
      ! in the element.
      call lalg_copyVectorInt (IrefTagsCoarse(1:),IrefTagsFine(1:),&
          NVT+NMT+rtriaCoarse%InelOfType(TRIA_NVEQUAD2D))
    
      ! Now to the refinement tags for the edges, this is more complicated.
      ! 
      !          3           
      !    +-12--+--8--+          +               
      !    |     |     |          |`.             
      !   11  4  9  3 10          6  7.           
      !    |     |     |          | 3  `.  2      
      !  4 +--3--+--6--+ 2      3 +--1---+.       
      !    |     |     |          |`.  1 | `.     
      !    4  1  2  2  5          9  2.  3   4.   
      !    |     |     |          | 4  `.| 2   `. 
      !    +--1--+--7--+          +--8---+---5---+
      !          1                       1        
      !
      ! For every coarse grid element, we have to look at the fine grid elements
      ! inside and transfer the refinement tags of the coarse grid edges to the
      ! fine grid ones 1,3,4,5,7,8,10,11,12 and the refinement tag of the element
      ! to the edges 2,3,6,9.
      !
      ! So loop through the coarse grid elements
      
      do iel = 1,nel
      
        ! How to transfer the tags now depends on whether we have
        ! a triangle or quad.
        nve = ubound(p_IedgesAtElementF,1)
        do while ((nve .gt. 1) .and. (p_IedgesAtElementF(nve,iel) .eq. 0))
          nve = nve - 1
        end do
        
        if (nve .eq. 4) then
          
          ! A coarse grid quad is decomponed into the four fine grid
          ! quads...
          iellocal(1) = iel
          iellocal(2) = rtriaCoarse%NEL+3*(iel-1)+1
          iellocal(3) = iellocal(2)+1
          iellocal(4) = iellocal(2)+2
        
          ! Get the tags on the coarse grid edges and the coarse grid element
          IrefTag(1) = IrefTagsCoarse(p_IedgesAtElementC(1,iel)+rtriaCoarse%NVT)
          IrefTag(2) = IrefTagsCoarse(p_IedgesAtElementC(2,iel)+rtriaCoarse%NVT)
          IrefTag(3) = IrefTagsCoarse(p_IedgesAtElementC(3,iel)+rtriaCoarse%NVT)
          IrefTag(4) = IrefTagsCoarse(p_IedgesAtElementC(4,iel)+rtriaCoarse%NVT)
          IrefTag(5) = IrefTagsCoarse(nvt+nmt+iel)
          
          ! Transfer the tags according to the above figure.
          IrefTagsFine(p_IedgesAtElementF(1,iellocal(1))+rtriaFine%NVT) = IrefTag(1)
          IrefTagsFine(p_IedgesAtElementF(2,iellocal(1))+rtriaFine%NVT) = IrefTag(5)
          IrefTagsFine(p_IedgesAtElementF(3,iellocal(1))+rtriaFine%NVT) = IrefTag(5)
          IrefTagsFine(p_IedgesAtElementF(4,iellocal(1))+rtriaFine%NVT) = IrefTag(4)
          
          IrefTagsFine(p_IedgesAtElementF(1,iellocal(2))+rtriaFine%NVT) = IrefTag(2)
          IrefTagsFine(p_IedgesAtElementF(2,iellocal(2))+rtriaFine%NVT) = IrefTag(5)
          IrefTagsFine(p_IedgesAtElementF(4,iellocal(2))+rtriaFine%NVT) = IrefTag(1)
          
          IrefTagsFine(p_IedgesAtElementF(1,iellocal(3))+rtriaFine%NVT) = IrefTag(3)
          IrefTagsFine(p_IedgesAtElementF(2,iellocal(3))+rtriaFine%NVT) = IrefTag(5)
          IrefTagsFine(p_IedgesAtElementF(4,iellocal(3))+rtriaFine%NVT) = IrefTag(2)
          
          IrefTagsFine(p_IedgesAtElementF(1,iellocal(4))+rtriaFine%NVT) = IrefTag(4)
          IrefTagsFine(p_IedgesAtElementF(4,iellocal(4))+rtriaFine%NVT) = IrefTag(3)
        else
          ! A coarse grid triangle is decomponed into the four fine grid
          ! quads...
          iellocal(1) = iel
          iellocal(2) = rtriaCoarse%NEL+3*(iel-1)+1
          iellocal(3) = iellocal(2)+1
          iellocal(4) = iellocal(2)+2

          ! Get the tags on the coarse grid edges and the coarse grid element
          IrefTag(1) = IrefTagsCoarse(p_IedgesAtElementC(1,iel)+rtriaCoarse%NVT)
          IrefTag(2) = IrefTagsCoarse(p_IedgesAtElementC(2,iel)+rtriaCoarse%NVT)
          IrefTag(3) = IrefTagsCoarse(p_IedgesAtElementC(3,iel)+rtriaCoarse%NVT)
          IrefTag(4) = IrefTagsCoarse(nvt+nmt+iel)
          
          ! Transfer the tags according to the above figure.
          IrefTagsFine(p_IedgesAtElementF(1,iellocal(1))+rtriaFine%NVT) = IrefTag(4)
          IrefTagsFine(p_IedgesAtElementF(2,iellocal(1))+rtriaFine%NVT) = IrefTag(4)
          IrefTagsFine(p_IedgesAtElementF(3,iellocal(1))+rtriaFine%NVT) = IrefTag(4)
          
          IrefTagsFine(p_IedgesAtElementF(1,iellocal(2))+rtriaFine%NVT) = IrefTag(2)
          IrefTagsFine(p_IedgesAtElementF(3,iellocal(2))+rtriaFine%NVT) = IrefTag(1)
          
          IrefTagsFine(p_IedgesAtElementF(1,iellocal(3))+rtriaFine%NVT) = IrefTag(3)
          IrefTagsFine(p_IedgesAtElementF(3,iellocal(3))+rtriaFine%NVT) = IrefTag(2)

          IrefTagsFine(p_IedgesAtElementF(1,iellocal(4))+rtriaFine%NVT) = IrefTag(1)
          IrefTagsFine(p_IedgesAtElementF(3,iellocal(4))+rtriaFine%NVT) = IrefTag(3)
        end if
        
        ! The four new elements receive the tag of the element
        IrefTagsFine(nvtfine+nmtfine+iellocal(1)) = IrefTagsCoarse(nvt+nmt+iel)
        IrefTagsFine(nvtfine+nmtfine+iellocal(2)) = IrefTagsCoarse(nvt+nmt+iel)
        IrefTagsFine(nvtfine+nmtfine+iellocal(3)) = IrefTagsCoarse(nvt+nmt+iel)
        IrefTagsFine(nvtfine+nmtfine+iellocal(4)) = IrefTagsCoarse(nvt+nmt+iel)
        
      end do
      
    end subroutine calcRefTags2D
  
  end subroutine tria_genRefTags2lv

  ! ***************************************************************************************

!<subroutine>

  subroutine tria_hangingNodeRefinement (rtriaSource,Ielements,rtriaDest)
  
!<description>
  ! Performs a local hanging-node refinement of the elements Ielements in the 
  ! triangulation rtriaSource.
  !
  ! Note: The implementation is still rudimentary. Currently, there are
  ! only 2D quad meshes allowed. 'Iterative' refinement (i.e. the refinement
  ! of elements which are not marked to be refined, due to a too high refinement
  ! level are not supported!
  ! Furthermore It is only possible to refine elements with hanging vertices
  ! if the subelements hanging on a hanging vertex are not refined
  ! at the same time.
!</description>
  
!<input>
  ! Source mesh to be refines.
  type(t_triangulation), intent(in) :: rtriaSource
  
  ! List of elements to be refined.
  integer, dimension(:), intent(in) :: Ielements
  
  ! OPTIONAL: Boundary structure that defines the domain.
  ! If not specified, information about boundary vertices (e.g. 
  ! parameter values of edge midpoints in 2D) are not initialised.
  !type(t_boundary), INTENT(in), OPTIONAL :: rboundary
!</input>

!<output>
  ! Refined mesh. Will be a standard mesh.
  type(t_triangulation), intent(out) :: rtriaDest
!</output>
  
!</subroutine>

    ! local variables
    integer :: i,nnve
    integer :: ivt,ivt1,ivt2,ive,nvtnew,ive2
    integer :: iel,nelnew,nelsum,ieldest,nquads
    integer :: imt,imtdest,nmtsum,imt2,iedge1,iedge2,iedge3,iedge4
    integer :: iel1,iel2,iel3
    integer, dimension(:), allocatable :: Iedges
    integer, dimension(:), allocatable :: IelementsSorted
    integer, dimension(:), allocatable :: IelementRef
    integer, dimension(:), allocatable :: IelementLocalId
    integer, dimension(:), allocatable :: IedgeHang
    integer, dimension(:), allocatable :: IedgeLocalId
    integer, dimension(:,:), pointer :: p_IelementsAtEdge
    integer, dimension(:,:), pointer :: p_IelementsAtEdgeDest
    integer, dimension(2) :: Isize
    integer, dimension(:,:), pointer :: p_IverticesAtElementSrc
    integer, dimension(:,:), pointer :: p_IneighboursAtElementSrc
    integer, dimension(:,:), pointer :: p_IneighboursAtElementDest
    integer, dimension(:,:), pointer :: p_IedgesAtElementSrc
    integer, dimension(:,:), pointer :: p_IedgesAtElementDest
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtElementDest
    real(dp), dimension(:,:), pointer :: p_DvertexCoordsSrc,p_DvertexCoordsDest
    integer, dimension(:), pointer :: p_InodalPropertySrc,p_InodalPropertyDst
    integer, dimension(:), pointer :: p_IelementsAtVertex,p_IelementsAtVertexIdx
    integer :: NVT
    
    nnve = rtriaSource%NNVE
    if (nnve .ne. TRIA_NVEQUAD2D) then
    
      call output_line (&
          '2-level refinement supports only quad meshes!', &
          OU_CLASS_ERROR,OU_MODE_STD,'tria_hangingNodeRefinement')
      call sys_halt()
      
    end if

    ! At the beginning, mark the elements to be refined
    allocate(IelementRef(rtriaSource%NEL))
    IelementRef(:) = 0
    do iel = 1,size(Ielements)
      IelementRef(Ielements(iel)) = 1
    end do
    
    call storage_getbase_int2d(rtriaSource%h_IelementsAtEdge,p_IelementsAtEdge)

    call storage_getbase_int2D(&
        rtriaSource%h_IedgesAtElement,p_IedgesAtElementSrc)

    call storage_getbase_int2D(&
        rtriaSource%h_IverticesAtEdge,p_IverticesAtEdge)

    call storage_getbase_int2D(&
        rtriaSource%h_IneighboursAtElement,p_IneighboursAtElementSrc)
    
    call storage_getbase_int(&
        rtriaSource%h_InodalProperty,p_InodalPropertySrc)

    NVT = rtriaSource%NVT

    ! Loop through all elements and mark the edges to be refined.
    allocate(IedgeHang(rtriaSource%NMT))
    IedgeHang(:) = 0
    do iel = 1,size(Ielements)
      do ive=1,ubound(p_IedgesAtElementSrc,1)
        if (p_IedgesAtElementSrc(ive,Ielements(iel)) .ne. 0) then
          imt = p_IedgesAtElementSrc(ive,Ielements(iel))

          ! Increase IedgeHang. Inner vertices get a 2 (as they are touched twice
          ! or when they have already a hanging vertex), 
          ! edges with hanging nodes get a 1.
          IedgeHang(imt) = IedgeHang(imt) + 1
          if (p_IneighboursAtElementSrc(ive,Ielements(iel)) .le. 0) then
            ! For boundary edges on refined elements, we manually increase
            ! IedgeHang to 2 as these do not produce hanging nodes!
            IedgeHang(imt) = 2
          end if
          
          if (p_InodalPropertySrc(imt+NVT) .lt. 0) then
            ! Mark edges with already hanging vertices with a 3;
            ! these edges must not produce new vertices.
            IedgeHang(imt) = 3
          end if
        end if
      end do
    end do
    
    ! We generate a 2-level-ordering like refinement. In our case that
    ! means:
    ! - Vertex numbers in the old (coarse) mesh get vertex numbers
    !   in the refined mesh
    ! - Refined edges procude new vertices which are attached to the
    !   coarse grid vertices
    ! - Refined elements produce new vertices which are attached to !
    !   this vertex set.
    ! - Edges that have already hanging vertices do not produce new ones.
    !
    ! In a first step, set up an array that assigns each refined
    ! edge a 'local' number -- so an edge will be the i-th refined edge.
    ! Not refined edges are dummy values in this array and not used later on.
    allocate(IedgeLocalId(rtriaSource%NMT))
    
    if (IedgeHang(1) .ne. 0) then
      IedgeLocalId(1) = 1
    else
      IedgeLocalId(1) = 0
    end if
    
    do imt = 2,rtriaSource%NMT
      if ((IedgeHang(imt) .ne. 0) .and. (IedgeHang(imt) .ne. 3)) then
        IedgeLocalId(imt) = IedgeLocalId(imt-1) + 1
      else
        IedgeLocalId(imt) = IedgeLocalId(imt-1)
      end if
    end do
    
    ! In the same way, assign each element a local id.
    allocate(IelementLocalId(rtriaSource%NEL))
    IelementLocalId(1) = IelementRef(1)
    
    do iel = 2,rtriaSource%NEL
      IelementLocalId(iel) = IelementLocalId(iel-1) + IelementRef(iel)
    end do

    ! Count the number elements we have now.
    ! Every refined element creates 3 new -- additional to itself.
    nelsum = size(Ielements)
    nelnew = nelsum * 3
    
    ! Count the number of vertices. Every refined element gives one vertex,
    ! every refined edge gets one vertex.
    nmtsum = IedgeLocalId(rtriaSource%NMT)
    nvtnew = nmtsum + nelsum
    
    ! Generate a list of all edges / elements (sorted) to be refined.
    allocate(Iedges(nmtsum))
    imtdest = 0
    do imt = 1,rtriaSource%NMT
      if ((IedgeHang(imt) .ne. 0)  .and. (IedgeHang(imt) .ne. 3)) then
        imtdest = imtdest + 1
        Iedges(imtdest) = imt
      end if
    end do

    allocate(IelementsSorted(nelsum))
    ieldest = 0
    do iel = 1,rtriaSource%NEL
      if (IelementRef(iel) .ne. 0) then
        ieldest = ieldest + 1
        IelementsSorted(ieldest) = iel
      end if
    end do
    
    ! Set up basic information
    rtriaDest%ndim = rtriaSource%ndim
    rtriaDest%NMT = 0
    rtriaDest%NNVE = rtriaSource%NNVE
    rtriaDest%NNEE = rtriaSource%NNEE
    rtriaDest%NBCT = rtriaSource%NBCT
    rtriaDest%NblindBCT = rtriaSource%NblindBCT
    rtriaDest%NEL = rtriaSource%NEL + nelnew
    rtriaDest%NVT = rtriaSource%NVT + nvtnew
    
    ! Vertex numbering:
    ! 1..rtriaCoarse%NVT : old vertices
    ! rtriaCoarse%NVT+1 .. rtriaCoarse%NVT + nvtnew : 
    !                      new vertices from edges
    ! rtriaCoarse%NVT+nvtnew+1 .. rtriaCoarse%NVT+nvtnew+nelnew :
    !                      new vertices from elements
    !    
    ! To understand the numbering, keep the following pictures in mind:
    !
    ! Coarse mesh:
    !
    !   4------21--------7-------17-------3
    !   |                |                |
    !   |                |                |
    !   |                |                |
    !  20       4       18        3       19
    !   |                |                |
    !   |                |                |
    !   |                |                |
    !   8------12--------9-------15-------6
    !   |                |                |
    !   |                |                |
    !   |                |                |
    !  13       1       11        2       14
    !   |                |                |
    !   |                |                |
    !   |                |                |
    !   1------10--------5-------16-------2    
    !
    ! Mesh with cell 1+4 refined. Edges / element midpoints are renumbered to a 
    ! consecutive order: 10,11,12,13,18,... -> E1,E2,... are the vertices 
    ! that stem from edges, 1,4 -> M1,M2 the vertices stemming from element 
    ! midpoints.
    !
    !   4------E6--------7----------------3
    !   |       |        |                |
    !   |   4   |   10   |                |
    !   |       |        |                |
    !  E7------M2-------E5        3       | 
    !   |       |        |                |
    !   |   8   |   9    |                |
    !   |       |        |                |
    !   8------E3--------9----------------6
    !   |       |        |                |
    !   |   7   |   6    |                |
    !   |       |        |                |
    !  E4------M1-------E2        2       | 
    !   |       |        |                |
    !   |   1   |   5    |                |
    !   |       |        |                |
    !   1------E1--------5----------------2
    !
    ! Mesh with cell 1+4 refined. Fine grid vertices renumbered.
    !
    !   4------15--------7----------------3
    !   |       |        |                |
    !   |   4   |   10   |                |
    !   |       |        |                |
    !  16------18-------14        3       | 
    !   |       |        |                |
    !   |   8   |   9    |                |
    !   |       |        |                |
    !   8------12--------9----------------6
    !   |       |        |                |
    !   |   7   |   6    |                |
    !   |       |        |                |
    !  13------17-------11        2       | 
    !   |       |        |                |
    !   |   1   |   5    |                |
    !   |       |        |                |
    !   1------10--------5----------------2
    !
    ! Allocate memory for IverticesAtElement.
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriaDest%NNVE,rtriaDest%NEL/)
    call storage_new ('tria_hangingNodeRefinement', 'KVERT', Isize, ST_INT, &
        rtriaDest%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)

    ! Get the pointer to the IverticesAtElement array
    call storage_getbase_int2D(&
        rtriaDest%h_IverticesAtElement,p_IverticesAtElementDest)
    call storage_getbase_int2D(&
        rtriaSource%h_IverticesAtElement,p_IverticesAtElementSrc)

    ! Create fine grid elements. Loop through the elements to be refined.
    do iel = 1,rtriaSource%NEL
    
      ! Is that element one to be refined?
      if (IelementRef(iel) .eq. 0) then
      
        ! Transfer 'Coarse grid' vertices
        do i=1,ubound(p_IverticesAtElementSrc,1)
          p_IverticesAtElementDest(i,iel) = p_IverticesAtElementSrc(i,iel)
        end do
      
      else
      
        ! This is a fine grid element...
      
        ! Determine number of subelements.
        iel1 = rtriaSource%NEL+3*(IelementLocalId(iel)-1)+1
        iel2 = iel1+1
        iel3 = iel1+2
        
        ! In nquads we count the number of the quad +NVT+NMT we process.
        ! That is the number of the midpoint of that element!
        nquads = rtriaSource%NVT + nmtsum + IelementLocalId(iel)
        
        ! Get the vertex numbers that stem from the four edges of the coarse grid element.
        ! If this element has already a hanging vertex on an edge,
        ! the vertex number is taken from the coarse grid.

        if (p_InodalPropertySrc(p_IedgesAtElementSrc(1,iel)+NVT) .lt. 0) then
          iedge1 = -p_InodalPropertySrc(p_IedgesAtElementSrc(1,iel)+NVT)
        else        
          iedge1 = IedgeLocalId(p_IedgesAtElementSrc (1,iel)+NVT)
        end if
        
        if (p_InodalPropertySrc(p_IedgesAtElementSrc(2,iel)+NVT) .lt. 0) then
          iedge2 = -p_InodalPropertySrc(p_IedgesAtElementSrc(2,iel)+NVT)
        else        
          iedge2 = IedgeLocalId(p_IedgesAtElementSrc (2,iel))
        end if
        
        if (p_InodalPropertySrc(p_IedgesAtElementSrc(3,iel)+NVT) .lt. 0) then
          iedge3 = -p_InodalPropertySrc(p_IedgesAtElementSrc(3,iel)+NVT)
        else        
          iedge3 = IedgeLocalId(p_IedgesAtElementSrc (3,iel))
        end if
        
        if (p_InodalPropertySrc(p_IedgesAtElementSrc(4,iel)+NVT) .lt. 0) then
          iedge4 = -p_InodalPropertySrc(p_IedgesAtElementSrc(4,iel)+NVT)
        else        
          iedge4 = IedgeLocalId(p_IedgesAtElementSrc (4,iel))
        end if
        
        
        ! To convert edges vertex numbers, we have to convert:
        !  old edge number 
        !  -> local id of the refined edge
        !  -> new vertex number (=nvt + local id of the refined edge)
        
        ! Step 1: Initialise IverticesOnElement for element IEL
        p_IverticesAtElementDest(1,iel) = p_IverticesAtElementSrc (1,iel)
        p_IverticesAtElementDest(2,iel) = iedge1+rtriaSource%NVT
        p_IverticesAtElementDest(3,iel) = nquads
        p_IverticesAtElementDest(4,iel) = iedge4+rtriaSource%NVT
        
        ! Step 2: Initialise IverticesOnElement for element IEL1
        p_IverticesAtElementDest(1,iel1) = p_IverticesAtElementSrc (2,iel)
        p_IverticesAtElementDest(2,iel1) = iedge2+rtriaSource%NVT
        p_IverticesAtElementDest(3,iel1) = nquads
        p_IverticesAtElementDest(4,iel1) = iedge1+rtriaSource%NVT
      
        ! Step 3: Initialise IverticesOnElement for element IEL2
        p_IverticesAtElementDest(1,iel2) = p_IverticesAtElementSrc (3,iel)
        p_IverticesAtElementDest(2,iel2) = iedge3+rtriaSource%NVT
        p_IverticesAtElementDest(3,iel2) = nquads
        p_IverticesAtElementDest(4,iel2) = iedge2+rtriaSource%NVT

        ! Step 4: Initialise IverticesOnElement for element IEL3
        p_IverticesAtElementDest(1,iel3) = p_IverticesAtElementSrc (4,iel)
        p_IverticesAtElementDest(2,iel3) = iedge4+rtriaSource%NVT
        p_IverticesAtElementDest(3,iel3) = nquads
        p_IverticesAtElementDest(4,iel3) = iedge3+rtriaSource%NVT
        
      end if
    end do
        
    ! Create the coordinate array.
    Isize = (/rtriaSource%ndim,rtriaDest%NVT/)
    call storage_new ('tria_hangingNodeRefinement', 'DCORVG',&
        Isize, ST_DOUBLE, &
        rtriaDest%h_DvertexCoords, ST_NEWBLOCK_ZERO)
        
    call storage_getbase_double2D(&
        rtriaSource%h_DvertexCoords,p_DvertexCoordsSrc)
    call storage_getbase_double2D(&
        rtriaDest%h_DvertexCoords,p_DvertexCoordsDest)
    
    ! 'Coarse grid' vertices
    do ivt = 1,rtriaSource%NVT
      do i=1,ubound(p_DvertexCoordsDest,1)
        p_DvertexCoordsDest(i,ivt) = p_DvertexCoordsSrc(i,ivt)
      end do
    end do
    
    ! New vertices by edges
    do imt = 1,nmtsum
      do i=1,ubound(p_DvertexCoordsDest,1)
        ivt1 = p_IverticesAtEdge(1,Iedges(imt))
        ivt2 = p_IverticesAtEdge(2,Iedges(imt))
        p_DvertexCoordsDest(i,rtriaSource%NVT+imt) = &
            0.5_DP * (p_DvertexCoordsSrc(i,ivt1) + p_DvertexCoordsSrc(i,ivt2))
      end do
    end do
    
    ! New vertices by elements
    elementloop: do iel = 1,nelsum
    
      p_DvertexCoordsDest(:,rtriaSource%NVT+nmtsum+iel) = 0.0_DP
      
      do i=1,ubound(p_IverticesAtElementSrc,1)
        ivt = p_IverticesAtElementSrc(i,IelementsSorted(iel))
        if (ivt .ne. 0) then
          p_DvertexCoordsDest(:,rtriaSource%NVT+nmtsum+iel) = &
              p_DvertexCoordsDest(:,rtriaSource%NVT+nmtsum+iel) + &
              p_DvertexCoordsSrc(:,ivt)
        else
          ! Triangle in a quad mesh e.g.
          ! Divide to get the mean.
          p_DvertexCoordsDest(:,rtriaSource%NVT+nmtsum+iel) = &
            p_DvertexCoordsDest(:,rtriaSource%NVT+nmtsum+iel) / real(i,dp)
          cycle elementloop
        end if
      end do

      ! Divide to get the mean. Note that i is nve+1 here, so subtract 1...
      p_DvertexCoordsDest(:,rtriaSource%NVT+nmtsum+iel) = &
        p_DvertexCoordsDest(:,rtriaSource%NVT+nmtsum+iel) / real(i-1,dp)
      
    end do elementloop
    
    ! Set up the nodal property array.
    ! The nodal property of old vertices can be copied, that of new
    ! vertices has to be set up manually
    call storage_new ('tria_hangingNodeRefinement', 'KNPR', &
        rtriaDest%NVT, ST_INT, &
        rtriaDest%h_InodalProperty, ST_NEWBLOCK_ZERO)

    call storage_getbase_int(&
        rtriaDest%h_InodalProperty,p_InodalPropertyDst)
    
    call lalg_copyVectorInt(p_InodalPropertySrc,p_InodalPropertyDst,&
        rtriaSource%NVT)
    
    ! Old edges -> New vertices
    do imt = 1,nmtsum
      p_InodalPropertyDst(NVT+imt) = &
        p_InodalPropertySrc(NVT+Iedges(imt))
    end do
    
    ! Old elements -> new vertices.
    ! they are in the domain, so they receive nodal property 0 -- what they 
    ! already have by initialisation.
    !
    ! Generate basic information about boundary vertices.
    call tria_genRawBoundary2D (rtriaDest)
    
    ! Now start to generate a standard mesh.
    call tria_initStandardMeshFromRaw(rtriaDest)!,rboundary)
    
    ! But that is not enough. The KNPR array as well as the IneighboursAtElement
    ! array are wrong up to now.
    !
    ! At first, we correct KNPR. It may be recalculated, so retrieve
    ! the pointer again.
    call storage_getbase_int(&
        rtriaDest%h_InodalProperty,p_InodalPropertyDst)

    call storage_getbase_int(&
      rtriaDest%h_IelementsAtVertex,p_IelementsAtVertex)
    call storage_getbase_int(&
      rtriaDest%h_IelementsAtVertexIdx,p_IelementsAtVertexIdx)
    call storage_getbase_int2d(&
      rtriaDest%h_IverticesAtElement,p_IverticesAtElementDest)
    call storage_getbase_int2d(&
      rtriaDest%h_IedgesAtElement,p_IedgesAtElementDest)

    ! For the hanging nodes, initialise KNPR appropriately.
    ! The hanging nodes can be identified by a "2" in the IedgeHang array.
    !
    ! We have to search the elements in the source mesh as we do not
    ! know the edge numbers in the target mesh.
    do iel = 1,rtriaSource%NEL
    
      ! Search the not refined elements for hanging vertices
      if (IelementRef(iel) .eq. 0) then
    
        ! Search the edges for hanging vertices
        do ive2 = 1,ubound(p_IedgesAtElementSrc,1)
         
          imt = p_IedgesAtElementSrc(ive2,iel)-rtriaSource%NVT
          
          if (IedgeHang(imt) .eq. 1) then
          
            ! Here`s a hanging vertex. Get the vertex number.
            ivt = rtriaSource%NVT + IedgeLocalId(imt)
            
            ! Get the two adjacent vertices on the parent edge.
            ivt1 = p_IverticesAtElementSrc(ive2,iel)
            ivt2 = p_IverticesAtElementSrc(mod(ive2,nnve)+1,iel)
            
            ! Get the edge number in the destination mesh. As element numbers
            ! of not refined elements coincide in both meshes, we can simply
            ! calculate the new edge number from the position in the element.
            imt2 = p_IedgesAtElementDest (ive2,iel)
            
            ! Put it to the nodal property array.
            p_InodalPropertyDst(imt2+rtriaDest%NVT) = -ivt
            
            ! Mark the vertex as hanging vertex
            p_InodalPropertyDst(ivt) = -imt2
            
            ! Find the two subedges of the big edge and mark them.
            ! Note that there are exactly 2 elements adjacent to that
            ! vertex by construction!
            iel1 = p_IelementsAtVertex(p_IelementsAtVertexIdx(ivt))
            iel2 = p_IelementsAtVertex(p_IelementsAtVertexIdx(ivt)+1)
            
            do ive = 1,ubound(p_IverticesAtElementDest,1)
              if ((p_IverticesAtElementDest(ive,iel1) .eq. ivt) .and. &
                  (p_IverticesAtElementDest(mod(ive,nnve)+1,iel1) .eq. ivt1)) then
                ! The edge 'after' the vertex is the neighbour.
                imt2 = p_IedgesAtElementDest(ive,iel1)
                p_InodalPropertyDst(imt2+NVT) = -IedgeLocalId(imt)
                exit
              end if
              
              if ((p_IverticesAtElementDest(ive,iel1) .eq. ivt2) .and. &
                  (p_IverticesAtElementDest(mod(ive,nnve)+1,iel1) .eq. ivt)) then
                ! The edge 'before' the vertex is the neighbour.
                imt2 = p_IedgesAtElementDest(ive,iel1)
                p_InodalPropertyDst(imt2+NVT) = -IedgeLocalId(imt)
                exit
              end if
            end do
            
            do ive = 1,ubound(p_IverticesAtElementDest,1)
              if ((p_IverticesAtElementDest(ive,iel2) .eq. ivt) .and. &
                  (p_IverticesAtElementDest(mod(ive,nnve)+1,iel2) .eq. ivt1)) then
                ! The edge 'after' the vertex is the neighbour.
                imt2 = p_IedgesAtElementDest(ive,iel2)
                p_InodalPropertyDst(imt2+NVT) = -IedgeLocalId(imt)
                exit
              end if
              
              if ((p_IverticesAtElementDest(ive,iel2) .eq. ivt2) .and. &
                  (p_IverticesAtElementDest(mod(ive,nnve)+1,iel2) .eq. ivt)) then
                ! The edge 'before' the vertex is the neighbour.
                imt2 = p_IedgesAtElementDest(ive,iel2)
                p_InodalPropertyDst(imt2+NVT) = -IedgeLocalId(imt)
                exit
              end if
            end do
            
          end if      
        end do
        
      end if
    
    end do
    
    ! The last thing: We have to correct the neighbourhood of the elements
    ! at hanging vertices. Up to now, the elements sharing an edge with
    ! a hanging vertex are not connected to a neighbourhood -- this we have 
    ! to change.
    call storage_getbase_int2d(&
      rtriaDest%h_IneighboursAtElement,p_IneighboursAtElementDest)
      
    call storage_getbase_int2d(rtriaDest%h_IelementsAtEdge,p_IelementsAtEdgeDest)
      
    ! Find the hanging vertices
    do ivt = 1,rtriaDest%NVT
      if (p_InodalPropertyDst(ivt) .lt. 0) then
        ! Get the coarse grid edge the vertex is the midpoint from
        imt = -p_InodalPropertyDst(ivt)
        
        ! Get the coarse grid element. The tria_initStandardMeshFromRaw
        ! recognised only this one on the edge, so by taking the first element
        ! adjacent to the edge, we have it.
        iel = p_IelementsAtEdgeDest(1,imt)
        
        ! Find the two subelements attached to the vertex. We can find
        ! them as they are the only elements adjacent to the vertex.
        iel1 = p_IelementsAtVertex(p_IelementsAtVertexIdx(ivt))
        iel2 = p_IelementsAtVertex(p_IelementsAtVertexIdx(ivt)+1)
        
        ! Save the first subelement as neighbour of the coarse grid element
        do ive = 1,nnve
          if (p_IedgesAtElementDest(ive,iel)-rtriaDest%NVT .eq. imt) then
            p_IneighboursAtElementDest(ive,iel) = iel1
            exit
          end if
        end do
        
        ! Save the coarse grid element as neighbour of the fine grid elements.
        ! Note that we have to search for the correct edge.
        ! For that purpose, search for the vertices. In ivt1/ivt2 we save the two
        ! endpoints of the edge in the coarse grid; these are adjacent to
        ! the hanging vertex on the fine grid. Beginning and end of the edge
        ! can be accessed via the ive calculated above.
        ivt1 = p_IverticesAtElementDest(ive,iel)
        ivt2 = p_IverticesAtElementDest(mod(ive,nnve)+1,iel)
        
        do ive = 1,nnve

          if (((p_IverticesAtElementDest(ive,iel1) .eq. ivt) .and. &
               (p_IverticesAtElementDest(mod(ive,nnve)+1,iel1) .eq. ivt1)) .or. &
              ((p_IverticesAtElementDest(ive,iel1) .eq. ivt2) .and. &
               (p_IverticesAtElementDest(mod(ive,nnve)+1,iel1) .eq. ivt))) then
            p_IneighboursAtElementDest(ive,iel1) = iel
          end if

          if (((p_IverticesAtElementDest(ive,iel2) .eq. ivt) .and. &
               (p_IverticesAtElementDest(mod(ive,nnve)+1,iel2) .eq. ivt1)) .or. &
              ((p_IverticesAtElementDest(ive,iel2) .eq. ivt2) .and. &
               (p_IverticesAtElementDest(mod(ive,nnve)+1,iel2) .eq. ivt))) then
            p_IneighboursAtElementDest(ive,iel2) = iel
          end if
        end do
        
      end if
    end do

  end subroutine tria_hangingNodeRefinement

  !====================================================================
  !
  !        ++       ++++
  !       + +       +   +  
  !         +       +    +
  !         +       +    +
  !         +       +    + 
  !         +       +   +  
  !       +++++     ++++
  ! tag@1D
  !====================================================================
  
!<subroutine>

  subroutine tria_createRawTria1D(rtriangulation, dleft, dright, nintervals)
  
!<description>
  ! This routine creates a 'raw' 1D triangulation with nintervals
  ! sub-intervals of same length.
!</description>

!<input>
  ! The left end of the interval. Must be < dright.
  real(DP), intent(in) :: dleft
  
  ! The right end of the interval. Must be > dleft.
  real(DP), intent(in) :: dright
  
  ! OPTIONAL: The number of sub-intervals to create. If given, nintervals
  ! must be > 0. If not given, one interval is created.
  integer, optional, intent(in) :: nintervals
!</input>

!<output>
  ! The triangulation.
  type(t_triangulation), intent(out) :: rtriangulation
!</output>

!</subroutine>

  ! Some local variables
  integer :: i,nintv
  real(DP) :: t, s
  real(DP), dimension(:,:), pointer :: p_Dcoords
  integer, dimension(:,:), pointer :: p_Iverts
  integer, dimension(:), pointer :: p_Idata
  integer, dimension(2) :: Isize
    
    ! Check parameters
    if (dleft .ge. dright) then
      call output_line ('dleft must be less than dright!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_createRawTria1D')
      call sys_halt()
    end if
    
    ! Set number of intervals
    nintv = 1
    if (present(nintervals)) then
      if (nintervals .gt. 0) nintv = nintervals
    end if
  
    ! Set the triangulation`s dimension
    rtriangulation%ndim = NDIM1D

    ! We have nintv+1 vertices
    rtriangulation%NVT = nintv+1
    
    rtriangulation%NNVE = 2

    ! Allocate vertices
    Isize = (/1, nintv+1/)
    call storage_new('tria_createRawTria1D', 'DCORVG',&
        Isize, ST_DOUBLE,&
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords,&
        p_Dcoords)
    
    ! Initialise vertices
    p_Dcoords(1,1) = dleft
    s = 1.0_DP / real(nintv, DP)
    do i=2, nintv
      t = real(i-1, DP) * s
      p_Dcoords(1,i) = (1.0_DP - t) * dleft + t * dright
    end do
    p_Dcoords(1,nintv+1) = dright
    
    ! And we have nintv elements
    rtriangulation%NEL = nintv
    rtriangulation%InelOfType(TRIA_NVELINE1D) = nintv

    ! Allocate elements
    Isize = (/2, nintv/)
    call storage_new('tria_createRawTria1D', 'KVERT', Isize, &
        ST_INT, rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement, p_Iverts)
    
    ! Initialise elements
    do i=1, nintv
      p_Iverts(1,i) = i
      p_Iverts(2,i) = i+1
    end do
    
    ! There are two boundary components 
    ! - the interval start and end point
    rtriangulation%NBCT = 2
    
    ! Allocate memory for boundary components
    call storage_new ('tria_createRawTria1D', 'KBCT', 3, ST_INT, &
        rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the boundary components
    call storage_getbase_int(rtriangulation%h_IboundaryCpIdx, p_Idata)
    p_Idata = (/ 1, 2, 3 /)

    ! There is one vertice per boundary component
    rtriangulation%NVBD = 2
    
    ! Allocate memory for boundary components
    call storage_new ('tria_createRawTria1D', 'KVBD', 2, ST_INT, &
        rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the boundary components
    call storage_getbase_int(rtriangulation%h_IverticesAtBoundary, p_Idata)
    p_Idata = (/ 1, nintv+1 /)

    ! Allocate memory for nodal property
    call storage_new ('tria_createRawTria1D', 'KNPR', nintv+1, ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(rtriangulation%h_InodalProperty,p_Idata)
    
    ! Set up nodal property
    p_Idata(1) = 1
    p_Idata(nintv+1) = 2

    ! That is it
      
  end subroutine tria_createRawTria1D

  !************************************************************************

!<subroutine>

  subroutine tria_readTriFile1D(rtriangulation, sfilename, bnoExtendedRaw)

!<description>
  ! This routine reads a .TRI file of a 1D triangulation into memory
  ! and creates a 'raw' triangulation (i.e. a triangulation that contains
  ! only basic information, see below). 
  !
  ! The triangulation structure rtriangulation is initialised with the data 
  ! from the file. The parameter sfilename gives the name of the .tri 
  ! file to read.
  !
  ! This reads only the very basic information that is needed to create
  ! the coarse grid. No information about boundary points, subdivisions
  ! or whatever is created.
  ! The following arrays / information tags will be initialised:
  !
  ! NEL,NVT,NMT,NBCT,NVBD,InelOfType,
  ! DvertexCoords, IverticesAtElement, InodalProperty, 
  ! IboundaryCpIdx, IverticesAtBoundary, DvertexParameterValue.
  ! The nodal property array InodalProperty contains only the data for the
  ! vertices of the triangulation.
  ! The arrays IverticesAtBoundary and DvertexParameterValue are sorted
  ! for the boundary component (according to IboundaryCpIdx) but not 
  ! for the parameter value.
  !
  ! If bnoExtendedRaw is not specified or set to .FALSE., an extended
  ! raw mesh is generated that provides also a proper numbering for edges.
  ! In that case, the following arrays are initialised as well:
  !
  ! IelementsAtVertexIdx,IelementsAtVertexIdx,IneighboursAtElement,
  ! IedgesAtElement,IfacesAtElement
  !
  ! The triangulation structure rtriangulation must be empty. All previous
  ! information in this structure (if there is any) is lost!
!</description>

!<input>
  ! The name of the .tri file to read.
  character(LEN=*), intent(in) :: sfilename
  
  ! OPTIONAL: Prevent creation of an extended raw mesh. If set to .false.,
  ! an 'extended raw' mesh will be created that provides a proper numbering
  ! for edges (standard). If set to '.true', the result will be a 'really raw'
  ! raw mesh with minimum information and no numbering for edges.
  logical, intent(in), optional :: bnoExtendedRaw
! </input>
  
!<output>
  ! Triangulation structure, to be filled with data
  type(t_triangulation), intent(out) :: rtriangulation
!</output>
  
!</subroutine>

    ! Input channel for reading
    integer :: iunit
    
    ! Open the file
    call io_openFileForReading(sfilename, iunit)

    ! We create a 1D triangulation here.
    rtriangulation%ndim = NDIM1D

    ! Read the basic mesh
    call tria_readRawTriangulation1D (iunit, rtriangulation)

    ! Create the basic boundary information
    call tria_genRawBoundary1D (rtriangulation)

    ! Extend the raw mesh by basic edge numbering,
    ! initialise an extended raw mesh.
    if (.not.present(bnoExtendedRaw)) then
      call tria_initExtendedRawMesh (rtriangulation)
    else
      if (.not.bnoExtendedRaw) then
        call tria_initExtendedRawMesh (rtriangulation)
      end if
    end if

    ! Close the file, finish
    close(iunit)

  end subroutine tria_readTriFile1D

  !************************************************************************

!<subroutine>

  subroutine tria_readRawTriangulation1D (iunit, rtriangulation)

!<description>  
  ! Auxiliary routine of tria_readTriFile1D.
  ! Reads basic information from a triangulation file into rtriangulation.
  ! That means, the following information arrays / tags are initialised:
  ! NEL,NVT,NMT,NBCT,NNEE,InelOfType,
  ! DvertexCoords, IverticesAtElement and InodalProperty.
  ! The data is read from the file without being changed!
!</description>
  
!<input>
  ! Unit number of the file to be read
  integer, intent(in) :: iunit
!</input>
  
!<inputoutput> 
  ! Triangulation to be initialised with basic data.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! Local variables
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    integer, dimension(:,:), pointer :: p_Idata2D
    integer, dimension(:), pointer :: p_Idata
    integer :: ivt, iel
    integer :: ive
    integer, dimension(2) :: Isize
    
    ! The first two lines in the file are comments.
    read(iunit,*)
    read(iunit,*)
    
    ! Read NEL,NVT,NMT,NVE,NBCT from the file
    ! and store this information in the structure.
    read (iunit,*) rtriangulation%NEL, rtriangulation%NVT,&
                   rtriangulation%NMT, rtriangulation%NNVE,&
                   rtriangulation%NBCT
       
    ! Check consistency: NNVE = 2
    if (rtriangulation%NNVE .ne. 2) then
      call output_line ('Triangulation structure is invalid: NNVE does not match 2!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_readRawTriangulation1D')
      call sys_halt()
    end if
    
    ! Comment: 'DCORVG'
    read (iunit,*)

    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM2D, NVT)
    Isize = (/NDIM1D,rtriangulation%NVT/)
    call storage_new ('tria_readRawTriangulation1D', 'DCORVG',&
        Isize, ST_DOUBLE,&
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_Ddata2D)
        
    ! Read the data from the file, store it in the array.
    ! read data into p_Ddata:
    ! read nvt x-coordinates into p_Ddata(1,ivt)
    read (iunit,*) (p_Ddata2D(NDIM1D,ivt), ivt=1,rtriangulation%NVT)

    ! Comment: 'KVERT'
    read (iunit,*)

    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(2, NEL)
    Isize = (/2,rtriangulation%NEL/)
    call storage_new ('tria_readRawTriangulation1D', 'KVERT', Isize,&
        ST_INT, rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    call storage_getbase_int2D(rtriangulation%h_IverticesAtElement, p_Idata2D)

    ! read ive=1 indices to 2 into p_Idata2D(ive,iel) where iel=1 to NEL
    read (iunit,*) ((p_Idata2D(ive,iel),ive=1,2), iel=1,rtriangulation%NEL)

    ! We have only 1D elements
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(TRIA_NVELINE1D) = rtriangulation%NEL
    
    ! Comment: 'KNPR'
    read (iunit,*)

    ! Allocate memory for InodalProperty 
    call storage_new ('tria_readRawTriangulation1D', 'KNPR', &
        rtriangulation%NVT, ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(rtriangulation%h_InodalProperty, p_Idata)

    ! Read the data
    read (iunit,*) (p_Idata(ivt), ivt=1,rtriangulation%NVT)

  end subroutine tria_readRawTriangulation1D

  !************************************************************************

!<subroutine>

  subroutine tria_genRawBoundary1D (rtriangulation)

!<description>  
  ! Auxiliary routine of tria_readTriFile1D.
  ! This routine initialises basic boundary arrays and cleans up 
  ! a basic triangulation. That means:
  ! -> NVBD is calculated
  ! -> IboundaryCpIdx is created and generated
  ! -> IverticesAtBoundary is created and generated.
  !    The vertices are ordered for the boundary component according
  !    to IboundaryCpIdx but not ordered for their parameter value.
!</description>
  
!<inputoutput>
  ! Triangulation to be initialised with basic data.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_InodalProperty
    integer :: ivbd,ivt
    integer :: ibct

    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty, p_InodalProperty)

    ! Each boundary component corresponds to a single boundary vertex
    rtriangulation%NVBD = rtriangulation%NBCT

    ! Allocate memory for IverticesAtBoundary.
    call storage_new ('tria_genRawBoundary1D', 'KVBD', &
        rtriangulation%NVBD, ST_INT, &
        rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    call storage_new ('tria_genRawBoundary1D', 'KBCT', &
        rtriangulation%NBCT+1, ST_INT, &
        rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    call storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    call storage_getbase_double2D (&
        rtriangulation%h_DvertexCoords,p_DvertexCoords)
        
    call storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! The first element in p_IboundaryCpIdx is (as the head) always =1.
    p_IboundaryCpIdx(1) = 1

    ! Perform a first loop over all vertices to check which are on the
    ! boundary. Count them. This way, we at first set up the boundary
    ! component index vector p_IboundaryCpIdx.
    ! In this first step, we save the number of vertices on each 
    ! boundary component in p_IboundaryCpIdx(2:NBCT+1)!
    do ivt=1,rtriangulation%NVT
      if (p_InodalProperty(ivt) .gt. 0) then
        ibct = p_InodalProperty(ivt)
        
        ! Increase the number of vertices in that boundary component by 1.
        ! The number of vertices on boundary component i is saved here
        ! at p_IboundaryCpIdx(i+1) for later.
        ! Note that the array was initialised with zero during the creation
        ! process!
        p_IboundaryCpIdx(ibct+1) = p_IboundaryCpIdx(ibct+1)+1
      end if
    end do
    
    ! Sum up the number of vertices on each boundary component to get the
    ! actual index vector.
    do ibct = 2, rtriangulation%NBCT+1
      p_IboundaryCpIdx(ibct) = p_IboundaryCpIdx(ibct)+p_IboundaryCpIdx(ibct-1)
    end do
    
    ! Shift the p_IboundaryCpIdx array by one position. That is a little trick in
    ! the use of p_IboundaryCpIdx!
    ! Imagine, we have 3 boundary components with 8,6 and 4 edges. 
    ! Before summing the entries up, p_IboundaryCpIdx may have looked like this:
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   8   6   4
    !
    ! Summing that up gave the actual indices:
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   9  15  19
    !
    ! which is already the correct index array.
    ! Shifting p_IboundaryCpIdx by 1 will give:
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   1   9  15
    
    p_IboundaryCpIdx(2:rtriangulation%NBCT+1) = p_IboundaryCpIdx(1:rtriangulation%NBCT)
   
    ! Then, we again loop through all vertices and collect those on the
    ! boundary. In that loop, we use p_IboundaryCpIdx(2:NBCT+1) as pointer and
    ! increase them for every point we find. The loop will behave like
    ! the first one, i.e. we will again find 8,6 and 4 vertices on the boundary components
    ! 1,2 and 3. Increasing the p_IboundaryCpIdx(2:NBCT+1) for boundary components
    ! 1..NBCT the same way as done in the first loop will therefore again lead to
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   9  15  19
    !
    ! Ok, let us catch the actual vertices.
    ! Check all vertices to find out, which vertices are on the boundary.
    do ivt=1,rtriangulation%NVT
      if (p_InodalProperty(ivt) .gt. 0) then
        ! id of the boundary component
        ibct = p_InodalProperty(ivt)
        
        ! set ivbd to the number of vertices on that boundary component
        ! thus ivbd holds the current number of vertices found for
        ! boundary component ibct and ivbd represents the current
        ! position in the p_IverticesAtBoundary array
        ivbd = p_IboundaryCpIdx(ibct+1)
        ! we have found a new point on that boundary component
        ! so increase the number of points by one
        p_IboundaryCpIdx(ibct+1) = ivbd+1
        
        ! Store the vertex as boundary vertex
        p_IverticesAtBoundary (ivbd) = ivt
      end if
    end do
    
  end subroutine tria_genRawBoundary1D

  !************************************************************************

!<subroutine>

  subroutine tria_genNeighboursAtElement1D(rtriangulation)

!<description>
  ! This routine generates the array IneighboursAtElement (KADJ). 
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IelementsAtVertexIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
  integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(2) :: Isize  
    integer :: iel1, iel2, ivi1, ivi2, ivt, ive
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      Isize = (/rtriangulation%NNVE,rtriangulation%NEL/)
      call storage_new ('tria_genNeighboursAtElement2D', 'KADJ', &
          Isize, ST_INT, &
          rtriangulation%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IneighboursAtElement, Isize)
      if (Isize(2) .ne. rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genNeighboursAtElement2D', &
            rtriangulation%NEL, rtriangulation%h_IneighboursAtElement, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Fill the array with 0. We overwrite only those positions <> 0.
    call storage_clear (rtriangulation%h_IneighboursAtElement)
    
    ! Get the array which is to be filled with data.
    call storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
        
    ! Get some data arrays about the vertices.
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
   
    ! Get the index array that tells us how many elements are adjacent to
    ! each vertex.
    call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
    call storage_getbase_int (rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)

    ! We have to look at all vertices to find neighbouring information.
    !
    ! Loop through all vertices.
    do ive = 1, rtriangulation%NVT
    
      ! Get both elements at this vertice
      ivi1 = p_IelementsAtVertexIdx(ive)
      ivi2 = p_IelementsAtVertexIdx(ive+1)
      
      ! If (ivi2 - ivi1) < 2, then we have no neighbour elements on
      ! this vertice
      if ((ivi2 - ivi1) .lt. 2) cycle
      
      ! Get indices of elements on this vertice
      iel1 = p_IelementsAtVertex(ivi1)
      iel2 = p_IelementsAtVertex(ivi1+1)
      
      ! Add iel2 as a neighbour to iel1
      do ivt = 1, rtriangulation%NNVE
        if (p_IverticesAtElement(ivt, iel1) .eq. ive) then
          p_IneighboursAtElement(ivt, iel1) = iel2
          exit
        end if
      end do
    
      ! Add iel1 as a neighbour to iel2
      do ivt = 1, rtriangulation%NNVE
        if (p_IverticesAtElement(ivt, iel2) .eq. ive) then
          p_IneighboursAtElement(ivt, iel2) = iel1
          exit
        end if
      end do

    end do ! ive
    
    ! That is it

  end subroutine tria_genNeighboursAtElement1D

  !************************************************************************

!<subroutine>

  subroutine tria_genElementVolume1D(rtriangulation)

!<description>
  ! This routine generates the element volume array DelementVolume (DAREA). 
  ! For this purpose, the following arrays are used:
  !    DvertexCoordinates, IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DelementVolume
    integer :: iel
    integer :: isize
    real(DP) :: dtotalVolume
    real(DP), dimension(TRIA_MAXNVE1D) :: Dpoints
    integer :: ive

    ! Is everything here we need?
    if (rtriangulation%h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line ('DvertexCoords not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume1D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement  not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume1D')
      call sys_halt()
    end if
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_DelementVolume .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementVolume1D', 'DAREA', &
          rtriangulation%NEL+1, ST_DOUBLE, &
          rtriangulation%h_DelementVolume, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_DelementVolume, isize)
      if (isize .ne. rtriangulation%NEL+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementVolume1D', &
            rtriangulation%NEL+1, rtriangulation%h_DelementVolume, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the arrays
    call storage_getbase_double2D (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double (rtriangulation%h_DelementVolume,&
        p_DelementVolume)
        
    dtotalVolume = 0.0_DP
        
    ! Calculate the element volume for all elements
    do iel=1,rtriangulation%NEL
    
      ! line element
      do ive=1,TRIA_NVELINE1D
        Dpoints(ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
      end do
      p_DelementVolume(iel) = abs(Dpoints(1) - Dpoints(2))
      
      dtotalVolume = dtotalVolume+p_DelementVolume(iel)
    end do
    
    ! Store the total volume in the last element of DelementVolume
    p_DelementVolume(rtriangulation%NEL+1) = dtotalVolume
    
  end subroutine tria_genElementVolume1D

  !====================================================================
  !
  !       ++++      ++++
  !           +     +   +  
  !          +      +    +
  !         +       +    +
  !        +        +    + 
  !       +         +   +  
  !       ++++      ++++
  ! tag@2D    
  !====================================================================
  
!<subroutine>

  subroutine tria_readTriFile2D(rtriangulation, sfilename, rboundary, &
                                bnoExtendedRaw)

!<description>
  ! This routine reads a .TRI file of a 2D triangulation into memory
  ! and creates a 'raw' triangulation (i.e. a triangulation that contains
  ! only basic information, see below). 
  !
  ! The triangulation structure rtriangulation is initialised with the data 
  ! from the file. The parameter sfilename gives the name of the .tri 
  ! file to read.
  !
  ! This reads only the very basic information that is needed to create
  ! the coarse grid. No information about boundary points, subdivisions
  ! or whatever is created.
  ! The following arrays / information tags will be initialised:
  !
  ! NEL,NVT,NMT,NBCT,NVBD,InelOfType,
  ! DvertexCoords, IverticesAtElement, InodalProperty, 
  ! IboundaryCpIdx, IverticesAtBoundary, DvertexParameterValue.
  ! The nodal property array InodalProperty contains only the data for the
  ! vertices of the triangulation.
  ! The arrays IverticesAtBoundary and DvertexParameterValue are sorted
  ! for the boundary component (according to IboundaryCpIdx) but not 
  ! for the parameter value.
  !
  ! If bnoExtendedRaw is not specified or set to .FALSE., an extended
  ! raw mesh is generated that provides also a proper numbering for edges.
  ! In that case, the following arrays are initialised as well:
  !
  ! IelementsAtVertexIdx,IelementsAtVertexIdx,IneighboursAtElement,
  ! IedgesAtElement,IfacesAtElement
  !
  ! The triangulation structure rtriangulation must be empty. All previous
  ! information in this structure (if there is any) is lost!
!</description>

!<input>
  ! The name of the .tri file to read.
  character(LEN=*), intent(in) :: sfilename

  ! OPTIONAL: An rboundary object specifying the underlying domain.
  ! If not specified, the routine assumes that the TRI file does not specify
  ! boundary parameter values, i.e. the point coordinates in the TRI file
  ! are all real coordinates. The array DvertexParameterValue is not
  ! generated in this case.
  type(t_boundary), intent(in), optional :: rboundary
  
  ! OPTIONAL: Prevent creation of an extended raw mesh. If set to .false.,
  ! an 'extended raw' mesh will be created that provides a proper numbering
  ! for edges (standard). If set to '.true', the result will be a 'really raw'
  ! raw mesh with minimum information and no numbering for edges.
  logical, intent(in), optional :: bnoExtendedRaw
! </input>
  
!<output>
  ! Triangulation structure, to be filled with data
  type(t_triangulation), intent(out) :: rtriangulation
!</output>
  
!</subroutine>

    ! Input channel for reading
    integer :: iunit
    
    ! Open the file
    call io_openFileForReading(sfilename, iunit)

    ! We create a 2D triangulation here.
    rtriangulation%ndim = NDIM2D

    ! Read the basic mesh
    call tria_readRawTriangulation2D (iunit, rtriangulation)

    ! Create the basic boundary information
    call tria_genRawBoundary2D (rtriangulation, rboundary)

    ! Extend the raw mesh by basic edge numbering,
    ! initialise an extended raw mesh.
    if (.not.present(bnoExtendedRaw)) then
      call tria_initExtendedRawMesh (rtriangulation)
    else
      if (.not.bnoExtendedRaw) then
        call tria_initExtendedRawMesh (rtriangulation)
      end if
    end if

    ! Close the file, finish
    close(iunit)

  end subroutine tria_readTriFile2D

  !************************************************************************

!<subroutine>

  subroutine tria_readRawTriangulation2D (iunit, rtriangulation)

!<description>  
  ! Auxiliary routine of tria_readTriFile2D.
  ! Reads basic information from a triangulation file into rtriangulation.
  ! That means, the following information arrays / tags are initialised:
  ! NEL,NVT,NMT,NBCT,NNEE,InelOfType,
  ! DvertexCoords, IverticesAtElement and InodalProperty.
  ! The data is read from the file without being changed!
!</description>
  
!<input>
  ! Unit number of the file to be read
  integer, intent(in) :: iunit
!</input>
  
!<inputoutput> 
  ! Triangulation to be initialised with basic data.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! Local variables
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    integer, dimension(:,:), pointer :: p_Idata2D
    integer, dimension(:), pointer :: p_Idata
    integer :: ivt, iel
    integer :: idim, ive, nve
    integer, dimension(2) :: Isize
    
    ! The first two lines in the file are comments.
    read(iunit,*)
    read(iunit,*)
    
    ! Read NEL,NVT,NMT,NVE,NBCT from the file
    ! and store this information in the structure.
    read (iunit,*) rtriangulation%NEL, rtriangulation%NVT,&
                   rtriangulation%NMT, rtriangulation%NNVE,&
                   rtriangulation%NBCT
        
    nve = rtriangulation%NNVE
    
    rtriangulation%NNEE  = rtriangulation%NNVE

    ! Comment: 'DCORVG'
    read (iunit,*)

    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM2D, NVT)
    Isize = (/NDIM2D,rtriangulation%NVT/)
    call storage_new ('tria_readRawTriangulation2D', 'DCORVG',&
        Isize, ST_DOUBLE,&
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_Ddata2D)
        
    ! Read the data from the file, store it in the array.
    ! read data into p_Ddata2D:
    ! first read nvt x-coordinates into p_Ddata2D(1,ivt)
    ! then read nvt  y-coordinates into p_Ddata2D(2,ivt)
    read (iunit,*) ((p_Ddata2D(idim,ivt),idim=1,NDIM2D), ivt=1,rtriangulation%NVT)

    ! Comment: 'KVERT'
    read (iunit,*)

    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/nve,rtriangulation%NEL/)
    call storage_new ('tria_readRawTriangulation2D', 'KVERT', Isize,&
        ST_INT, rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    call storage_getbase_int2D(rtriangulation%h_IverticesAtElement, p_Idata2D)

    ! read ive=1 indices to nve into p_Idata2D(ive,iel) where iel=1 to NEL
    read (iunit,*) ((p_Idata2D(ive,iel),ive=1,nve), iel=1,rtriangulation%NEL)

    ! Loop through the elements and determine how many elements
    ! of each element type we have.
    rtriangulation%InelOfType(:) = 0
    do iel=1,rtriangulation%NEL
      ! start at the last index of element iel down to the first
      do ive=nve,1,-1
        if (p_Idata2D(ive,iel) .ne. 0) then
          rtriangulation%InelOfType(ive) = rtriangulation%InelOfType(ive)+1
          exit
        end if
      end do
    end do

    ! Comment: 'KNPR'
    read (iunit,*)

    ! Allocate memory for InodalProperty 
    call storage_new ('tria_readRawTriangulation2D', 'KNPR', &
        rtriangulation%NVT, ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(rtriangulation%h_InodalProperty, p_Idata)

    ! Read the data
    read (iunit,*) (p_Idata(ivt), ivt=1,rtriangulation%NVT)

  end subroutine tria_readRawTriangulation2D

  !************************************************************************

!<subroutine>

  subroutine tria_genRawBoundary2D (rtriangulation, rboundary)

!<description>  
  ! Auxiliary routine of tria_readTriFile2D.
  ! This routine initialises basic boundary arrays and cleans up 
  ! a basic triangulation. That means:
  ! -> NVBD is calculated
  ! -> IboundaryCpIdx is created and generated
  ! -> IverticesAtBoundary is created and generated.
  !    The vertices are ordered for the boundary component according
  !    to IboundaryCpIdx but not ordered for their parameter value.
  ! -> If rboundary is specified,
  !    DvertexParameterValue is created and generated.
  !    The parameter values are ordered for the boundary component 
  !    according to IboundaryCpIdx but not ordered for the parameter value.
  ! -> If rboundary is specified,
  !    the parameter values of the boundary vertices are extracted
  !    from the first coordinate in DvertexCoords and put into
  !    the DvertexParameterValue array.
  ! -> If rboundary is specified,
  !    based on the parameter value of each boundary vertex, the 
  !    DvertexCoords array receives the actual coordinates of 
  !    the boundary vertices.
  !    If not specified, the routine assumes that DvertexCoords
  !    already contains the real point coordinates.
!</description>
  
!<input>
  ! OPTIONAL: The parametrisation that specifies the coordinates of the 
  ! boundary points.
  ! If specified, DvertexParameterValue is generated from DvertexCoords
  ! and the coordinates of boundary vertices are (re-)generated
  ! by DvertexParameterValue.
  type(t_boundary), intent(in), optional :: rboundary
!</input>
  
!<inputoutput>
  ! Triangulation to be initialised with basic data.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer :: ivbd,ivt
    integer :: ibct
    integer, dimension(:), pointer :: p_InodalProperty

    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty, p_InodalProperty)

    ! Calculate NVBD by simply counting how many elements
    ! in p_InodalProperty are <> 0.
    rtriangulation%NVBD = 0
    ivbd = 0
    do ivt=1,rtriangulation%NVT
      !IF (p_InodalProperty(ivt) .NE. 0) rtriangulation%NVBD = rtriangulation%NVBD+1
      if (p_InodalProperty(ivt) .gt. 0) ivbd = ivbd+1
    end do
    rtriangulation%NVBD = ivbd

    ! Allocate memory for IverticesAtBoundary.
    call storage_new ('tria_genRawBoundary2D', 'KVBD', &
        rtriangulation%NVBD, ST_INT, &
        rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    call storage_new ('tria_genRawBoundary2D', 'KBCT', &
        rtriangulation%NBCT+1, ST_INT, &
        rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    call storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    call storage_getbase_double2D (&
        rtriangulation%h_DvertexCoords,p_DvertexCoords)
        
    call storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! The first element in p_IboundaryCpIdx is (as the head) always =1.
    p_IboundaryCpIdx(1) = 1

    ! Perform a first loop over all vertices to check which are on the
    ! boundary. Count them. This way, we at first set up the boundary
    ! component index vector p_IboundaryCpIdx.
    ! In this first step, we save the number of vertices on each 
    ! boundary component in p_IboundaryCpIdx(2:NBCT+1)!
    do ivt=1,rtriangulation%NVT
      if (p_InodalProperty(ivt) .gt. 0) then
        ibct = p_InodalProperty(ivt)
        
        ! Increase the number of vertices in that boundary component by 1.
        ! The number of vertices on boundary component i is saved here
        ! at p_IboundaryCpIdx(i+1) for later.
        ! Note that the array was initialised with zero during the creation
        ! process!
        p_IboundaryCpIdx(ibct+1) = p_IboundaryCpIdx(ibct+1)+1
      end if
    end do
    
    ! Sum up the number of vertices on each boundary component to get the
    ! actual index vector.
    do ibct = 2, rtriangulation%NBCT+1
      p_IboundaryCpIdx(ibct) = p_IboundaryCpIdx(ibct)+p_IboundaryCpIdx(ibct-1)
    end do
    
    ! Shift the p_IboundaryCpIdx array by one position. That is a little trick in
    ! the use of p_IboundaryCpIdx!
    ! Imagine, we have 3 boundary components with 8,6 and 4 edges. 
    ! Before summing the entries up, p_IboundaryCpIdx may have looked like this:
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   8   6   4
    !
    ! Summing that up gave the actual indices:
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   9  15  19
    !
    ! which is already the correct index array.
    ! Shifting p_IboundaryCpIdx by 1 will give:
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   1   9  15
    
    p_IboundaryCpIdx(2:rtriangulation%NBCT+1) = p_IboundaryCpIdx(1:rtriangulation%NBCT)
    
    ! Then, we again loop through all vertices and collect those on the
    ! boundary. In that loop, we use p_IboundaryCpIdx(2:NBCT+1) as pointer and
    ! increase them for every point we find. The loop will behave like
    ! the first one, i.e. we will again find 8,6 and 4 vertices on the boundary components
    ! 1,2 and 3. Increasing the p_IboundaryCpIdx(2:NBCT+1) for boundary components
    ! 1..NBCT the same way as done in the first loop will therefore again lead to
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   9  15  19
    !
    ! Ok, let us catch the actual vertices.
    !
    ! The loop must be slightly modified if rboundary is not present!
    if (present(rboundary)) then

      ! Allocate memory for  and DvertexParameterValue
      call storage_new ('tria_genRawBoundary2D', &
          'DVBDP', rtriangulation%NVBD, &
          ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
      
      ! Get the array where to store boundary parameter values.
      call storage_getbase_double (&
          rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
          
      ! Check all vertices to find out, which vertices are on the boundary.
      do ivt=1,rtriangulation%NVT
        if (p_InodalProperty(ivt) .gt. 0) then
          ibct = p_InodalProperty(ivt)
          
          ! Create a new point on that boundary component 
          ! and get the number, the point will have.
          ! Note that the array was initialised with zero during the creation
          ! process!
          ivbd = p_IboundaryCpIdx(ibct+1)
          p_IboundaryCpIdx(ibct+1) = ivbd+1
          
          ! Store the vertex as boundary vertex
          p_IverticesAtBoundary (ivbd) = ivt
          
          ! Store the parameter value; it is saved in DvertexCoords(1,.)
          p_DvertexParameterValue (ivbd) = p_DvertexCoords(1,ivt)
          
          ! Replace the coordinates in DvertexCoords by those
          ! given by the parametrisation.
          call boundary_getCoords(rboundary, ibct, p_DvertexParameterValue (ivbd), &
              p_DvertexCoords(1,ivt), p_DvertexCoords(2,ivt))
          
        end if
      end do
      
    else
    
      ! No parametrisation available, the array with boundary parameter values 
      ! is not generaterd.
      !
      ! Check all vertices to find out, which vertices are on the boundary.
      do ivt=1,rtriangulation%NVT
        if (p_InodalProperty(ivt) .gt. 0) then
          ! id of the boundary component
          ibct = p_InodalProperty(ivt)
          
          ! set ivbd to the number of vertices on that boundary component
          ! thus ivbd holds the current number of vertices found for
          ! boundary component ibct and ivbd represents the current
          ! position in the p_IverticesAtBoundary array
          ivbd = p_IboundaryCpIdx(ibct+1)

          ! we have found a new point on that boundary component
          ! so increase the number of points by one
          p_IboundaryCpIdx(ibct+1) = ivbd+1
          
          ! Store the vertex as boundary vertex
          p_IverticesAtBoundary (ivbd) = ivt
        end if
      end do
      
    end if
    
  end subroutine tria_genRawBoundary2D

  ! ***************************************************************************

!<subroutine>    

  subroutine tria_rawGridToTri (rtriangulation)
  
!<description>
  ! This routine converts a 2D 'raw' mesh into a triangular 2D mesh.
  ! All elements are converted to triangles.
  ! Warning: This routine can only applied to a 'raw' mesh, e.g. a mesh which 
  !  comes directly from a .TRI file. It is not advisable to apply this routine to
  !  a 'standard' mesh that contains already further mesh information 
  !  (adjacencies,...), as any information about possible two-level ordering
  !  is lost!
  !  However, when applied to a 'standard' mesh, the mesh information can be
  !  updated using tria_initStandardMeshFromRaw to form a valid standard mesh
  !  again.
!</description>

!<inputoutput>
  ! The triangulation structure to be converted.
  ! Is replaced by a triangular mesh.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>

    integer :: icount
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: h_IverticesAtElementTri
    integer, dimension(:,:), pointer :: p_IverticesAtElementTri
    integer, dimension(2) :: Isize
   
    ! For this routine we currently assume that there are only triangles
    ! and quads in the triangulation. Might be a matter of change in 
    ! the future...
   
    ! Get the points-at-element array
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    
    icount = rtriangulation%InelOfType(TRIA_NVEQUAD2D)

    ! Check if quadrilaterals exists at all
    if (icount .eq. 0) return
    
    ! Create a new p_IverticesAtElement array for the triangular mesh.
    Isize = (/TRIA_NVETRI2D,icount+rtriangulation%NEL/)
    call storage_new ('tria_quadToTri', 'KVERTTRI', Isize, ST_INT, &
        h_IverticesAtElementTri,ST_NEWBLOCK_NOINIT)
    call storage_getbase_int2d (h_IverticesAtElementTri,p_IverticesAtElementTri)
    
    ! Convert the array
    call convert_QuadToTria (&
        rtriangulation%NEL, icount,&
        p_IverticesAtElement, p_IverticesAtElementTri)
    
    ! Replace the old array by the new, release the old one.
    call storage_free (rtriangulation%h_IverticesAtElement)
    rtriangulation%h_IverticesAtElement = h_IverticesAtElementTri
    
    ! Finally, set up NEL and InelOfType.
    ! Every quad got two triangles, so the number of elements increases by the number
    ! of quads!
    rtriangulation%NEL = rtriangulation%NEL + icount
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(TRIA_NVETRI2D) = rtriangulation%NEL
    rtriangulation%nnve = 3
    rtriangulation%NNEE = 3
    
    ! If the mesh is an extended raw mesh, regenerate extended raw mesh
    ! information.
    if (rtriangulation%h_IelementsAtVertexIdx .ne. ST_NOHANDLE) then
      ! Remove the old extended raw mesh arrays. Recreate everything.
      call tria_resetToRaw(rtriangulation,.false.)
      call tria_initExtendedRawMesh (rtriangulation)    
    end if
    
    ! That is it.

  contains

  !<subroutine>

    subroutine convert_QuadToTria (nel, nelquad, IverticesAtElement, Kvert_triang)

  !<description>
    ! Purpose: Convert mixed quad/tri mesh to triangular mesh
    !
    ! This routine creates a triangular KVERT structure from a 
    ! quadrilateral KVERT structure.
  !</description>

  !<input>
    ! Number of elements in the old grid
    integer,intent(in) :: nel

    ! Number of quads in the old grid
    integer,intent(in) :: nelquad

    ! KVERT structure of the old mesh
    ! array [1..4,1..nel] of integer
    integer, dimension(:,:), intent(in) :: IverticesAtElement
  !</input>

  !<output>
    ! KVERT structure of the tri mesh
    ! array [1..4,1..nel+nelquad] of integer
    integer, dimension(:,:), intent(out)  :: Kvert_triang
  !</output>

  !</subroutine>
  
      ! local variables
      integer :: i,j
        
      j = nel
      do i=1,nel
        ! Copy the first three entries of each IverticesAtElement subarray
        ! to the first half of Kvert_triang. They form NEL quads.
        Kvert_triang(1,i) = IverticesAtElement(1,i)
        Kvert_triang(2,i) = IverticesAtElement(2,i)
        Kvert_triang(3,i) = IverticesAtElement(3,i)
        
        ! For every quad we find, we produce a second triangle with triangle
        ! number NEL+1,...
        if (IverticesAtElement(4,i) .ne. 0) then
        
          ! Get the next free element number behind the first set of triangles
          j = j+1
          
          ! The second triangle in each quad consists of vertices 1,3,4.
          Kvert_triang(1,j) = IverticesAtElement(1,i)
          Kvert_triang(2,j) = IverticesAtElement(3,i)
          Kvert_triang(3,j) = IverticesAtElement(4,i)
        
        end if
      end do
      
      ! That is it.
      
    end subroutine convert_QuadToTria
    
  end subroutine tria_rawGridToTri

  !************************************************************************

!<subroutine>

  subroutine tria_genElementsAtVertex1D2D(rtriangulation)

!<description>
  ! This routine generates NNelAtVertex and the array IelementsAtVertex.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement.
  ! If necessary, new memory is allocated.
  !
  ! This routine works for 1D and 2D meshes.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer :: ive,nnve
    integer :: ivt,isize,isize2
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_Iaux1
    integer :: haux1,iel

    ! Do we have (enough) memory for the index array?
    if (rtriangulation%h_IelementsAtVertexIdx .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtVertex1D2D', 'IelementsAtVertexIdx', &
          rtriangulation%NVT+1, ST_INT, &
          rtriangulation%h_IelementsAtVertexIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IelementsAtVertexIdx, isize)
      if (isize .ne. rtriangulation%NVT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtVertex1D2D', &
            rtriangulation%NVT+1, rtriangulation%h_IelementsAtVertexIdx, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the index array.
    call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
        
    ! Get some data arrays about the vertices.
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)

    ! Fill the index array with zero.
    call storage_clear (rtriangulation%h_IelementsAtVertexIdx)
    
    nnve = rtriangulation%NNVE

    ! We create the index array in two steps. In the first step,
    ! we loop over the elements to find out, how many elements
    ! meet at each vertex. We store this information in the index
    ! array at position 2..NVT+1 (i.e. shifted by one) for later.
    
    do iel = 1,rtriangulation%NEL
      do ive = 1,nnve
        ivt = p_IverticesAtElement(ive,iel)
        ! Cancel that element if we reached the end. Might happen if there
        ! are triangles in a quad mesh e.g.
        if (ivt .eq. 0) exit

        p_IelementsAtVertexIdx(ivt+1)=p_IelementsAtVertexIdx(ivt+1)+1
      end do
    end do
    
    ! Set the first element in p_IverticesAtElement to 1. Then, sum up
    ! all the length information to get the index array.
    ! Simultaneously calculate NNelAtVertex.
    p_IelementsAtVertexIdx(1) = 1
    rtriangulation%NNelAtVertex = 0
    
    do ivt = 2,rtriangulation%nvt+1
      rtriangulation%NNelAtVertex = &
          max(rtriangulation%NNelAtVertex,p_IelementsAtVertexIdx(ivt))
          
      p_IelementsAtVertexIdx(ivt) = &
          p_IelementsAtVertexIdx(ivt) + p_IelementsAtVertexIdx(ivt-1)
    end do
    
    isize = p_IelementsAtVertexIdx(rtriangulation%NVT+1)-1

    ! isize contains now the length of the array where we store the adjacency
    ! information (IelementsAtVertex).
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IelementsAtVertex .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtVertex1D2D', 'IelementsAtVertex', &
          isize, ST_INT, &
          rtriangulation%h_IelementsAtVertex, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IelementsAtVertex, isize2)
      if (isize .ne. isize2) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtVertex1D2D', &
            isize, rtriangulation%h_IelementsAtVertex, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the array.
    call storage_getbase_int (rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)

    ! Duplicate the p_IelementsAtVertexIdx array. We use that as pointer and index
    ! when new elements at a vertex are found.
    haux1 = ST_NOHANDLE
    call storage_copy (rtriangulation%h_IelementsAtVertexIdx,haux1)
    call storage_getbase_int (haux1,p_Iaux1)

    ! Now, we perform a second loop over the elements and the vertices.
    ! This time, we fetch the element number and write it to the    
    ! p_IelementsAtVertex array.
    ! p_Iaux1 counts how many positions in p_IelementsAtVertex are occupied.
    do iel = 1,rtriangulation%NEL
      do ive = 1,nnve
      
        ivt = p_IverticesAtElement(ive,iel)
        ! Cancel that element if we reached the end. Might happen if there
        ! are triangles in a quad mesh e.g.
        if (ivt .eq. 0) exit
        
        ! Remember the element number and increase the pointer in the
        ! elements-adjacent-to-that-vertex list.
        p_IelementsAtVertex( p_Iaux1(ivt) ) = iel
        p_Iaux1(ivt) = p_Iaux1(ivt)+1
      end do
    end do
    
    ! Release the auxiliary array, that is it.
    call storage_free (haux1)

  end subroutine tria_genElementsAtVertex1D2D

  !************************************************************************

!<subroutine>

  subroutine tria_genNeighboursAtElement2D(rtriangulation)

!<description>
  ! This routine generates the array IneighboursAtElement (KADJ). 
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IelementsAtVertexIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer :: ive,nve,nnve
    integer, dimension(2) :: Isize
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
    integer :: iel
    integer :: ivt,ivtneighbour
    
    integer :: haux1, haux2
    integer, dimension(:,:), pointer :: p_IedgeAtVertex
    integer :: iidxEdge, iedge, iedgeneighbour
    integer, dimension(:), pointer :: p_IedgeIdx

    nnve = rtriangulation%NNVE

    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      Isize = (/nnve,rtriangulation%NEL/)
      call storage_new ('tria_genNeighboursAtElement2D', 'KADJ', &
          Isize, ST_INT, &
          rtriangulation%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IneighboursAtElement, Isize)
      if (Isize(2) .ne. rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genNeighboursAtElement2D', &
            rtriangulation%NEL, rtriangulation%h_IneighboursAtElement, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Fill the array with 0. We overwrite only those positions <> 0.
    call storage_clear (rtriangulation%h_IneighboursAtElement)
    
    ! Get the array which is to be filled with data.
    call storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
        
    ! Get some data arrays about the vertices.
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
        
    ! In the following, we create an array IedgeAtVertex that saves
    ! information about the edges adjacent to each element in 
    ! counterclockwise sense. The array has the same length as IelementsAtVertex 
    ! and saves information about which edges are adjacent to a vertex:
    !  IedgeAtVertex(1) = vertex
    !  IedgeAtVertex(2) = next neighbour of the vertex in counterclockwise sense
    !  IedgeAtVertex(3) = number of the element from which we look at that edge
    !  IedgeAtVertex(4) = local number of the vertex in the element
    !
    ! The array IedgeIdx remembers which positions in the IedgeAtVertex
    ! array are occupied. Each index starts from `one behind possible
    ! free memory locations` and is decreased for every edge found.
    ! Note that for every element, there is exactly one edge following a vertex
    ! in counterclockwise sense!
    !
    !     +---1---+      +---+---1
    !     | 1 | 4 |      |   | 1 |
    !     2---X---4  or  +---2---X
    !     | 2 | 3 |      |   | 2 |
    !     +---3---+      +---+---+

    ! Get the index array that tells us how many elements are adjacent to
    ! each vertex. Make a copy of that array.
    haux2 = ST_NOHANDLE
    call storage_copy (rtriangulation%h_IelementsAtVertexIdx,haux2)
    call storage_getbase_int (haux2,p_IedgeIdx)
    
    call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
    
    ! Actually, we need the entries 2..NVT+1 of Iaux2
    p_IedgeIdx => p_IedgeIdx(2:)

    ! Create the auxiliary array for edge information.
    Isize = (/4,p_IedgeIdx(size(p_IedgeIdx))-1/)
    call storage_new ('tria_genNeighboursAtElement2D', 'edgeAtVertex', &
        Isize, ST_INT, haux1, ST_NEWBLOCK_ZERO)
    call storage_getbase_int2d (haux1,p_IedgeAtVertex)
    
    ! Loop through the elements to calculate edge information.
    do iel=1,rtriangulation%NEL
     
      ! Get the number of vertices of that element.
      nve = nnve
      do while (p_IverticesAtElement(nve,iel) .eq. 0)
        nve = nve-1
      end do
    
      ! Loop through the vertices
      do ive=1,nve
        ! Number of current vertex?
        ivt = p_IverticesAtElement(ive,iel)
        
        ! What`s the neighbour?
        ivtneighbour = p_IverticesAtElement(mod(ive,nve)+1,iel)
        
        ! Get the next free entry in the array with the edge information
        ! p_IelementsAtVertexIdx
        iidxEdge = p_IedgeIdx(ivt)-1
        
        ! Save information about the edge we look at from our current element.
        !
        !   +---1 IVTNEIGHBOUR
        !   | 1 |
        !   +---X IVT
        !
        p_IedgeAtVertex(1,iidxEdge) = ivt
        p_IedgeAtVertex(2,iidxEdge) = ivtneighbour
        p_IedgeAtVertex(3,iidxEdge) = iel
        p_IedgeAtVertex(4,iidxEdge) = ive
        
        ! Remember that we occupied one entry in p_IedgeAtVertex
        p_IedgeIdx(ivt) = iidxEdge
      
      end do ! ive
      
    end do ! iel
    
    ! Next, we have to look at all edges to find neighbouring information.
    !
    ! Loop through all edges adjacent to all vertices.
    do iedge = 1,ubound(p_IedgeAtVertex,2)
    
      ! Get the edge information
      ivt          = p_IedgeAtVertex(1,iedge) 
      ivtneighbour = p_IedgeAtVertex(2,iedge) 
      iel          = p_IedgeAtVertex(3,iedge) 
    
      ! Now, loop through all edges adjacent to the neighbour vertex
      ! ivtneighbour to find possible adjacent elements at the current
      ! edge.
      ! p_IedgeIdx(ivtneighbour) is the number of the first entry
      ! in p_IedgeAtVertex that is occupied!
      do iedgeneighbour = p_IedgeIdx(ivtneighbour), &
                          p_IelementsAtVertexIdx(ivtneighbour+1)-1
        
        ! When the vertices of that edge coincide (i.e. if the counterclockwise
        ! neighbour of the neighbour vertex is the current vertex)
        ! we found an adjacency!
        !
        !   +---+
        !   | 2 |
        !  1X<->X2
        !   | 1 |
        !   +---+
        
        if (p_IedgeAtVertex(2,iedgeneighbour) .eq. ivt) then
        
          ! Save the adjacency information to our current element.
          ! We do not save it for the neighbour element (although we have it here)
          ! since we arrive at the neighbour element, too -- which would mean
          ! to store every information twice...
          p_IneighboursAtElement(p_IedgeAtVertex(4,iedge),iel) = &
            p_IedgeAtVertex(3,iedgeneighbour)
        
        end if
    
      end do ! iedgeneighbour
    
    end do ! iedge
    
    ! Release memory, finish.
    call storage_free (haux1)
    call storage_free (haux2)

  end subroutine tria_genNeighboursAtElement2D

  !************************************************************************

!<subroutine>

  subroutine tria_genEdgesAtElement2D(rtriangulation)

!<description>
  ! This routine calculates the edge numbering. That means,
  ! it generates information about the edges adjacent to
  ! each element IedgesAtElement (KMID) and calculates the correct NMT.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer :: ive, iveneighbour
    integer :: iel, ielneighbour
    integer :: iedge
    integer, dimension(2) :: Isize

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call storage_getsize (rtriangulation%h_IneighboursAtElement, Isize)
      call storage_new ('tria_genEdgesAtElement2D', 'KMID', &
          Isize, ST_INT, &
          rtriangulation%h_IedgesAtElement, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IedgesAtElement, Isize)
      if (Isize(2) .ne. rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genEdgesAtElement2D', &
            rtriangulation%NEL, rtriangulation%h_IedgesAtElement, &
            ST_NEWBLOCK_NOINIT, .false.)
        Isize(2) = rtriangulation%NEL
      end if
    end if
    
    ! Fill IedgesAtElement with 0. That is important in case some
    ! elements in the array are not tackled when searching for edges
    ! (e.g. in meshes where triangles and quads are mixed).
    call storage_clear (rtriangulation%h_IedgesAtElement)

    call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    
    ! iedge counts the edges and specifies the last given edge number.
    iedge = 0
    
    ! Loop through all elements
    do iel = 1,Isize(2)
    
      ! Loop through all edges on each element
      do ive = 1,Isize(1)
      
        ! Stop if we handled all edges; this is important if there are triangles
        ! in a quad mesh e.g.
        if (p_IverticesAtElement(ive,iel) .eq. 0) exit
        
        ! Check the neighbour element.
        ! If the neightbour element has number =0 (no neighbour) or a number
        ! greater than IEL, we found the edge the first time and give it a number.
        if ((p_IneighboursAtElement(ive,iel) .eq. 0) .or. &
            (p_IneighboursAtElement(ive,iel) .gt. iel)) then
        
          iedge = iedge + 1
          
          ! Add NVT to iedge to get the edge number
          p_IedgesAtElement(ive,iel) = iedge
        
        else
        
          ! Otherweise, we had that edge already. Look into the neighbour element
          ! (which definitely exists, p_IneighboursAtElement cannot be =0 there!)
          ! to find the edge number.
          ielneighbour = p_IneighboursAtElement(ive,iel)
          
          do iveneighbour = 1,Isize(1)
            if (p_IneighboursAtElement(iveneighbour,ielneighbour) .eq. iel) then
              p_IedgesAtElement(ive,iel) = &
                  p_IedgesAtElement(iveneighbour,ielneighbour)
              exit
            end if
          end do
        
        end if
        
      end do
    
    end do
    
    ! Save the correct NMT.
    rtriangulation%NMT = iedge
    
  end subroutine tria_genEdgesAtElement2D

  !************************************************************************

!<subroutine>

  subroutine tria_genElementVolume2D(rtriangulation)

!<description>
  ! This routine generates the element volume array DelementVolume (DAREA). 
  ! For this purpose, the following arrays are used:
  !    DvertexCoordinates, IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DelementVolume
    integer :: iel
    integer :: isize
    real(DP) :: dtotalVolume
    real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: Dpoints
    integer :: ive

    ! Is everything here we need?
    if (rtriangulation%h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line ('DvertexCoords not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement  not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      call sys_halt()
    end if
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_DelementVolume .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementVolume2D', 'DAREA', &
          rtriangulation%NEL+1, ST_DOUBLE, &
          rtriangulation%h_DelementVolume, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_DelementVolume, isize)
      if (isize .ne. rtriangulation%NEL+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementVolume2D', &
            rtriangulation%NEL+1, rtriangulation%h_DelementVolume, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the arrays
    call storage_getbase_double2D (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double (rtriangulation%h_DelementVolume,&
        p_DelementVolume)
        
    dtotalVolume = 0.0_DP
        
    ! Currently, we support triangules and quads.
    if (ubound(p_IverticesAtElement,1) .eq. TRIA_NVETRI2D) then

      ! Calculate the element volume for all elements
      do iel=1,rtriangulation%NEL
        ! triangular element
        do ive=1,TRIA_NVETRI2D
          Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
          Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        end do
        p_DelementVolume(iel) = gaux_getArea_tria2D(Dpoints)
        
        dtotalVolume = dtotalVolume+p_DelementVolume(iel)
      end do
    
    else

      ! Calculate the element volume for all elements
      do iel=1,rtriangulation%NEL
      
        if (p_IverticesAtElement(4,iel) .eq. 0) then
          ! triangular element
          do ive=1,TRIA_NVETRI2D
            Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
            Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
          end do
          p_DelementVolume(iel) = gaux_getArea_tria2D(Dpoints)
        else
          ! quad element
          do ive=1,TRIA_NVEQUAD2D
            Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
            Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
          end do
          p_DelementVolume(iel) = gaux_getArea_quad2D(Dpoints)
        end if
        
        dtotalVolume = dtotalVolume+p_DelementVolume(iel)
      end do
      
    end if
    
    ! Store the total volume in the last element of DelementVolume
    p_DelementVolume(rtriangulation%NEL+1) = dtotalVolume
    
  end subroutine tria_genElementVolume2D

  !************************************************************************

!<subroutine>

  subroutine tria_sortBoundaryVertices1D2D (rtriangulation)

!<description>  
  ! This routine sorts the arrays IverticesAtBoundary and 
  ! DvertexParameterValue such that they are ordered for increasing
  ! parameter values of the vertices.
!</description>
  
!<inputoutput>
  ! The triangulation.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer :: istart,iend
    integer :: ibct,ivbd
    integer :: hresort
    integer, dimension(:), pointer :: p_Iresort

    if (rtriangulation%h_DvertexParameterValue .eq. ST_NOHANDLE) then
      ! We cannot sort the boundary vertices without parameter values!
      return
    end if

    ! Get pointers to the arrays
    call storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    call storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
        
    call storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Sort the p_DvertexCoords (sub)-arrays for increasing
    ! parameter value. 
    call storage_new ('tria_sortBoundaryVertices1D2D', &
        'resort', rtriangulation%NVBD, &
        ST_INT, hresort, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (hresort,p_Iresort)
    
    ! Fill p_Iresort with 1,2,3,...
    do ivbd=1,rtriangulation%NVBD
      p_Iresort(ivbd) = ivbd
    end do
    
    ! For each boundary component on the physical boundary, call the sorting routine 
    ! and calculate the mapping how the entries are sorted.
    do ibct = 1,rtriangulation%NBCT
      istart = p_IboundaryCpIdx(ibct)
      iend = p_IboundaryCpIdx(ibct+1)-1
      
      ! Sort the sub-array of boundary component ibct.
      ! Remember, how the parameter values are resorted when sorting the array.
      call sort_dp(p_DvertexParameterValue(istart:iend),SORT_QUICK,&
          p_Iresort(istart:iend))
      
    end do

    ! Resort the vertices according to the calculated mapping.
    ! Afterwards, the vertices are ordered according to the increasing
    ! parameter value.
    do ivbd=1,size(p_IverticesAtBoundary)
      p_Iresort(ivbd) = p_IverticesAtBoundary(p_Iresort(ivbd))
    end do
    
    do ivbd=1,size(p_IverticesAtBoundary)
      p_IverticesAtBoundary(ivbd) = p_Iresort(ivbd)
    end do
    
    call storage_free (hresort)

  end subroutine tria_sortBoundaryVertices1D2D
  
  !************************************************************************

!<subroutine>

  subroutine tria_genElementsAtBoundary1D2D(rtriangulation)

!<description>
  ! This routine generates information about the elements at the boundary
  ! IelementsAtBoundary (KEBD). 
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement, IverticesAtBoundary,
  !    IboundaryCpIdx, IelementsAtVertexIdx, IelementsAtVertex.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer :: ive, ibct, ivbd, iadjElement, nnve
    integer :: iel
    integer :: isize, ivt

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary1D2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary1D2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary1D2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IboundaryCpIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary1D2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtVertexIdx .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtVertexIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary1D2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtVertex .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtVertex not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary1D2D')
      call sys_halt()
    end if

    ! Get the arrays.
    ! get the vertices at the element
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    ! get the neighbours at the element
    call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    ! get the vertices on the boundary
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    ! get the index array for the boundary component
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    ! get the elements at a vertex
    call storage_getbase_int (rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)
    ! get the elements at vertex index array
    call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IelementsAtBoundary .eq. ST_NOHANDLE) then
      ! We have as many elements on the boundary as vertices!
      ! Initialise the array with zero to give proper values for elements on a
      ! possible 'blind' boundary component -- as we do not calculate elements
      ! adjacent to vertices on the 'blind' BC!
      call storage_new ('tria_genElementsAtBoundary1D2D', 'KEBD', &
          rtriangulation%NVBD, ST_INT, &
          rtriangulation%h_IelementsAtBoundary, ST_NEWBLOCK_ZERO)
    else
      call storage_getsize (rtriangulation%h_IelementsAtBoundary, isize)
      if (isize .ne. rtriangulation%NVBD) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtBoundary1D2D', &
            rtriangulation%NVBD, rtriangulation%h_IelementsAtBoundary, &
            ST_NEWBLOCK_ZERO, .false.)
      end if
    end if
    
    call storage_getbase_int (rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)

    nnve = rtriangulation%NNVE

    ! Loop through all boundary components
    do ibct = 1,rtriangulation%NBCT+rtriangulation%NblindBCT
    
      ! On each boundary component, loop through all vertices
      do ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
      
        !   +---+---+---+
        !   |   |   |   |
        !   +---+---+---+
        !   |   | 1 | 2 |
        !   +---+---X---+
        !          ivt
      
        ! Get the current boundary vertex
        ivt = p_IverticesAtBoundary(ivbd)
      
        ! Loop through all elements adjacent to the current vertex
        do iadjElement = p_IelementsAtVertexIdx(ivt),p_IelementsAtVertexIdx(ivt+1)-1
        
          ! Get the element number
          iel = p_IelementsAtVertex(iadjElement)
        
          ! Find the local number of the vertex in the element
          do ive=1,nnve
            if (p_IverticesAtElement (ive,iel) .eq. ivt) exit
          end do
          
          ! Test if the element has a neighbour at the edge that is starting
          ! with our current vertex
          if (p_IneighboursAtElement(ive,iel) .eq. 0) then
          
            ! Yes, that is the boundary edge we are searching for!
            ! It starts with ivt and is present in the boundary component
            ! we are currently processing.
            ! So we can save the boundary element number and proceed with the next 
            ! boundary vertex.
            
            p_IelementsAtBoundary (ivbd) = p_IelementsAtVertex(iadjElement)
            
            exit
           
          end if
        
        end do
      
      end do
    
    end do
    
  end subroutine tria_genElementsAtBoundary1D2D

  !************************************************************************

!<subroutine>

  subroutine tria_genBoundaryVertexPos1D2D(rtriangulation)

!<description>
  ! This routine generates the array IboundaryVertexPos which is used
  ! to search for the position of a boundary vertex.
  ! For this purpose, the following arrays are used:
  !    IverticesAtBoundary, IboundaryCpIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:,:), pointer :: p_IboundaryVertexPos
    integer :: ivbd, ibct
    integer, dimension(2) :: Isize

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genBoundaryVertexPos1D2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IboundaryCpIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genBoundaryVertexPos1D2D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IboundaryVertexPos .eq. ST_NOHANDLE) then
      ! We have as many elements on the boundary as vertices!
      Isize = (/2,rtriangulation%NVBD/)
      call storage_new ('tria_genBoundaryVertexPos1D2D', 'KVBDI', &
          Isize, ST_INT, &
          rtriangulation%h_IboundaryVertexPos, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IboundaryVertexPos, Isize)
      if (Isize(2) .ne. rtriangulation%NVBD) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genBoundaryVertexPos1D2D', &
            rtriangulation%NVBD, rtriangulation%h_IboundaryVertexPos, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    call storage_getbase_int2D (rtriangulation%h_IboundaryVertexPos,&
        p_IboundaryVertexPos)

    ! Fill the array.
    ! Store the IverticesAtBoundary in the first entry and the index in the 2nd.
    do ivbd = 1,rtriangulation%NVBD
      p_IboundaryVertexPos(1,ivbd) = p_IverticesAtBoundary(ivbd)
      p_IboundaryVertexPos(2,ivbd) = ivbd
    end do
    
    ! Sort the array -- inside of each boundary component.
    ! Use the vertex number as key.
    ! We only sort the first NBCT BC`s. The boundary component NBCT+1 (if it exists)
    ! collects all vertices that are not on the physical boundary but stem from
    ! the creation of subdomains. These vertices have no parameter value.
    do ibct = 1,rtriangulation%NBCT
      call arraySort_sortByIndex_int (&
          p_IboundaryVertexPos(:,p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),1)
    end do
    
  end subroutine tria_genBoundaryVertexPos1D2D

  !************************************************************************

!<subroutine>

  subroutine tria_genElementsAtEdge2D(rtriangulation)

!<description>
  ! This routine generates information about the elements adjacent to each 
  ! edge IelementsAtEdge (KMID) and NNelAtEdge. 
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement.
  ! NMT must be set up correctly.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IelementsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer :: ive
    integer :: iel
    integer :: iedge
    integer, dimension(2) :: Isize

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      call sys_halt()
    end if

    if (rtriangulation%NMT .eq. 0) then
      call output_line ('Edge information (NMT) not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    
    ! in 3D this is a list because there no fixed number
    ! of elements at an edge
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IelementsAtEdge .eq. ST_NOHANDLE) then
      Isize = (/2,rtriangulation%NMT/)
      call storage_new ('tria_genElementsAtEdge2D', 'KMID', &
          Isize, ST_INT, &
          rtriangulation%h_IelementsAtEdge, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IelementsAtEdge, Isize)
      if (Isize(2) .ne. rtriangulation%NMT) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtEdge2D', &
            rtriangulation%NMT, rtriangulation%h_IelementsAtEdge, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! We have at most 2 elements per edge.
    rtriangulation%NNelAtEdge = 2
    
    call storage_getbase_int2D (rtriangulation%h_IelementsAtEdge,&
        p_IelementsAtEdge)
    
    ! Loop through all elements and all edges on the elements
    ! all elements
    do iel = 1,ubound(p_IedgesAtElement,2)
      ! loop 1-4
      do ive = 1,ubound(p_IedgesAtElement,1)

        ! Is there a neighbour?
        if (p_IneighboursAtElement(ive,iel) .eq. 0) then
          ! No. Store the element as only adjacent one to that edge.
          !
          ! Do not do anything if we are looking at a nonexisting edge here!
          ! (i.e. edge 4 of a tri element in a quad mesh -- tri elements
          !  do not have 4 edges :-) )
          if (p_IedgesAtElement(ive,iel) .ne. 0) then

            iedge = p_IedgesAtElement (ive,iel)
            p_IelementsAtEdge(1,iedge) = iel
            p_IelementsAtEdge(2,iedge) = 0
            
          end if
          
        elseif (p_IneighboursAtElement(ive,iel) .lt. iel) then
        
          ! There is a neighbour and it has a smaller number than the current element --
          ! so we have not had that edge! Store the two adjacent elements.
        
          iedge = p_IedgesAtElement (ive,iel)
          p_IelementsAtEdge(1,iedge) = iel
          p_IelementsAtEdge(2,iedge) = p_IneighboursAtElement(ive,iel)
        
        end if
      
      end do
      
    end do
    
  end subroutine tria_genElementsAtEdge2D

  !************************************************************************

!<subroutine>

  subroutine tria_genVerticesAtEdge2D(rtriangulation)

!<description>
  ! This routine generates information about the vertices adjacent to 
  ! each edge IverticesAtEdge.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer :: ive,nnve
    integer :: ivtneighbour
    integer :: iel
    integer :: iedge
    integer, dimension(2) :: Isize

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) then
      Isize = (/2,rtriangulation%NMT/)
      call storage_new ('tria_genVerticesAtEdge2D', 'KMID', &
          Isize, ST_INT, &
          rtriangulation%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IverticesAtEdge, Isize)
      if (Isize(2) .ne. rtriangulation%NMT) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genVerticesAtEdge2D', &
            rtriangulation%NMT, rtriangulation%h_IverticesAtEdge, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    call storage_getbase_int2D (rtriangulation%h_IverticesAtEdge,&
        p_IverticesAtEdge)
    
    nnve = ubound(p_IedgesAtElement,1)
    
    ! Loop through all elements and all edges on the elements
    do iel = 1,ubound(p_IedgesAtElement,2)
      
      do ive = 1,ubound(p_IedgesAtElement,1)
      
        ! Stop if we handled all edges; this is important if there are triangles
        ! in a quad mesh e.g.
        if (p_IverticesAtElement(ive,iel) .eq. 0) exit
        
        ! Is there a neighbour which number is less than iel? Or even =0?
        ! If yes, we did not tackle the edge.
          
        if (p_IneighboursAtElement(ive,iel) .lt. iel) then

          ! Save the numbers of the adjacent vertices.
          iedge = p_IedgesAtElement (ive,iel)
          
          p_IverticesAtEdge(1,iedge) = p_IverticesAtElement(ive,iel)
          
          ! Also save the neighbour. Note that we have to check the number of the
          ! neighbour against 0 because it may be that the p_IverticesAtElement(:,.)
          ! array is not completely filled -- e.g. if there are triangles in a quad
          ! mesh!
          ivtneighbour = p_IverticesAtElement(mod(ive,nnve)+1,iel)
          if (ivtneighbour .eq. 0) ivtneighbour = p_IverticesAtElement(1,iel)
          p_IverticesAtEdge(2,iedge) = ivtneighbour
        
        end if
      
      end do
      
    end do
    
  end subroutine tria_genVerticesAtEdge2D

  !************************************************************************

!<subroutine>

  subroutine tria_genEdgeNodalProperty2D(rtriangulation)

!<description>
  ! This routine generates the nodal property tags for all edges 
  ! InodalProperty(NVT+1:NVT+NMT) (KNPR). 
  ! For this purpose, the following arrays are used:
  !    InodalProperty(1:NVT), IverticesAtElement, 
  !    IedgesAtElement, IneighboursAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer :: ive
    integer :: iel
    integer :: isize
    integer :: ivt1,ivt2,NVT

    ! Is everything here we need?
    if (rtriangulation%h_InodalProperty .eq. ST_NOHANDLE) then
      call output_line ('InodalPropertys not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      call sys_halt()
    end if
    
    ! Do we have (enough) memory for that array?
    call storage_getsize (rtriangulation%h_InodalProperty, isize)
    if (isize .lt. rtriangulation%NVT+rtriangulation%NMT) then
      ! If the size is wrong, reallocate memory.
      ! Copy the old content as we must not destroy the old nodal property
      ! tags of the vertices.
      ! Fill new entries with 0. This will set the nodal property
      ! array of new vertices stemming from element midpoints to
      ! 'inner nodes'.
      call storage_realloc ('tria_genEdgeNodalProperty2D', &
          rtriangulation%NVT+rtriangulation%NMT, &
          rtriangulation%h_InodalProperty, &
          ST_NEWBLOCK_ZERO, .true.)
    end if
    
    ! Get the arrays.
    call storage_getbase_int (rtriangulation%h_InodalProperty,&
        p_InodalProperty)
    call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
        
    ! Initialise the nodal property with 0 by default.
    call lalg_clearVectorInt (&
        p_InodalProperty(rtriangulation%NVT+1:rtriangulation%NVT+rtriangulation%NMT))
        
    ! Get NVT
    NVT = rtriangulation%NVT
    
    ! Loop through all elements and all edges on the elements
    do iel = 1,ubound(p_IedgesAtElement,2)
    
      do ive = 1,ubound(p_IedgesAtElement,1)
      
        ! Stop if we handled all edges; this is important if there are triangles
        ! in a quad mesh e.g.
        if (p_IedgesAtElement(ive,iel) .eq. 0) exit
        
        ! The edge nodal property is initialised with 0 by default -- inner edge.
        ! Is there a neighbour? If yes, we have an inner edge. If not, this is 
        ! a boundary edge. 
        if (p_IneighboursAtElement(ive,iel) .eq. 0) then
        
          ! Check the two vertices adjacent to that edge. If both are on the
          ! same BC, the nodal property of the vertex is the chosen one.
          ! if they are different, that edge is an edge on the 'blind' boundary
          ! that occurs if a mesh is a subdomain of a larger mesh.
          ! In that case, the nodal property must be NBCT+1 to assign the edge
          ! to the 'blind' part of the boundary.
          ivt1 = p_IverticesAtElement(ive,iel)
          ivt2 = p_IverticesAtElement(mod(ive,ubound(p_IedgesAtElement,1))+1,iel)
          
          ! In case, ivt2=0, there is e.g. a triangle in a quad mesh and we hit a
          ! non-assigned position in IverticesAtElement. So take the next 
          ! assigned position -- which is of course the first position in the array
          ! corresponding to the first vertex of the cell.
          if (ivt2 .eq. 0) ivt2 = p_IverticesAtElement(1,iel)
          
          if (p_InodalProperty(ivt1) .eq. p_InodalProperty(ivt2)) then
            ! Get the number of the boundary component from the vertex preceeding
            ! the edge and store it as information for the edge. 
            p_InodalProperty(p_IedgesAtElement(ive,iel)+NVT) = &
                p_InodalProperty(ivt1)
          else
            ! 'blind' edge
            p_InodalProperty(p_IedgesAtElement(ive,iel)+NVT) = rtriangulation%NBCT+1
          end if
        
        end if
        
      end do
      
    end do

  end subroutine tria_genEdgeNodalProperty2D

  !************************************************************************

!<subroutine>

  subroutine tria_genEdgesAtBoundary2D(rtriangulation)

!<description>
  ! This routine generates information about the edges at the boundary
  ! IedgesAtBoundary (KMBD). Initialises NMBD.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IelementsAtBoundary, IverticesAtBoundary, 
  !    IedgesAtElement, IboundaryCpIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer :: ive, ibct, ivbd, nnve
    integer :: iel
    integer :: isize, ivt

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IboundaryCpIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IedgesAtBoundary .eq. ST_NOHANDLE) then
      ! We have as many elements on the boundary as vertices!
      call storage_new ('tria_genEdgesAtBoundary2D', 'KMBD', &
          rtriangulation%NVBD, ST_INT, &
          rtriangulation%h_IedgesAtBoundary, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IedgesAtBoundary, isize)
      if (isize .ne. rtriangulation%NVBD) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genEdgesAtBoundary2D', &
            rtriangulation%NVBD, rtriangulation%h_IedgesAtBoundary, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)

    nnve = rtriangulation%NNVE

    ! Loop through all boundary components
    do ibct = 1,rtriangulation%NBCT+rtriangulation%NblindBCT
    
      ! On each boundary component, loop through all elements
      do ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
      
        !   +---+---+
        !   |   |   |
        !   +---+---+
        !   |   |IEL|
        !   +---X===+
        !      ivt
      
        ! Get the current boundary vertex and element number
        ivt = p_IverticesAtBoundary(ivbd)
        iel = p_IelementsAtBoundary(ivbd)
        
        ! Find the local number of the vertex in the element
        do ive=1,nnve
          if (p_IverticesAtElement (ive,iel) .eq. ivt) exit
        end do
        
        ! Save the edge following the vertex on that element.
        p_IedgesAtBoundary(ivbd) = p_IedgesAtElement(ive,iel)
      
      end do
    
    end do
    
    ! We have as many edges on the boundary as vertices.
    rtriangulation%NMBD = rtriangulation%NVBD
    
  end subroutine tria_genEdgesAtBoundary2D

  !************************************************************************

!<subroutine>

  subroutine tria_genEdgeParameterValue2D(rtriangulation,rboundary)

!<description>
  ! This routine generates the parameter value of edge midpoints
  ! on the boundary DedgeParameterValue (DMBDP).
  ! For this purpose, the following arrays are used:
  !    DvertexParameterValue, IverticesAtBoundary, IedgesAtBoundary,
  !    IboundaryCpIdx, InodalProperty.
  ! If necessary, new memory is allocated.
!</description>

!<input>
  ! Boundary structure that defines the parametrisation of the boundary.
  type(t_boundary), intent(in) :: rboundary
!</input>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_InodalProperty
    
    integer :: ibct, ivbd, hvertAtBd
    integer :: isize,NVT
    real(DP) :: dpar1,dpar2,dmaxPar

    ! Is everything here we need?
    if (rtriangulation%h_DvertexParameterValue .eq. ST_NOHANDLE) then
      call output_line ('DvertexParameterValue not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IboundaryCpIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)

    call storage_getbase_int (rtriangulation%h_InodalProperty,&
        p_InodalProperty)
    
    ! Allocate an auxiliary array containing a copy of the parameter values 
    ! of the vertices.
    hvertAtBd = ST_NOHANDLE
    call storage_copy (rtriangulation%h_DvertexParameterValue,hvertAtBd)
    call storage_getbase_double (hvertAtBd,p_DvertexParameterValue)
    
    ! Convert the parameter values of the vertices from 0-1 into length
    ! parametrisation. We need this later to get the correct parameter
    ! values of the edge mitpoints.
    !
    ! Note that parameter values -1 are not converted in this routine!
    ! That allows to later detect points that are not on the physical boundary,
    ! as these have parameter value -1.
    do ibct=1,rtriangulation%NBCT
      call boundary_convertParameterList (rboundary,ibct,&
          p_DvertexParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
          p_DvertexParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
          BDR_PAR_01,BDR_PAR_LENGTH)
    end do
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_DedgeParameterValue .eq. ST_NOHANDLE) then
      ! We have as many elements on the boundary as vertices!
      call storage_new ('tria_genEdgeParameterValue2D', 'KMID', &
          rtriangulation%NVBD, ST_DOUBLE, &
          rtriangulation%h_DedgeParameterValue, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_DedgeParameterValue, isize)
      if (isize .ne. rtriangulation%NVBD) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genEdgeParameterValue2D', &
            rtriangulation%NVBD, rtriangulation%h_DedgeParameterValue, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    
    call storage_getbase_double (rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)

    NVT = rtriangulation%NVT

    ! Loop through all boundary components
    do ibct = 1,rtriangulation%NBCT
    
      ! Check if the BC is empty:
      if (p_IboundaryCpIdx(ibct) .lt. p_IboundaryCpIdx(ibct+1)) then
      
        if (ibct .le. rtriangulation%NBCT) then
          ! On the physical boundary, get the maximum parameter value of that BC.
          dmaxPar = boundary_dgetMaxParVal(rboundary, ibct, BDR_PAR_LENGTH)
        end if
      
        ! On each boundary component, loop through all vertices
        do ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-2
        
          ! Check if the edge is really on the boundary. If yes, calculate
          ! its parameter value by taking the mean of the parameter values
          ! of the two endpoints. If the edge belongs to the 'blind'
          ! boundary (happens on subdomains, where the domain boundary os not
          ! the physical boundary), 
          if (p_InodalProperty(p_IedgesAtBoundary(ivbd)+NVT) .le. &
              rtriangulation%NBCT) then
        
            ! Get the parameter value of the current vertex and its neighbour.
            dpar1 = p_DvertexParameterValue(ivbd)
            dpar2 = p_DvertexParameterValue(ivbd+1)
            
            ! The edge parameter value is the mean.
            p_DedgeParameterValue(ivbd) = 0.5_DP*(dpar1+dpar2)
          
          else
            
            ! Otherwise, use -1.0 as parameter value.
            p_DedgeParameterValue(ivbd) = -1.0_DP
          
          end if
        
        end do
      
        ! The 'last' vertex is a special case as there is no real neighbour.
        ! The parameter value of the 'neighbour' is the maximum parameter
        ! value of the boundary component + the parameter value of the
        ! first vertex of the boundary component.
        
        ivbd = p_IboundaryCpIdx(ibct+1)-1
        
        if ((ibct .le. rtriangulation%NBCT) .and. &
            (p_InodalProperty(p_IedgesAtBoundary(ivbd)+NVT) .le. &
             rtriangulation%NBCT)) then
      
          ! Get the parameter value of the current vertex and its neighbour.
          dpar1 = p_DvertexParameterValue(ivbd)
          dpar2 = p_DvertexParameterValue(p_IboundaryCpIdx(ibct)) + dmaxPar
          
          ! The edge parameter value is the mean.
          p_DedgeParameterValue(ivbd) = 0.5_DP*(dpar1+dpar2)
        
        else
          
          ! Otherwise, use -1.0 as parameter value.
          p_DedgeParameterValue(ivbd) = -1.0_DP
        
        end if

        ! On the real boundary...
        if (ibct .le. rtriangulation%NBCT) then
          ! Convert the parameter values of the edge midpoints back from
          ! length parametrisation to 0-1 parametrisation. 
          ! This automatically 'rounds down' parameter values that are > dmaxPar!   
          call boundary_convertParameterList (rboundary,ibct,&
            p_DedgeParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
            p_DedgeParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
            BDR_PAR_LENGTH,BDR_PAR_01)
        end if
        
      end if
        
    end do

    ! Release temporary array, finish.
    call storage_free(hvertAtBd)
    
  end subroutine tria_genEdgeParameterValue2D

  !************************************************************************

!<subroutine>

  subroutine tria_genBoundaryEdgePos2D(rtriangulation)

!<description>
  ! This routine generates the array IboundaryEdgePos which is used
  ! to search for the position of a boundary edges.
  ! For this purpose, the following arrays are used:
  !    IedgesAtBoundary, IboundaryCpIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:,:), pointer :: p_IboundaryEdgePos
    integer :: ivbd, ibct
    integer, dimension(2) :: Isize

    ! Is everything here we need?
    if (rtriangulation%h_IedgesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genBoundaryEdgePos2D')
      call sys_halt()
    end if

    if (rtriangulation%h_IboundaryCpIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genBoundaryEdgePos2D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IboundaryEdgePos .eq. ST_NOHANDLE) then
      ! We have as many elements on the boundary as vertices!
      Isize = (/2,rtriangulation%NMBD/)
      call storage_new ('tria_genBoundaryVertexPos2D', 'KEBDI', &
          Isize, ST_INT, &
          rtriangulation%h_IboundaryEdgePos, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IboundaryEdgePos, Isize)
      if (Isize(2) .ne. rtriangulation%NMBD) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genBoundaryEdgePos2D', &
            rtriangulation%NMBD, rtriangulation%h_IboundaryEdgePos, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    call storage_getbase_int2D (rtriangulation%h_IboundaryEdgePos,&
        p_IboundaryEdgePos)

    ! Fill the array.
    ! Store the IverticesAtBoundary in the first entry and the index in the 2nd.
    do ivbd = 1,rtriangulation%NMBD
      p_IboundaryEdgePos(1,ivbd) = p_IedgesAtBoundary(ivbd)
      p_IboundaryEdgePos(2,ivbd) = ivbd
    end do
    
    ! We only sort the first NBCT BC`s. The boundary component NBCT+1 (if it exists)
    ! collects all vertices that are not on the physical boundary but stem from
    ! the creation of subdomains. These vertices have no parameter value.
    do ibct = 1,rtriangulation%NBCT
      call arraySort_sortByIndex_int (&
          p_IboundaryEdgePos(:,p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),1)
    end do
    
  end subroutine tria_genBoundaryEdgePos2D

  !************************************************************************  

!<subroutine>  

  subroutine tria_genEdgesAtVertex2D (rtriangulation)

!<description>
    ! The routine produces an array pair that
    ! helps you find all edges that are attached to
    ! a particular vertex
!</description>
  
!<inputoutput>
    ! The triangulation structure to be updated.
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    
    integer, dimension(:), pointer :: p_IedgesAtVertexIdx
    integer, dimension(:), pointer :: p_IedgesAtVertex
    
    integer :: iee, iGlobal, isize, index
    integer :: ivt
    ! edgesatelement, dann verticesatedge und fertig
    
    ! Is everything here we need?
    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtVertex2D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtEdge not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtVertex2D')
      call sys_halt()
    end if
    
    ! Get the array out of the triangulation structure
    call storage_getbase_int2d(rtriangulation%h_IedgesAtElement, &
        p_IedgesAtElement)
    
    call storage_getbase_int2d(rtriangulation%h_IverticesAtEdge, &
        p_IverticesAtEdge)
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IedgesAtVertexIdx .eq. ST_NOHANDLE) then
      isize = rtriangulation%NVT+1
      call storage_new ('tria_genEdgesAtVertex2D', 'IedgesAtVertexIdx', &
          isize, ST_INT, &
          rtriangulation%h_IedgesAtVertexIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IedgesAtVertexIdx, isize)
      if (isize .ne. rtriangulation%NVT+1) then
        ! If the size is wrong, reallocate memory.
        isize = rtriangulation%NVT+1
        call storage_realloc ('tria_genEdgesAtVertex2D', &
            isize, rtriangulation%h_IedgesAtVertexIdx, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    call storage_getbase_int(rtriangulation%h_IedgesAtVertexIdx, &
        p_IedgesAtVertexIdx)
    
    ! initialize the p_IedgesAtVertexIdx array
    p_IedgesAtVertexIdx(1) = 1
    p_IedgesAtVertexIdx(2:rtriangulation%NVT+1) = 0
    
    ! loop over all edges
    do iee=1,rtriangulation%NMT
      
      ! get the global vertex index 
      iGlobal=p_IverticesAtEdge(1,iee)
      ! increase the edge count for this vertex
      p_IedgesAtVertexIdx(iGlobal+1) = p_IedgesAtVertexIdx(iGlobal+1) + 1
      
      ! get the global vertex index 
      iGlobal=p_IverticesAtEdge(2,iee)
      ! increase the edge count for this vertex
      p_IedgesAtVertexIdx(iGlobal+1) = p_IedgesAtVertexIdx(iGlobal+1) + 1
      
    end do ! iee
    
    ! sum up the entries to get the index array
    do iee=2,rtriangulation%NVT+1
      p_IedgesAtVertexIdx(iee) = p_IedgesAtVertexIdx(iee) + & 
          p_IedgesAtVertexIdx(iee-1)
    end do
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IedgesAtVertex .eq. ST_NOHANDLE) then
      isize = p_IedgesAtVertexIdx(rtriangulation%NVT+1)-1
      call storage_new ('tria_genEdgesAtVertex2D', 'IedgesAtVertex', &
          isize, ST_INT, &
          rtriangulation%h_IedgesAtVertex, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IedgesAtVertex, isize)
      if (isize .ne. p_IedgesAtVertexIdx(rtriangulation%NVT+1)-1) then
        isize = p_IedgesAtVertexIdx(rtriangulation%NVT+1)-1
        call storage_realloc ('tria_genEdgesAtVertex2D', &
            isize, rtriangulation%h_IedgesAtVertex, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    call storage_getbase_int(rtriangulation%h_IedgesAtVertex, &
        p_IedgesAtVertex)  
    
    ! Shift the array positions.
    !
    ! Do not use an array operation here, may cause a stack overflow
    ! (stupid Intel compiler!)
    !
    ! p_IedgesAtVertexIdx(2:rtriangulation%NVT+1) = &
    !     p_IedgesAtVertexIdx(1:rtriangulation%NVT)    
    do ivt = rtriangulation%NVT,1,-1
      p_IedgesAtVertexIdx(ivt+1) = p_IedgesAtVertexIdx(ivt)
    end do
    
    ! loop over all edges
    do iee=1,rtriangulation%NMT
      
      ! get the global vertex index
      iGlobal=p_IverticesAtEdge(1,iee)
      
      ! get the index into the p_IedgesAtVertex array
      index = p_IedgesAtVertexIdx(iGlobal+1)
      
      ! write the edge number into the array
      p_IedgesAtVertex(index) = iee
      
      ! increase the edge count by one
      p_IedgesAtVertexIdx(iGlobal+1) = p_IedgesAtVertexIdx(iGlobal+1) + 1
      
      ! repeat procedure for the 2nd vertex of the edge
      ! get the global vertex index
      iGlobal=p_IverticesAtEdge(2,iee)
      ! get the index into the p_IedgesAtVertex array
      index = p_IedgesAtVertexIdx(iGlobal+1)
      
      ! write the edge number into the array
      p_IedgesAtVertex(index) = iee
      
      ! increase the edge count by one
      p_IedgesAtVertexIdx(iGlobal+1) = p_IedgesAtVertexIdx(iGlobal+1) + 1
      
    end do ! iee
    
  end subroutine tria_genEdgesAtVertex2D

  !====================================================================
  !
  !       ++++      ++++
  !           +     +   +  
  !           +     +    +
  !        +++      +    +
  !           +     +    + 
  !           +     +   +  
  !       ++++      ++++
  ! tag@3D
  !====================================================================
  
!<subroutine>


  subroutine tria_readTriFile3D(rtriangulation, sfilename, rboundary, &
                                bnoExtendedRaw)

!<description>
  ! This routine reads a .TRI file of a 3D triangulation into memory
  ! and creates a 'raw' triangulation (i.e. a triangulation that contains
  ! only basic information, see below). 
  !
  ! The triangulation structure rtriangulation is initialised with the data 
  ! from the file. The parameter sfilename gives the name of the .tri 
  ! file to read.
  !
  ! This reads only the very basic information that is needed to create
  ! the coarse grid. No information about boundary points, subdivisions
  ! or whatever is created.
  ! The following arrays / information tags will be initialised:
  !
  ! NEL,NVT,NMT,NBCT,NVBD,InelOfType,
  ! DvertexCoords, IverticesAtElement, InodalProperty, 
  ! IboundaryCpIdx, IverticesAtBoundary, DvertexParameterValue.
  ! The nodal property array InodalProperty contains only the data for the
  ! vertices of the triangulation.
  ! The arrays IverticesAtBoundary and DvertexParameterValue are sorted
  ! for the boundary component (according to IboundaryCpIdx) but not 
  ! for the parameter value.
  !
  ! If bnoExtendedRaw is not specified or set to .FALSE., an extended
  ! raw mesh is generated that provides also a proper numbering for edges.
  ! In that case, the following arrays are initialised as well:
  !
  ! IelementsAtVertexIdx,IelementsAtVertexIdx,IneighboursAtElement,
  ! IedgesAtElement,IfacesAtElement
  !
  ! The triangulation structure rtriangulation must be empty. All previous
  ! information in this structure (if there is any) is lost!
!</description>

!<input>
  ! The name of the .tri file to read.
  character(len=*), intent(in) :: sfilename

  ! OPTIONAL: An rboundary object specifying the underlying domain.
  ! If not specified, the routine assumes that the TRI file does not specify
  ! boundary parameter values, i.e. the point coordinates in the TRI file
  ! are all real coordinates. The array DvertexParameterValue is not
  ! generated in this case.
  type(t_boundary), intent(in), optional :: rboundary

  ! OPTIONAL: Prevent creation of an extended raw mesh. If set to .false.,
  ! an 'extended raw' mesh will be created that provides a proper numbering
  ! for edges and faces (standard). If set to '.true', the result will be a 
  ! 'really raw' raw mesh with minimum information and no numbering for edges
  ! and faces.
  logical, intent(in), optional :: bnoExtendedRaw
! </input>
  
!<output>
  ! Triangulation structure, to be filled with data
  type(t_triangulation), intent(out) :: rtriangulation
!</output>
  
!</subroutine>

    ! input channel for reading
    integer :: iunit
    
    ! open the file
    call io_openfileforreading(sfilename, iunit)

    ! we create a 3D triangulation here.
    rtriangulation%ndim = NDIM3D

    ! read the basic mesh
    call tria_readRawTriangulation3D (iunit, rtriangulation)

    ! create the basic boundary information
    call tria_genRawBoundary3D (rtriangulation, rboundary)

    ! Extend the raw mesh by basic edge numbering,
    ! initialise an extended raw mesh.
    if (.not. present(bnoExtendedRaw)) then
      call tria_initExtendedRawMesh (rtriangulation)
    else
      if (.not. bnoExtendedRaw) then
        call tria_initExtendedRawMesh (rtriangulation)
      end if
    end if

    ! close the file, finish
    close(iunit)
  
  end subroutine tria_readTriFile3D

  !************************************************************************  

!<subroutine>  
  
  subroutine tria_readRawTriangulation3D(iunit, rtriangulation)
  
!<description>  
  ! Auxiliary routine of tria_readTriFile3D.
  ! Reads basic information from a triangulation file into rtriangulation.
  ! That means, the following information arrays / tags are initialised:
  ! NEL,NVT,NMT,NBCT,InelOfType,
  ! DvertexCoords, IverticesAtElement and InodalProperty.
  ! The data is read from the file without being changed! haha
!</description>

  
!<input>
  ! Unit number of the file to be read
  integer, intent(in) :: iunit
!</input>
  
!<output>
  ! Triangulation structure, to be filled with data
  type(t_triangulation), intent(inout) :: rtriangulation
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    integer, dimension(:,:), pointer :: p_Idata2D
    integer, dimension(:), pointer :: p_Idata    
    integer, dimension(2) :: Isize
    integer :: idim,ivt,ive,iel
    
    
    ! The first two lines in the file are comments.
    read(iunit,*)
    read(iunit,*)

    ! Read NEL,NVT,NBCT,NVE,NEE,NAE from the file
    ! and store this information in the structure.
    read(iunit,*) rtriangulation%NEL, rtriangulation%NVT,&
                  rtriangulation%NBCT, rtriangulation%NNVE,&
                  rtriangulation%NNEE, rtriangulation%NNAE
    
    ! Vertices per face. That is simple: Only tetrahedral elements
    ! have exactly three vertices per face. All other elements 
    ! have three and four vertices per face.
    if (rtriangulation%NNAE .eq. 4) then
      rtriangulation%NNVA = 3
    else
      rtriangulation%NNVA = 4
    end if
    
    ! skip Comment: 'DCORVG'
    read (iunit,*) 
    
    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM3D, NVT)
    Isize = (/NDIM3D, rtriangulation%NVT/)   
    call storage_new('tria_readRawTriangulation3D', 'DCORVG',&
        Isize, ST_DOUBLE,&
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
    
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_Ddata2D)
    
    ! Read the data from the file, store it in the array.
    ! read data into p_Ddata2D :
    ! first read nvt x-coordinates into p_Ddata2D(1,ivt)
    ! then read nvt  y-coordinates into p_Ddata2D(2,ivt)
    ! then read nvt  z-coordinates into p_Ddata2D(3,ivt)
    read (iunit,*) ((p_Ddata2D(idim,ivt),idim=1,NDIM3D), ivt=1,rtriangulation%NVT)
    
    ! skip Comment: 'KVERT'
    read (iunit,*)
    
    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriangulation%NNVE,rtriangulation%NEL/)
    call storage_new('tria_readRawTriangulation3D', 'KVERT', Isize,&
        ST_INT, rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
    
    ! Get the pointer to the IverticesAtElement array and read the array
    call storage_getbase_int2D(rtriangulation%h_IverticesAtElement, p_Idata2D)
    
    ! read ive=1 indices to nve into p_Idata2D(ive,iel) where iel=1 to NEL     
    read (iunit,*) ((p_Idata2D(ive,iel),ive=1,rtriangulation%NNVE),&
                                        iel=1,rtriangulation%NEL)
    
    ! Loop through the elements and determine how many elements
    ! of each element type we have.
    rtriangulation%InelOfType(:) = 0
    do iel=1,rtriangulation%NEL
      ! start at the last index of element iel down to the first
      do ive=rtriangulation%NNVE,1,-1
        if (p_Idata2D(ive,iel) .ne. 0) then
          rtriangulation%InelOfType(ive) = rtriangulation%InelOfType(ive)+1
          exit
        end if
      end do
    end do

    ! skip Comment: 'KNPR'
    read (iunit,*)    
    
    ! Allocate memory for InodalProperty
    call storage_new('tria_readRawTriangulation3D', 'KNPR',&
        rtriangulation%NVT, ST_INT, &
        rtriangulation%h_InodalProperty,  ST_NEWBLOCK_ZERO)  
    
    ! get a pointer to the memory
    call storage_getbase_int(rtriangulation%h_InodalProperty, p_Idata)
    
    ! Read the data   
    read (iunit,*) (p_idata(ivt),ivt=1,rtriangulation%NVT)   
        
  end subroutine tria_readRawTriangulation3D

  !************************************************************************  
  
!<subroutine>  
  
  subroutine tria_genRawBoundary3D(rtriangulation, rboundary)
  
!<description>  
  ! Auxiliary routine of tria_readTriFile3D.
  ! This routine initialises basic boundary arrays and cleans up 
  ! a basic triangulation. That means:
  ! -> NVBD is calculated
  ! -> IboundaryCpIdx is created and generated
  ! -> IverticesAtBoundary is created and generated.
  !    The vertices are ordered for the boundary component according
  !    to IboundaryCpIdx but not ordered for their parameter value.
  ! -> If rboundary is specified,
  !    DvertexParameterValue is created and generated.
  !    The parameter values are ordered for the boundary component 
  !    according to IboundaryCpIdx but not ordered for the parameter value.
  ! -> If rboundary is specified,
  !    the parameter values of the boundary vertices are extracted
  !    from the first coordinate in DvertexCoords and put into
  !    the DvertexParameterValue array.
  ! -> If rboundary is specified,
  !    based on the parameter value of each boundary vertex, the 
  !    DvertexCoords array receives the actual coordinates of 
  !    the boundary vertices.
  !    If not specified, the routine assumes that DvertexCoords
  !    already contains the real point coordinates.
!</description>

  
!<input>
  ! OPTIONAL: The parametrisation that specifies the coordinates of the 
  ! boundary points.
  ! If specified, DvertexParameterValue is generated from DvertexCoords
  ! and the coordinates of boundary vertices are (re-)generated
  ! by DvertexParameterValue.
  type(t_boundary), intent(in), optional :: rboundary
!</input>
  
!<inputoutput>
  ! Triangulation to be initialised with basic data.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IverticesAtBoundary
    integer :: ivbd,ivt
    integer :: ibct
    integer, dimension(:), pointer :: p_InodalProperty
    
    ! get a pointer to the Inodalproperty array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty, p_InodalProperty)
    
    ! calculate the number of vertices on the boundary (NVBD)
    ! by counting the number of non-zero elements in p_InodalProperty
    rtriangulation%NVBD = 0
    
    ! initialize with zero
    ivbd = 0
    ibct = 0
    
    ! count number of elements on the boundary
    do ivt = 1,rtriangulation%NVT
      if(p_InodalProperty(ivt) .ne. 0) ivbd = ivbd + 1
    end do
    
    ! assign number of vertices on the boundary
    rtriangulation%NVBD = ivbd
    
    ! Allocate memory for IverticesAtBoundary.
    call storage_new('tri_genRawBoundary3D',&
        'KVBD', rtriangulation%NVBD, &
        ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
    
    ! allocate memory for the boundary compnent index vector and
    ! init with zeros
    call storage_new('tri_genRawBoundary3D', &
        'KBCT', rtriangulation%NBCT+1, &
        ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! get pointers to the arrays just created
    call storage_getbase_int(&
        rtriangulation%h_IverticesAtBoundary, p_IverticesAtBoundary)
    
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords, p_DvertexCoords)
    
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
    
    ! the first element in p_IboundaryCpIdx is always 1
    p_IboundaryCpIdx(1) = 1
    
    ! assign the indices of the boundary vertices
    ! first save the number of vertices in each boundary component in
    ! p_IboundaryCpIdx(2:NBCT+1)
    do ivt = 1,rtriangulation%NVT
      if(p_InodalProperty(ivt) .ne. 0) then
        ibct = p_InodalProperty(ivt)
        p_iboundaryCpIdx(ibct+1) = p_iboundaryCpIdx(ibct+1) + 1
      end if
    end do
    
    ! now create the actual index array
    do ibct = 2, rtriangulation%NBCT+1
      p_iboundaryCpIdx(ibct) = p_iboundaryCpIdx(ibct)+p_iboundaryCpIdx(ibct-1)
    end do
    
    
    ! shift indices ah ok !... we increase it again later... haha
    p_IboundaryCpIdx(2:rtriangulation%NBCT+1) = p_IboundaryCpIdx(1:rtriangulation%NBCT)
    
    ! assign the vertices at boundary and the component index array
    do ivt=1, rtriangulation%NVT   
      ! if the vertex is not an inner vertex
      if(p_InodalProperty(ivt) .ne. 0) then
        ! get the id of the boundary component
        ibct = p_InodalProperty(ivt)
        
        ! set ivbd to the number of vertices on that boundary component
        ! thus ivbd holds the current number of vertices found for
        ! boundary component ibct and ivbd represents the current
        ! position in the p_IverticesAtBoundary array
        ivbd = p_IboundaryCpIdx(ibct+1)
        
        ! we have found a new point on that boundary component
        ! so increate the number of points by one
        p_IboundaryCpIdx(ibct+1) = ivbd + 1
        ! store the vertex as boundary vertex
        p_IverticesAtBoundary(ivbd) = ivt
        
      end if
    end do
    
  end subroutine tria_genRawBoundary3D
  
  !************************************************************************

!<subroutine>

    subroutine tria_genElementsAtVertex3D(rtriangulation)    
    
!<description>
  ! This routine generates NnelAtElements and the array IelementsAtVertex.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>
    
!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
    
!</subroutine>
    
    ! local variables
    integer , dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:), pointer :: p_Iaux
    integer :: isize,isize2,haux
    integer :: ive,ivt,iel
    
    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtVertex3D')
      call sys_halt()
    end if

    ! get the array.
    call storage_getbase_int2d(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IelementsAtVertexIdx .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtVertex3D',&
          'IelementsAtVertexIdx', rtriangulation%NVT+1, ST_INT, &
          rtriangulation%h_IelementsAtVertexIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IelementsAtVertexIdx, isize)
      if (isize .ne. rtriangulation%NVT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtVertex3D', &
            rtriangulation%NVT+1,&
            rtriangulation%h_IelementsAtVertexIdx, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Fill the index array with zero.
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)
    call lalg_clearVectorInt (p_IelementsAtVertexIdx)
    
    ! first we calculate the number of elements at each vertex simply
    ! by counting; thus, loop over all elements
    do iel = 1, rtriangulation%NEL

      ! loop over all vertices at the element
      do ive = 1, rtriangulation%NNVE

        ! ivt is the ive-th vertex at element iel
        ivt = p_IverticesAtElement(ive,iel)

        ! check if ivt is 'empty'
        if(ivt .eq. 0) exit
        
        ! increase the number of elements by one
        p_IelementsAtVertexIdx(ivt+1) = p_IelementsAtVertexIdx(ivt+1) + 1
        
      end do ! end ive
    end do ! end iel
    
    ! set the first index to 1
    p_IelementsAtVertexIdx(1) = 1
    rtriangulation%NNelAtVertex = 0
    
    ! In the next step we sum up the number of elements at two
    ! successive vertices to create the index array thus at the
    ! penultimate position of p_IelementsAtVertexIdx we find the
    ! length of p_IelementsAtVertex. Calculate NNelAtVertex.
    do ivt = 2, rtriangulation%NVT+1
      rtriangulation%NNelAtVertex = max(rtriangulation%NNelAtVertex, &
                                        p_IelementsAtVertexIdx(ivt))

      p_IelementsAtVertexIdx(ivt) = p_IelementsAtVertexIdx(ivt) + &
                                    p_IelementsAtVertexIdx(ivt-1)
    end do
    
    ! set the size
    isize = p_IelementsAtVertexIdx(rtriangulation%NVT+1)-1
       
    ! Isize contains now the length of the array where we store the
    ! adjacency information (IelementsAtVertex).  Do we have (enough)
    ! memory for that array?
    if (rtriangulation%h_IelementsAtVertex .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtVertex3D', 'IelementsAtVertex', &
          isize, ST_INT, rtriangulation%h_IelementsAtVertex, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IelementsAtVertex, isize2)
      if (isize .ne. isize2) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtVertex3D', isize,&
            rtriangulation%h_IelementsAtVertex, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get the pointer to the array
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)
      
    ! Duplicate the p_IelementsAtVertexIdx array. We use that as
    ! pointer and index if new elements at a vertex are found.
    haux = ST_NOHANDLE
    call storage_copy (rtriangulation%h_IelementsAtVertexIdx, haux)
    call storage_getbase_int (haux, p_Iaux)
    
    ! loop over all elements
    do iel = 1, rtriangulation%NEL
      ! loop over all vertices of the element
      do ive = 1, rtriangulation%NNVE
        
        ! ivt is the ive-th vertex at element iel
        ivt = p_IverticesAtElement(ive,iel)
        
        ! check if ivt is 'empty'
        if( ivt .eq. 0) exit
        
        ! store the adjacency information at position p_Iaux1(ivt)        
        p_IelementsAtVertex( p_Iaux(ivt) ) = iel

        ! increase the position of the next element in p_Iaux1(ivt)
        p_Iaux(ivt) = p_Iaux(ivt) + 1
        
      end do ! end iel
    end do ! end ive
    
    call storage_free(haux)
    
  end subroutine tria_genElementsAtVertex3D

  !************************************************************************   
  
!<subroutine>  
  
  subroutine tria_genNeighboursAtElement3D(rtriangulation)

!<description>
  ! This routine generates the array IneighboursAtElement.
  ! For this purpose, a list of element connectors is generated.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

!</subroutine>
    
    ! local variables
    integer :: j,iel,iElements  
    integer, dimension(2) :: Isize
    
    ! a pointer to the array this routine is supposed to build
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    
    ! the list of connectors
    type(t_connector3D), dimension(:), pointer :: p_IConnectList
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      Isize = (/rtriangulation%NNAE, rtriangulation%NEL/)
      call storage_new('tria_genNeighboursAtElement3D','KADJ',&
          Isize, ST_INT,&
          rtriangulation%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IneighboursAtElement, Isize)
      if (Isize(2) .ne. rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genNeighboursAtElement3D', &
            rtriangulation%NEL, rtriangulation%h_IneighboursAtElement, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get a pointer to the memory just allocated
    call storage_getbase_int2d(&
        rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    
    ! fill vector with zeros
    call lalg_clearVectorInt2D(p_IneighboursAtElement)   
    
    ! compute number of items for mixed triangulations
    iElements = rtriangulation%InelOfType(TRIA_NVETET3D)  * TRIA_NAETET3D  +&
                rtriangulation%InelOfType(TRIA_NVEPYR3D)  * TRIA_NAEPYR3D  +&
                rtriangulation%InelOfType(TRIA_NVEPRIS3D) * TRIA_NAEPRIS3D +&
                rtriangulation%InelOfType(TRIA_NVEHEXA3D) * TRIA_NAEHEXA3D

    ! first build the connector list
    allocate(p_IConnectList(iElements))
    call tria_buildConnectorList(p_IConnectList, rtriangulation)
        
    ! ConnectorList is build, now sort it
    call tria_sortElements3DInt(p_IConnectList, iElements)
    call tria_sortElements3D(p_IConnectList, iElements)
    
    ! assign the neighbours at elements
    ! traverse the connector list
    do iel = 2, iElements
     
      ! check for equivalent connectors... that means:
      ! check if all 4 vertices that define the face are equal.
      ! For mixed triangulations the fourth vertex may be zero but
      ! this is the case for both items of the connector list
      j = 0
      do while(p_IConnectList(iel-1)%I_conData(j+1) .eq. &
               p_IConnectList(iel)%I_conData(j+1) )
        ! increment counter
        j = j+1 
      end do
      
      ! assign information
      if(j .eq. 4) then
        ! iel is a neighbour of iel-1 at the p_IConnectList(iel-1)%I_conData(6) face
        p_IneighboursAtElement(p_IConnectList(iel-1)%I_conData(6), &
                               p_IConnectList(iel-1)%I_conData(5)) = &
                               p_IConnectList(iel)%I_conData(5)
        
        ! iel-1 is a neighbour of iel at the p_IConnectList(iel)%I_conData(6) face
        p_IneighboursAtElement(p_IConnectList(iel)%I_conData(6), &
                               p_IConnectList(iel)%I_conData(5)) = &
                               p_IConnectList(iel-1)%I_conData(5)
      end if
      
    end do
    
    
    ! free list of connectors
    deallocate(p_IConnectList)

  end subroutine tria_genNeighboursAtElement3D

  !************************************************************************   
  
!<subroutine>  

  subroutine tria_genEdgesAtElement3D(rtriangulation)    

!<description>
  ! This routine creates the information which edges belong to an 
  ! element and stores them in IedgesAtElement. Furthermore the 
  ! edge numbering is calculated. That means, information about 
  ! the edges adjacent to each element IedgesAtElement (KMID) is
  ! generated and the correct value of NMT is calculated.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement,
  !    IelementsAtVertexIdx, IelementsAtVertex.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! local variables
  
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(2) :: Isize
    integer :: iloc1,iloc2,ivt1,ivt2,iedge,iel,ied,iSCElement

    ! list of local edges
    integer, dimension(2,TRIA_NNETET3D), parameter :: IedgesTet =&
             reshape((/1,2, 2,3, 3,1, 1,4, 2,4, 3,4/),&
                     (/2,TRIA_NNETET3D/))
    integer, dimension(2,TRIA_NNEPYR3D), parameter :: IedgesPyr =&
             reshape((/1,2, 2,3, 3,4, 4,1, 1,5, 2,5, 3,5, 4,5/),&
                     (/2,TRIA_NNEPYR3D/))
    integer, dimension(2,TRIA_NNEPRIS3D), parameter :: IedgesPri =&
             reshape((/1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4/),&
                     (/2,TRIA_NNEPRIS3D/))
    integer, dimension(2,TRIA_NNEHEXA3D), parameter :: IedgesHex =&
             reshape((/1,2, 2,3, 3,4, 4,1, 1,5, 2,6,&
                       3,7, 4,8, 5,6, 6,7, 7,8, 8,5/),&
                     (/2,TRIA_NNEHEXA3D/))
    
    ! check if the arrays are present      
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtVertexIdx .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtVertexIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtVertex .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtVertex not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement3D')
      call sys_halt()
    end if
    
    ! get the arrays.
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    
    call storage_getbase_int2D(&
        rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)   
    
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)

    call storage_getbase_int(&
        rtriangulation%h_IelementsAtVertexIdx, p_IelementsAtVertexIdx)

    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      Isize=(/rtriangulation%NNEE, rtriangulation%NEL/)
      call storage_new('tria_genEdgesAtElement3D', 'KMID', Isize,&
          ST_INT, rtriangulation%h_IedgesAtElement, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IedgesAtElement, Isize)
      if (Isize(2) .ne. rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genEdgesAtElement3D', &
            rtriangulation%NEL, rtriangulation%h_IedgesAtElement, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! get the pointer to the memory and fill the array with zeros
    call storage_getbase_int2D(&
        rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call lalg_clearVector(p_IedgesAtElement)
    
    ! iedge counts the edges and specifies the last given edge number.
    iedge = 0
    
    ! loop over all elements
    do iel = 1, rtriangulation%NEL
      
      ! What type of element are we?
      select case(tria_getNVE(p_IverticesAtElement, iel))
      case (TRIA_NVETET3D)
        !=========================================================  
        ! loop over all edges of the tetrahedron
        do ied = 1, TRIA_NNETET3D
          
          ! get the local vertex indices of the current edge
          iloc1 = IedgesTet(1,ied)
          iloc2 = IedgesTet(2,ied)
          
          ! get the corresponding global vertex numbers
          ivt1 = p_IverticesAtElement(iloc1,iel)
          ivt2 = p_IverticesAtElement(iloc2,iel)
          
          ! get the element with the lowest element number that contains 
          ! the edge consisting of the vertices ivt1 and ivt2
          iSCElement = findSmallestCommonElement(p_IelementsAtVertexIdx,&
                                                 p_IelementsAtVertex, ivt1, ivt2, iel)
          
          ! if the smallest common element index greater or equal
          ! to the current iel the edge does not yet have an index
          ! so assign an index to the edge                   
          if(iSCElement .ge. iel) then
            
            ! increment edge number
            iedge = iedge + 1
            
            ! assign the edge number
            p_IedgesAtElement(ied,iel) = iedge
            
          else
            ! the smallest common element index is less than the current
            ! iel the edge already has a number, so search for edge 
            ! (iloc1,iloc2) in the smallest common element
            p_IedgesAtElement(ied,iel) = findEdgeInSmallestCommonElement(&
                                            p_IverticesAtElement, p_IedgesAtElement,&
                                            ivt1, ivt2, iSCElement)
          end if
        end do


      case (TRIA_NVEPYR3D)
        !=========================================================  
        ! loop over all edges of the pyramid
        do ied = 1, TRIA_NNEPYR3D
          
          ! get the local vertex indices of the current edge
          iloc1 = IedgesPyr(1,ied)
          iloc2 = IedgesPyr(2,ied)
          
          ! get the corresponding global vertex numbers
          ivt1 = p_IverticesAtElement(iloc1,iel)
          ivt2 = p_IverticesAtElement(iloc2,iel)
          
          ! get the element with the lowest element number that contains 
          ! the edge consisting of the vertices ivt1 and ivt2
          iSCElement = findSmallestCommonElement(p_IelementsAtVertexIdx,&
                                                 p_IelementsAtVertex, ivt1, ivt2, iel)
          
          ! if the smallest common element index greater or equal to
          ! the current iel the edge does not yet have an index
          ! so assign an index to the edge                   
          if(iSCElement .ge. iel) then
            
            ! increment edge number
            iedge = iedge + 1
            
            ! assign the edge number
            p_IedgesAtElement(ied,iel) = iedge
            
          else
            ! the smallest common element index is less than the current
            ! iel the edge already has a number, so search for edge 
            ! (iloc1,iloc2) in the smallest common element
            p_IedgesAtElement(ied,iel) = findEdgeInSmallestCommonElement(&
                                            p_IverticesAtElement, p_IedgesAtElement,&
                                            ivt1, ivt2, iSCElement)
          end if
        end do


      case (TRIA_NVEPRIS3D)
        !=========================================================  
        ! loop over all edges of the prism
        do ied = 1, TRIA_NNEPRIS3D
          
          ! get the local vertex indices of the current edge
          iloc1 = IedgesPri(1,ied)
          iloc2 = IedgesPri(2,ied)
          
          ! get the corresponding global vertex numbers
          ivt1 = p_IverticesAtElement(iloc1,iel)
          ivt2 = p_IverticesAtElement(iloc2,iel)
          
          ! get the element with the lowest element number that contains 
          ! the edge consisting of the vertices ivt1 and ivt2
          iSCElement = findSmallestCommonElement(p_IelementsAtVertexIdx,&
                                                 p_IelementsAtVertex, ivt1, ivt2, iel)
          
          ! if the smallest common element index greater or equal to
          ! the current iel the edge does not yet have an index
          ! so assign an index to the edge                   
          if(iSCElement .ge. iel) then
            
            ! increment edge number
            iedge = iedge + 1
            
            ! assign the edge number
            p_IedgesAtElement(ied,iel) = iedge
            
          else
            ! the smallest common element index is less than the current
            ! iel the edge already has a number, so search for edge 
            ! (iloc1,iloc2) in the smallest common element
            p_IedgesAtElement(ied,iel) = findEdgeInSmallestCommonElement(&
                                            p_IverticesAtElement, p_IedgesAtElement,&
                                            ivt1, ivt2, iSCElement)
          end if
        end do


      case (TRIA_NVEHEXA3D)
        !=========================================================  
        ! loop over all edges of the hexahedron
        do ied = 1, TRIA_NNEHEXA3D
          
          ! get the local vertex indices of the current edge
          iloc1 = IedgesHex(1,ied)
          iloc2 = IedgesHex(2,ied)
          
          ! get the corresponding global vertex numbers
          ivt1 = p_IverticesAtElement(iloc1,iel)
          ivt2 = p_IverticesAtElement(iloc2,iel)
          
          ! get the element with the lowest element number that contains 
          ! the edge consisting of the vertices ivt1 and ivt2
          iSCElement = findSmallestCommonElement(p_IelementsAtVertexIdx,&
                                                 p_IelementsAtVertex, ivt1, ivt2, iel)
          
          ! if the smallest common element index greater or equal to
          ! the current iel the edge does not yet have an index
          ! so assign an index to the edge                   
          if(iSCElement .ge. iel) then
            
            ! increment edge number
            iedge = iedge + 1
            
            ! assign the edge number
            p_IedgesAtElement(ied,iel) = iedge
            
          else
            ! the smallest common element index is less than the current
            ! iel the edge already has a number, so search for edge 
            ! (iloc1,iloc2) in the smallest common element
            p_IedgesAtElement(ied,iel) = findEdgeInSmallestCommonElement(&
                                            p_IverticesAtElement, p_IedgesAtElement,&
                                            ivt1, ivt2, iSCElement)
          end if
        end do

      case default
        call output_line('Unsupported type of element shape',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement3D')
        call sys_halt()
      end select
    end do ! end iel
    
    ! Save the correct NMT.
    rtriangulation%NMT = iedge

  contains

    ! ---------------------------------------------------------------

    function findSmallestCommonElement(IelementsAtVertexIdx, IelementsAtVertex,&
                                       ivt1, ivt2, iel) result(iSCElement)

      integer, dimension(:), intent(in) :: IelementsAtVertexIdx
      integer, dimension(:), intent(in) :: IelementsAtVertex
      integer, intent(in) :: ivt1,ivt2,iel

      integer :: iSCElement

      ! local variables
      integer :: iel1,iel2,i,j
      
      ! Initialise smallest common element
      iSCElement = iel

      ! loop over all the elements attached to vertex ivt1
      do i = IelementsAtVertexIdx(ivt1), IelementsAtVertexIdx(ivt1+1)-1
        
        ! current element at vertex ivt1
        iel1 = IelementsAtVertex(i)
        
        ! loop over all elements attached to vertex ivt2
        do j = IelementsAtVertexIdx(ivt2), IelementsAtVertexIdx(ivt2+1)-1
          
          ! get the current element at vertex iVGlobal2
          iel2 = IelementsAtVertex(j) 
          
          ! check the vertices share element iel1
          if( (iel2 .eq. iel1) .and.&
              (iel2 .lt. iSCElement) ) then
            iSCElement = iel2      
          end if

        end do ! end j
      end do ! end i     

    end function findSmallestCommonElement

    ! ---------------------------------------------------------------

    function findEdgeInSmallestCommonElement(IverticesAtElement, IedgesAtElement,&
                                             ivt1, ivt2, iel) result(iedge)

      integer, dimension(:,:), intent(in) :: IverticesAtElement
      integer, dimension(:,:), intent(in) :: IedgesAtElement
      integer, intent(in) :: ivt1,ivt2,iel

      integer :: iedge

      ! local variables
      integer :: iVertexAtEdge1, iVertexAtEdge2,i

      ! What type of element are we?
      select case(tria_getNVE(IverticesAtElement, iel))
      case(TRIA_NVETET3D)
        ! Loop over all edges of the tetrahedron
        do i = 1, TRIA_NNETET3D
          
          ! get the indices of the current edge of iel
          iVertexAtEdge1 = IVerticesAtElement(IedgesTet(1,i), iel)
          iVertexAtEdge2 = IVerticesAtElement(IedgesTet(2,i), iel)
          
          if( (ivt1 .eq. iVertexAtEdge1) .and. & 
              (ivt2 .eq. iVertexAtEdge2) .or. &
              (ivt1 .eq. iVertexAtEdge2) .and. &
              (ivt2 .eq. iVertexAtEdge1)) then
            
            ! assign the edge number
            iedge = IedgesAtElement(i,iel)
            exit
          end if
        end do

      case (TRIA_NVEPYR3D)
        ! Loop over all edges of the pyramid
        do i = 1, TRIA_NNEPYR3D
          
          ! get the indices of the current edge of iel
          iVertexAtEdge1 = IVerticesAtElement(IedgesPyr(1,i), iel)
          iVertexAtEdge2 = IVerticesAtElement(IedgesPyr(2,i), iel)
          
          if( (ivt1 .eq. iVertexAtEdge1) .and. & 
              (ivt2 .eq. iVertexAtEdge2) .or. &
              (ivt1 .eq. iVertexAtEdge2) .and. &
              (ivt2 .eq. iVertexAtEdge1)) then
            
            ! assign the edge number
            iedge = IedgesAtElement(i,iel)
            exit
          end if
        end do
        
      case (TRIA_NVEPRIS3D)
        ! Loop over all edges of the prism
        do i = 1, TRIA_NNEPRIS3D
          
          ! get the indices of the current edge of iel
          iVertexAtEdge1 = IVerticesAtElement(IedgesPri(1,i), iel)
          iVertexAtEdge2 = IVerticesAtElement(IedgesPri(2,i), iel)
          
          if( (ivt1 .eq. iVertexAtEdge1) .and. & 
              (ivt2 .eq. iVertexAtEdge2) .or. &
              (ivt1 .eq. iVertexAtEdge2) .and. &
              (ivt2 .eq. iVertexAtEdge1)) then
            
            ! assign the edge number
            iedge = IedgesAtElement(i,iel)
            exit
          end if
        end do

      case (TRIA_NVEHEXA3D)
        ! Loop over all edges of the hexahedron
        do i = 1, TRIA_NNEHEXA3D
          
          ! get the indices of the current edge of iel
          iVertexAtEdge1 = IVerticesAtElement(IedgesHex(1,i), iel)
          iVertexAtEdge2 = IVerticesAtElement(IedgesHex(2,i), iel)
          
          if( (ivt1 .eq. iVertexAtEdge1) .and. & 
              (ivt2 .eq. iVertexAtEdge2) .or. &
              (ivt1 .eq. iVertexAtEdge2) .and. &
              (ivt2 .eq. iVertexAtEdge1)) then
            
            ! assign the edge number
            iedge = IedgesAtElement(i,iel)
            exit
          end if
        end do

      case default
        call output_line('Unsupported type of element shape',&
                         OU_CLASS_ERROR,OU_MODE_STD,'findEdgeInSmallestCommonElement')
        call sys_halt()
      end select

    end function findEdgeInSmallestCommonElement
      
  end subroutine tria_genEdgesAtElement3D

  !************************************************************************   
  
!<subroutine>  

  subroutine tria_genFacesAtElement3D(rtriangulation)

!<description>
  ! This routine builds the FacesAtElement array and
  ! assigns global face numbers.
  ! For this purpose, the following arrays are used:
  !    IneighboursAtElement, IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>
  
!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
  
!</subroutine>  

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IfacesAtElement
    integer, dimension(2) :: Isize    
    integer :: iel,jel,iface,jface,ifaceGlobal

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtElement3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtElement3D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)

    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IfacesAtElement .eq. ST_NOHANDLE) then
      Isize = (/rtriangulation%NNAE, rtriangulation%NEL/)
      call storage_new('tria_genFacesAtElement3D', 'IfacesAtElement',&
          Isize, ST_INT, rtriangulation%h_IfacesAtElement, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IfacesAtElement, Isize)
      if (Isize(2) .ne. rtriangulation%NEL) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genFacesAtElement3D', &
            rtriangulation%NEL, rtriangulation%h_IfacesAtElement, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get the pointer and fill the array with zero
    call storage_getbase_int2D(&
        rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    call lalg_clearVector(p_IfacesAtElement)
    
    ! initialize the global face number
    ifaceGlobal = 0

    ! loop over all elements
    do iel = 1, rtriangulation%NEL
      
      ! loop over all local faces
      do iface = 1, tria_getNAE(p_IverticesAtElement, iel)
        
        ! check if a face number was already assigned
        if((p_IneighboursAtElement(iface,iel) .eq. 0) .or. &
            (p_IneighboursAtElement(iface,iel) >  iel)) then
          
          ! a face number was not yet assigned
          ! thus increment face number
          ifaceGlobal = ifaceGlobal + 1
          
          ! assign the global face number
          p_IfacesAtElement(iface,iel) = ifaceGlobal
          
        else
          
          ! a face number was already assigned
          
          ! get the element number of the neighbour
          jel = p_IneighboursAtElement(iface,iel)
          
          ! initialise the local face number of the neighbouring element
          jface = 1
          
          ! determine the local face number of the neighbouring element
          do while(iel .ne. p_IneighboursAtElement(jface,jel))
            jface = jface+1
          end do
          
          ! assign the global face number
          p_IfacesAtElement(iface,iel) = p_IfacesAtElement(jface,jel)
          
        end if
        
      end do ! end iface
    end do ! end iel
    
    ! number of faces in total
    rtriangulation%NAT = ifaceGlobal

  end subroutine tria_genFacesAtElement3D

  !************************************************************************   

!<subroutine>  

  subroutine tria_genElementsAtEdge3D(rtriangulation)
  
!<description>
  ! This routine generates information about the elements adjacent
  ! to each edge IelementsAtEdge (KMID) and NNelAtEdge. 
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement,
  !    IedgesAtElement.
  ! NMT must be set up correctly.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
    
!</subroutine>  

    ! local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IelementsAtEdgeIdx3D
    integer, dimension(:), pointer :: p_IelementsAtEdge3D
    integer, dimension(:), pointer :: p_Iaux1
    integer :: iel,iedge,iglobalEdge,isize,haux1
    
    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
      call sys_halt()
    end if
    
    if (rtriangulation%NMT .eq. 0) then
      call output_line ('Edge information (NMT) not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
      call sys_halt()
    end if
    
    ! Get the arrays.
    call storage_getbase_int2D(&
        rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)

    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IelementsAtEdgeIdx3D .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtEdge3D',&
          'IelementsAtVertexIdx', rtriangulation%NMT+1, ST_INT, &
          rtriangulation%h_IelementsAtEdgeIdx3D, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IelementsAtEdgeIdx3D, isize)
      if (isize .ne. rtriangulation%NMT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtEdge3D', &
            rtriangulation%NMT+1, rtriangulation%h_IelementsAtEdgeIdx3D, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the array and fill the index array with zero.
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtEdgeIdx3D, p_IelementsAtEdgeIdx3D)
    call lalg_clearVector(p_IelementsAtEdgeIdx3D)

    ! first we calculate the number of elements at each edge simply by counting
    
    ! loop over all elements
    do iel = 1, rtriangulation%NEL
      ! loop over all local edges
      do iedge = 1, rtriangulation%NNEE
        
        ! iglobalEdge is the iedge-th edge at element iel
        iglobalEdge = p_IedgesAtElement(iedge,iel)
        
        ! increase the number of elements by one
        if (iglobalEdge > 0)&
            p_IelementsAtEdgeIdx3D(iglobalEdge+1) =&
            p_IelementsAtEdgeIdx3D(iglobalEdge+1) + 1
        
      end do ! end iedge      
    end do ! end iel
    
    ! set the first index to 1
    p_IelementsAtEdgeIdx3D(1) = 1
    rtriangulation%NNelAtEdge = 0
    
    ! In the next step we sum up the number of elements at two
    ! successive edges to create the index array, thus at
    ! the penultimate position of p_IelementsAtEdgeIdx3D
    ! we find the length of p_IelementsAtEdge3D.
    !
    ! Simultaneously calculate NNelAtEdge.
    do iedge = 2, rtriangulation%NMT+1
      rtriangulation%NNelAtEdge = max(rtriangulation%NNelAtEdge,&
                                      p_IelementsAtEdgeIdx3D(iedge))
      
      p_IelementsAtEdgeIdx3D(iedge) = p_IelementsAtEdgeIdx3D(iedge)+&
                                      p_IelementsAtEdgeIdx3D(iedge-1)
    end do
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IelementsAtEdge3D .eq. ST_NOHANDLE) then
      isize = p_IelementsAtEdgeIdx3D(rtriangulation%NMT+1)-1
      call storage_new ('tria_genElementsAtEdge3D', 'IelementsAtEdge3D', &
          isize, ST_INT,&
          rtriangulation%h_IelementsAtEdge3D, ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rtriangulation%h_IelementsAtEdge3D, isize)
      if (isize .ne. p_IelementsAtEdgeIdx3D(rtriangulation%NMT+1)-1) then
        ! If the size is wrong, reallocate memory.
        isize = p_IelementsAtEdgeIdx3D(rtriangulation%NMT+1)-1
        call storage_realloc ('tria_genElementsAtEdge3D', &
            isize, rtriangulation%h_IelementsAtEdge3D, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get the array pointer  
    call storage_getbase_int(rtriangulation%h_IelementsAtEdge3D,&
        p_IelementsAtEdge3D)
    
    ! Duplicate the p_IelementsAtVertexIdx array. We use that as 
    ! pointer and index when new elements at a vertex are found.
    haux1 = ST_NOHANDLE
    call storage_copy (rtriangulation%h_IelementsAtEdgeIdx3D, haux1)
    call storage_getbase_int (haux1, p_Iaux1)  
    
    do iel = 1, rtriangulation%NEL
      ! loop over all local edges
      do iedge = 1, rtriangulation%NNEE
        
        ! iglobalEdge is the iedge-th edge at element iel
        iglobalEdge = p_IedgesAtElement(iedge,iel)
        
        if (iglobalEdge > 0) then
          ! store the adjacency information at position p_Iaux1(ivt)        
          p_IelementsAtEdge3D( p_Iaux1(iglobalEdge) ) = iel
          ! increase the position of the next element in p_Iaux1(ivt)
          p_Iaux1(iglobalEdge) = p_Iaux1(iglobalEdge) + 1
        end if

      end do ! end iedge
    end do ! end iel
    
    call storage_free(haux1)

  end subroutine tria_genElementsAtEdge3D

  !************************************************************************   

!<subroutine>  

  subroutine tria_genVerticesAtEdge3D(rtriangulation)  

!<description>
  ! This routine assigns the Vertices at Edge array.
  ! For this purpose, the following arrays are used:
  !    IedgesAtElement, IverticesAtElement,
  !    IelementsAtEdgeIdx3D, IelementsAtEdge3D.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
    
!</subroutine>  
  
    ! local variables
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IelementsAtEdgeIdx3D
    integer, dimension(:), pointer :: p_IelementsAtEdge3D
    integer, dimension(2) :: Isize
    integer :: iel,ilocEdge,iedge,iglobalEdge
    
    ! list of local edges
    integer, dimension(2,TRIA_NNETET3D), parameter :: IedgesTet =&
             reshape((/1,2, 2,3, 3,1, 1,4, 2,4, 3,4/),&
                     (/2,TRIA_NNETET3D/))
    integer, dimension(2,TRIA_NNEPYR3D), parameter :: IedgesPyr =&
             reshape((/1,2, 2,3, 3,4, 4,1, 1,5, 2,5, 3,5, 4,5/),&
                     (/2,TRIA_NNEPYR3D/))
    integer, dimension(2,TRIA_NNEPRIS3D), parameter :: IedgesPri =&
             reshape((/1,2, 2,3, 3,1, 1,4, 2,5, 3,6, 4,5, 5,6, 6,4/),&
                     (/2,TRIA_NNEPRIS3D/))
    integer, dimension(2,TRIA_NNEHEXA3D), parameter :: IedgesHex =&
             reshape((/1,2, 2,3, 3,4, 4,1, 1,5, 2,6,&
                       3,7, 4,8, 5,6, 6,7, 7,8, 8,5/),&
                     (/2,TRIA_NNEHEXA3D/))
       
    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtEdgeIdx3D .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtEdgeIdx3D not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtEdge3D .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtEdge3D not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge3D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtEdgeIdx3D, p_IelementsAtEdgeIdx3D)
    call storage_getbase_int(&
        rtriangulation%h_IelementsAtEdge3D, p_IelementsAtEdge3D)
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) then
      Isize = (/2, rtriangulation%NMT/)
      call storage_new('tria_genElementsAtEdge3D', 'IverticesAtEdge',&
          Isize, ST_INT,  rtriangulation%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IverticesAtEdge, Isize)
      if (Isize(2) .ne. rtriangulation%NMT) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtEdge3D', &
            rtriangulation%NMT, rtriangulation%h_IverticesAtEdge, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get the pointer to the memory
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtEdge, p_IverticesAtEdge)
        
    ! loop over all edges (with global edge numbers)
    do iedge = 1, rtriangulation%NMT
      
      ! get an element that is attached to that edge
      iel = p_IelementsAtEdge3D(p_IelementsAtEdgeIdx3D(iedge))
      
      select case(tria_getNVE(p_IverticesAtElement, iel))
      case (TRIA_NVETET3D)
        !========================================================= 
        ! loop over all local edges of the tetrahedron
        do ilocEdge = 1, TRIA_NNETET3D
          ! get the global edge number
          iglobalEdge = p_IedgesAtElement(ilocEdge,iel)
          
          ! check if this edge`s global number equal to the current edge
          if(iedge .eq. iglobalEdge) then
            ! get the global indices of the vertices attached to that
            ! edge and assign the vertex numbers to the edge array
            p_IverticesAtEdge(1,iedge) = p_IverticesAtElement(IedgesTet(1,ilocEdge),iel)
            p_IverticesAtEdge(2,iedge) = p_IverticesAtElement(IedgesTet(2,ilocEdge),iel)
            
            ! That is it
            exit
          end if
        end do


      case (TRIA_NVEPYR3D)
        !========================================================= 
        ! loop over all local edges of the pyramid
        do ilocEdge = 1, TRIA_NNEPYR3D
          ! get the global edge number
          iglobalEdge = p_IedgesAtElement(ilocEdge,iel)

          ! check if this edge`s global number equal to the current edge
          if(iedge .eq. iglobalEdge) then
            ! get the global indices of the vertices attached to that
            ! edge and assign the vertex numbers to the edge array
            p_IverticesAtEdge(1,iedge) = p_IverticesAtElement(IedgesPyr(1,ilocEdge),iel)
            p_IverticesAtEdge(2,iedge) = p_IverticesAtElement(IedgesPyr(2,ilocEdge),iel)
            
            ! That is it
            exit
          end if
        end do


      case (TRIA_NVEPRIS3D)
        !========================================================= 
        ! loop over all local edges of the prism
        do ilocEdge = 1, TRIA_NNEPRIS3D
          ! get the global edge number
          iglobalEdge = p_IedgesAtElement(ilocEdge,iel)

          ! check if this edge`s global number equal to the current edge
          if(iedge .eq. iglobalEdge) then
            ! get the global indices of the vertices attached to that
            ! edge and assign the vertex numbers to the edge array
            p_IverticesAtEdge(1,iedge) = p_IverticesAtElement(IedgesPri(1,ilocEdge),iel)
            p_IverticesAtEdge(2,iedge) = p_IverticesAtElement(IedgesPri(2,ilocEdge),iel)
            
            ! That is it
            exit
          end if
        end do


      case (TRIA_NVEHEXA3D)
        !========================================================= 
        ! loop over all local edges of the hexahedron
        do ilocEdge = 1, TRIA_NNEHEXA3D
          ! get the global edge number
          iglobalEdge = p_IedgesAtElement(ilocEdge,iel)
          
          ! check if this edge`s global number equal to the current edge
          if(iedge .eq. iglobalEdge) then
            ! get the global indices of the vertices attached to that
            ! edge and assign the vertex numbers to the edge array
            p_IverticesAtEdge(1,iedge) = p_IverticesAtElement(IedgesHex(1,ilocEdge),iel)
            p_IverticesAtEdge(2,iedge) = p_IverticesAtElement(IedgesHex(2,ilocEdge),iel)
            
            ! That is it
            exit
          end if
        end do
        

      case DEFAULT
        call output_line('Unsupported type of element shape',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge3D')
        call sys_halt()
      end select
    end do ! end iedge  
  
  end subroutine tria_genVerticesAtEdge3D

  !************************************************************************  

!<subroutine>

  subroutine tria_genEdgeNodalProperty3D(rtriangulation)

!<description>
  ! This routine generates the nodalproperty for the nodes 
  ! that will be formed by the current set of edges.
  ! For this purpose, the following arrays are used:
  !    InodalProperty(1:NMT), IverticesAtEdge,
  !    IboundaryCpEdgeIdx,IedgesAtBoundary.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:), pointer :: p_IboundaryCpEdgesIdx
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer :: isize,ibct,imt,IcpIdx1,IcpIdx2
    
    ! Is everything here we need?
    if (rtriangulation%h_InodalProperty .eq. ST_NOHANDLE) then
      call output_line ('InodalProperty not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtEdge not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IedgesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IboundaryCpEdgesIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpEdgesIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty3D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int(&
        rtriangulation%h_IedgesAtBoundary, p_IedgesAtBoundary)
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpEdgesIdx, p_IboundaryCpEdgesIdx)
    call storage_getbase_int2d(&
        rtriangulation%h_IverticesAtEdge, p_IverticesAtEdge)

    ! Check the size of the p_InodalProperty array its size should be
    ! NVT+NMT+NAT, because in the nodal property array we want to
    ! store the properties of NVT vertices, NMT edges and NAT faces
    call storage_getsize (rtriangulation%h_InodalProperty, isize)
    if (isize .ne. rtriangulation%NVT + &
                   rtriangulation%NMT + &
                   rtriangulation%NAT) then
      ! If the size is wrong, reallocate memory.
      call storage_realloc ('tria_genEdgeNodalProperty3D', &
          rtriangulation%NVT+rtriangulation%NMT+rtriangulation%NAT,&
          rtriangulation%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)
    end if
    
    ! get the array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty, p_InodalProperty)

    ! Loop over all edges to see if the current edge imt is on the
    ! boundary, if it is on the boundary we write the number of the
    ! boundary component to position NVT+imt in p_InodalProperty
    do imt = 1, rtriangulation%NMT
      
      ! get the boundary component index of the vertices that build
      ! the edge and calculate the boundary component number
      ibct = p_InodalProperty(p_IverticesAtEdge(1,imt))
      
      if (ibct .eq. 0) then

        ! mark this edge as inner edge
        p_InodalProperty(rtriangulation%NVT+imt) = 0

      else
        
        IcpIdx1 = p_IboundaryCpEdgesIdx(ibct)
        IcpIdx2 = p_IboundaryCpEdgesIdx(ibct+1)-1
        
        ! check if the edge is located on the boundary
        if (tria_BinSearch(p_IedgesAtBoundary,&
            imt, IcpIdx1, IcpIdx2) .eq. 1) then
          
          ! the edge is located on the boundary; thus, 
          ! write that number to the p_InodalProperty array
          p_InodalProperty(rtriangulation%NVT+imt) = ibct

        else
          
          ! the edge is an inner edge assign a zero
          p_InodalProperty(rtriangulation%NVT+imt) = 0
          
        end if
      end if

    end do ! end imt
           
  end subroutine tria_genEdgeNodalProperty3D

  !************************************************************************  

!<subroutine>

  subroutine tria_genFaceNodalProperty3D(rtriangulation)

!<description>
  ! This routine generates the nodalproperty for the nodes 
  ! that will be formed by the current set of faces.
  ! For this purpose, the following arrays are used:
  
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>
      
    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtFace
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:), pointer :: p_IboundaryCpFacesIdx    
    integer, dimension(:), pointer :: p_IfacesAtBoundary
    integer :: isize,ibct,iface,IcpIdx1,IcpIdx2

    ! Is everything here we need?
    if (rtriangulation%h_InodalProperty .eq. ST_NOHANDLE) then
      call output_line ('InodalProperty not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFaceNodalProperty3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IverticesAtFace .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFaceNodalProperty3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IfacesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IfacesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFaceNodalProperty3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IboundaryCpFacesIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpFacesIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFaceNodalProperty3D')
      call sys_halt()
    end if
    
    ! Get the arrays.
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtBoundary, p_IfacesAtBoundary)
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpFacesIdx, p_IboundaryCpFacesIdx)
    call storage_getbase_int2d(&
        rtriangulation%h_IverticesAtFace, p_IverticesAtFace)
      
    ! Check the size of the p_InodalProperty array its size should be
    ! NVT+NMT+NAT, because in the nodal property array we want to
    ! store the properties of NVT vertices, NMT edges and NAT faces
    call storage_getsize (rtriangulation%h_InodalProperty, isize)
    if (isize .ne. rtriangulation%NVT + &
                   rtriangulation%NMT + &
                   rtriangulation%NAT) then
      ! If the size is wrong, reallocate memory.
      call storage_realloc ('tria_genFaceNodalProperty3D', &
          rtriangulation%NVT+rtriangulation%NMT+rtriangulation%NAT,&
          rtriangulation%h_InodalProperty, ST_NEWBLOCK_NOINIT, .true.)
    end if

    ! get the array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty, p_InodalProperty)
    
    ! Loop over all faces and check if the current face iface is on
    ! the boundary, if it is on the boundary we write the number of
    ! the boundary component to the position NVT+NMT+iface in
    ! p_InodalProperty
    do iface = 1, rtriangulation%NAT
    
      ! get the boundary component index of the vertices that build
      ! the face and calculate the boundary component number
      ibct = p_InodalProperty(p_IverticesAtFace(1,iface))
      
      if(ibct .eq. 0) then
        
        ! mark this face as inner face
        p_InodalProperty(rtriangulation%NVT+rtriangulation%NMT+iface) = 0

      else
        
        IcpIdx1 = p_IboundaryCpFacesIdx(ibct)
        IcpIdx2 = p_IboundaryCpFacesIdx(ibct+1)-1
        
        ! check if the face is located on the boundary
        if(tria_BinSearch(p_IfacesAtBoundary,&
            iface, IcpIdx1, IcpIdx2) .eq. 1) then

          ! the face is located on the boundary; thus,
          ! write that number to the p_InodalProperty array
          p_InodalProperty(rtriangulation%NVT+rtriangulation%NMT+iface) = ibct

        else
          
          ! the face is an inner face so assign zero
          p_InodalProperty(rtriangulation%NVT+rtriangulation%NMT+iface) = 0

        end if
      end if
      
    end do ! end iface
    
  end subroutine tria_genFaceNodalProperty3D

   !************************************************************************  

!<subroutine>  
   
  subroutine tria_genEdgesAtBoundary3D(rtriangulation)

!<description>
  ! This routine generates the edgesAtBoundary array.
  ! For this purpose, the following arrays are used:
  !    IverticesAtEdge, InodalProperty, IfacesAtEdge,
  !    IfacesAtEdgeIdx, IfacesAtBoundary,
  !    
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:), pointer :: p_IfacesAtEdgeIdx
    integer, dimension(:), pointer :: p_IfacesAtEdge
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IfacesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpEdgesIdx
    integer, dimension(:), pointer :: p_IbdyComponents
    integer, dimension(:), pointer :: p_IboundaryCpFacesIdx
    integer, dimension(:), pointer :: p_Iaux
    integer :: haux, h_IbdyComponents,ibct,ibdyFace,imt,inotFound
    integer :: iface,nface,ifaceIndex,iBdyComp,isize
 
    ! Is everything here we need?
    if (rtriangulation%h_InodalProperty .eq. ST_NOHANDLE) then
      call output_line ('InodalPropertys not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtEdge .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtEdge not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IfacesAtBoundary .eq. ST_NOHANDLE) then
      call output_line ('IfacesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IboundaryCpFacesIdx .eq. ST_NOHANDLE) then
      call output_line ('IboundaryCpFacesIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IfacesAtEdgeIdx .eq. ST_NOHANDLE) then
      call output_line ('IfacesAtEdgeIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IfacesAtEdge .eq. ST_NOHANDLE) then
      call output_line ('IfacesAtEdge not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3D')
      call sys_halt()
    end if
    
    ! Get the arrays.
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty, p_InodalProperty)
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtBoundary, p_IfacesAtBoundary)
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtEdgeIdx, p_IfacesAtEdgeIdx)
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtEdge, p_IfacesAtEdge)
    call storage_getbase_int2d(&
        rtriangulation%h_IverticesAtEdge, p_IverticesAtEdge)
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpFacesIdx, p_IboundaryCpFacesIdx)


    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IboundaryCpEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new ('tria_genEdgesAtBoundary3D', 'IboundaryCpEdgesIdx', &
          rtriangulation%NBCT+rtriangulation%NblindBCT+1, ST_INT, &
          rtriangulation%h_IboundaryCpEdgesIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IboundaryCpEdgesIdx, isize)
      if (isize .ne. rtriangulation%NBCT+rtriangulation%NblindBCT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genEdgesAtBoundary3D', &
            rtriangulation%NBCT+rtriangulation%NblindBCT+1,&
            rtriangulation%h_IboundaryCpEdgesIdx, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! Get the array.
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpEdgesIdx, p_IboundaryCpEdgesIdx)
    
    ! create an auxilliary array p_Iaux of size NAT, whereby
    ! p_Iaux(i) = 0 if face i is an inner face and it is
    ! greater 0 if face i is not located at the boundary
    haux = ST_NOHANDLE
    call storage_new ('tria_genEdgesAtBoundary3D', 'Iaux', &
        rtriangulation%NAT, ST_INT, haux, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(haux, p_Iaux)

    ! another auxilliary array which is used as a pointer to
    ! the free positions for the different boundary components
    h_IbdyComponents = ST_NOHANDLE
    call storage_copy(&
        rtriangulation%h_IboundaryCpFacesIdx, h_IbdyComponents)
    call storage_getbase_int (h_IbdyComponents, p_IbdyComponents)
    
    
    ! the first index is one; the remaining indices are initialised with zeros
    call lalg_clearVector(p_IboundaryCpEdgesIdx)
    p_IboundaryCpEdgesIdx(1) = 1
    
    ! nface points to the position in the facesAtBoundary array
    nface = 1
    
    ! Loop over all faces and count the number of edges on the boundary
    do iface = 1, rtriangulation%NAT
      
      ! check if the maximum number of faces on the boundary is reached
      if(nface > rtriangulation%NABD) exit
      
      inotFound = 0
      ! check for this face in every boundary component
      do ibct = 1, rtriangulation%NBCT+rtriangulation%NblindBCT
        
        ! check if face iface a boundary face
        ibdyFace = p_IbdyComponents(ibct)

        ! if we are at or beyond the end of that 
        ! boundary component then skip it
        if(ibdyFace .ge. p_IboundaryCpFacesIdx(ibct+1)) cycle
        
        if(p_IfacesAtBoundary(ibdyFace) .eq. iface) then

          ! we have identified a boundary face
          p_Iaux(iface) = 1

          ! increase the number of boundary faces found
          nface = nface + 1

          ! increase the pointer to the boundary comp index
          p_IbdyComponents(ibct) = p_IbdyComponents(ibct) + 1
          exit ! break out of the ibct loop 

        else

          ! increase the not found counter
          inotFound = inotFound + 1

        end if
      end do ! end do ibct
      
      ! if not found
      if(inotFound .eq. rtriangulation%NBCT+rtriangulation%NblindBCT) then
        p_Iaux(iface) = 0
      end if
      
    end do ! end iface
    
    ! initialise the number of edges at the boundary
    rtriangulation%NMBD = 0
    
    ! Loop over all edges
    do imt = 1, rtriangulation%NMT
      
      ! check all faces connected to the current edge
      do iface = p_IfacesAtEdgeIdx(imt), p_IfacesAtEdgeIdx(imt+1)-1
        
        ifaceIndex = p_IfacesAtEdge(iface)
        if(p_Iaux(ifaceIndex) .eq. 1) then

          ! the edge is connected to a boundary face so it is a boundary
          ! edge; thus, increase the number of edges on the boundary
          rtriangulation%NMBD = rtriangulation%NMBD + 1
          
          ! get the boundary component number of a vertex of the edge
          iBdyComp = p_InodalProperty(p_IverticesAtEdge(1,imt))
          
          ! increase the number of edges by one
          p_IboundaryCpEdgesIdx(iBdyComp+1) = p_IboundaryCpEdgesIdx(iBdyComp+1) + 1

          ! break out of this loop
          exit
          
        end if
        
      end do ! end iface
    end do ! end ive
    
    ! Sum up the entries to compute the actual index array
    do ibct = 2, rtriangulation%NBCT+rtriangulation%NblindBCT+1
      p_IboundaryCpEdgesIdx(ibct) = p_IboundaryCpEdgesIdx(ibct) + &
                                    p_IboundaryCpEdgesIdx(ibct-1)
    end do ! end ive
    
    ! shift the indices...
    p_IboundaryCpEdgesIdx(2:rtriangulation%NBCT+rtriangulation%NblindBCT+1) = &
        p_IboundaryCpEdgesIdx(1:rtriangulation%NBCT+rtriangulation%NblindBCT)


    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IedgesAtBoundary .eq. ST_NOHANDLE) then
      call storage_new ('tria_genEdgesAtBoundary3D', 'IedgesAtBoundary', &
          rtriangulation%NMBD, ST_INT, &
          rtriangulation%h_IedgesAtBoundary, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IedgesAtBoundary, isize)
      if (isize .ne. rtriangulation%NMBD) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genEdgesAtBoundary3D', rtriangulation%NMBD,&
            rtriangulation%h_IedgesAtBoundary, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the array.
    call storage_getbase_int(&
        rtriangulation%h_IedgesAtBoundary, p_IedgesAtBoundary)
    
    ! Loop over all edges
    do imt = 1, rtriangulation%NMT
      
      ! check all faces connected to the current edge
      do iface = p_IfacesAtEdgeIdx(imt), p_IfacesAtEdgeIdx(imt+1)-1
      
        ifaceIndex = p_IfacesAtEdge(iface)
        if(p_Iaux(ifaceIndex) .eq. 1) then
          
          ! get the corresponding boundary component for the boundary edge
          iBdyComp = p_InodalProperty(p_IverticesAtEdge(1,imt))
          
          ! store the edge number 
          p_IedgesAtBoundary(p_IboundaryCpEdgesIdx(iBdyComp+1)) = imt
          
          ! increase the pointer in the edgesAtBoundary array
          p_IboundaryCpEdgesIdx(iBdyComp+1) = p_IboundaryCpEdgesIdx(iBdyComp+1) + 1

          ! break out of this loop
          exit  

        end if

      end do ! end iface    
    end do ! end iface
    
    ! free memory
    call storage_free(haux)
    call storage_free(h_IbdyComponents)
    
  end subroutine tria_genEdgesAtBoundary3D

  !************************************************************************  

!<subroutine> 

  subroutine tria_genFacesAtBoundary3D(rtriangulation)

!<description>
  ! This routine generates the facesAtBoundary array.
  ! For this purpose, the following arrays are used:
  !    IverticesAtFace, IelementsAtFace, InodalProperty.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IelementsAtFace
    integer, dimension(:,:), pointer :: p_IverticesAtFace
    integer, dimension(:), pointer :: p_InodalProperty
    integer, dimension(:), pointer :: p_IfacesAtBoundary
    integer, dimension(:), pointer :: p_IboundaryCpFacesIdx    
    integer :: iface, ibct, ifbd, isize
    
    ! Is everything here we need?
    if (rtriangulation%h_InodalProperty .eq. ST_NOHANDLE) then
      call output_line ('InodalPropertys not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtBoundary3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtFace .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtBoundary3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IverticesAtFace .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtBoundary3D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty, p_InodalProperty)
    call storage_getbase_int2D(&
        rtriangulation%h_IelementsAtFace, p_IelementsAtFace)
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtFace, p_IverticesAtFace)
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IboundaryCpFacesIdx .eq. ST_NOHANDLE) then
      call storage_new ('tria_genFacesAtBoundary3D', 'IboundaryCpFacesIdx', &
          rtriangulation%NBCT+rtriangulation%NblindBCT+1, ST_INT, &
          rtriangulation%h_IboundaryCpFacesIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IboundaryCpFacesIdx, isize)
      if (isize .ne. rtriangulation%NBCT+rtriangulation%NblindBCT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genFacesAtBoundary3D', &
            rtriangulation%NBCT+rtriangulation%NblindBCT+1,&
            rtriangulation%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase_int(&
        rtriangulation%h_IboundaryCpFacesIdx, p_IboundaryCpFacesIdx)
    
    ! the first index is one; the remaining indices are
    ! initialized with zeros and updated step-by-step
    call lalg_clearVectorInt(p_IboundaryCpFacesIdx)
    p_IboundaryCpFacesIdx(1) = 1
    
    ! loop over all faces to count the
    ! number of faces on the boundary 
    do iface = 1, rtriangulation%NAT
    
      ! Check if we are a boundary face, that is, the second
      ! element adjacent to the current face is empty
      if (p_IelementsAtFace(2,iface) .eq. 0) then

        ! Get the number of the boundary component
        ibct = p_InodalProperty(p_IverticesAtFace(1,iface))
        
        ! increase the number of entries for the corresponding 
        ! boundary component. Note that the array was initialised
        ! with zero during the creation process!
        p_IboundaryCpFacesIdx(ibct+1) = p_IboundaryCpFacesIdx(ibct+1) + 1
      end if
    end do

    ! Sum up the entries to compute the actual index array
    do ibct = 2, rtriangulation%NBCT+rtriangulation%NblindBCT+1
      p_IboundaryCpFacesIdx(ibct) = p_IboundaryCpFacesIdx(ibct) + &
                                    p_IboundaryCpFacesIdx(ibct-1)
    end do

    ! assign the number of faces on the boundary
    rtriangulation%NABD = p_IboundaryCpFacesIdx(rtriangulation%NBCT + &
                                                rtriangulation%NblindBCT+1) - 1
    
    ! shift the indices
    p_IboundaryCpFacesIdx(2:rtriangulation%NBCT+rtriangulation%NblindBCT+1) = &
        p_IboundaryCpFacesIdx(1:rtriangulation%NBCT+rtriangulation%NblindBCT)

    ! allocate memory for the p_IfacesAtBoundary array
    if(rtriangulation%h_IfacesAtBoundary .eq. ST_NOHANDLE) then
      call storage_new ('tria_genFacesAtBoundary3D', 'IfacesAtBoundary', &
          rtriangulation%NABD, ST_INT,&
          rtriangulation%h_IfacesAtBoundary, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IfacesAtBoundary, isize)
      if (isize .ne. rtriangulation%NABD) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genFacesAtBoundary3D', &
            rtriangulation%NABD, rtriangulation%h_IfacesAtBoundary, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtBoundary, p_IfacesAtBoundary)
    
    ! loop over all faces and fill the array p_IfacesAtBoundary
    do iface = 1, rtriangulation%NAT
      
      ! Check if we are a boundary face, that is, the second
      ! element adjacent to the current face is empty
      if (p_IelementsAtFace(2,iface) .eq. 0) then
        
        ! Get the number of the boundary component
        ibct = p_InodalProperty(p_IverticesAtFace(1,iface))
                
        ! Get the next free number of the face on the boundary
        ifbd = p_IboundaryCpFacesIdx(ibct+1)
        
        ! we have found a new face on that boundary component
        ! so increase the number of faces by one
        p_IboundaryCpFacesIdx(ibct+1) = ifbd+1
        
        ! store the face as boundary face
        p_IfacesAtBoundary(ifbd) = iface
      end if
    end do

  end subroutine tria_genFacesAtBoundary3D

  !************************************************************************  

!<subroutine>    

  subroutine tria_genFacesAtVertex3D(rtriangulation)

!<description>
  ! This routine builds the FacesAtVertex list.
  ! For this purpose, the following arrays are used:
  !    IverticesAtFace
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtFace
    integer, dimension(:), pointer :: p_IfacesAtVertexIdx
    integer, dimension(:), pointer :: p_IfacesAtVertex
    integer, dimension(:), pointer :: p_Iaux1
    integer :: iface, ivt, ive, haux1, isize

    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtFace .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtVertex3D')
      call sys_halt()
    end if
    
    ! Get the array.
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtFace, p_IverticesAtFace)
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IfacesAtVertexIdx .eq. ST_NOHANDLE) then
      call storage_new ('tria_genFacesAtVertex3D',&
          'IfacesAtVertexIdx', rtriangulation%NVT+1, ST_INT, &
          rtriangulation%h_IfacesAtVertexIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IfacesAtVertexIdx, isize)
      if (isize .ne. rtriangulation%NVT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genFacesAtVertex3D', &
            rtriangulation%NVT+1, rtriangulation%h_IfacesAtVertexIdx, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if

    ! get the pointer and fill the index array with zero.
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtVertexIdx, p_IfacesAtVertexIdx)
    call lalg_clearVector(p_IfacesAtVertexIdx)
    
    ! set the first value to 1
    p_IfacesAtVertexIdx(1) = 1
    
    ! create the index array
    do iface = 1, rtriangulation%NAT
      do ive = 1, 4
        
        ! get the global vertex number
        ivt = p_IverticesAtFace(ive,iface)
        
        ! increment the facecount for this vertex
        if (ivt > 0) p_IfacesAtVertexIdx(ivt+1) = &
                     p_IfacesAtVertexIdx(ivt+1) + 1

      end do ! end ive
    end do ! end iface
    
    
    ! create the actual index array  
    do iface = 2, rtriangulation%NVT+1
      p_IfacesAtVertexIdx(iface) = p_IfacesAtVertexIdx(iface) + &
                                   p_IfacesAtVertexIdx(iface-1);
    end do ! end iface
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IfacesAtVertex .eq. ST_NOHANDLE) then
      call storage_new ('tria_genFacesAtVertex3D', 'IfacesAtVertex', &
          p_IfacesAtVertexIdx(rtriangulation%NVT+1)-1, ST_INT, &
          rtriangulation%h_IfacesAtVertex, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IfacesAtVertex, isize)
      if (isize .ne. p_IfacesAtVertexIdx(rtriangulation%NVT+1)-1) then
        ! If the size is wrong, reallocate memory.
        isize = p_IfacesAtVertexIdx(rtriangulation%NVT+1)-1
        call storage_realloc ('tria_genFacesAtVertex3D', &
            isize, rtriangulation%h_IfacesAtVertex, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get the pointer  
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtVertex, p_IfacesAtVertex)
    
    ! build the auxilliary array  
    haux1 = ST_NOHANDLE
    call storage_copy (rtriangulation%h_IfacesAtVertexIdx, haux1)
    call storage_getbase_int (haux1, p_Iaux1) 
    
    ! assign the connectivity info
    do iface = 1, rtriangulation%NAT
      do ive = 1, 4

        ! get the global vertex number
        ivt = p_IverticesAtFace(ive,iface)
        
        if (ivt > 0) then
          ! store the adjacency information at position p_Iaux1(ivt)        
          p_IfacesAtVertex( p_Iaux1(ivt) ) = iface
          
          ! increase the position of the next element in p_Iaux1(ivt)
          p_Iaux1(ivt) = p_Iaux1(ivt) + 1    
        end if
        
      end do ! end ive
    end do ! end iface
    
    call storage_free(haux1)   
    
  end subroutine tria_genFacesAtVertex3D

  !************************************************************************  

!<subroutine>

  subroutine tria_genFacesAtEdge3D (rtriangulation) 
  
!<description>
  ! This routine builds the FacesAtEdge array.
  ! For this purpose, the following arrays are used:
  !    IedgesAtFace
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>

!</subroutine>
    
    ! local parameters
    integer, dimension(:,:), pointer :: p_IedgesAtFace
    integer, dimension(:), pointer :: p_IfacesAtEdgeIdx
    integer, dimension(:), pointer :: p_IfacesAtEdge
    integer, dimension(:), pointer :: p_Iaux1
    integer :: haux1, iface, iglobalEdge, iedge, isize
    
    ! Is everything here we need?
    if (rtriangulation%h_IedgesAtFace .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtEdge3D')
      call sys_halt()
    end if
    
    ! Get the array.
    call storage_getbase_int2d(&
        rtriangulation%h_IedgesAtFace, p_IedgesAtFace)
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IfacesAtEdgeIdx .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtEdge3D', 'IfacesAtEdgeIdx', &
          rtriangulation%NMT+1, ST_INT, &
          rtriangulation%h_IfacesAtEdgeIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IfacesAtEdgeIdx, isize)
      if (isize .ne. rtriangulation%NMT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genFacesAtEdge3D', rtriangulation%NMT+1,&
            rtriangulation%h_IfacesAtEdgeIdx, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the array
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtEdgeIdx, p_IfacesAtEdgeIdx)
    
    ! Initialialise the index arrax
    call lalg_clearVector(p_IfacesAtEdgeIdx)
    p_IfacesAtEdgeIdx(1) = 1
    
    ! Loop over all faces and count the number of faces for each edge
    do iface = 1, rtriangulation%NAT
      
      ! Loop over all edges at the current face
      do iedge = 1, 4

        ! Get the global edge number
        iglobalEdge = p_IedgesAtFace(iedge,iface)

        ! increase the number of faces at these edges by one
        if (iglobalEdge > 0)&
            p_IfacesAtEdgeIdx(iglobalEdge+1) = &
            p_IfacesAtEdgeIdx(iglobalEdge+1) + 1 
        
      end do ! end iedge
    end do ! end iface
    
    ! Sum up the number of faces at each edge to get the actual index vector
    do iedge = 2, rtriangulation%NMT+1
      p_IfacesAtEdgeIdx(iedge) = p_IfacesAtEdgeIdx(iedge) +&
                                 p_IfacesAtEdgeIdx(iedge-1);
    end do ! end iedge
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IfacesAtEdge .eq. ST_NOHANDLE) then
      isize = p_IfacesAtEdgeIdx(rtriangulation%NMT+1)-1
      call storage_new ('tria_genFacesAtEdge3D', 'IfacesAtEdge', &
          p_IfacesAtEdgeIdx(rtriangulation%NMT+1)-1, ST_INT, &
          rtriangulation%h_IfacesAtEdge, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IfacesAtEdge, isize)
      if (isize .ne. p_IfacesAtEdgeIdx(rtriangulation%NMT+1)-1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genFacesAtEdge3D',&
            p_IfacesAtEdgeIdx(rtriangulation%NMT+1)-1,&
            rtriangulation%h_IfacesAtEdge, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the array
    call storage_getbase_int(&
        rtriangulation%h_IfacesAtEdge, p_IfacesAtEdge)
    
    haux1 = ST_NOHANDLE
    call storage_copy (rtriangulation%h_IfacesAtEdgeIdx, haux1)
    call storage_getbase_int (haux1, p_Iaux1) 
    
    ! Loop over all faces and assign the connectivity info
    do iface = 1, rtriangulation%NAT
      
      ! Loop over all edges at the current face
      do iedge = 1, 4

        ! Get the global edge number
        iglobalEdge = p_IedgesAtFace(iedge,iface)
        
        if (iglobalEdge > 0) then
          ! store the adjacency information at position p_Iaux1(ivt)        
          p_IfacesAtEdge( p_Iaux1(iglobalEdge) ) = iface
          ! increase the position of the next element in p_Iaux1(ivt)
          p_Iaux1(iglobalEdge) = p_Iaux1(iglobalEdge) + 1    
        end if
        
      end do ! end iedge
    end do ! end iface
    
    call storage_free(haux1)   
    
  end subroutine tria_genFacesAtEdge3D

  !************************************************************************  

!<subroutine>
  
  subroutine tria_genEdgesAtFace3D(rtriangulation)  

!<description>
  ! This routine builds the EdgesAtFace array.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IfacesAtElement,
  !    IedgesAtElement, IelementsAtFace.
  ! If necessary, new memory is allocated.
!</description>
  
!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IfacesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IelementsAtFace
    integer, dimension(:,:), pointer :: p_IedgesAtFace
    integer :: iface, iel, ilocalFace, iglobalFace
    integer, dimension(2) :: Isize
    
    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IfacesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IfacesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IedgesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtFace .eq. ST_NOHANDLE) then
      call output_line ('IelementsAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IelementsAtFace, p_IelementsAtFace)
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IedgesAtFace .eq. ST_NOHANDLE) then
      Isize = (/4,rtriangulation%NAT/)
      call storage_new('tria_genEdgesAtFace3D', 'IedgesAtFace',&
          Isize, ST_INT,&
          rtriangulation%h_IedgesAtFace, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IedgesAtFace, Isize)
      if (Isize(2) .ne. rtriangulation%NAT) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genEdgesAtFace3D', rtriangulation%NEL,&
            rtriangulation%h_IedgesAtFace, ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    
    ! Get the array.
    call storage_getbase_int2D(&
        rtriangulation%h_IedgesAtFace, p_IedgesAtFace)
    
    ! loop over all faces
    do iface = 1, rtriangulation%NAT
      
      ! get the first face this element is conncected to
      iel = p_IelementsAtFace(1,iface)
      
      ! determine which local face is iface 
      do ilocalFace = 1, rtriangulation%NNAE
        iglobalFace = p_IfacesAtElement(ilocalFace,iel)
        if(iglobalFace .eq. iface) exit
      end do

      ! What type of element are we?
      select case(tria_getNVE(p_IverticesAtElement, iel))
      case (TRIA_NVETET3D)
        !=========================================================  
        select case (ilocalFace)
        case (1)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(4,iface)=0

        case (2)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(4,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(4,iface)=0

        case (3)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(4,iface)=0

        case (4)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(4,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(4,iface)=0

        case default
          call output_line('Invalid local face number',&
                           OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
          call sys_halt()
        end select


      case (TRIA_NVEPYR3D)
        !=========================================================  
        select case (ilocalFace)
        case (1)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(4,iel)

        case (2)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(4,iface)=0

        case (3)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(7,iel)
          p_IedgesAtFace(4,iface)=0

        case (4)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(7,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(8,iel)
          p_IedgesAtFace(4,iface)=0

        case (5)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(4,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(8,iel)
          p_IedgesAtFace(4,iface)=0

        case default
          call output_line('Invalid local face number',&
                           OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
          call sys_halt()
        end select

        
      case (TRIA_NVEPRIS3D)
        !=========================================================  
        select case (ilocalFace)
        case (1)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(4,iface)=0
        
        case (2)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(4,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(7,iel)

        case (3)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(8,iel)

        case (4)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(4,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(9,iel)

        case (5)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(7,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(8,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(9,iel)
          p_IedgesAtFace(4,iface)=0

        case default
          call output_line('Invalid local face number',&
                           OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
          call sys_halt()
        end select
          

      case (TRIA_NVEHEXA3D)
        !=========================================================  
        select case (ilocalFace)
        case (1)
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(4,iel)

        case (2)    
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(1,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(9,iel)
          
        case (3)    
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(2,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(6,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(7,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(10,iel)
          
        case (4)    
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(3,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(7,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(8,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(11,iel)
          
        case (5)    
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(4,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(8,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(5,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(12,iel)
          
        case (6)    
          ! assign the edges
          p_IedgesAtFace(1,iface)=p_IedgesAtElement(9,iel)
          p_IedgesAtFace(2,iface)=p_IedgesAtElement(10,iel)
          p_IedgesAtFace(3,iface)=p_IedgesAtElement(11,iel)
          p_IedgesAtFace(4,iface)=p_IedgesAtElement(12,iel)
          
        case default
          call output_line('Invalid local face number',&
                           OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported type of element shape',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtFace3D')
        call sys_halt()
      end select
      
    end do ! end iface
    
  end subroutine tria_genEdgesAtFace3D

  !************************************************************************  
  
!<subroutine>  

  subroutine tria_genElementsAtFace3D(rtriangulation) 

!<description>
  ! This routine builds the ElementsAtFace array.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement,
  !    IfacesAtElement.
  ! If necessary, new memory is allocated.
!</description>
  
!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>

!</subroutine>
  
    ! local variables
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IfacesAtElement
    integer, dimension(:,:), pointer :: p_IelementsAtFace
    integer :: iface,iel,jel,ifaceNumber
    integer, dimension(2) :: Isize
    

    ! Is everything here we need?
    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtFace3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtFace3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IfacesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IfacesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtFace3D')
      call sys_halt()
    end if

    ! Get the arrays.
    call storage_getbase_int2D(&
        rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IfacesAtElement, p_IfacesAtElement)

    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IelementsAtFace .eq. ST_NOHANDLE) then
      Isize = (/2,rtriangulation%NAT/)
      call storage_new('tria_genElementsAtFace3D', 'IelementsAtFace', &
          Isize, ST_INT, &
          rtriangulation%h_IelementsAtFace, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IelementsAtFace, Isize)
      if (Isize(2) .ne. rtriangulation%NAT) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtFace3D', &
            rtriangulation%NAT, rtriangulation%h_IelementsAtFace,&
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get the pointer
    call storage_getbase_int2D(&
        rtriangulation%h_IelementsAtFace, p_IelementsAtFace)
    
    
    ! loop over all elements
    do iel = 1, rtriangulation%NEL
      ! loop over all faces of this element
      do iface = 1, tria_getNAE(p_IverticesAtElement, iel)
        
        ! get the number of the neighbouring element
        jel = p_IneighboursAtElement(iface,iel)

        ! if there is no neighbour at this element
        if(jel .eq. 0) then
          
          ifaceNumber = p_IfacesAtElement(iface,iel)
          
          p_IelementsAtFace(1,ifaceNumber) = iel
          p_IelementsAtFace(2,ifaceNumber) = 0
          
        elseif (jel < iel) then
          
          ! There is a neighbour and it has a smaller number than the
          ! current element. That is, we have not had that face!
          ! Store the two adjacent elements.
          
          ifaceNumber = p_IfacesAtElement(iface,iel)
          
          p_IelementsAtFace(1,ifaceNumber) = iel
          p_IelementsAtFace(2,ifaceNumber) = jel
         
        end if
        
      end do ! end iface
    end do ! end iel
    
  end subroutine tria_genElementsAtFace3D

  !************************************************************************  

!<subroutine>  

  subroutine tria_genVerticesAtFace3D(rtriangulation)

!<description>
  ! This routine builds the VerticesAtFace array.
  ! For this purpose, the following arrays are used:
  !    IverticesAtElement, IneighboursAtElement,
  !    IfacesAtElement.
  ! If necessary, new memory is allocated.
!</description>
  
!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
  
!</subroutine>    

    ! local variables
    integer, dimension(:,:), pointer  :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IfacesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtFace
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    integer, dimension(2) :: Isize
    integer :: iel,jel,iface,jface,ive,ivt,ifaceGlobal,ifaceNumber

    ! list of local face numbers
    integer, dimension(4,TRIA_NAETET3D), parameter :: IfacesTet =&
             reshape((/1,2,3,0, 1,4,2,0, 2,4,3,0,&
                       3,4,1,0/), (/4,TRIA_NAETET3D/))
    integer, dimension(4,TRIA_NAEPYR3D), parameter :: IfacesPyr =&
             reshape((/1,2,3,4, 1,5,2,0, 2,5,3,0,&
                       3,5,4,0, 4,5,1,0/), (/4,TRIA_NAEPYR3D/))
    integer, dimension(4,TRIA_NAEPRIS3D), parameter :: IfacesPri =&
             reshape((/1,2,3,0, 1,4,5,2, 2,5,6,3,&
                       3,6,4,1, 4,6,5,0/), (/4,TRIA_NAEPRIS3D/))
    integer, dimension(4,TRIA_NAEHEXA3D), parameter :: IfacesHex =&
             reshape((/1,2,3,4, 1,5,6,2, 2,6,7,3,&
                       3,7,8,4, 1,4,8,5, 5,8,7,6/), (/4,TRIA_NAEHEXA3D/))
    
    ! Is everything here we need?
    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtFace3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IneighboursAtElement .eq. ST_NOHANDLE) then
      call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtFace3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IfacesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IfacesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtFace3D')
      call sys_halt()
    end if
    
    ! Get the arrays.
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    call storage_getbase_int2D(&
        rtriangulation%h_IfacesAtElement, p_IfacesAtElement)
    
    ! Do we have (enough) memory for that array?
    if(rtriangulation%h_IverticesAtFace .eq. ST_NOHANDLE) then
      Isize = (/4, rtriangulation%NAT/)
      call storage_new('tria_genVerticesAtFace3D', 'IverticesAtFace', &
          Isize, ST_INT, &
          rtriangulation%h_IverticesAtFace, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rtriangulation%h_IverticesAtFace, Isize)
      if (Isize(2) .ne. rtriangulation%NAT) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genVerticesAtFace3D', &
            rtriangulation%NAT, rtriangulation%h_IverticesAtFace, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! get the pointer
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtFace, p_IverticesAtFace)
    
    ! initialize the global face number
    ifaceGlobal = 0
    
    ! loop over all elements
    do iel = 1, rtriangulation%NEL

      ! What kind of element shape are we?
      select case(tria_getNVE(p_IverticesAtElement, iel))
      case (TRIA_NVETET3D)
        !=========================================================  
        ! loop over all local faces of the tetrahedron
        do iface = 1, TRIA_NAETET3D
        
          ! check if a face number was already assigned
          if((p_IneighboursAtElement(iface,iel) .eq. 0) .or. &
              (p_IneighboursAtElement(iface,iel) >  iel)) then
          
            ! a face number was not yet assigned
            ! thus increment face number
            ifaceGlobal = ifaceGlobal + 1
            
            ! assign the vertices at this face
            do ive = 1, 3
              p_IverticesAtFace(ive,ifaceGlobal) = &
                  p_IverticesAtElement(IfacesTet(ive,iface),iel)
            end do
            p_IverticesAtFace(4,  ifaceGlobal) = 0
            
          else
            
            ! a face number was already assigned
            
            ! get the element number of the neighbour
            jel = p_IneighboursAtElement(iface,iel)
            
            ! initialise the local face number of the neighbouring element
            jface = 1
            
            ! determine the local face number of the neighbouring element
            do while(iel .ne. p_IneighboursAtElement(jface,jel))
              jface = jface+1
            end do
            
            ! index of the face in the array is obtained by nmt and nvt
            ifaceNumber = p_IfacesAtElement(jface,jel)
            
            ! assign the vertices at this face
            do ive = 1, 3
            p_IverticesAtFace(ive,ifaceNumber) = &
                p_IverticesAtElement(IfacesTet(ive,iface),iel)
          end do
            p_IverticesAtFace(4,  ifaceNumber) = 0

          end if
        end do ! end iface


      case (TRIA_NVEPYR3D)
        !=========================================================  
        ! loop over all local faces of the pyramid
        do iface = 1, TRIA_NAEPYR3D
        
          ! check if a face number was already assigned
          if((p_IneighboursAtElement(iface,iel) .eq. 0) .or. &
             (p_IneighboursAtElement(iface,iel) >  iel)) then
          
            ! a face number was not yet assigned
            ! thus increment face number
            ifaceGlobal = ifaceGlobal + 1
            
            ! assign the vertices at this face
            do ive = 1, 4
              ivt = IfacesPyr(ive,iface)
              if (ivt .eq. 0) then
                p_IverticesAtFace(ive,ifaceGlobal) = 0
              else
                p_IverticesAtFace(ive,ifaceGlobal) = &
                    p_IverticesAtElement(ivt,iel)
              end if
            end do
            
          else
            
            ! a face number was already assigned
            
            ! get the element number of the neighbour
            jel = p_IneighboursAtElement(iface,iel)
            
            ! initialise the local face number of the neighbouring element
            jface = 1
            
            ! determine the local face number of the neighbouring element
            do while(iel .ne. p_IneighboursAtElement(jface,jel))
              jface = jface+1
            end do
            
            ! index of the face in the array is obtained by nmt and nvt
            ifaceNumber = p_IfacesAtElement(jface,jel)
            
            ! assign the vertices at this face
            do ive = 1, 4
              ivt = IfacesPyr(ive,iface)
              if (ivt .eq. 0) then
                p_IverticesAtFace(ive,ifaceNumber) = 0
              else
                p_IverticesAtFace(ive,ifaceNumber) = &
                    p_IverticesAtElement(ivt,iel)
              end if
            end do
            
          end if
        end do ! end iface


      case (TRIA_NVEPRIS3D)
        !=========================================================  
        ! loop over all local faces of the prism
        do iface = 1, TRIA_NAEPRIS3D
        
          ! check if a face number was already assigned
          if((p_IneighboursAtElement(iface,iel) .eq. 0) .or. &
             (p_IneighboursAtElement(iface,iel) >  iel)) then
          
            ! a face number was not yet assigned
            ! thus increment face number
            ifaceGlobal = ifaceGlobal + 1
            
            ! assign the vertices at this face
            do ive = 1, 4
              ivt = IfacesPri(ive,iface)
              if (ivt .eq. 0) then
                p_IverticesAtFace(ive,ifaceGlobal) = 0
              else
                p_IverticesAtFace(ive,ifaceGlobal) = &
                    p_IverticesAtElement(ivt,iel)
              end if
            end do
            
          else
            
            ! a face number was already assigned
            
            ! get the element number of the neighbour
            jel = p_IneighboursAtElement(iface,iel)
            
            ! initialise the local face number of the neighbouring element
            jface = 1
            
            ! determine the local face number of the neighbouring element
            do while(iel .ne. p_IneighboursAtElement(jface,jel))
              jface = jface+1
            end do
            
            ! index of the face in the array is obtained by nmt and nvt
            ifaceNumber = p_IfacesAtElement(jface,jel)
            
            ! assign the vertices at this face
            do ive = 1, 4
              ivt = IfacesPri(ive,iface)
              if (ivt .eq. 0) then
                p_IverticesAtFace(ive,ifaceNumber) = 0
              else
                p_IverticesAtFace(ive,ifaceNumber) = &
                    p_IverticesAtElement(ivt,iel)
              end if
            end do
            
          end if
        end do ! end iface


      case (TRIA_NVEHEXA3D)
        !=========================================================  
        ! loop over all local faces of the hexahedron
        do iface = 1, TRIA_NAEHEXA3D
        
          ! check if a face number was already assigned
          if((p_IneighboursAtElement(iface,iel) .eq. 0) .or. &
             (p_IneighboursAtElement(iface,iel) >  iel)) then
          
            ! a face number was not yet assigned
            ! thus increment face number
            ifaceGlobal = ifaceGlobal + 1
            
            ! assign the vertices at this face
            do ive = 1, 4
              p_IverticesAtFace(ive,ifaceGlobal) = &
                  p_IverticesAtElement(IfacesHex(ive,iface),iel)
            end do
            
          else
            
            ! a face number was already assigned
            
            ! get the element number of the neighbour
            jel = p_IneighboursAtElement(iface,iel)
            
            ! initialise the local face number of the neighbouring element
            jface = 1
            
            ! determine the local face number of the neighbouring element
            do while(iel .ne. p_IneighboursAtElement(jface,jel))
              jface = jface+1
            end do
            
            ! index of the face in the array is obtained by nmt and nvt
            ifaceNumber = p_IfacesAtElement(jface,jel)
            
            ! assign the vertices at this face
            do ive = 1, 4
              p_IverticesAtFace(ive,ifaceNumber) = &
                  p_IverticesAtElement(IfacesHex(ive,iface),iel)
            end do
            
          end if
        end do ! end iface


      case default
        call output_line('Unsupported type of element shape',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtFace3D')
        call sys_halt()
      end select
      
    end do ! end iel

  end subroutine tria_genVerticesAtFace3D

  ! ***************************************************************************

!<subroutine>

  subroutine tria_genElementVolume3D(rtriangulation)

!<description>
  ! This routine generates the element volume array DelementVolume (DAREA). 
  ! For this purpose, the following arrays are used:
  !    DvertexCoordinates, IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DelementVolume
    integer :: iel
    integer :: isize
    real(DP) :: dtotalVolume
    real(DP), dimension(NDIM3D,TRIA_MAXNVE3D) :: Dpoints
    integer :: ive

    ! Is everything here we need?
    if (rtriangulation%h_DvertexCoords .eq. ST_NOHANDLE) then
      call output_line ('DvertexCoords not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume3D')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtElement .eq. ST_NOHANDLE) then
      call output_line ('IverticesAtElement  not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume3D')
      call sys_halt()
    end if
    
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_DelementVolume .eq. ST_NOHANDLE) then
      call storage_new ('tria_genElementVolume3D', 'DAREA', &
          rtriangulation%NEL+1, ST_DOUBLE, &
          rtriangulation%h_DelementVolume, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_DelementVolume, isize)
      if (isize .ne. rtriangulation%NEL+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementVolume3D', &
            rtriangulation%NEL+1, rtriangulation%h_DelementVolume, &
            ST_NEWBLOCK_NOINIT, .false.)
      end if
    end if
    
    ! Get the arrays
    call storage_getbase_double2D (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double (rtriangulation%h_DelementVolume,&
        p_DelementVolume)
        
    dtotalVolume = 0.0_DP
        
    ! Calculate the element volume for all elements
    do iel=1,rtriangulation%NEL
      ! triangular element
      do ive=1,TRIA_NVEHEXA3D
        Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
        Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        Dpoints(3,ive) = p_DvertexCoords(3,p_IverticesAtElement(ive,iel))
      end do
      p_DelementVolume(iel) = gaux_getVolume_hexa3D(Dpoints)
      
      dtotalVolume = dtotalVolume+p_DelementVolume(iel)
    end do
    
    ! Store the total volume in the last element of DelementVolume
    p_DelementVolume(rtriangulation%NEL+1) = dtotalVolume
    
  end subroutine tria_genElementVolume3D
  
  !************************************************************************

!<subroutine>  

  pure subroutine tria_getVerticesAtFaceDirect(iface, nva, IverticesAtElement,&
                                               bpositive, IverticesAtFace)

!<description>
  ! This routine calculates the vertices on face iface of element iel.
  ! The vertex numbers are returned in mathematically positive or negative 
  ! sense when 'looking' at the vertices from the center of the element.
!</description>

!<input>
  ! Local face number (1..NAE=6 usually).
  integer, intent(in) :: iface
  
  ! Maximum number of vertices per face. 3 for 3D tetraheral meshes,
  ! 4 for 3D hexahedral meshes. 
  integer, intent(in) :: nva

  ! Array containing the vertices on the current element
  integer, dimension(:), intent(in) :: IverticesAtElement
  
  ! TRUE=return the vertices in mathematically positive sense.
  ! FALSE=return the vertices in mathematically negative sense.
  logical, intent(in) :: bpositive
!</input>
  
!<output>
  ! Array receiving the vertices on local face iface.
  integer, dimension(:), intent(out) :: IverticesAtFace
!</output>
  
!</subroutine>

    ! local variables
    integer :: i

    ! List of vertices on each of the faces -- for hexahedral meshes. 
    ! This array defines the ordering of the vertices. 
    ! We have the following local vertex vertices:
    !
    !           8----------------7
    !          /|               /|
    !         / |              / |
    !        /  |             /  |
    !       /   |            /   |
    !      5----------------6    |
    !      |    |           |    |
    !      |    4-----------|----3
    !      |   /            |   / 
    !      |  /             |  /  
    !      | /              | /   
    !      |/               |/ 
    !      1----------------2  
    ! 
    ! The faces are ordered as follows; the smallest local vertex number is
    ! always vertex 1 of the face. The ordering of the vertices is always
    ! in mathematically positive sense when looking 'from the center'
    ! of the element to the face:
    !                 
    !                                                                    -7 
    !                                                                    /| 
    !                                                                   / | 
    !                                                                  /  | 
    !                            /                /                   /   | 
    !                           5----------------6                  -6    | 
    !       |                |  |                |                   | 3  | 
    !       4----------------3  |                |                   |   -3 
    !      /                /   |                |                   |   /  
    !     /       1        /    |       2        |                   |  /   
    !    /                /     |                |                   | /    
    !  |/               |/      |/               |/                  |/     
    !  1----------------2       1----------------2                  -2      
    !
    !   8----------------7         8-                    8----------------7
    !  /|               /|        /|                    /|               /|
    !   |                |       / |                   /      6         /  
    !   |       4        |      /  |                  /                /   
    !   |                |     /   |                 /                /    
    !   |                |    5-   |                5----------------6     
    !   |                |    |  5 |                |                |     
    !   4----------------3    |    4-                                         
    !  /                /     |   /                                           
    !                         |  /                                            
    !                         | /                                             
    !                         |/                                              
    !                         1-
    !
    ! Note that the real element may be inverted (face 6 below face 1). In that
    ! case, all meanings of positive/negative orientation are changed to the
    ! opposite. Fortunately, this does usually not harm any algorithm working
    ! with the cell connectivity.
          
    integer, dimension(3,TRIA_NAETET3D), parameter :: IverticesTet =&
             reshape((/1,2,3, 1,4,2, 2,4,3, 1,4,3/), (/3,TRIA_NAETET3D/))

    integer, dimension(4,TRIA_NAEPYR3D), parameter :: IverticesPyr =&
             reshape((/1,2,3,4, 1,5,2,0, 2,5,3,0,&
                       3,4,5,0, 1,5,4,0/), (/4,TRIA_NAEPYR3D/))

    integer, dimension(4,TRIA_NAEPRIS3D), parameter :: IverticesPri =&
             reshape((/1,2,3,0, 1,4,5,2, 2,3,6,5,&
                       1,4,6,3, 4,5,6,0/), (/4,TRIA_NAEPRIS3D/))
    
    integer, dimension(4,TRIA_NAEHEXA3D), parameter :: IverticesHexa =&
             reshape((/1,2,3,4, 1,5,6,2, 2,6,7,3,&
                       3,7,8,4, 1,4,8,5, 5,8,7,6/), (/4,TRIA_NAEHEXA3D/))
    
    ! What type of element are we
    select case(ubound(IverticesAtElement,1))
    case(TRIA_NVETET3D)
      !=========================================================  
      ! Tetrahedron element: Get the vertices on the face
      if (bpositive) then
        do i=1,nva
          IverticesAtFace(i) = IverticesAtElement(IverticesTet(i,iface))
        end do
      else
        do i=1,nva
          IverticesAtFace(nva-i+1) = IverticesAtElement(IverticesTet(i,iface))
        end do
      end if

    case(TRIA_NVEHEXA3D)
      !=========================================================  
      ! Hexahedron element: Get the vertices on the face
      if (bpositive) then
        do i=1,nva
          IverticesAtFace(i) = IverticesAtElement(IverticesHexa(i,iface))
        end do
      else
        do i=1,nva
          IverticesAtFace(nva-i+1) = IverticesAtElement(IverticesHexa(i,iface))
        end do
      end if

    case default
      IverticesAtFace = 0
    end select

  end subroutine tria_getVerticesAtFaceDirect

  !====================================================================
  ! AUXILIARY ROUTINES FOR CONNECTOR LISTS
  ! tag@aux
  !====================================================================
  
!<subroutine>      

  subroutine tria_buildConnectorList(IConnectList, rtriangulation)

!<description>
  ! This routine builds the connector list used 
  ! in the neighbours at elements routine
!</description>

!<input>
  type(t_triangulation), intent(in) :: rtriangulation
!</input>
    
!<output>    
  ! the list of connectors this routine is supposed to build
  type(t_connector3D), dimension(:), intent(out) :: IConnectList    
!</output>  
  
!</subroutine>    

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: iel,k,nfaces

    ! function body
    
    ! Get some data arrays about the vertices.
    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
    
    ! initialise the number of faces
    nfaces = 0
    
    ! loop through all elements
    do iel = 1, rtriangulation%NEL
      select case(tria_getNVE(p_IverticesAtElement, iel))
      case (TRIA_NVETET3D)
        ! build connectors for each tetrahedron

        !=========================================================  
        ! first face
        nfaces = nfaces+1
        
        do k=1,3
          IConnectList(nfaces)%I_conData(k) = p_IverticesAtElement(k,iel)
        end do
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 1
        
        !=========================================================
        ! second face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(4) = 0

        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 2
        
        !=========================================================  
        ! third face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 3
        
        !=========================================================  
        ! fourth face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 4
        
        !========================================================= 


      case (TRIA_NVEPYR3D)
        ! build connectors for each pyramid

        !=========================================================  
        ! first face
        nfaces = nfaces+1

        do k=1,4
          IConnectList(nfaces)%I_conData(k) = p_IverticesAtElement(k,iel)
        end do
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 1

        !=========================================================
        ! second face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(5,iel)
        IConnectList(nfaces)%I_conData(4) = 0

        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 2
        
        !=========================================================  
        ! third face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(5,iel)
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 3
        
        !=========================================================  
        ! fourth face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(5,iel)
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 4

        !=========================================================  
        ! fifth face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(5,iel)
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 5

        !=========================================================  


      case (TRIA_NVEPRIS3D)
        ! build connectors for each prism

        !=========================================================  
        ! first face
        nfaces = nfaces+1
        
        do k=1,3
          IConnectList(nfaces)%I_conData(k) = p_IverticesAtElement(k,iel)
        end do
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 1

        !=========================================================  
        ! fifth face
        nfaces = nfaces+1
        
        do k=4,6
          IConnectList(nfaces)%I_conData(k-3) = p_IverticesAtElement(k,iel)
        end do
        IConnectList(nfaces)%I_conData(4) = 0
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 5

        !=========================================================
        ! second face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(5,iel)
        IConnectList(nfaces)%I_conData(4) = p_IverticesAtElement(4,iel)
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 2
        
        !=========================================================
        ! third face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(6,iel)
        IConnectList(nfaces)%I_conData(4) = p_IverticesAtElement(5,iel)
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 3
        
        !=========================================================
        ! fourth face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(4) = p_IverticesAtElement(6,iel)
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 4

        !=========================================================  


      case (TRIA_NVEHEXA3D)
        ! build connectors for each hexahedron
        
        !=========================================================  
        ! first face
        nfaces = nfaces+1
        
        do k=1,4
          IConnectList(nfaces)%I_conData(k) = p_IverticesAtElement(k,iel)
        end do
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 1

        !=========================================================
        ! sixth face
        nfaces = nfaces+1
        
        do k=5,8
          IConnectList(nfaces)%I_conData(k-4) = p_IverticesAtElement(k,iel)
        end do
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 6

        !=========================================================
        ! second face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(5,iel)
        IConnectList(nfaces)%I_conData(4) = p_IverticesAtElement(6,iel)
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 2
        
        !=========================================================  
        ! fourth face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(7,iel)
        IConnectList(nfaces)%I_conData(4) = p_IverticesAtElement(8,iel)
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 4
        
        !=========================================================  
        ! third face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(2,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(3,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(6,iel)
        IConnectList(nfaces)%I_conData(4) = p_IverticesAtElement(7,iel)
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 3
        
        !=========================================================  
        ! fifth face
        nfaces = nfaces+1
        
        IConnectList(nfaces)%I_conData(1) = p_IverticesAtElement(1,iel)
        IConnectList(nfaces)%I_conData(2) = p_IverticesAtElement(4,iel)
        IConnectList(nfaces)%I_conData(3) = p_IverticesAtElement(8,iel)
        IConnectList(nfaces)%I_conData(4) = p_IverticesAtElement(5,iel)
        
        ! save the number of the element this face was found from
        IConnectList(nfaces)%I_conData(5) = iel
        
        ! assign the local face number
        IConnectList(nfaces)%I_conData(6) = 5
        
        !=========================================================  

      case default
        call output_line('Unsupported type of element shape',&
                         OU_CLASS_ERROR,OU_MODE_STD,'tria_buildConnectorList')
        call sys_halt()
      end select
      
    end do
    
  end subroutine tria_buildConnectorList

  !************************************************************************   

!<subroutine>

  subroutine tria_sortElements3D(IConnectList, iElements)
  
!<description>
  ! This subroutine establishes the lexicographic
  ! ordering on the list of connectors in 3D
!</description>

!<input>    
  integer, intent(in) :: iElements
!</input>  

!<inputoutput>        
  type(t_connector3D), dimension(:), intent(inout) :: IConnectList
!</inputoutput>

!</subroutine>

    ! local
    integer :: j
    
    do j = TRIA_NCONNECT3D, 1, -1
      call tria_mergesort(IConnectList, 1, iElements, j)
    end do
    
  end subroutine tria_sortElements3D

  !************************************************************************   
  
!<subroutine>

  subroutine tria_sortElements3DInt(IConnectList, iElements)

!<description>
  ! This subroutine establishes the sorted numbering 
  ! on the list of connectors in 3D
!</description>
    
  ! parameter values

!<input>    
  integer, intent(in) :: iElements
!</input>

!<inputoutput>        
  type(t_connector3D), dimension(:), intent(inout) :: IConnectList
!</inputoutput>          
        
!</subroutine>
    
  ! local variables
    integer :: i

    ! create a sorted numbering in all connectors    
    do i = 1, iElements
      call sort(IConnectList(i)%I_conData(1:4))
    end do
    
  contains

    ! ---------------------------------------------------------------

    pure subroutine sort(Idata)
      integer, intent(inout), dimension(4) :: Idata

      if (Idata(2) < Idata(1)) call swap(Idata(2), Idata(1))
      if (Idata(3) < Idata(2)) call swap(Idata(3), Idata(2))
      if (Idata(4) < Idata(3)) call swap(Idata(4), Idata(3))
      if (Idata(2) < Idata(1)) call swap(Idata(2), Idata(1))
      if (Idata(3) < Idata(2)) call swap(Idata(3), Idata(2))
      if (Idata(2) < Idata(1)) call swap(Idata(2), Idata(1))
    end subroutine sort

    ! ---------------------------------------------------------------

    elemental pure subroutine swap(a,b)
      integer, intent(inout) :: a,b

      ! local variables
      integer :: c

      c = a
      a = b
      b = c
    end subroutine swap
    
  end subroutine tria_sortElements3DInt

  !************************************************************************   

!<subroutine>      

  recursive subroutine tria_mergesort(IConnectList, l, r, pos)
    
!<description>
  ! This routine sorts a connector list it is used as an 
  ! auxilliary routine during the Neighbours at elements routine
!</description>
    
!<input>
  ! the array positions l...r will be sorted
  ! the sorting key is element 'pos' of the connector
  integer, intent(in) :: l,r,pos
!</input>
  
!<inputoutput>
  ! the list of connectors    
  type(t_connector3D), dimension(:), intent(inout) :: IConnectList
!</inputoutput>

!</subroutine>      
    
    ! local variables
    integer :: m
    
    if(l < r) then
      
      m = l + (r-l)/2
      
      call tria_mergesort(IConnectList, l,   m, pos)
      call tria_mergesort(IConnectList, m+1, r, pos)
      call tria_merge(IConnectList, l, m, r, pos)
    
    end if
    
  end subroutine tria_mergesort

  !************************************************************************      

!<subroutine> 

  subroutine tria_merge(IConnectList, l, m, r, pos)
    
!<description>
  ! 
  ! standard auxilliary routine in the mergesort algorithm
  ! 
!</description>
    
!<input> 
  ! the array positions l...r will be sorted
  ! the sorting key is element 'pos' of the connector
  integer, intent(in) :: l,r,m,pos
!</input>

!<inputoutput>
  ! the list of connectors     
  type(t_connector3D), dimension(:), intent(inout) :: IConnectList
!</inputoutput>

!</subroutine>     

    ! local variables
    integer :: i,j,n1,n2,k
    type(t_connector3D), dimension(:), pointer :: p_L, p_R
    
    ! function body
    
    
    ! init counters
    n1 = m - l + 1
    
    n2 = r - m 
    
    k = l    
    
    ! allocate memory for merging
    allocate(p_L(n1))
    allocate(p_R(n2))
    
    ! fill left array
    do i = 1, n1
      p_L(i) = IConnectList(l+i-1)
    end do
    
    ! fill right array
    do j = 1, n2
      p_R(j) = IConnectList(m+j)
    end do
    
    i = 1
    j = 1
    
    ! merge 
    do
      if( (i > n1 ) .or. (j > n2) ) exit
      
      ! if the current element of the left array is smaller
      ! copy it to p_ConnectorList
      ! else
      ! copy the element from the right array
      if(p_L(i)%I_conData(pos) .le. p_R(j)%I_conData(pos)) then
        IConnectList(k) = p_L(i)
        i = i + 1
        k = k + 1
      else
        IConnectList(k) = p_R(j)                 
        j = j + 1
        k = k + 1
      end if
      
    end do
    
    ! copy the remaining entries of p_L (if present)
    do
      if(i > n1) exit
      
      IConnectList(k) = p_L(i)
      ! increment counters
      k = k + 1
      i = i + 1
      
    end do
    
    ! copy the remaining entries of p_R (if present)
    do
      if(j > n2) exit
      
      IConnectList(k) = p_R(j)
      ! increment counters
      k = k + 1
      j = j + 1
      
    end do
    
    ! done merging
    
    ! free p_L and p_R
    deallocate(p_L)
    deallocate(p_R)
    
  end subroutine tria_merge

  !************************************************************************     

!<function>

  function tria_BinSearch(p_Iarray, Ivalue, Ilbound, Iubound)

!<description>
  ! This function performs binary searching for integer arrays
!</description>
  
!<input>
    ! The integer array to search in
    integer, dimension(:),intent(in) :: p_Iarray

    ! The integer to search for
    integer, intent(in) :: Ivalue

    ! The lower and uppers bounds used in the array
    integer, intent(in) :: Ilbound, Iubound
!</input>

!<result>
    integer :: tria_BinSearch
!</result>
!</function>
    
    ! local variables
    integer :: Imid, Ilboundloc, Iuboundloc
    
    Ilboundloc = Ilbound
    Iuboundloc = Iubound
    
    ! standard binary search scheme...
    do while( Ilboundloc .le. Iuboundloc )
    
      Imid = (Ilboundloc + Iuboundloc) / 2
    
      if(p_Iarray(Imid) > Ivalue) then
        Iuboundloc = Imid-1
      elseif(p_Iarray(Imid) < Ivalue) then
        Ilboundloc = Imid+1
      else ! found Ivalue
        tria_BinSearch = 1
        return
      end if
    
    end do ! end while
    
    ! Ivalue was not found...
    tria_BinSearch = 0
  
  end function tria_BinSearch
  
  !************************************************************************      

!<subroutine> 

  subroutine tria_calcBoundingBox(rtriangulation,DboundingBoxMin,DboundingBoxMax)
    
!<description>
  ! Calculates the X/Y/Z coordinates of a bounding box surrounding the mesh.
!</description>
    
!<input> 
  ! Underlying triangulation
  type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<output>
  ! Minimum X/Y/Z coordinate of a bounding box around the mesh.
  real(DP), dimension(:), intent(out) :: DboundingBoxMin

  ! Maximum X/Y/Z coordinate of a bounding box around the mesh.
  real(DP), dimension(:), intent(out) :: DboundingBoxMax
!</output>

!</subroutine>     

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer :: ipt,idim
    
    ! Get the coordinate array
    call storage_getbase_double2d(rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    DboundingBoxMin(1:rtriangulation%ndim) = p_DvertexCoords(1:rtriangulation%ndim,1)
    DboundingBoxMax(1:rtriangulation%ndim) = p_DvertexCoords(1:rtriangulation%ndim,1)
    
    ! Find the minimum x/y/z-coordinate in each direction
    do ipt = 2,rtriangulation%nvt
      do idim = 1,rtriangulation%ndim
        DboundingBoxMin(idim) = min(DboundingBoxMin(idim),p_DvertexCoords(idim,ipt))
        DboundingBoxMax(idim) = max(DboundingBoxMax(idim),p_DvertexCoords(idim,ipt))
      end do
    end do

  end subroutine
  
end module triangulation

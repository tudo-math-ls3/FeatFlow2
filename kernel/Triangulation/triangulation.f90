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
!# 1.) tria_wrp_tria2Structure
!#     -> Wrapper. Create a FEAT 2.0 triangulation structure from a
!#        FEAT 1.x triangulation structure.
!#
!# 2.) tria_readTriFile2D
!#     -> Reads a .TRI file and creates a 'raw' mesh with only basic 
!#        information.
!#
!# 3.) tria_initStandardMeshFromRaw
!#     -> Generates all standard arrays for a mesh, i.e. converts a 'raw' mesh
!#        (as set up by tria_readTriFile2D e.g.) to a standard mesh.
!#
!# 4.) tria_refine2LevelOrdering
!#     -> Refines a mesh according to the 2-level ordering algorithm.
!#        Creates a 'raw' fine mesh from a 'standard' coarse mesh.
!#
!# 5.) tria_compress2LevelOrdHierarchy
!#     -> Can be used to compress a mesh hierarchy created with the 2-level
!#        ordering. Saves memory. Shares vertex coordinates between fine
!#        and coarse mesh.
!#
!# 6.) tria_quickRefine2LevelOrdering
!#     -> Refines a mesh multiple times according to the 2-level ordering 
!#        algorithm. Creates a 'raw' fine mesh from a 'raw' or 'standard' 
!#        coarse mesh.
!#
!# 7.) tria_done
!#     -> Cleans up a triangulation structure, releases memory from the heap.
!#
!# 8.) tria_rawGridToTri
!#     -> Converts a raw mesh into a triangular mesh.
!#
!# 9.) tria_duplicate
!#     -> Creates a duplicate / backup of a triangulation.
!#        Some information may be shared between two triangulation structures.
!#
!# 10.) tria_restore
!#      -> Restores a triangulation previously backed up with tria_backup.
!#
!# 11.) tria_searchBoundaryNode
!#      -> Search for the position of a boundary vertex / edge on the boundary.
!#
!# 12.) tria_getPointsOnEdge
!#      -> For all edges in a triangulation, calculate the coordinates 
!#         of a number of points on each edge
!#
!# Auxiliary routines:
!#
!# 1.) tria_readRawTriangulation2D
!#     -> Reads basic data from a TRI file.
!#
!# 2.) tria_genRawBoundary2D
!#     -> Generates basic information about boundary vertices.
!#
!# 3.) tria_sortBoundaryVertices2D
!#     -> Sorts the boundary vertices for increasing parameter values.
!#
!# 4.) tria_genElementsAtVertex2D
!#     -> Generates the IelementsAtVertex array for a 2D triangulation
!#
!# 5.) tria_genNeighboursAtElement2D
!#     -> Generates the IneighboursAtElement array for a 2D triangulation
!#
!# 6.) tria_genEdgesAtElement2D
!#     -> Generates the IedgesAtElement array for a 2D triangulation
!#
!# 7.) tria_genElementsAtEdge2D
!#     -> Generates the IelementsAtEdge array for a 2D triangulation
!#
!# 8.) tria_genVerticesAtEdge2D
!#     -> Generates the IverticesAtEdge array for a 2D triangulation
!#
!# 9.) tria_genEdgeNodalProperty2D
!#     -> Generates the InodalProperty array part for all edges
!#
!# 10.) tria_genElementVolume2D
!#      -> Generates the DelementVolume array for a 2D triangulation
!#
!# 11.) tria_genElementsAtBoundary2D
!#      -> Generates the IelementsAtBoundary array for a 2D triangulation
!#
!# 12.) tria_genEdgesAtBoundary2D
!#      -> Generates the IedgesAtBoundary array for a 2D triangulation
!#
!# 13.) tria_genEdgeParameterValue2D
!#      -> Generates the DedgeParameterValue array for a 2D triangulatioon
!# 
!# 14.) tria_genBoundaryVertexPos2D
!#      -> Generates the IboundaryVertexPos2D array for a 2D triangulatioon
!#
!# 15.) tria_genBoundaryEdgePos2D
!#      -> Generates the IboundaryEdgePos2D array for a 2D triangulatioon
!#
!#
!#  FAQ - Some explainations
!# --------------------------
!# 1.) When reading the refinement routine, there's written something about
!#     'raw' meshes and 'standard' meshes. What does that mean?
!#
!#     A 'raw' mesh is a very basic form of a mesh with most information
!#     missing. E.g. there is no information about adjacencies.
!#     Such meshes come e.g. from .TRI files when reading them with
!#     tria_readTriFile2D or similar routines. These basic meshes can normally
!#     not be used for any computations; all the missing informationhas first
!#     to be 'extracted' or 'generated' based on them. They can only be used
!#     for very low-level modifications; e.g. the routine tria_rawGridToTri
!#     allowes to convert a quad mesh in a triangular mesh, which would
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
!# 2.) And what does that mean if I want to read e.g. a 2D mesh and to refine it?
!#
!#     Well, to read a mesh and prepare it to be used as single mesh, use:
!#
!#       CALL tria_readTriFile2D (rtriangulation, 'somemesh.tri', rboundary)
!#       CALL tria_initStandardMeshFromRaw (rtria,rboundary)
!#
!#     If you want to refine a mesh let's say 4 times after reading it, use:
!#
!#       CALL tria_readTriFile2D (rtriangulation, 'somemesh.tri', rboundary)
!#       CALL tria_quickRefine2LevelOrdering(4,rtriangulation,rboundary)
!#       CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)
!# 
!#     If you want to generate level 1..4 for a multiple-grid structure, use
!#
!#       TYPE(t_triangulation), DIMENSION(4) :: rtriangulation
!#
!#       CALL tria_readTriFile2D (rtriangulation(1), 'somemesh.tri', rboundary)
!#       CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)
!#       DO ilev=2,4
!#         CALL tria_refine2LevelOrdering (rtriangulation(ilev-1),&
!#             rtriangulation(ilev), rboundary)
!#         CALL tria_initStandardMeshFromRaw (rtriangulation(ilev),rboundary)
!#       END DO
!#
!# 3.) What is the tria_compress2LevelOrdHierarchy for?
!#
!#     This routine can be used after the refinement process to save memory.
!#     As the vertex coordinates of all points in the coarse mesh are the
!#     same in the fine mesh, there's redundant data! The
!#     tria_compress2LevelOrdHierarchy now removes this redundant data by
!#     releasing the vertex coordinate array on the coarse mesh and replacing
!#     it by that of the fine mesh. The routine can be used as follows:
!#
!#     Generate the triangulation:
!#
!#       TYPE(t_triangulation), DIMENSION(4) :: rtriangulation
!#
!#       CALL tria_readTriFile2D (rtriangulation(1), 'somemesh.tri', rboundary)
!#       CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)
!#       DO ilev=2,4
!#         CALL tria_refine2LevelOrdering (rtriangulation(ilev-1),&
!#             rtriangulation(ilev), rboundary)
!#         CALL tria_initStandardMeshFromRaw (rtriangulation(ilev),rboundary)
!#       END DO
!#
!#    Remove redundant data
!#
!#       DO ilev=4-1,1
!#         CALL tria_compress2LevelOrdHierarchy (rtriangulation(ilev+1),&
!#              rtriangulation(ilev))
!#       END DO
!#
!#    One should note that afterwards, each change of vertex coordinates on the
!#    fine grid directly affects the vertex coordinates on all coarse grid
!#    and vice versa. But that's not a bug, it's a feature :-)
!# 
!# </purpose>
!##############################################################################

MODULE triangulation

  USE fsystem
  USE genoutput
  USE storage
  USE basicgeometry
  USE linearalgebra
  USE boundary
  USE geometryaux
  USE sort

  IMPLICIT NONE
  
  INCLUDE 'stria.inc'

!<constants>

!<constantblock description="Triangulation constants">
  
  ! Maximum number of vertices per element supported by the triangulation
  ! module. Currently, this is 8 which is the case for 3D hexahedral meshes.
  INTEGER, PARAMETER :: TRIA_MAXNVE   = 8

  ! Maximum number of corner-vertices in each 2D element.
  ! We set this to 4 to allow triangle and quadrilateral element shapes.
  ! May be set to higher vales in the future for the support of
  ! isoparametric elements!
  ! This is the old NNVE.
  INTEGER, PARAMETER :: TRIA_MAXNVE2D = 4

  ! Maximum number of edges in each element.
  ! We set this to 4 to allow triangle and quadrilateral element shapes.
  ! This is the old NNVE, too.
  INTEGER, PARAMETER :: TRIA_MAXNME2D = TRIA_MAXNVE2D
  
  ! Number of vertices per element for triangular element shapes in 2D.
  INTEGER, PARAMETER :: TRIA_NVETRI2D  = 3

  ! Number of vertices per element for quadrilateral element shapes in 2D.
  INTEGER, PARAMETER :: TRIA_NVEQUAD2D = 4
  
!</constantblock>

!<constantblock description="KIND values for triangulation data">
  
  ! kind value for indexing the vertices in a triangulation
  INTEGER, PARAMETER :: PREC_VERTEXIDX  = I32

  ! Alternative name for PREC_VERTEXIDX
  INTEGER, PARAMETER :: PREC_POINTIDX   = PREC_VERTEXIDX

  ! kind value for indexing the edges in a triangulation
  INTEGER, PARAMETER :: PREC_EDGEIDX    = I32

  ! kind value for indexing the elements in a triangulation
  INTEGER, PARAMETER :: PREC_ELEMENTIDX = I32

!</constantblock>
  
!<constantblock description="Duplication flags. Specifies which information is shared \
!                            between triangulation structures">

  INTEGER(I32), PARAMETER :: TR_SHARE_DVERTEXCOORDS          = 2** 0  ! DCORVG 
  INTEGER(I32), PARAMETER :: TR_SHARE_DFREEVERTEXCOORDINATES = 2** 1  ! DCORMG
  INTEGER(I32), PARAMETER :: TR_SHARE_IVERTICESATELEMENT     = 2** 2  ! KVERT 
  INTEGER(I32), PARAMETER :: TR_SHARE_IEDGESATELEMENT        = 2** 3  ! KMID  
  INTEGER(I32), PARAMETER :: TR_SHARE_INEIGHBOURSATELEMENT   = 2** 4  ! KADJ  
  INTEGER(I32), PARAMETER :: TR_SHARE_IELEMENTSATVERTEX      = 2** 5  ! KVEL  
  INTEGER(I32), PARAMETER :: TR_SHARE_IELEMENTSATEDGE        = 2** 6  ! KMEL  
  INTEGER(I32), PARAMETER :: TR_SHARE_INODALPROPERTY         = 2** 7  ! KNPR  
  INTEGER(I32), PARAMETER :: TR_SHARE_KMM                    = 2** 8  ! KMM   
  INTEGER(I32), PARAMETER :: TR_SHARE_IVERTICESATBOUNDARY    = 2** 9  ! KVBD  
  INTEGER(I32), PARAMETER :: TR_SHARE_IELEMENTSATBOUNDARY    = 2**10  ! KEBD  
  INTEGER(I32), PARAMETER :: TR_SHARE_IBOUNDARYCPIDX         = 2**11  ! KBCT  
  INTEGER(I32), PARAMETER :: TR_SHARE_DVERTEXPARAMETERVALUE  = 2**12  ! DVBDP 
  INTEGER(I32), PARAMETER :: TR_SHARE_DEDGEPARAMETERVALUE    = 2**13  ! DMBDP 
  INTEGER(I32), PARAMETER :: TR_SHARE_IEDGESATBOUNDARY       = 2**14  ! KMBD  
  INTEGER(I32), PARAMETER :: TR_SHARE_IVERTICESATEDGE        = 2**16  ! KEAN  
  INTEGER(I32), PARAMETER :: TR_SHARE_IBOUNDARYVERTEXPOS     = 2**17  ! KVBDI 
  INTEGER(I32), PARAMETER :: TR_SHARE_IBOUNDARYEDGEPOS       = 2**18  ! KMBDI 
  INTEGER(I32), PARAMETER :: TR_SHARE_DELEMENTAREA           = 2**19  ! DAREA 
  
  ! Share everything
  INTEGER(I32), PARAMETER :: TR_SHARE_ALL = NOT(0_I32)

!</constantblock>
  
!</constants>


!<types>

!<typeblock>
  
  ! Each element consists of at most TRIA_MAXNVE2D points.
  ! Each point has a number, which is usually an integer value.
  
  TYPE t_elementCorners2D
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D) :: Icorners
  END TYPE

  ! Each element contains at most TRIA_MAXNME2D edges.
  ! Each edge has a number, which is usually an integer value.
  
  TYPE t_elementEdges2D
    INTEGER(PREC_EDGEIDX), DIMENSION(TRIA_MAXNME2D) :: Iedges
  END TYPE
  
  ! Each element contains at most NMAXEDGES neighbour elements,
  ! each meeting the element in an edge.
  ! Each neighbour element has a number, which is usually an integer value.
  
  TYPE t_elementNeighbours2D
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNME2D) :: Ineighbours
  END TYPE
  
!<typeblock>
  
!</typeblock>
  
  ! The basic triangulation structure for a triangulation.
  TYPE t_triangulation
  
    ! The 'old' triangulation structure for compatibility.
    ! This was just an array of length SZTRIA emulating a structure,
    ! where the entries could be accessed by means of the Oxxxx constants.
    ! !!!THIS WILL BE DELETED IN THE FINAL VERSION!!!
    INTEGER, DIMENSION(SZTRIA) :: Itria
    
    ! Duplication flag. Bitfield. Used by TRIDUP/TRIRST/TRIDEL to 
    ! mark which information of a triangulation structure 
    ! coincides with those of another triangulation structure and
    ! must not be released from memory when deleting a copy of a
    ! triangulation structure.
    ! When a bit is set to 1, the corresponding array is
    ! maintained by another triangulation structure and must not
    ! be deleted by TRIDEL. When the bit is 0, the array is a real
    ! copy of another array and must be deleted in TRIDEL.
    ! Bit  0: DCORVG is a copy of another structure
    ! Bit  1: DCORMG is a copy of another structure
    ! Bit  2: KVERT  is a copy of another structure
    ! Bit  3: KMID   is a copy of another structure
    ! Bit  4: KADJ   is a copy of another structure
    ! Bit  5: KVEL   is a copy of another structure
    ! Bit  6: KMEL   is a copy of another structure
    ! Bit  7: KNPR   is a copy of another structure
    ! Bit  8: KMM    is a copy of another structure
    ! Bit  9: KVBD   is a copy of another structure
    ! Bit 10: KEBD   is a copy of another structure
    ! Bit 11: KBCT   is a copy of another structure
    ! Bit 12: DVBDP  is a copy of another structure
    ! Bit 13: DMBDP  is a copy of another structure
    ! Bit 14: KMBD   is a copy of another structure
    ! Bit 16: KEAN   is a copy of another structure
    ! Bit 17: KVBDI  is a copy of another structure
    ! Bit 18: KMBDI  is a copy of another structure
    ! Bit 19: DAREA  is a copy of another structure
    INTEGER(I32)             :: iduplicationFlag
  
    ! Dimension of the triangulation.
    ! NDIM2D=2D triangulation, NDIM2D=3D triangulation, 0=not initialised
    INTEGER                  :: ndim = 0
  
    ! Number of points in the domain corresponding to corners of elements;
    ! coincides with SIZE(RcornerCoordinates)
    INTEGER(PREC_VERTEXIDX)   :: NVT = 0
    
    ! Number of edges in the domain belonging to elements
    INTEGER(PREC_EDGEIDX)    :: NMT = 0
    
    ! Number of elements in the domain; 
    ! corresponding to SIZE(RverticesOnElement)
    INTEGER(PREC_ELEMENTIDX) :: NEL = 0
    
    ! Number of boundary components
    INTEGER             :: NBCT = 0
    
    ! Number of vertices on the boundary.
    ! For 2D domains, this coincides with the number of edges on 
    ! the boundary: Every vertex on the boundary has an edge following
    ! the vertex in mathematical positive sense.
    INTEGER             :: NVBD = 0
    
    ! Number of edges on the boundary; coincides with NVBD
    ! for 2D domains.
    INTEGER             :: NMBD = 0
    
    ! Number of elements with a defined number of vertices per element.
    ! InelOfType(TRIA_NVETRI2D)  = number of triangles in the mesh (2D).
    ! InelOfType(TRIA_NVEQUAD2D) = number of quads in the mesh (2D).
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE) :: InelOfType = 0
  
    ! Number of vertices per edge; normally = 0.
    ! If a regular distribution of vertices on edges is given,
    ! NVPED saves the number of vertices on each edge;
    ! e.g. 1 if midpoints on edges exist in p_DfreeVertexCoordinates.
    INTEGER             :: nverticesPerEdge = 0
    
    ! Total number of vertices on edges; normally = 0.
    ! Total number of vertices on all edges, realized in p_DfreeVertexCoordinates. 
    ! E.g. if midpoints on edges exist, there is NVEDT=NMT.
    INTEGER(PREC_VERTEXIDX)   :: nVerticesOnAllEdges = 0
    
    ! Number of inner-element vertices; normally = 0.
    ! If a regular distribution of vertices in the inner of 
    ! each element is given, NIELV saves the number of vertices 
    ! in the inner of each element; e.g. 1 if element-midpoints 
    ! exist in p_DfreeVertexCoordinates.
    INTEGER             :: nverticesInEachElement = 0

    ! Total number of vertices in elements; normally = 0.
    ! Total number of vertices on all elements, realized in p_DfreeVertexCoordinates. 
    ! E.g. if element-midpoints exist in DCORMG, there is NIEVT=NEL.
    INTEGER(PREC_VERTEXIDX)   :: nverticesInAllElements = 0
    
    ! Number of additional vertices; normally = 0.
    ! Can be set <> 0 if there are any additional vertices 
    ! realized in p_DfreeVertexCoordinates, that don't belong to a regular 
    ! distribution of vertices in corners, on edges or on elements.
    INTEGER(PREC_VERTEXIDX)    :: nadditionalVertices = 0
  
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..ndim,1..NVT] of double
    ! with
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    ! fo 2D meshes and
    !   p_DvertexCoords(1,.) = X-coordinate.
    !   p_DvertexCoords(2,.) = Y-coordinate.
    !   p_DvertexCoords(3,.) = Y-coordinate.
    ! for 3D meshes.
    ! This is a handle to the old DCORVG-array.
    !
    ! Note that the array may be longer than NVT in general!
    ! (May happen in case of a mesh hierarchy generated by a 2-level
    ! refinement, where the coordinates of the points on the
    ! coarser levels are contained in te coordinates of the
    ! finer levels.)
    ! In such a case, only the first NVT n-tuples in this array are valid!
    INTEGER        :: h_DvertexCoords = ST_NOHANDLE
    
    ! Vertices Adjacent to an Element.
    ! Handle to h_IverticesAtElement=array [1..TRIA_MAXNVE2D,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! This is a handle to the old KVERT array.
    INTEGER        :: h_IverticesAtElement = ST_NOHANDLE

    ! Edges Adjacent to an Element.
    ! Handle to 
    !       p_RedgesAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element the node numbers of the edges following the
    ! corner vertices in mathematically positive sense.
    ! This is the old KMID array.
    ! To be able to distinguish a number of an edge from a vertex number, 
    ! edges are numbered in the range NVT+1..NVT+NMT. 
    INTEGER        :: h_IedgesAtElement = ST_NOHANDLE
    
    ! Neighbour Elements Adjacent to an Element.
    ! Handle to 
    !       p_RneighboursAtElement = array [1..TRIA_MAXNME2D,1..NEL] of integer
    ! For each element, the numbers of adjacent elements
    ! in mathematically positive sense, meeting the element in an edge.
    ! p_RneighbourElement(IEL)\%Ineighbours(.) describes the elements adjacent 
    ! to IEL along the edges (p_RedgesOnElement(IEL)\%Iedges(.)-NVT).
    ! This is the old KADJ array.
    INTEGER        :: h_IneighboursAtElement = ST_NOHANDLE
    
    ! Elements Adjacent to an Edge. 
    ! Handle to 
    !       p_IelementsAtEdge = array [1..2,1..NMT] of integer.
    ! The numbers of the two elements adjacent to an edge IMT in 2D. 
    ! For boundary edges, p_IelementsOnEdge(2,IMT) is set to 0.
    ! This is the old KMEL array.
    INTEGER        :: h_IelementsAtEdge = ST_NOHANDLE

    ! Vertices Adjacent to an Edge. 
    ! Handle to 
    !       p_IverticesAtEdge = array [1..2,1..NMT]
    ! The numbers of the two vertices adjacent to an edge IMT. 
    INTEGER        :: h_IverticesAtEdge = ST_NOHANDLE
    
    ! Nodal property array. 
    ! Handle to 
    !       p_InodalProperty=array [1..NVT+NMT] of integer.
    ! p_InodalProperty(i) defines for each vertex i=(1..NVT) 
    ! and each edge i=(NVT+1..NVT+NMT) its function inside of the
    ! geometry. Generally said, the range of the p_InodalProperty-array 
    ! characterizes the type of the node (=vertex/edge):
    ! = 0    : The vertex/edge is an inner vertex/edge
    ! > 0    : The vertex/edge is a boundary vertex/edge on the real
    !           boundary. KNPR(.) defines the number of the boundary
    !           component.
    ! < 0,
    ! >= -NEL: The vertex/edge is an invalid node to element -KNPR(.)
    !          (-> hanging node, not implemented for now)
    ! This is the old KNPR-array, slightly modified for edges!
    INTEGER         :: h_InodalProperty = ST_NOHANDLE
    
    ! 2D triangulation: Array with area of each element.
    ! 3D triangulation: Array with volume of each element.
    ! Handle to 
    !       p_DelementArea = array [1..NEL+1] of double.
    ! p_DelementArea [NEL+1] gives the total area/voloume of the domain.
    INTEGER         :: h_DelementVolume = ST_NOHANDLE
    
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
    INTEGER        :: h_IelementsAtVertexIdx = ST_NOHANDLE
    
    ! Array containing the Elements Adjacent to a Vertex.
    ! Handle to 
    !       p_IelementsAtVertex = array(1..*) of integer
    ! p_IelementsAtVertex ( p_IelementsAtVertexIdx(IVT)..p_IelementsAtVertexIdx(IVT+1)-1 )
    ! contains the number of the adjacent element in a vertex.
    ! This replaces the old KVEL array.
    INTEGER        :: h_IelementsAtVertex = ST_NOHANDLE
    
    ! Boundary component index vector of length NBCT+1 for 
    ! p_IverticesAtBoundary / p_IedgesAtBoundary / ... arrays.
    ! Handle to 
    !      p_IboundaryCpIdx = array [1..NBCT+1] of integer.
    ! For a (real) boundary component i all corner nodes 
    ! for that boundary component are saved in 
    !   p_IverticesAtBoundary ( p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1 ).
    ! All boundary edges of this boundary component are saved in 
    !   p_IedgesAtBoundary ( p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1 ).
    ! p_IboundaryCpIdx(NBCT+1) points to NVBD+1 for easier access to the last
    ! boundary component.
    ! This is the old KBCT array.
    INTEGER          :: h_IboundaryCpIdx = ST_NOHANDLE

    ! Vertices on boundary. 
    ! Handle to 
    !       p_IverticesAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all vertices on the (real) boundary
    ! in mathematically positive sense.
    ! The boundary vertices of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KVBD array.
    INTEGER          :: h_IverticesAtBoundary = ST_NOHANDLE

    ! Edges Adjacent to the boundary. 
    ! Handle to 
    !       p_IedgesAtBoundary = array [1..NVBD] of integer.
    ! This array contains a list of all edges on the (real) boundary
    ! in mathematically positive sense.
    ! The boundary edges of boundary component i are saved at
    !        p_IboundaryCpIdx(i)..p_IboundaryCpIdx(i+1)-1.
    ! This is the old KMBD array.
    INTEGER          :: h_IedgesAtBoundary = ST_NOHANDLE

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
    INTEGER          :: h_IelementsAtBoundary = ST_NOHANDLE
    
    ! Parameter values of vertices on the boundary.
    ! Handle to 
    !       p_DvertexParameterValue = array [1..NVBD] of real
    ! p_DvertexParameterValue(i) contains the parameter value of 
    ! boundary vertex i, which corresponds to the corner vertex 
    ! p_IverticesAtBoundary(I).
    INTEGER          :: h_DvertexParameterValue = ST_NOHANDLE

    ! Parameter values of edge midpoints on the boundary.
    ! Handle to 
    !       p_DedgeParameterValue = array [1..NMBD] of real
    ! p_DedgeParameterValue(i) contains the parameter value of
    ! the midpoint of boundary edge i, which corresponds to the 
    ! edge p_IedgesAtBoundary(I).
    INTEGER          :: h_DedgeParameterValue = ST_NOHANDLE
    
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
    INTEGER          :: h_IboundaryVertexPos = ST_NOHANDLE

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
    INTEGER          :: h_IboundaryEdgePos = ST_NOHANDLE
    
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
    ! don't belong to any of the previous two groups.
    !
    ! The numbers of regularly distributed vertices on edges/elements can
    ! be calculated directly with an appropriate formula. E.g. let's assume,
    ! we have n regularly distributed vertices on each edge (excluding the
    ! starting/ending point). Then the corresponding vertices on 
    ! egde E (=1..NMT) have the numbers:
    !       (NVT-1)+(NMT-1) + (E-1)*n + 1
    !    .. (NVT-1)+(NMT-1) + (E-1)*n + n
    ! the formula for regular distributed vertices on elements is similar.
    INTEGER           :: h_DfreeVertexCoordinates = ST_NOHANDLE
    
  END TYPE

!</typeblock>
  
!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_wrp_tria2Structure (TRIA, rtriangulation)
  
  USE afcutil
  
!<description>
  ! Wrapper routine. Accepts an 'old' triangulation structure array of CC2D
  ! and converts it completely to a triangulation structure. All 'old'
  ! information in the triangulation structure is overwritten.
  ! If rtriangulation contains an old structure, the existing structure 
  ! is recreated/updated.
!</description>
  
!<input>
  ! The old triangulation structure array that should be converted to
  ! the new triangulation structure.
  INTEGER, DIMENSION(SZTRIA), INTENT(IN) :: TRIA
!</input>
  
!<inputoutput>
  ! The triangulation structure which will be overwritten by the information
  ! in TRIA.
  TYPE(t_triangulation), INTENT(INOUT)      :: rtriangulation
!</inputoutput>
  
!</subroutine>

  ! Do we have a structure?
  IF (rtriangulation%ndim .NE. 0) THEN
    ! Release the old structure without removing it from the heap.
    CALL tria_done (rtriangulation)
  END IF

  ! We take the ownership of the triangulation structure - for compatibility
  ! to old FEAT 1.x:
  rtriangulation%Itria = TRIA
  
  ! Copy static entries
  rtriangulation%NVT                      = TRIA(ONVT  )
  rtriangulation%NMT                      = TRIA(ONMT  )
  rtriangulation%NEL                      = TRIA(ONEL  )
  rtriangulation%NBCT                     = TRIA(ONBCT )
  rtriangulation%NVBD                     = TRIA(ONVBD )
  rtriangulation%NMBD                     = TRIA(ONVBD ) ! NMBD=NVBD !!!
  rtriangulation%nverticesPerEdge         = TRIA(ONVPED)
  rtriangulation%nVerticesOnAllEdges      = TRIA(ONVEDT)
  rtriangulation%nverticesInEachElement   = TRIA(ONVEDT)
  rtriangulation%nverticesInAllElements   = TRIA(ONIEVT)
  rtriangulation%nadditionalVertices      = TRIA(ONANT )
  
  ! Set ndim <> 0; this declares the structure as 'initialised'.
  rtriangulation%ndim = NDIM2D
  
  ! The duplication flag stays at 0 - as for now, all
  ! arrays are created from the feat arrays as new arrays, and
  ! so they are not a copy of another array.
  
  rtriangulation%iduplicationFlag = 0
  
  ! *******************************************************
  ! Copy DCORVG, create p_RcornerCoordinates.
  
  CALL copy_featarray_double2d ('DCORVG',2,INT(rtriangulation%NVT),TRIA(OLCORVG),&
                                rtriangulation%h_DvertexCoords)
  
  ! *******************************************************
  ! Copy KVERT, create p_RverticesAtElement.
  
  CALL copy_featarray_int2d ('KVERT',4,INT(rtriangulation%NEL),TRIA(OLVERT),&
                             rtriangulation%h_IverticesAtElement)

  ! *******************************************************
  ! Copy KMID, create p_RedgesAtElement.
  
  CALL copy_featarray_int2d ('KMID',4,INT(rtriangulation%NEL),TRIA(OLMID),&
                             rtriangulation%h_IedgesAtElement)

  ! *******************************************************
  ! Copy KADJ, create p_RneighboursAtElement.
  
  CALL copy_featarray_int2d ('KADJ',4,INT(rtriangulation%NEL),TRIA(OLADJ),&
                             rtriangulation%h_IneighboursAtElement)

  ! *******************************************************
  ! Copy KMEL, create p_IelementsAtEdge.
  
  CALL copy_featarray_int2d ('KMEL',2,INT(rtriangulation%NMT),TRIA(OLMEL),&
                             rtriangulation%h_IelementsAtEdge)

  ! *******************************************************
  ! Copy KEAN, create p_IverticesAtEdge.
  
  CALL copy_featarray_int2d ('KEAN',2,INT(rtriangulation%NMT),TRIA(OLEAN),&
                             rtriangulation%h_IverticesAtEdge)

  ! *******************************************************
  ! Copy KNPR, create p_InodalProperty.
  
  CALL copy_featarray_int1d ('KNPR',INT(rtriangulation%NVT+rtriangulation%NMT),&
                             TRIA(OLNPR),rtriangulation%h_InodalProperty)
                             
  ! *******************************************************
  ! Copy KAREA, create hpDelementArea.
  
  CALL copy_featarray_double1d ('KAREA',INT(rtriangulation%NEL+1),&
                                TRIA(OLAREA),rtriangulation%h_DelementVolume)
                             
  ! *******************************************************
  ! Initialise the new KVEL, create p_IelementsAtVertexIdx/p_IelementsAtVertex.
                             
  CALL translate_KVEL (TRIA(ONVEL),INT(rtriangulation%NVT),TRIA(OLVEL), &
                       rtriangulation%h_IelementsAtVertex,&
                       rtriangulation%h_IelementsAtVertexIdx)

  ! *******************************************************
  ! Copy KBCT, create p_IboundaryCpIdx.
  
  CALL copy_featarray_int1d ('KBCT',(rtriangulation%NBCT+1),&
                             TRIA(OLBCT),rtriangulation%h_IboundaryCpIdx)
  
  ! *******************************************************
  ! Copy KVBD, create p_IverticesAtBoundary.
  
  CALL copy_featarray_int1d ('KVBD',INT(rtriangulation%NVBD),&
                             TRIA(OLVBD),rtriangulation%h_IverticesAtBoundary)

  ! *******************************************************
  ! Copy KMBD, create p_IedgesAtBoundary.
  
  CALL copy_featarray_int1d ('KMBD',(rtriangulation%NMBD),&
                             TRIA(OLMBD),rtriangulation%h_IedgesAtBoundary)

  ! *******************************************************
  ! Copy KEBD, create p_IelementsAtBoundary.
  
  CALL copy_featarray_int1d ('KEBD',INT(rtriangulation%NMBD),&
                             TRIA(OLEBD),rtriangulation%h_IelementsAtBoundary)

  ! *******************************************************
  ! Copy DVBDP, create p_DvertexParameterValue.
  
  CALL copy_featarray_double1d ('DVBDP',INT(rtriangulation%NVBD),&
                                TRIA(OLVBDP),rtriangulation%h_DvertexParameterValue)

  ! *******************************************************
  ! Copy DMBDP, create p_DedgeParameterValue.
  
  CALL copy_featarray_double1d ('DMBDP',INT(rtriangulation%NMBD),&
                                TRIA(OLMBDP),rtriangulation%h_DedgeParameterValue)


  ! *******************************************************
  ! Copy DCORMG, create p_DfreeVertexCoordinates.
  
  CALL copy_featarray_double2d ('DCORMG',2,INT(rtriangulation%NMT),TRIA(OLCORMG),&
                                rtriangulation%h_DfreeVertexCoordinates)
  

  CONTAINS
  
    ! Now the COPY-sub-subroutines used above.
  
    ! *************************************************************************  
    ! Copies the FEAT array with feat-handle ifeathandle to the FEAT2.0
    ! array identified by ihandle. If the size associated to ihandle is
    ! wrong, the array is reallocated.
    ! The FEAT array is assumed to have the shape 
    !   array [1..idim1,1..idim2] of integer
    
    SUBROUTINE copy_featarray_int2d (name,idim1,idim2,ifeathandle,ihandle)
    
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: idim1, idim2, ifeathandle
    INTEGER, INTENT(INOUT) :: ihandle
    
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_array
    !INTEGER, DIMENSION(:), POINTER :: p_array2
    
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(I32) :: j,i,kpos
    
    INCLUDE 'cmem.inc'
    
    ! Clear the array if no FEAT handle assigned.
    
    IF ((ifeathandle .EQ. 0) .AND. (ihandle .NE. ST_NOHANDLE)) THEN
      CALL storage_free (ihandle)
      RETURN
    END IF
    
    Isize = (/idim1,idim2/)
    
    ! Do we have to reallocate?
    
    IF (ihandle .NE. ST_NOHANDLE) THEN
    
      CALL storage_getbase_int2D (ihandle,p_array)
    
      IF (SIZE(p_array,2) .NE. idim2) THEN
        ! Size wrong. Deallocate the old array, create a new one.
        CALL storage_free (ihandle)
        CALL storage_new2D ('tria_wrp_tria2Structure', name, &
                            Isize, ST_INT, ihandle, &
                            ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the array
      CALL storage_new2D ('tria_wrp_tria2Structure', name, &
                          Isize, ST_INT, ihandle, &
                          ST_NEWBLOCK_NOINIT)
    END IF
    CALL storage_getbase_int2D (ihandle,p_array)

    ! Copy the FEAT array
    !p_array2 => feat_htpint(idim1*idim2,ifeathandle)
    kpos = L(ifeathandle)
    !p_array2 => KWORK(L(ifeathandle):L(ifeathandle)+idim1*idim2-1)

    DO j=0,idim2-1
      DO i=0,idim1-1
        p_array(i+1,j+1) = KWORK(kpos+idim1*j+i)
      END DO
    END DO
    
    END SUBROUTINE

    ! *************************************************************************
    ! Copies the FEAT array with feat-handle ifeathandle to the FEAT2.0
    ! array identified by ihandle. If the size associated to ihandle is
    ! wrong, the array is reallocated.
    ! The FEAT array is assumed to have the shape 
    !   array [1..idim1] of integer
    
    SUBROUTINE copy_featarray_int1d (name,idim1,ifeathandle,ihandle)
    
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: idim1,ifeathandle
    INTEGER, INTENT(INOUT) :: ihandle
    
    INTEGER(I32), DIMENSION(:), POINTER :: p_array
    !INTEGER, DIMENSION(:), POINTER :: p_array2
    INTEGER(I32) :: i,kpos
    
    INCLUDE 'cmem.inc'
    
    ! Clear the array if no FEAT handle assigned.
    
    IF ((ifeathandle .EQ. 0) .AND. (ihandle .NE. ST_NOHANDLE)) THEN
      CALL storage_free (ihandle)
      RETURN
    END IF
    
    ! Do we have to reallocate?

    IF (ihandle .NE. ST_NOHANDLE) THEN
    
      CALL storage_getbase_int (ihandle,p_array)
    
      IF (SIZE(p_array) .NE. idim1) THEN
        ! Size wrong. Deallocate the old array, create a new one.
        CALL storage_free (ihandle)
        CALL storage_new ('tria_wrp_tria2Structure', name, &
                          INT(idim1,I32), ST_INT, ihandle, &
                          ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the array
      CALL storage_new ('tria_wrp_tria2Structure', name, &
                        INT(idim1,I32),ST_INT, ihandle, &
                        ST_NEWBLOCK_NOINIT)
    END IF
    CALL storage_getbase_int (ihandle,p_array)

    ! Copy the FEAT array
    !p_array2 => feat_htpint(idim1,ifeathandle)
    kpos = L(ifeathandle)
    !p_array2 => KWORK(L(ifeathandle):L(ifeathandle)+idim1-1)
    
    DO i=0,idim1-1
      p_array(i+1) = KWORK(kpos+i) !p_array2(i)
    END DO
    
    END SUBROUTINE

    ! *************************************************************************
    ! Copies the FEAT array with feat-handle ifeathandle to the FEAT2.0
    ! array identified by ihandle. If the size associated to ihandle is
    ! wrong, the array is reallocated.
    ! The FEAT array is assumed to have the shape 
    !   array [1..idim1,1..idim2] of double
    
    SUBROUTINE copy_featarray_double2d (name,idim1,idim2,ifeathandle,ihandle)
    
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: idim1, idim2, ifeathandle
    INTEGER, INTENT(INOUT) :: ihandle
    
    REAL(DP), DIMENSION(:,:), POINTER :: p_array
    !DOUBLE PRECISION, DIMENSION(:), POINTER :: p_array2
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(I32) :: j,i,kpos
    
    INCLUDE 'cmem.inc'
    
    ! Clear the array if no FEAT handle assigned.
    
    IF ((ifeathandle .EQ. 0) .AND. (ihandle .NE. ST_NOHANDLE)) THEN
      CALL storage_free (ihandle)
      RETURN
    END IF
    
    Isize = (/idim1,idim2/)
    
    ! Do we have to reallocate?

    IF (ihandle .NE. ST_NOHANDLE) THEN
    
      CALL storage_getbase_double2D (ihandle,p_array)
    
      IF (SIZE(p_array,2) .NE. idim2) THEN
        ! Size wrong. Deallocate the old array, create a new one.
        CALL storage_free (ihandle)
        CALL storage_new2D ('tria_wrp_tria2Structure', name, &
                            ISize, ST_DOUBLE, ihandle, &
                            ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the array
      CALL storage_new2D ('tria_wrp_tria2Structure', name, &
                          ISize, ST_DOUBLE, ihandle, &
                          ST_NEWBLOCK_NOINIT)
    END IF
    CALL storage_getbase_double2D (ihandle,p_array)

    ! Copy the FEAT array
    !p_array2 => feat_htpdouble(idim1*idim2,ifeathandle)
    kpos = L(ifeathandle)
    !p_array2 => DWORK(L(ifeathandle):L(ifeathandle)+idim1*idim2-1)
    
    DO j=0,idim2-1
      DO i=0,idim1-1
        p_array(i+1,j+1) = DWORK(kpos+idim1*j+i)
      END DO
    END DO
    
    END SUBROUTINE

    ! *************************************************************************
    ! Copies the FEAT array with feat-handle ifeathandle to the FEAT2.0
    ! array identified by ihandle. If the size associated to ihandle is
    ! wrong, the array is reallocated.
    ! The FEAT array is assumed to have the shape 
    !   array [1..idim1] of integer
    
    SUBROUTINE copy_featarray_double1d (name,idim1,ifeathandle,ihandle)
    
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: idim1,ifeathandle
    INTEGER, INTENT(INOUT) :: ihandle
    
    REAL(DP), DIMENSION(:), POINTER :: p_array
    !DOUBLE PRECISION, DIMENSION(:), POINTER :: p_array2
    
    INTEGER :: i,kpos
    
    INCLUDE 'cmem.inc'
    
    ! Clear the array if no FEAT handle assigned.
    
    IF ((ifeathandle .EQ. 0) .AND. (ihandle .NE. ST_NOHANDLE)) THEN
      CALL storage_free (ihandle)
      RETURN
    END IF
    
    ! Do we have to reallocate?

    IF (ihandle .NE. ST_NOHANDLE) THEN
    
      CALL storage_getbase_double (ihandle,p_array)
    
      IF (SIZE(p_array) .NE. idim1) THEN
        ! Size wrong. Deallocate the old array, create a new one.
        CALL storage_free (ihandle)
        CALL storage_new ('tria_wrp_tria2Structure', name, &
                          INT(idim1,I32), ST_DOUBLE, ihandle, &
                          ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the array
      CALL storage_new ('tria_wrp_tria2Structure', name, &
                        INT(idim1,I32),ST_DOUBLE, ihandle, &
                        ST_NEWBLOCK_NOINIT)
    END IF
    CALL storage_getbase_double (ihandle,p_array)

    ! Copy the FEAT array
    !p_array2 => feat_htpdouble(idim1,ifeathandle)
    kpos = L(ifeathandle)
    !p_array2 => DWORK(L(ifeathandle):L(ifeathandle)+idim1-1)
    DO i=0,idim1-1
      p_array(i+1) = DWORK(kpos+i)
    END DO
    
    END SUBROUTINE

    ! *************************************************************************
    ! Builds a new KADJ structure, translates the old KVEL.
    ! ihandle represents the translated array, ihandleidx corresponds
    ! to the index array inside the new KVEL.
    
    SUBROUTINE translate_KVEL (NVEL,NVT,ifeathandle,ihandle,ihandleidx)
    
    INTEGER, INTENT(IN) :: NVEL,NVT, ifeathandle
    INTEGER, INTENT(INOUT) :: ihandle,ihandleidx
    
    INTEGER(I32), DIMENSION(:), POINTER :: p_array, p_arrayidx
    !INTEGER, DIMENSION(:), POINTER :: p_kadj
    INTEGER :: i,j, nentries, kpos
    
    INCLUDE 'cmem.inc'
    
    ! Clear the array if no FEAT handle assigned.
    
    IF (ifeathandle .EQ. 0) THEN
      IF (ihandle .NE. ST_NOHANDLE) THEN
        CALL storage_free (ihandle)
      END IF
      IF (ihandleidx .NE. ST_NOHANDLE) THEN
        CALL storage_free (ihandleidx)
      END IF
      RETURN
    END IF
    
    ! Do we have to reallocate?

    ! Get the KADJ array - as 1D representation of length NVEL*NVT
    !p_kadj => feat_htpint(NVEL*NVT,ifeathandle)
    kpos = L(ifeathandle)
    
    ! Count the number of entries in the translated KADJ array.
    ! This is the number of nonzero entries in KADJ - and at least 1
    ! (if there's only one cell...)
    nentries = 0
    DO i=1,NVEL*NVT
      IF (KWORK(kpos+i-1) .NE. 0) nentries = nentries+1
      !IF (p_kadj(i) .NE. 0) nentries = nentries+1
    END DO
    nentries = MAX(1,nentries)
    
    ! Check the existance of the translated array
    
    IF (ihandle .NE. ST_NOHANDLE) THEN
    
      CALL storage_getbase_int (ihandle,p_array)
    
      IF (SIZE(p_array) .NE. nentries) THEN
        ! Size wrong. Deallocate the old array, create a new one.
        CALL storage_free (ihandle)
        IF (ihandleidx .NE. 0) CALL storage_free (ihandleidx)
        CALL storage_new ('tria_wrp_tria2Structure', 'KADJ', &
                          INT(nentries,I32), ST_INT, ihandle, &
                          ST_NEWBLOCK_NOINIT)
        CALL storage_new ('tria_wrp_tria2Structure', 'KADJIDX', &
                          INT(NVT+1,I32), ST_INT, ihandleidx, &
                          ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the new array
      IF (ihandleidx .NE. 0) CALL storage_free (ihandleidx)
      CALL storage_new ('tria_wrp_tria2Structure', 'KADJ', &
                        INT(nentries,I32), ST_INT, ihandle, &
                        ST_NEWBLOCK_NOINIT)
      CALL storage_new ('tria_wrp_tria2Structure', 'KADJIDX', &
                        INT(NVT+1,I32), ST_INT, ihandleidx, &
                        ST_NEWBLOCK_NOINIT)
    END IF
    CALL storage_getbase_int (ihandle,p_array)
    CALL storage_getbase_int (ihandleidx,p_arrayidx)

    ! *************************************************************************
    ! Copy the FEAT array, set up the translated KADJ.
    ! h_IelementsAtVertex receives the entries of KADJ.
    ! h_IelementsAtVertexIdx receives indices in h_IelementsAtVertex;
    ! more precisely, h_IelementsAtVertexIdx(i) points to the first
    ! element adjacent to vertex i in h_IelementsAtVertex.
    ! h_IelementsAtVertex( h_IelementsAtVertexIdx(i) .. h_IelementsAtVertexIdx(i+1)-1 )
    ! is then a list of all adjacent elements to vertex i.
    nentries = 0
    DO i=0,NVT-1
      
      ! Build the index array
      p_arrayidx(i+1) = nentries+1
      
      ! Copy the entries
      DO j=0,NVEL-1
        !IF (p_kadj(i*NVEL+j) .NE. 0) THEN
        IF (KWORK(kpos+i*NVEL+j) .NE. 0) THEN
          nentries = nentries+1
          p_array(nentries) = KWORK(kpos+i*NVEL+j) ! p_kadj(i*NVEL+j)
        ELSE
          EXIT
        END IF
      END DO
    END DO
    
    ! Set up last element in the index array.
    p_arrayidx(NVT+1) = nentries+1
    
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_duplicate (rtriangulation,rbackupTriangulation,&
                             iduplicationFlag,bupdate)
  
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
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
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
  INTEGER(I32), INTENT(IN)          :: iduplicationFlag
  
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
  !    from rtriangulation to rbackupTriangulation. It is OR'ed
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
  LOGICAL, INTENT(IN), OPTIONAL     :: bupdate
!</input>

!<inputoutput>
  ! The "destination" discretisation structure receives the information.
  ! Depending on iduplicationFlag, the arrays are duplicated or shared
  ! between rtriangulation and rbackupTriangulation
  TYPE(t_triangulation), INTENT(INOUT), TARGET :: rbackupTriangulation
!</inputoutput>
  
!</subroutine>
  
    ! local variables
    INTEGER(I32) :: idupFlag
    
    LOGICAL :: bupd
    
    bupd = .FALSE.
    IF (PRESENT(bupdate)) bupd = bupdate

    IF (.NOT. bupd) THEN
      ! Release any old data.
      CALL tria_done (rbackupTriangulation)
      
      rbackupTriangulation%ndim                   = rtriangulation%ndim                  
      rbackupTriangulation%NVT                    = rtriangulation%NVT                   
      rbackupTriangulation%NMT                    = rtriangulation%NMT                   
      rbackupTriangulation%NEL                    = rtriangulation%NEL                   
      rbackupTriangulation%NBCT                   = rtriangulation%NBCT                  
      rbackupTriangulation%NVBD                   = rtriangulation%NVBD                  
      rbackupTriangulation%NMBD                   = rtriangulation%NMBD                  
      rbackupTriangulation%nverticesPerEdge       = rtriangulation%nverticesPerEdge      
      rbackupTriangulation%nVerticesOnAllEdges    = rtriangulation%nVerticesOnAllEdges   
      rbackupTriangulation%nverticesInEachElement = rtriangulation%nverticesInEachElement
      rbackupTriangulation%nverticesInAllElements = rtriangulation%nverticesInAllElements
      rbackupTriangulation%nadditionalVertices    = rtriangulation%nadditionalVertices   
      rbackupTriangulation%InelOfType(:)          = rtriangulation%InelOfType(:)
      
      ! Decide on IDPFLG which arrays to copy
      rbackupTriangulation%iduplicationFlag = iduplicationFlag
      idupFlag = iduplicationFlag
      
    ELSE

      ! Create a bitfield what to copy by ORing iduplicationFlag with what
      ! we have in rbackupTriangulation. That way, only arrays that exist as
      ! real duplicates are copied from rtriangulation to rbackupTriangulation.
      
      idupFlag = IOR(iduplicationFlag,rbackupTriangulation%iduplicationFlag)
    
    END IF
    
    ! Call checkAndCopy for all the arrays. this will either copy the handle
    ! or allocate new memory and copy the content of the array.
    
    ! Bit  0: DCORVG 
    CALL checkAndCopy(idupflag, TR_SHARE_DVERTEXCOORDS,&
          rtriangulation%h_DvertexCoords, &
          rbackupTriangulation%h_DvertexCoords)

    ! Bit  2: KVERT  
    CALL checkAndCopy(idupflag, TR_SHARE_IVERTICESATELEMENT,&
          rtriangulation%h_IverticesAtElement, &
          rbackupTriangulation%h_IverticesAtElement)

    ! Bit  3: KMID   
    CALL checkAndCopy(idupflag, TR_SHARE_IEDGESATELEMENT,&
          rtriangulation%h_IedgesAtElement, &
          rbackupTriangulation%h_IedgesAtElement)
    
    ! Bit  4: KADJ   
    CALL checkAndCopy(idupflag, TR_SHARE_INEIGHBOURSATELEMENT,&
          rtriangulation%h_IneighboursAtElement, &
          rbackupTriangulation%h_IneighboursAtElement)

    ! Bit  6: KMEL   
    CALL checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge, &
          rbackupTriangulation%h_IelementsAtEdge)

    ! Bit 16: KEAN   
    CALL checkAndCopy(idupflag,TR_SHARE_IVERTICESATEDGE,&
          rtriangulation%h_IverticesAtEdge, &
          rbackupTriangulation%h_IverticesAtEdge)

    ! Bit  7: KNPR   
    CALL checkAndCopy(idupflag, TR_SHARE_INODALPROPERTY,&
          rtriangulation%h_InodalProperty, &
          rbackupTriangulation%h_InodalProperty)

    ! Bit 19: DAREA  
    CALL checkAndCopy(idupflag,TR_SHARE_DELEMENTAREA,&
          rtriangulation%h_DelementVolume, &
          rbackupTriangulation%h_DelementVolume)

    ! Bit  5: KVEL   
    CALL checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertexIdx, &
          rbackupTriangulation%h_IelementsAtVertexIdx)
    CALL checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertex, &
          rbackupTriangulation%h_IelementsAtVertex)

    ! Bit 11: KBCT   
    CALL checkAndCopy(idupflag,TR_SHARE_IBOUNDARYCPIDX,&
          rtriangulation%h_IboundaryCpIdx, &
          rbackupTriangulation%h_IboundaryCpIdx)

    ! Bit  9: KVBD   
    CALL checkAndCopy(idupflag, TR_SHARE_IVERTICESATBOUNDARY,&
          rtriangulation%h_IverticesAtBoundary, &
          rbackupTriangulation%h_IverticesAtBoundary)

    ! Bit 14: KMBD   
    CALL checkAndCopy(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IedgesAtBoundary, &
          rbackupTriangulation%h_IedgesAtBoundary)

    ! Bit 10: KEBD   
    CALL checkAndCopy(idupflag,TR_SHARE_IELEMENTSATBOUNDARY,&
          rtriangulation%h_IelementsAtBoundary, &
          rbackupTriangulation%h_IelementsAtBoundary)

    ! Bit 12: DVBDP  
    CALL checkAndCopy(idupflag,TR_SHARE_DVERTEXPARAMETERVALUE,&
          rtriangulation%h_DvertexParameterValue, &
          rbackupTriangulation%h_DvertexParameterValue)

    ! Bit 13: DMBDP  
    CALL checkAndCopy(idupflag,TR_SHARE_DEDGEPARAMETERVALUE,&
          rtriangulation%h_DedgeParameterValue, &
          rbackupTriangulation%h_DedgeParameterValue)

    ! Bit 17: KVBDI  
    CALL checkAndCopy(idupflag,TR_SHARE_IBOUNDARYVERTEXPOS,&
          rtriangulation%h_IboundaryVertexPos, &
          rbackupTriangulation%h_IboundaryVertexPos)

    ! Bit 18: KMBDI  
    CALL checkAndCopy(idupflag,TR_SHARE_IBOUNDARYEDGEPOS,&
          rtriangulation%h_IboundaryEdgePos, &
          rbackupTriangulation%h_IboundaryEdgePos)
    
    ! Bit  1: DCORMG 
    CALL checkAndCopy(idupflag, TR_SHARE_DFREEVERTEXCOORDINATES,&
          rtriangulation%h_DfreeVertexCoordinates, &
          rbackupTriangulation%h_DfreeVertexCoordinates)
    
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
      
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_restore (rbackupTriangulation,rtriangulation)
  
!<description>
  ! This routine restores data of a triangulation structure. All
  ! information arrays not shared between rbackupTriangulation and another 
  ! triangulation structure is copied (back) into the rtriangulation.
!</description>

!<input>
  ! Backup of a triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rbackupTriangulation
!</input>

!<inputoutput>
  ! Destination triangulation.
  ! All arrays where a duplicates exist in rtriangulation are copied
  ! to rtriangulation, overwriting the old information arrays.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER(I32) :: idupFlag
  
    idupFlag = rtriangulation%iduplicationFlag
      
    ! Call checkAndCopy for all the arrays. this will either copy the handle
    ! or copy the content of the array.
    
    ! Bit  0: DCORVG 
    CALL checkAndCopy(idupflag, TR_SHARE_DVERTEXCOORDS,&
          rtriangulation%h_DvertexCoords, &
          rbackupTriangulation%h_DvertexCoords)

    ! Bit  2: KVERT  
    CALL checkAndCopy(idupflag, TR_SHARE_IVERTICESATELEMENT,&
          rtriangulation%h_IverticesAtElement, &
          rbackupTriangulation%h_IverticesAtElement)

    ! Bit  3: KMID   
    CALL checkAndCopy(idupflag, TR_SHARE_IEDGESATELEMENT,&
          rtriangulation%h_IedgesAtElement, &
          rbackupTriangulation%h_IedgesAtElement)
    
    ! Bit  4: KADJ   
    CALL checkAndCopy(idupflag, TR_SHARE_INEIGHBOURSATELEMENT,&
          rtriangulation%h_IneighboursAtElement, &
          rbackupTriangulation%h_IneighboursAtElement)

    ! Bit  6: KMEL   
    CALL checkAndCopy(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge, &
          rbackupTriangulation%h_IelementsAtEdge)

    ! Bit 16: KEAN   
    CALL checkAndCopy(idupflag,TR_SHARE_IVERTICESATEDGE,&
          rtriangulation%h_IverticesAtEdge, &
          rbackupTriangulation%h_IverticesAtEdge)

    ! Bit  7: KNPR   
    CALL checkAndCopy(idupflag, TR_SHARE_INODALPROPERTY,&
          rtriangulation%h_InodalProperty, &
          rbackupTriangulation%h_InodalProperty)

    ! Bit 19: DAREA  
    CALL checkAndCopy(idupflag,TR_SHARE_DELEMENTAREA,&
          rtriangulation%h_DelementVolume, &
          rbackupTriangulation%h_DelementVolume)

    ! Bit  5: KVEL   
    CALL checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertexIdx, &
          rbackupTriangulation%h_IelementsAtVertexIdx)
    CALL checkAndCopy(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertex, &
          rbackupTriangulation%h_IelementsAtVertex)

    ! Bit 11: KBCT   
    CALL checkAndCopy(idupflag,TR_SHARE_IBOUNDARYCPIDX,&
          rtriangulation%h_IboundaryCpIdx, &
          rbackupTriangulation%h_IboundaryCpIdx)

    ! Bit  9: KVBD   
    CALL checkAndCopy(idupflag, TR_SHARE_IVERTICESATBOUNDARY,&
          rtriangulation%h_IverticesAtBoundary, &
          rbackupTriangulation%h_IverticesAtBoundary)

    ! Bit 14: KMBD   
    CALL checkAndCopy(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IedgesAtBoundary, &
          rbackupTriangulation%h_IedgesAtBoundary)

    ! Bit 10: KEBD   
    CALL checkAndCopy(idupflag,TR_SHARE_IELEMENTSATBOUNDARY,&
          rtriangulation%h_IelementsAtBoundary, &
          rbackupTriangulation%h_IelementsAtBoundary)

    ! Bit 12: DVBDP  
    CALL checkAndCopy(idupflag,TR_SHARE_DVERTEXPARAMETERVALUE,&
          rtriangulation%h_DvertexParameterValue, &
          rbackupTriangulation%h_DvertexParameterValue)

    ! Bit 13: DMBDP  
    CALL checkAndCopy(idupflag,TR_SHARE_DEDGEPARAMETERVALUE,&
          rtriangulation%h_DedgeParameterValue, &
          rbackupTriangulation%h_DedgeParameterValue)

    ! Bit 17: KVBDI  
    CALL checkAndCopy(idupflag,TR_SHARE_IBOUNDARYVERTEXPOS,&
          rtriangulation%h_IboundaryVertexPos, &
          rbackupTriangulation%h_IboundaryVertexPos)

    ! Bit 18: KMBDI  
    CALL checkAndCopy(idupflag,TR_SHARE_IBOUNDARYEDGEPOS,&
          rtriangulation%h_IboundaryEdgePos, &
          rbackupTriangulation%h_IboundaryEdgePos)
    
    ! Bit  1: DCORMG 
    CALL checkAndCopy(idupflag, TR_SHARE_DFREEVERTEXCOORDINATES,&
          rtriangulation%h_DfreeVertexCoordinates, &
          rbackupTriangulation%h_DfreeVertexCoordinates)
    
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
      
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_done (rtriangulation)
  
!<description>
  ! This routine cleans up a triangulation structure.
  ! All memory allocated by handles in the structure is released from the heap.
  ! The routine does not clean up old FEAT 1.x legacy arrays in the
  ! rtriangulation%Itria substructure!
!</description>

!<inputoutput>
  ! The triangulation structure to be cleaned up.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    INTEGER(I32) :: idupflag
    
    IF (rtriangulation%ndim .EQ. 0) RETURN
    
    idupflag = rtriangulation%iduplicationFlag
    
    ! Bit  8: KMM    is a copy of another structure
    ! ... does not exist!?!

    ! Just release all allocated handles.
    ! Take care of which handles are duplicates from other structures - 
    ! these must not be released, as we are not the owner of them!
    
    ! Bit  0: DCORVG is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_DVERTEXCOORDS,&
          rtriangulation%h_DvertexCoords)

    ! Bit  2: KVERT  is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_IVERTICESATELEMENT,&
          rtriangulation%h_IverticesAtElement)

    ! Bit  3: KMID   is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_IEDGESATELEMENT,&
          rtriangulation%h_IedgesAtElement)
    
    ! Bit  4: KADJ   is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_INEIGHBOURSATELEMENT,&
          rtriangulation%h_IneighboursAtElement)

    ! Bit  6: KMEL   is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge)

    ! Bit 16: KEAN   is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_IVERTICESATEDGE,&
          rtriangulation%h_IverticesAtEdge)

    ! Bit  7: KNPR   is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_INODALPROPERTY,&
          rtriangulation%h_InodalProperty)

    ! Bit 19: DAREA  is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_DELEMENTAREA,&
          rtriangulation%h_DelementVolume)

    ! Bit  5: KVEL   is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertexIdx)
    CALL checkAndRelease(idupflag, TR_SHARE_IELEMENTSATVERTEX,&
          rtriangulation%h_IelementsAtVertex)

    ! Bit 11: KBCT   is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_IBOUNDARYCPIDX,&
          rtriangulation%h_IboundaryCpIdx)

    ! Bit  9: KVBD   is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_IVERTICESATBOUNDARY,&
          rtriangulation%h_IverticesAtBoundary)

    ! Bit 14: KMBD   is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_IEDGESATBOUNDARY,&
          rtriangulation%h_IedgesAtBoundary)

    ! Bit 10: KEBD   is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_IELEMENTSATBOUNDARY,&
          rtriangulation%h_IelementsAtBoundary)

    ! Bit 12: DVBDP  is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_DVERTEXPARAMETERVALUE,&
          rtriangulation%h_DvertexParameterValue)

    ! Bit 13: DMBDP  is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_DEDGEPARAMETERVALUE,&
          rtriangulation%h_DedgeParameterValue)

    ! Bit 17: KVBDI  is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_IBOUNDARYVERTEXPOS,&
          rtriangulation%h_IboundaryVertexPos)

    ! Bit 18: KMBDI  is a copy of another structure
    CALL checkAndRelease(idupflag,TR_SHARE_IBOUNDARYEDGEPOS,&
          rtriangulation%h_IboundaryEdgePos)
    
    ! Bit  1: DCORMG is a copy of another structure
    CALL checkAndRelease(idupflag, TR_SHARE_INODALPROPERTY,&
          rtriangulation%h_DfreeVertexCoordinates)
    
    ! Clean up the rest of the structure

    rtriangulation%iduplicationFlag = 0
    rtriangulation%ndim = 0
    rtriangulation%NVT = 0
    rtriangulation%NMT = 0
    rtriangulation%NEL = 0
    rtriangulation%NBCT = 0
    rtriangulation%NVBD = 0
    rtriangulation%NMBD = 0
    rtriangulation%nverticesPerEdge = 0
    rtriangulation%nVerticesOnAllEdges = 0
    rtriangulation%nverticesInEachElement = 0
    rtriangulation%nverticesInAllElements = 0
    rtriangulation%nadditionalVertices = 0
    rtriangulation%InelOfType(:) = 0

    ! That's it...

  CONTAINS
  
    ! **********************************************************
    ! Release handle ihandle if bitfield ibitfield in idubFlag is not set.
    ! Otherwise, ihandle is set to ST_NOHANDLE.    
    SUBROUTINE checkAndRelease (idupFlag,ibitfield,ihandle)
    
    INTEGER(I32), INTENT(IN) :: ibitfield
    INTEGER(I32), INTENT(IN) :: idupFlag
    INTEGER, INTENT(INOUT) :: ihandle
    
      IF (IAND(idupFlag,ibitfield) .NE. ibitfield) THEN
        IF (ihandle .NE. ST_NOHANDLE) CALL storage_free(ihandle)
      ELSE
        ihandle = ST_NOHANDLE
      END IF
      
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>    

  SUBROUTINE tria_rawGridToTri (rtriangulation)
  
!<description>
  ! This routine converts a 2D 'raw' mesh into a triangular 2D mesh. All elements are 
  ! converted to triangles.
  ! Warning: This routine can only applied to a 'raw' mesh, e.g. a mesh which 
  !  comes directly from a .TRI file. It's not advisable to apply this routine to
  !  a 'standard' mesh that contains already further mesh information 
  !  (adjacencies,...), as any information about possible two-level ordering
  !  is lost!
  !  However, when applied to a 'standard' mesh, the mesh information can be
  !  updated using tria_initStandardMeshFromRaw to form a valid standard mesh
  !  again.
!</description>

!<inputoutput>
  ! The triangulation structure to be converted. Is replaced by a triangular
  ! mesh.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>

    INTEGER(PREC_ELEMENTIDX) :: icount
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER :: h_IverticesAtElementTri
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementTri
    INTEGER(I32), DIMENSION(2) :: Isize
   
    ! For this routine we currently assume that there are only triangles
    ! and quads in the triangulation. Might be a matter of change in 
    ! the future...
   
    ! Get the points-at-element array
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
    icount = rtriangulation%InelOfType(TRIA_NVEQUAD2D)
    
    ! Create a new p_IverticesAtElement array for the triangular mesh.
    Isize = (/TRIA_NVETRI2D,icount+rtriangulation%NEL/)
    CALL storage_new2D ('tria_quadToTri', 'KVERTTRI', Isize, ST_INT, &
                        h_IverticesAtElementTri,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int2d (h_IverticesAtElementTri,p_IverticesAtElementTri)
    
    ! Convert the array
    CALL quadToTriang_aux1 (&
        rtriangulation%NEL, icount,&
        p_IverticesAtElement, p_IverticesAtElementTri)
    
    ! Replace the old array by the new, release the old one.
    CALL storage_free (rtriangulation%h_IverticesAtElement)
    rtriangulation%h_IverticesAtElement = h_IverticesAtElementTri
    
    ! Finally, set up NEL and InelOfType.
    ! Every quad got two triangles, so the number of elements increases by the number
    ! of quads!
    rtriangulation%NEL = rtriangulation%NEL + icount
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(TRIA_NVETRI2D) = rtriangulation%NEL
    
    ! That's it.

  CONTAINS

    SUBROUTINE quadToTriang_aux1 (nel, nelquad, IverticesAtElement, Kvert_triang)

  !<description>
    ! Purpose: Convert mixed quad/tri mesh to triangular mesh
    !
    ! This routine creates a triangular KVERT structure from a 
    ! quadrilateral KVERT structure.
  !</description>

  !<input>
    ! Number of elements in the old grid
    INTEGER(PREC_ELEMENTIDX),INTENT(IN) :: nel

    ! Number of quads in the old grid
    INTEGER,INTENT(IN) :: nelquad

    ! KVERT structure of the old mesh
    ! array [1..4,1..nel] of integer
    INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement
  !</input>

  !<output>
    ! KVERT structure of the tri mesh
    ! array [1..4,1..nel+nelquad] of integer
    INTEGER(I32), DIMENSION(:,:), INTENT(OUT)  :: Kvert_triang
  !</output>
  
      ! local variables
      INTEGER(PREC_ELEMENTIDX) :: i,j
        
      j = nel
      DO i=1,nel
        ! Copy the first three entries of each IverticesAtElement subarray
        ! to the first half of Kvert_triang. They form NEL quads.
        Kvert_triang(1,i) = IverticesAtElement(1,i)
        Kvert_triang(2,i) = IverticesAtElement(2,i)
        Kvert_triang(3,i) = IverticesAtElement(3,i)
        
        ! For every quad we find, we produce a second triangle with triangle
        ! number NEL+1,...
        IF (IverticesAtElement(4,i) .NE. 0) THEN
        
          ! Get the next free element number behind the first set of triangles
          j = j+1
          
          ! The second triangle in each quad consists of vertices 1,3,4.
          Kvert_triang(1,j) = IverticesAtElement(1,i)
          Kvert_triang(2,j) = IverticesAtElement(3,i)
          Kvert_triang(3,j) = IverticesAtElement(4,i)
        
        END IF
      END DO
      
      ! That's it.
      
    END SUBROUTINE  
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_readTriFile2D(rtriangulation, sfilename, rboundary)

!<description>
  ! This routine reads a .TRI file of a 2D triangulation into memory
  ! and creates a 'raw' triangulation (i.e. a triangulation that contains
  ! only basic information, see below). 
  !
  ! The triangulation structure rtriangulation is initialised with the data 
  ! from the file. The parameter sfilename gives the name of the .prm 
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
  ! The triangulation structure rtriangulation must be empty. All previous
  ! information in this structure (if there is any) is lost!
!</description>

!<input>
  ! The name of the .prm file to read.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename

  ! OPTIONAL: An rboundary object specifying the underlying domain.
  ! If not specified, the routine assumes that the TRI file does not specify
  ! boundary parameter values, i.e. the point coordinates in the TRI file
  ! are all real coordinates. The array DvertexParameterValue is not
  ! generated in this case.
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary
! </input>
  
!<output>
  ! Triangulation structure, to be filled with data
  TYPE(t_triangulation), INTENT(OUT) :: rtriangulation
!</output>
  
!</subroutine>

    ! Input channel for reading
    INTEGER :: iunit
    
    ! Open the file
    CALL io_openFileForReading(sfilename, iunit)

    ! We create a 2D triangulation here.
    rtriangulation%ndim = NDIM2D

    ! Read the basic mesh
    CALL tria_readRawTriangulation2D (iunit,rtriangulation)

    ! Create the basic boundary information
    CALL tria_genRawBoundary2D (rtriangulation,rboundary)

    ! Close the file, finish
    CLOSE(iunit)

  END SUBROUTINE
  
!************************************************************************

!<subroutine>

  SUBROUTINE tria_readRawTriangulation2D (iunit,rtriangulation)

!<description>  
  ! Auxiliary routine of tria_readTriFile2D.
  ! Reads basic information from a triangulation file into rtriangulation.
  ! That means, the following information arrays / tags are initialised:
  ! NEL,NVT,NMT,NBCT,InelOfType,
  ! DvertexCoords, IverticesAtElement and InodalProperty.
  ! The data is read from the file without being changed!
!</description>
  
!<input>
  ! Unit number of the file to be read
  INTEGER, INTENT(IN) :: iunit
!</input>
  
!<inputoutput> 
  ! Triangulation to be initialised with basic data.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>

    ! Local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata2D
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_Idata2D
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata
    INTEGER(I32) :: ivt, iel
    INTEGER :: idim, ive, nve
    INTEGER(I32), DIMENSION(2) :: Isize
    
    ! The first two lines in the file are comments.
    READ(iunit,*)
    READ(iunit,*)
    
    ! Read NEL,NVT,NMT,NVE,NBCT from the file
    ! and store this information in the structure.
    READ (iunit,*) rtriangulation%NEL,rtriangulation%NVT,rtriangulation%NMT,&
        NVE,rtriangulation%NBCT    

    ! Comment: 'DCORVG'
    READ (iunit,*)

    ! Allocate memory for the basic arrays on the heap
    Isize = (/NDIM2D,INT(rtriangulation%NVT,I32)/)
    CALL storage_new2D ('tria_read_tri2D', 'DCORVG', Isize, ST_DOUBLE, &
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    CALL storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
        
    ! Read the data from the file, store it in the array.
    READ (iunit,*) ((p_Ddata2D(idim,ivt),idim=1,NDIM2D), ivt=1,rtriangulation%NVT)

    ! Comment: 'KVERT'
    READ (iunit,*)

    ! Allocate memory for IverticesAtElement
    Isize = (/nve,INT(rtriangulation%NEL,I32)/)
    CALL storage_new2D ('tria_read_tri2D', 'KVERT', Isize, ST_INT, &
        rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    CALL storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)

    READ (iunit,*) ((p_Idata2D(ive,iel),ive=1,nve), iel=1,rtriangulation%NEL)

    ! Loop through the elements and determine how many elements
    ! of each element type we have.
    rtriangulation%InelOfType(:) = 0
    DO iel=1,rtriangulation%nel
      DO ive=nve,1,-1
        IF (p_Idata2D(ive,iel) .NE. 0) THEN
          rtriangulation%InelOfType(ive) = rtriangulation%InelOfType(ive)+1
          EXIT
        END IF
      END DO
    END DO

    ! Comment: 'KNPR'
    READ (iunit,*)

    ! Allocate memory for InodalProperty 
    CALL storage_new ('tria_read_tri2D', 'KNPR', &
        INT(rtriangulation%NVT,I32), ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    CALL storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_Idata)

    ! Read the data
    READ (iunit,*) (p_Idata(ivt),ivt=1,rtriangulation%NVT)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genRawBoundary2D (rtriangulation,rboundary)

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
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary
!</input>
  
!<inputoutput>
  ! Triangulation to be initialised with basic data.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_VERTEXIDX) :: ivbd,ivt
    INTEGER :: ibct
    INTEGER(I32), DIMENSION(:), POINTER :: p_InodalProperty

    ! Get the pointer to the InodalProperty array
    CALL storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_InodalProperty)

    ! Calculate NVBD by simply counting how many elements
    ! in p_InodalProperty are <> 0.
    rtriangulation%NVBD = 0
    DO ivt=1,rtriangulation%NVT
      IF (p_InodalProperty(ivt) .NE. 0) rtriangulation%NVBD = rtriangulation%NVBD+1
    END DO

    ! Allocate memory for IverticesAtBoundary.
    CALL storage_new ('tria_generateBasicBoundary', &
        'KVBD', INT(rtriangulation%NVBD,I32), &
        ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    CALL storage_new ('tria_generateBasicBoundary', &
        'KBCT', INT(rtriangulation%NBCT+1,I32), &
        ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    CALL storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    CALL storage_getbase_double2D (&
        rtriangulation%h_DvertexCoords,p_DvertexCoords)
        
    CALL storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! The first element in p_IboundaryCpIdx is (as the head) always =1.
    p_IboundaryCpIdx(1) = 1

    ! Perform a first loop over all vertices to check which are on the
    ! boundary. Count them. This way, we at first set up the boundary
    ! component index vector p_IboundaryCpIdx.
    ! In this first step, we save the number of vertices on each 
    ! boundary component in p_IboundaryCpIdx(2:NBCT+1)!
    DO ivt=1,rtriangulation%NVT
      IF (p_InodalProperty(ivt) .NE. 0) THEN
        ibct = p_InodalProperty(ivt)
        
        ! Increase the number of vertices in that boundary component by 1.
        ! The number of vertices on boundary component i is saved here
        ! at p_IboundaryCpIdx(i+1) for later.
        ! Note that the array was initialised with zero during the creation
        ! process!
        p_IboundaryCpIdx(ibct+1) = p_IboundaryCpIdx(ibct+1)+1
      END IF
    END DO
    
    ! Sum up the number of vertices on each boundary component to get the
    ! actual index vector.
    DO ibct = 2,rtriangulation%NBCT+1
      p_IboundaryCpIdx(ibct) = p_IboundaryCpIdx(ibct)+p_IboundaryCpIdx(ibct-1)
    END DO
    
    ! Shift the p_IboundaryCpIdx array by one position. That's a little trick in
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
    ! the first one, i.e. we'll again find 8,6 and 4 vertices on the boundary components
    ! 1,2 and 3. Increasing the p_IboundaryCpIdx(2:NBCT+1) for boundary components
    ! 1..NBCT the same way as done in the first loop will therefore again lead to
    !
    !         i            1   2   3   4
    ! p_IboundaryCpIdx(i)  1   9  15  19
    !
    ! Ok, let's catch the actual vertices.
    !
    ! The loop must be slightly modified if rboundary is not present!
    IF (PRESENT(rboundary)) THEN

      ! Allocate memory for  and DvertexParameterValue
      CALL storage_new ('tria_generateBasicBoundary', &
          'DVBDP', INT(rtriangulation%NVBD,I32), &
          ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
      
      ! Get the array where to store boundary parameter values.
      CALL storage_getbase_double (&
          rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
          
      ! Check all vertices to find out, which vertices are on the boundary.
      DO ivt=1,rtriangulation%NVT
        IF (p_InodalProperty(ivt) .NE. 0) THEN
          ibct = p_InodalProperty(ivt)
          
          ! Create a new point on that boundary component 
          ! and get the number, the point will have.
          ! Note that the array was initialised with zero during the creation
          ! process!
          ivbd = p_IboundaryCpIdx(ibct+1)
          p_IboundaryCpIdx(ibct+1) = ivbd+1
          
          ! Store the vertex as boundary vertex
          p_IverticesAtBoundary (ivbd) = ivt
          
          ! Store the parameter value; it's saved in DvertexCoords(1,.)
          p_DvertexParameterValue (ivbd) = p_DvertexCoords(1,ivt)
          
          ! Replace the coordinates in DvertexCoords by those
          ! given by the parametrisation.
          CALL boundary_getCoords(rboundary, ibct, p_DvertexParameterValue (ivbd), &
              p_DvertexCoords(1,ivt), p_DvertexCoords(2,ivt))
          
        END IF
      END DO
      
    ELSE
    
      ! No parametrisation available, the array with boundary parameter values 
      ! is not generaterd.
      !
      ! Check all vertices to find out, which vertices are on the boundary.
      DO ivt=1,rtriangulation%NVT
        IF (p_InodalProperty(ivt) .NE. 0) THEN
          ibct = p_InodalProperty(ivt)
          
          ! Create a new point on that boundary component 
          ! and get the number, the point will have.
          ! Note that the array was initialised with zero during the creation
          ! process!
          ivbd = p_IboundaryCpIdx(ibct+1)
          p_IboundaryCpIdx(ibct+1) = ivbd+1
          
          ! Store the vertex as boundary vertex
          p_IverticesAtBoundary (ivbd) = ivt
        END IF
      END DO
      
    END IF
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_initStandardMeshFromRaw(rtriangulation,rboundary)

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
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary
!</input>

!<inputoutput>
  ! Triangulation structure to be initialised.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    SELECT CASE (rtriangulation%ndim)
    CASE (NDIM1D)
    
    CASE (NDIM2D)
      ! Generate all standard arrays for 2D meshes.
      CALL tria_genElementsAtVertex2D    (rtriangulation)
      CALL tria_genNeighboursAtElement2D (rtriangulation)
      CALL tria_genEdgesAtElement2D      (rtriangulation)
      CALL tria_genElementsAtEdge2D      (rtriangulation)
      CALL tria_genVerticesAtEdge2D      (rtriangulation)
      CALL tria_genEdgeNodalProperty2D   (rtriangulation)
      CALL tria_genElementVolume2D       (rtriangulation)

      CALL tria_sortBoundaryVertices2D   (rtriangulation)
      CALL tria_genElementsAtBoundary2D  (rtriangulation)
      CALL tria_genEdgesAtBoundary2D     (rtriangulation)
      IF (PRESENT(rboundary)) THEN
        CALL tria_genEdgeParameterValue2D  (rtriangulation,rboundary)
      END IF
      CALL tria_genBoundaryVertexPos2D   (rtriangulation)
      CALL tria_genBoundaryEdgePos2D     (rtriangulation)
      
    CASE (NDIM3D)
    
    CASE DEFAULT
      CALL output_line ('Triangulation structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_generateStandardMeshFromRaw')
      CALL sys_halt()
    END SELECT

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_sortBoundaryVertices2D (rtriangulation)

!<description>  
  ! This routine sorts the arrays IverticesAtBoundary and 
  ! DvertexParameterValue such that they are ordered for increasing
  ! parameter values of the vertices.
!</description>
  
!<inputoutput>
  ! The triangulation.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_VERTEXIDX) :: istart,iend
    INTEGER :: ibct,ivbd
    INTEGER :: hresort
    INTEGER(I32), DIMENSION(:), POINTER :: p_Iresort

    IF (rtriangulation%h_DvertexParameterValue .EQ. ST_NOHANDLE) THEN
      ! We cannot sort the boundary vertices without parameter values!
      RETURN
    END IF

    ! Get pointers to the arrays
    CALL storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    CALL storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
        
    CALL storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Sort the p_DvertexCoords (sub)-arrays for increasing
    ! parameter value. 
    CALL storage_new ('tria_generateBasicBoundary', &
        'resort', INT(rtriangulation%NVBD,I32), &
        ST_INT, hresort, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int (hresort,p_Iresort)
    
    ! Fill p_Iresort with 1,2,3,...
    DO ivbd=1,rtriangulation%NVBD
      p_Iresort(ivbd) = ivbd
    END DO
    
    ! For each boundary component, call the sorting routine and calculate the
    ! mapping how the entries are sorted.
    DO ibct = 1,rtriangulation%NBCT
      istart = p_IboundaryCpIdx(ibct)
      iend = p_IboundaryCpIdx(ibct+1)-1
      
      ! Sort the sub-array of boundary component ibct.
      ! Remember, how the parameter values are resorted when sorting the array.
      CALL sort_dp(p_DvertexParameterValue(istart:iend),SORT_QUICK,&
          p_Iresort(istart:iend))
      
    END DO

    ! Resort the vertices according to the calculated mapping.
    ! Afterwards, the vertices are ordered according to the increasing
    ! parameter value.
    DO ivbd=1,SIZE(p_IverticesAtBoundary)
      p_Iresort(ivbd) = p_IverticesAtBoundary(p_Iresort(ivbd))
    END DO
    
    DO ivbd=1,SIZE(p_IverticesAtBoundary)
      p_IverticesAtBoundary(ivbd) = p_Iresort(ivbd)
    END DO
    
    CALL storage_free (hresort)

  END SUBROUTINE

!************************************************************************

!<function>

  INTEGER FUNCTION tria_getNNVE(rtriangulation)

!<description>
  ! Determines the maximum number of vertices per element in
  ! rtriangulation.
!</description>

!<inputoutput>
  ! The triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
!</inputoutput>
  
!<result>
  ! Maximum number of vertices per element.
!</result>
  
!</function>

    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_Idata2D
    
    ! NNVE is given by the first dimension of KVERT
    CALL storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)
        
    tria_getNNVE = UBOUND(p_Idata2D,1)

  END FUNCTION

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genElementsAtVertex2D(rtriangulation)

!<description>
  ! This routine generates the array IelementsAtVertex.
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER :: ive,nnve
    INTEGER(PREC_VERTEXIDX) :: ivt,isize
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    
    INTEGER(PREC_ELEMENTIDX) :: iel
    
    INTEGER :: haux1
    INTEGER(I32), DIMENSION(:), POINTER :: p_Iaux1

    ! Do we have (enough) memory for the index array?
    IF (rtriangulation%h_IelementsAtVertexIdx .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('tria_genElementsAtVertex2D', 'IelementsAtVertexIdx', &
          INT(rtriangulation%NVT+1,I32), ST_INT, &
          rtriangulation%h_IelementsAtVertexIdx, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize (rtriangulation%h_IelementsAtVertexIdx, isize)
      IF (isize .NE. rtriangulation%NVT+1) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genElementsAtVertex2D', &
            INT(rtriangulation%NVT+1,I32), rtriangulation%h_IelementsAtVertexIdx, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Get the index array.
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
        
    ! Get some data arrays about the vertices.
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)

    ! Fill the index array with zero.
    CALL storage_clear (rtriangulation%h_IelementsAtVertexIdx)
    
    nnve = tria_getNNVE(rtriangulation)

    ! We create the index array in two steps. In the first step,
    ! we loop over the elements to find out, how many elements
    ! meet at each vertex. We store this information in the index
    ! array at position 2..NVT+1 (i.e. shifted by one) for later.
    
    DO iel = 1,rtriangulation%NEL
      DO ive = 1,nnve
        ivt = p_IverticesAtElement(ive,iel)
        ! Cancel that element if we reached the end. Might happen if there
        ! are triangles in a quad mesh e.g.
        IF (ivt .EQ. 0) EXIT

        p_IelementsAtVertexIdx(ivt+1)=p_IelementsAtVertexIdx(ivt+1)+1
      END DO
    END DO
    
    ! Set the first element in p_IverticesAtElement to 1. Then, sum up
    ! all the length information to get the index array.
    p_IelementsAtVertexIdx(1) = 1
    DO ivt = 2,rtriangulation%nvt+1
      p_IelementsAtVertexIdx(ivt) = &
          p_IelementsAtVertexIdx(ivt) + p_IelementsAtVertexIdx(ivt-1)
    END DO
    
    isize = p_IelementsAtVertexIdx(rtriangulation%NVT+1)-1

    ! isize contains now the length of the array where we store the adjacency
    ! information (IelementsAtVertex).
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IelementsAtVertex .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('tria_genElementsAtVertex2D', 'IelementsAtVertex', &
          INT(isize,I32), ST_INT, &
          rtriangulation%h_IelementsAtVertex, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize (rtriangulation%h_IelementsAtVertex, isize)
      IF (isize .NE. isize) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genElementsAtVertex2D', &
            INT(isize,I32), rtriangulation%h_IelementsAtVertex, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Get the array.
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)

    ! Duplicate the p_IelementsAtVertexIdx array. We use that as pointer and index
    ! when new elements at a vertex are found.
    haux1 = ST_NOHANDLE
    CALL storage_copy (rtriangulation%h_IelementsAtVertexIdx,haux1)
    CALL storage_getbase_int (haux1,p_Iaux1)

    ! Now, we perform a second loop over the elements and the vertices.
    ! This time, we fetch the element number and write it to the    
    ! p_IelementsAtVertex array.
    ! p_Iaux1 counts how many positions in p_IelementsAtVertex are occupied.
    DO iel = 1,rtriangulation%NEL
      DO ive = 1,nnve
      
        ivt = p_IverticesAtElement(ive,iel)
        ! Cancel that element if we reached the end. Might happen if there
        ! are triangles in a quad mesh e.g.
        IF (ivt .EQ. 0) EXIT
        
        ! Remembner the element number and increase the pointer in the
        ! elements-adjacent-to-that-vertex list.
        p_IelementsAtVertex( p_Iaux1(ivt) ) = iel
        p_Iaux1(ivt) = p_Iaux1(ivt)+1
      END DO
    END DO
    
    ! Release the auxiliary array, that's it.
    CALL storage_free (haux1)

  END SUBROUTINE
    
!************************************************************************

!<subroutine>

  SUBROUTINE tria_genNeighboursAtElement2D(rtriangulation)

!<description>
  ! This routine generates the array IneighboursAtElement (KADJ). 
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IelementsAtVertexIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER :: ive,nve,nnve
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_VERTEXIDX) :: ivt,ivtneighbour
    
    INTEGER :: haux1, haux2
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgeAtVertex
    INTEGER(PREC_EDGEIDX) :: iidxEdge, iedge, iedgeneighbour
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IedgeIdx

    nnve = tria_getNNVE(rtriangulation)

    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      Isize = (/nnve,INT(rtriangulation%NEL,I32)/)
      CALL storage_new2D ('tria_genNeighboursAtElement2D', 'KADJ', &
          Isize, ST_INT, &
          rtriangulation%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize2D (rtriangulation%h_IneighboursAtElement, Isize)
      IF (Isize(2) .NE. rtriangulation%NEL) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genNeighboursAtElement2D', &
            rtriangulation%NEL, rtriangulation%h_IneighboursAtElement, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Fill the array with 0. We overwrite only those positions <> 0.
    CALL storage_clear (rtriangulation%h_IneighboursAtElement)
    
    ! Get the array which is to be filled with data.
    CALL storage_getbase_int2d (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
        
    ! Get some data arrays about the vertices.
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
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
    ! array are occupied. Each index starts from 'one behind possible
    ! free memory locations' and is decreased for every edge found.
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
    CALL storage_copy (rtriangulation%h_IelementsAtVertexIdx,haux2)
    CALL storage_getbase_int (haux2,p_IedgeIdx)
    
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
    
    ! Actually, we need the entries 2..NVT+1 of Iaux2
    p_IedgeIdx => p_IedgeIdx(2:)

    ! Create the auxiliary array for edge information.
    Isize = (/4,p_IedgeIdx(SIZE(p_IedgeIdx))-1/)
    CALL storage_new2D ('tria_genNeighboursAtElement2D', 'edgeAtVertex', &
        Isize, ST_INT, haux1, ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int2d (haux1,p_IedgeAtVertex)
    
    ! Loop through the elements to calculate edge information.
    DO iel=1,rtriangulation%NEL
     
      ! Get the number of vertices of that element.
      nve = nnve
      DO WHILE (p_IverticesAtElement(nve,iel) .EQ. 0)
        nve = nve-1
      END DO
    
      ! Loop through the vertices
      DO ive=1,nve
        ! Number of current vertex?
        ivt = p_IverticesAtElement(ive,iel)
        
        ! What's the neighbour?
        ivtneighbour = p_IverticesAtElement(MOD(ive,nve)+1,iel)
        
        ! Get the next free entry in the array with the edge information
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
      
      END DO ! ive
      
    END DO ! iel
    
    ! Next, we have to look at all edges to find neighbouring information.
    !
    ! Loop through all edges adjacent to all vertices.
    DO iedge = 1,UBOUND(p_IedgeAtVertex,2)
    
      ! Get the edge information
      ivt          = p_IedgeAtVertex(1,iedge) 
      ivtneighbour = p_IedgeAtVertex(2,iedge) 
      iel          = p_IedgeAtVertex(3,iedge) 
    
      ! Now, loop through all edges adjacent to the neighbour vertex
      ! ivtneighbour to find possible adjacent elements at the current
      ! edge.
      ! p_IedgeIdx(ivtneighbour) is the number of the first entry
      ! in p_IedgeAtVertex that's occupied!
      DO iedgeneighbour = p_IedgeIdx(ivtneighbour), &
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
        
        IF (p_IedgeAtVertex(2,iedgeneighbour) .EQ. ivt) THEN
        
          ! Save the adjacency information to our current element.
          ! We don't save it for the neighbour element (although we have it here)
          ! since we arrive at the neighbour element, too -- which would mean
          ! to store every information twice...
          p_IneighboursAtElement(p_IedgeAtVertex(4,iedge),iel) = &
            p_IedgeAtVertex(3,iedgeneighbour)
        
        END IF
    
      END DO ! iedgeneighbour
    
    END DO ! iedge
    
    ! Release memory, finish.
    CALL storage_free (haux1)
    CALL storage_free (haux2)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genElementsAtEdge2D(rtriangulation)

!<description>
  ! This routine generates information about the elements adjacent to each 
  ! edge IelementsAtEdge (KMID). 
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IneighboursAtElement.
  ! NMT must be set up correctly.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IelementsAtEdge
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_VERTEXIDX) :: NVT
    INTEGER(PREC_EDGEIDX) :: iedge
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%NMT .EQ. 0) THEN
      CALL output_line ('Edge information (NMT) not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IelementsAtEdge .EQ. ST_NOHANDLE) THEN
      Isize = (/2,rtriangulation%NMT/)
      CALL storage_new2D ('tria_genElementsAtEdge2D', 'KMID', &
          Isize, ST_INT, &
          rtriangulation%h_IelementsAtEdge, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize2D (rtriangulation%h_IelementsAtEdge, Isize)
      IF (Isize(2) .NE. rtriangulation%NEL) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genElementsAtEdge2D', &
            rtriangulation%NEL, rtriangulation%h_IelementsAtEdge, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int2D (rtriangulation%h_IelementsAtEdge,p_IelementsAtEdge)
    
    NVT = rtriangulation%NVT
    
    ! Loop through all elements and all edges on the elements
    DO iel = 1,UBOUND(p_IedgesAtElement,2)
      
      DO ive = 1,UBOUND(p_IedgesAtElement,1)

        ! Is there a neighbour?
        IF (p_IneighboursAtElement(ive,iel) .EQ. 0) THEN
          ! No. Store the element as only adjacent one to that edge.
          !
          ! Don't do anything if we are looking at a nonexisting edge here!
          ! (i.e. edge 4 of a tri element in a quad mesh -- tri elements
          !  don't have 4 edges :-) )
          IF (p_IedgesAtElement(ive,iel) .NE. 0) THEN

            iedge = p_IedgesAtElement (ive,iel)-NVT
            p_IelementsAtEdge(1,iedge) = iel
            p_IelementsAtEdge(2,iedge) = 0
            
          END IF
          
        ELSE IF (p_IneighboursAtElement(ive,iel) .LT. iel) THEN
        
          ! There is a neighbour and it has a smaller number than the current element --
          ! so we haven't had that edge! Store the two adjacent elements.
        
          iedge = p_IedgesAtElement (ive,iel)-NVT
          p_IelementsAtEdge(1,iedge) = iel
          p_IelementsAtEdge(2,iedge) = p_IneighboursAtElement(ive,iel)
        
        END IF
      
      END DO
      
    END DO
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genVerticesAtEdge2D(rtriangulation)

!<description>
  ! This routine generates information about the vertices adjacent to 
  ! each edge IverticesAtEdge.
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IneighboursAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive,nnve
    INTEGER(PREC_VERTEXIDX) :: ivtneighbour
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_EDGEIDX) :: iedge
    INTEGER(PREC_VERTEXIDX) :: NVT
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IverticesAtEdge .EQ. ST_NOHANDLE) THEN
      Isize = (/2,rtriangulation%NMT/)
      CALL storage_new2D ('tria_genVerticesAtEdge2D', 'KMID', &
          Isize, ST_INT, &
          rtriangulation%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize2D (rtriangulation%h_IverticesAtEdge, Isize)
      IF (Isize(2) .NE. rtriangulation%NEL) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genVerticesAtEdge2D', &
            rtriangulation%NEL, rtriangulation%h_IverticesAtEdge, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    
    nnve = UBOUND(p_IedgesAtElement,1)
    NVT = rtriangulation%NVT
    
    ! Loop through all elements and all edges on the elements
    DO iel = 1,UBOUND(p_IedgesAtElement,2)
      
      DO ive = 1,UBOUND(p_IedgesAtElement,1)
      
        ! Stop if we handled all edges; this is important if there are triangles
        ! in a quad mesh e.g.
        IF (p_IverticesAtElement(ive,iel) .EQ. 0) EXIT
        
        ! Is there a neighbour which number is less than iel? Or even =0?
        ! If yes, we didn't tackle the edge.
          
        IF (p_IneighboursAtElement(ive,iel) .LT. iel) THEN

          ! Save the numbers of the adjacent vertices.
          iedge = p_IedgesAtElement (ive,iel)-NVT
          
          p_IverticesAtEdge(1,iedge) = p_IverticesAtElement(ive,iel)
          
          ! Also save the neighbour. Note that we have to check the number of the
          ! neighbour against 0 because it may be that the p_IverticesAtElement(:,.)
          ! array is not completely filled -- e.g. if there are triangles in a quad
          ! mesh!
          ivtneighbour = p_IverticesAtElement(MOD(ive,nnve)+1,iel)
          IF (ivtneighbour .EQ. 0) ivtneighbour = p_IverticesAtElement(1,iel)
          p_IverticesAtEdge(2,iedge) = ivtneighbour
        
        END IF
      
      END DO
      
    END DO
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genEdgesAtElement2D(rtriangulation)

!<description>
  ! This routine calculates the edge numbering. That means,
  ! it generates information about the edges adjacent to
  ! each element IedgesAtElement (KMID) and calculates the correct NMT.
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IneighboursAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive, iveneighbour
    INTEGER(PREC_ELEMENTIDX) :: iel, ielneighbour
    INTEGER(PREC_EDGEIDX) :: iedge
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL storage_getsize2D (rtriangulation%h_IneighboursAtElement, Isize)
      CALL storage_new2D ('tria_genEdgesAtElement2D', 'KMID', &
          Isize, ST_INT, &
          rtriangulation%h_IedgesAtElement, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize2D (rtriangulation%h_IedgesAtElement, Isize)
      IF (Isize(2) .NE. rtriangulation%NEL) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genEdgesAtElement2D', &
            rtriangulation%NEL, rtriangulation%h_IedgesAtElement, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Fill IedgesAtElement with 0. That's important in case some
    ! elements in the array are not tackled when searching for edges
    ! (e.g. in meshes where triangles and quads are mixed).
    CALL storage_clear (rtriangulation%h_IedgesAtElement)

    CALL storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    
    ! iedge counts the edges and specifies the last given edge number.
    ! The numbering starts with NVT, the first edge gets the number NVT+1.
    iedge = rtriangulation%NVT
    
    ! Loop through all elements
    DO iel = 1,Isize(2)
    
      ! Loop through all edges on each element
      DO ive = 1,Isize(1)
      
        ! Stop if we handled all edges; this is important if there are triangles
        ! in a quad mesh e.g.
        IF (p_IverticesAtElement(ive,iel) .EQ. 0) EXIT
        
        ! Check the neighbour element.
        ! If the neightbour element has number =0 (no neighbour) or a number
        ! greater than IEL, we found the edge the first time and give it a number.
        IF ((p_IneighboursAtElement(ive,iel) .EQ. 0) .OR. &
            (p_IneighboursAtElement(ive,iel) .GT. iel)) THEN
        
          iedge = iedge + 1
          
          ! Add NVT to iedge to get the edge number
          p_IedgesAtElement(ive,iel) = iedge
        
        ELSE
        
          ! Otherweise, we had that edge already. Look into the neighbour element
          ! (which definitely exists, p_IneighboursAtElement cannot be =0 there!)
          ! to find the edge number.
          ielneighbour = p_IneighboursAtElement(ive,iel)
          
          DO iveneighbour = 1,Isize(1)
            IF (p_IneighboursAtElement(iveneighbour,ielneighbour) .EQ. iel) THEN
              p_IedgesAtElement(ive,iel) = &
                  p_IedgesAtElement(iveneighbour,ielneighbour)
              EXIT
            END IF
          END DO
        
        END IF
        
      END DO
    
    END DO
    
    ! Save the correct NMT.
    rtriangulation%NMT = iedge-rtriangulation%NVT
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genEdgeNodalProperty2D(rtriangulation)

!<description>
  ! This routine generates the nodal property tags for all edges 
  ! InodalProperty(NVT+1:NVT+NMT) (KNPR). 
  ! For this purpose, the following arrays are used:
  ! InodalProperty(1:NVT), IedgesAtElement, IneighboursAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(I32), DIMENSION(:), POINTER :: p_InodalProperty
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_VERTEXIDX) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_InodalProperty .EQ. ST_NOHANDLE) THEN
      CALL output_line ('InodalPropertys not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      CALL sys_halt()
    END IF
    
    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      CALL sys_halt()
    END IF
    
    ! Do we have (enough) memory for that array?
    CALL storage_getsize (rtriangulation%h_InodalProperty, isize)
    IF (Isize .NE. rtriangulation%NVT+rtriangulation%NMT) THEN
      ! If the size is wrong, reallocate memory.
      ! Copy the old content as we mustn't destroy the old nodal property
      ! tags of the vertices.
      CALL storage_realloc ('tria_genEdgeNodalProperty2D', &
          rtriangulation%NVT+rtriangulation%NMT, &
          rtriangulation%h_InodalProperty, &
          ST_NEWBLOCK_NOINIT, .TRUE.)
    END IF
    
    ! Get the arrays.
    CALL storage_getbase_int (rtriangulation%h_InodalProperty,p_InodalProperty)
    CALL storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IverticesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
    
    ! Loop through all elements and all edges on the elements
    DO iel = 1,UBOUND(p_IedgesAtElement,2)
    
      DO ive = 1,UBOUND(p_IedgesAtElement,1)
      
        ! Stop if we handled all edges; this is important if there are triangles
        ! in a quad mesh e.g.
        IF (p_IedgesAtElement(ive,iel) .EQ. 0) EXIT
        
        ! Is there a neighbour? If yes, we have an inner edge. If not, this is 
        ! a boundary edge.
        ! Note that by checking "neighbour > iel", we simultaneously check two things:
        ! 1.) Is there a neighbour at all? (neighbour <> 0)
        ! 2.) Has the neighbour a greater number? If not, we treated that
        !     edge already before and don't have to tackle it again.
        IF (p_IneighboursAtElement(ive,iel) .GT. iel) THEN
        
          ! Get the number of the boundary component from the vertex preceeding
          ! the edge and store it as information for the edge.
          p_InodalProperty(p_IedgesAtElement(ive,iel)) = &
              p_InodalProperty(p_IverticesAtElement(ive,iel))
        
        END IF
      
      END DO
      
    END DO

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genElementVolume2D(rtriangulation)

!<description>
  ! This routine generates the element volume array DelementVolume (DAREA). 
  ! For this purpose, the following arrays are used:
  ! DvertexCoordinates, IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
    REAL(DP), DIMENSION(:), POINTER :: p_DelementVolume
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_ELEMENTIDX) :: isize
    REAL(DP) :: dtotalVolume
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE2D) :: Dpoints
    INTEGER :: ive

    ! Is everything here we need?
    IF (rtriangulation%h_DvertexCoords .EQ. ST_NOHANDLE) THEN
      CALL output_line ('h_DvertexCoords not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement  not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      CALL sys_halt()
    END IF
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_DelementVolume .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('tria_genElementVolume2D', 'DAREA', &
          INT(rtriangulation%NEL+1,I32), ST_DOUBLE, &
          rtriangulation%h_DelementVolume, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize (rtriangulation%h_DelementVolume, isize)
      IF (isize .NE. rtriangulation%NEL+1) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genElementVolume2D', &
            INT(rtriangulation%NEL+1,I32), rtriangulation%h_DelementVolume, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    ! Get the arrays
    CALL storage_getbase_double2D (rtriangulation%h_DvertexCoords,&
        p_DvertexCoords)
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    CALL storage_getbase_double (rtriangulation%h_DelementVolume,&
        p_DelementVolume)
        
    dtotalVolume = 0.0_DP
        
    ! Currently, we support triangules and quads.
    IF (UBOUND(p_IverticesAtElement,1) .EQ. TRIA_NVETRI2D) THEN

      ! Calculate the element volume for all elements
      DO iel=1,rtriangulation%NEL
        ! triangular element
        DO ive=1,TRIA_NVETRI2D
          Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
          Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        END DO
        p_DelementVolume(iel) = gaux_getArea_tria2D(Dpoints)
        
        dtotalVolume = dtotalVolume+p_DelementVolume(iel)
      END DO
    
    ELSE

      ! Calculate the element volume for all elements
      DO iel=1,rtriangulation%NEL
      
        IF (p_IverticesAtElement(4,iel) .EQ. 0) THEN
          ! triangular element
          DO ive=1,TRIA_NVETRI2D
            Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
            Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
          END DO
          p_DelementVolume(iel) = gaux_getArea_tria2D(Dpoints)
        ELSE
          ! quad element
          DO ive=1,TRIA_NVEQUAD2D
            Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
            Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
          END DO
          p_DelementVolume(iel) = gaux_getArea_quad2D(Dpoints)
        END IF
        
        dtotalVolume = dtotalVolume+p_DelementVolume(iel)
      END DO
      
    END IF
    
    ! Store the total volume in the last element of DelementVolume
    p_DelementVolume(rtriangulation%NEL+1) = dtotalVolume
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genElementsAtBoundary2D(rtriangulation)

!<description>
  ! This routine generates information about the elements at the boundary
  ! IelementsAtBoundary (KEBD). 
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IneighboursAtElement, IverticesAtBoundary,
  ! IboundaryCpIdx, IelementsAtVertexIdx, IelementsAtVertex.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtBoundary
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER :: ive, ibct, ivbd, iadjElement, nnve
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_VERTEXIDX) :: isize, ivt

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IelementsAtVertexIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IelementsAtVertexIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IelementsAtVertex .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IelementsAtVertex not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
    CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,p_IelementsAtVertexIdx)
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IelementsAtBoundary .EQ. ST_NOHANDLE) THEN
      ! We have as many elements on the boudary as vertices!
      CALL storage_new ('tria_genElementsAtBoundary2D', 'KEBD', &
          rtriangulation%NVBD, ST_INT, &
          rtriangulation%h_IelementsAtBoundary, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize (rtriangulation%h_IelementsAtBoundary, isize)
      IF (isize .NE. rtriangulation%NVBD) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genElementsAtBoundary2D', &
            rtriangulation%NVBD, rtriangulation%h_IelementsAtBoundary, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int (rtriangulation%h_IelementsAtBoundary,p_IelementsAtBoundary)

    nnve = tria_getNNVE(rtriangulation)

    ! Loop through all boundary components
    DO ibct = 1,rtriangulation%NBCT
    
      ! On each boundary component, loop through all vertices
      DO ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
      
        !   +---+---+---+
        !   |   |   |   |
        !   +---+---+---+
        !   |   | 1 | 2 |
        !   +---+---X---+
        !          ivt
      
        ! Get the current boundary vertex
        ivt = p_IverticesAtBoundary(ivbd)
      
        ! Loop through all elements adjacent to the current vertex
        DO iadjElement = p_IelementsAtVertexIdx(ivt),p_IelementsAtVertexIdx(ivt+1)-1
        
          ! Get the element number
          iel = p_IelementsAtVertex(iadjElement)
        
          ! Find the local number of the vertex in the element
          DO ive=1,nnve
            IF (p_IverticesAtElement (ive,iel) .EQ. ivt) EXIT
          END DO
          
          ! Test if the element has a neighbour at the edge that's starting
          ! with our current vertex
          IF (p_IneighboursAtElement(ive,iel) .EQ. 0) THEN
          
            ! Yes, that's the boundary edge we are searching for!
            ! It starts with ivt and is present in the boundary component
            ! we are currently processing.
            ! So we can save the boundary element number and proceed with the next 
            ! boundary vertex.
            
            p_IelementsAtBoundary (ivbd) = p_IelementsAtVertex(iadjElement)
            
            EXIT
           
          END IF
        
        END DO
      
      END DO
    
    END DO
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genEdgesAtBoundary2D(rtriangulation)

!<description>
  ! This routine generates information about the edges at the boundary
  ! IedgesAtBoundary (KMBD). Initialises NMBD.
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IelementsAtBoundary, IverticesAtBoundary, 
  ! IedgesAtElement, IboundaryCpIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtBoundary
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IedgesAtBoundary
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER :: ive, ibct, ivbd, nnve
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_VERTEXIDX) :: isize, ivt

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IelementsAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IelementsAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IelementsAtBoundary,p_IelementsAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IedgesAtBoundary .EQ. ST_NOHANDLE) THEN
      ! We have as many elements on the boudary as vertices!
      CALL storage_new ('tria_genEdgesAtBoundary2D', 'KMBD', &
          rtriangulation%NVBD, ST_INT, &
          rtriangulation%h_IedgesAtBoundary, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize (rtriangulation%h_IelementsAtBoundary, isize)
      IF (isize .NE. rtriangulation%NVBD) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genEdgesAtBoundary2D', &
            rtriangulation%NVBD, rtriangulation%h_IedgesAtBoundary, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int (rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)

    nnve = tria_getNNVE(rtriangulation)

    ! Loop through all boundary components
    DO ibct = 1,rtriangulation%NBCT
    
      ! On each boundary component, loop through all elements
      DO ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-1
      
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
        DO ive=1,nnve
          IF (p_IverticesAtElement (ive,iel) .EQ. ivt) EXIT
        END DO
        
        ! Save the edge following the vertex on that element.
        p_IedgesAtBoundary(ivbd) = p_IedgesAtElement(ive,iel)
      
      END DO
    
    END DO
    
    ! We have as many edges on the boundary as vertices.
    rtriangulation%NMBD = rtriangulation%NVBD
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genBoundaryVertexPos2D(rtriangulation)

!<description>
  ! This routine generates the array IboundaryVertexPos which is used
  ! to search for the position of a boundary vertex.
  ! For this purpose, the following arrays are used:
  ! IverticesAtBoundary, IboundaryCpIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IboundaryVertexPos
    INTEGER :: ivbd, ibct
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genBoundaryVertexPos2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IboundaryVertexPos .EQ. ST_NOHANDLE) THEN
      ! We have as many elements on the boudary as vertices!
      Isize = (/2_I32,INT(rtriangulation%NVBD,I32)/)
      CALL storage_new2D ('tria_genBoundaryVertexPos2D', 'KVBDI', &
          Isize, ST_INT, &
          rtriangulation%h_IboundaryVertexPos, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize2D (rtriangulation%h_IboundaryVertexPos, Isize)
      IF (Isize(2) .NE. rtriangulation%NVBD) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genBoundaryVertexPos2D', &
            rtriangulation%NVBD, rtriangulation%h_IboundaryVertexPos, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int2D (rtriangulation%h_IboundaryVertexPos,p_IboundaryVertexPos)

    ! Fill the array.
    ! Store the IverticesAtBoundary in the first entry and the index in the 2nd.
    DO ivbd = 1,rtriangulation%NVBD
      p_IboundaryVertexPos(1,ivbd) = p_IverticesAtBoundary(ivbd)
      p_IboundaryVertexPos(2,ivbd) = ivbd
    END DO
    
    ! Sort the array -- inside of each boundary component.
    ! Use the vertex number as key.
    DO ibct = 1,rtriangulation%NBCT
      CALL arraySort_sortByIndex_int (&
          p_IboundaryVertexPos(:,p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),1)
    END DO
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genBoundaryEdgePos2D(rtriangulation)

!<description>
  ! This routine generates the array IboundaryEdgePos which is used
  ! to search for the position of a boundary edges.
  ! For this purpose, the following arrays are used:
  ! IedgesAtBoundary, IboundaryCpIdx.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IedgesAtBoundary
    INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IboundaryEdgePos
    INTEGER :: ivbd, ibct
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IedgesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genBoundaryEdgePos2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genBoundaryEdgePos2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int (rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IboundaryEdgePos .EQ. ST_NOHANDLE) THEN
      ! We have as many elements on the boudary as vertices!
      Isize = (/2_I32,INT(rtriangulation%NMBD,I32)/)
      CALL storage_new2D ('tria_genBoundaryVertexPos2D', 'KEBDI', &
          Isize, ST_INT, &
          rtriangulation%h_IboundaryEdgePos, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize2D (rtriangulation%h_IboundaryEdgePos, Isize)
      IF (Isize(2) .NE. rtriangulation%NMBD) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genBoundaryEdgePos2D', &
            rtriangulation%NMBD, rtriangulation%h_IboundaryEdgePos, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int2D (rtriangulation%h_IboundaryEdgePos,p_IboundaryEdgePos)

    ! Fill the array.
    ! Store the IverticesAtBoundary in the first entry and the index in the 2nd.
    DO ivbd = 1,rtriangulation%NMBD
      p_IboundaryEdgePos(1,ivbd) = p_IedgesAtBoundary(ivbd)
      p_IboundaryEdgePos(2,ivbd) = ivbd
    END DO
    
    ! Sort the array -- inside of each boundary component.
    ! Use the vertex number as key.
    DO ibct = 1,rtriangulation%NBCT
      CALL arraySort_sortByIndex_int (&
          p_IboundaryEdgePos(:,p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),1)
    END DO
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genEdgeParameterValue2D(rtriangulation,rboundary)

!<description>
  ! This routine generates the parameter value of edge midpoints
  ! on the boundary DedgeParameterValue (DMBDP).
  ! For this purpose, the following arrays are used:
  ! DvertexParameterValue, IverticesAtBoundary, IboundaryCpIdx.
  ! If necessary, new memory is allocated.
!</description>

!<input>
  ! Boundary structure that defines the parametrisation of the boundary.
  TYPE(t_boundary), INTENT(IN) :: rboundary
!</input>

!<inputoutput>
  ! The triangulation structure to be updated.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    REAL(DP), DIMENSION(:), POINTER :: p_DedgeParameterValue
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    
    INTEGER :: ibct, ivbd, hvertAtBd
    INTEGER(PREC_VERTEXIDX) :: isize
    REAL(DP) :: dpar1,dpar2,dmaxPar

    ! Is everything here we need?
    IF (rtriangulation%h_DvertexParameterValue .EQ. ST_NOHANDLE) THEN
      CALL output_line ('DvertexParameterValue not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      CALL sys_halt()
    END IF

    ! Get the arrays.
    CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
        
    ! Allocate an auxiliary array containing a copy of the parameter values 
    ! of the vertices.
    hvertAtBd = ST_NOHANDLE
    CALL storage_copy (rtriangulation%h_DvertexParameterValue,hvertAtBd)
    CALL storage_getbase_double (hvertAtBd,p_DvertexParameterValue)
    
    ! Convert the parameter values of the vertices from 0-1 into length
    ! parametrisation. We need this later to get the correct parameter
    ! values of the edge mitpoints.
    DO ibct=1,rtriangulation%NBCT
      CALL boundary_convertParameterList (rboundary,ibct,&
        p_DvertexParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
        p_DvertexParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
        BDR_PAR_01,BDR_PAR_LENGTH)
    END DO
    
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_DedgeParameterValue .EQ. ST_NOHANDLE) THEN
      ! We have as many elements on the boudary as vertices!
      CALL storage_new ('tria_genEdgeParameterValue2D', 'KMID', &
          rtriangulation%NVBD, ST_DOUBLE, &
          rtriangulation%h_DedgeParameterValue, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize (rtriangulation%h_DedgeParameterValue, isize)
      IF (isize .NE. rtriangulation%NVBD) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genEdgeParameterValue2D', &
            rtriangulation%NVBD, rtriangulation%h_DedgeParameterValue, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_double (rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)

    ! Loop through all boundary components
    DO ibct = 1,rtriangulation%NBCT
    
      ! Get the maximum parameter value of that BC.
      dmaxPar = boundary_dgetMaxParVal(rboundary, ibct, BDR_PAR_LENGTH)
    
      ! On each boundary component, loop through all vertices
      DO ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-2
      
        ! Get the parameter value of the current vertex and its neighbour.
        dpar1 = p_DvertexParameterValue(ivbd)
        dpar2 = p_DvertexParameterValue(ivbd+1)
        
        ! The edge parameter value is the mean.
        p_DedgeParameterValue(ivbd) = 0.5_DP*(dpar1+dpar2)
      
      END DO
      
      ! The 'last' vertex is a special case as there is no real neighbour.
      ! The parameter value of the 'neighbour' is the maximum parameter
      ! value of the boundary component + the parameter value of the
      ! first vertex of the boundary component.
      !
      ! Note that ivbd points to p_IboundaryCpIdx(ibct+1)-1 by
      ! Fortran standard of DO-loops!
      dpar1 = p_DvertexParameterValue(ivbd)
      dpar2 = p_DvertexParameterValue(p_IboundaryCpIdx(ibct)) + dmaxPar
      p_DedgeParameterValue(ivbd) = 0.5_DP*(dpar1+dpar2)

      ! Convert the parameter values of the edge midpoints back from
      ! length parametrisation to 0-1 parametrisation. 
      ! This automatically 'rounds down' parameter values that are > dmaxPar!   
      CALL boundary_convertParameterList (rboundary,ibct,&
        p_DedgeParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
        p_DedgeParameterValue(p_IboundaryCpIdx(ibct):p_IboundaryCpIdx(ibct+1)-1),&
        BDR_PAR_LENGTH,BDR_PAR_01)
        
    END DO

    ! Release temporary array, finish.
    CALL storage_free(hvertAtBd)
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_quickRefine2LevelOrdering(nfine,rtriangulation,rboundary)

!<description>
  ! This routine refines the given mesh rsourceTriangulation according to
  ! the 2-level ordering algorithm nfine times. The refined mesh is saved in 
  ! rdestTriangulation.
  !
  ! rtriangulation can be either a 'raw' mesh (as read from a .TRI
  ! file) or a 'standard' mesh. The mesh is overwritten by the refined
  ! one. The resulting mesh will be a 'raw' mesh), i.e. the caller must 
  ! add further information to it with routines like 
  ! tria_initStandardMeshFromRaw!
!</description>

!<input>
  ! Number of refinements that should be applied to rtriangulation. > 0.
  ! For a value <= 0, nothing will happen.
  INTEGER, INTENT(IN) :: nfine

  ! OPTIONAL: Boundary structure that defines the parametrisation of the boundary.
  ! If specified, the coordinates of the new boundary vertices are
  ! recomputed according to the analytic boundary.
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary
!</input>

!<inputoutput>
  ! The triangulation to be refined; can be a 'raw' or a 'standard' mesh.
  ! Is overwritten by the refined mesh.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>
 
    INTEGER :: ifine

    ! Refine nfine times:
    DO ifine = 1,nfine
    
      ! Create missing arrays in the source mesh.
      ! Create only those arrays we need to make the refinement process
      ! as fast as possible.
      IF (rtriangulation%h_IelementsAtVertex .EQ. ST_NOHANDLE) &
        CALL tria_genElementsAtVertex2D    (rtriangulation)
      IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) &
        CALL tria_genNeighboursAtElement2D (rtriangulation)
      IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) &
        CALL tria_genEdgesAtElement2D      (rtriangulation)
      IF (rtriangulation%h_IverticesAtEdge .EQ. ST_NOHANDLE) &
        CALL tria_genVerticesAtEdge2D      (rtriangulation)
      IF (rtriangulation%h_InodalProperty .EQ. ST_NOHANDLE) &
        CALL tria_genEdgeNodalProperty2D   (rtriangulation)

      CALL tria_sortBoundaryVertices2D   (rtriangulation)
      IF (rtriangulation%h_IelementsAtBoundary .EQ. ST_NOHANDLE) &
        CALL tria_genElementsAtBoundary2D  (rtriangulation)
      IF (rtriangulation%h_IedgesAtBoundary .EQ. ST_NOHANDLE) &
        CALL tria_genEdgesAtBoundary2D     (rtriangulation)
      IF (PRESENT(rboundary) .AND. &
          (rtriangulation%h_DedgeParameterValue .EQ. ST_NOHANDLE)) THEN
        CALL tria_genEdgeParameterValue2D  (rtriangulation,rboundary)
      END IF
     
      ! Refine the mesh, replace the source mesh.
      CALL tria_refine2LevelOrdering(rtriangulation,rboundary=rboundary)
          
    END DO
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_refine2LevelOrdering(&
      rsourceTriangulation,rdestTriangulation,rboundary)

!<description>
  ! This routine refines the given mesh rsourceTriangulation according to
  ! the 2-level ordering algorithm. The refined mesh is saved in 
  ! rdestTriangulation.
  !
  ! rsourceTriangulation must be a 'standard' mesh. The resulting refined
  ! mesh in rdestTriangulation will be a 'raw' mesh (like as being read 
  ! from a .TRI file), i.e. the caller must add further information 
  ! to it with routines like tria_initStandardMeshFromRaw!
!</description>

!<input>
  ! OPTIONAL: Boundary structure that defines the parametrisation of the boundary.
  ! If specified, the coordinates of the new boundary vertices are
  ! recomputed according to the analytic boundary.
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary
!</input>

!<inputoutput>
  ! The source triangulation to be refined; must be a 'standard' mesh.
  TYPE(t_triangulation), INTENT(INOUT) :: rsourceTriangulation
!</inputoutput>

!<output>
  ! OPTIONAL: Destination triangulation structure that receives the refined mesh.
  ! The refined mesh will be a 'raw' mesh as being read from a .TRI file e.g..
  ! If not specified, the source triangulation rsourceTriangulation is replaced
  ! by the refined mesh.
  TYPE(t_triangulation), INTENT(OUT), OPTIONAL :: rdestTriangulation
!</output>
  
!</subroutine>
 
    TYPE(t_triangulation) :: rdestTria
    
    ! Call the correct submethod depending on the dimension.
    SELECT CASE (rsourceTriangulation%ndim)
    CASE (NDIM1D)
    
    CASE (NDIM2D)
      ! Refine the basic mesh
      CALL tria_refineMesh2lv2D(rsourceTriangulation,rdestTria)
      
      ! Refine the boundary
      CALL tria_refineBdry2lv2D(rsourceTriangulation,rdestTria,rboundary)
      
    CASE (NDIM3D)
    
    CASE DEFAULT
      CALL output_line ('Triangulation structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_generateStandardMeshFromRaw')
      CALL sys_halt()
    END SELECT
    
    ! Either copy rdestTria to rdestTriangulation or overwrite the source
    ! triangulation.
    IF (PRESENT(rdestTriangulation)) THEN
      rdestTriangulation = rdestTria
    ELSE
      CALL tria_done (rsourceTriangulation)
      rsourceTriangulation = rdestTria
    END IF

  CONTAINS

    ! ---------------------------------------------------------------
  
    SUBROUTINE tria_refineMesh2lv2D(rsourceTriangulation,rdestTriangulation)

    ! This routine refines the given 2D mesh rsourceTriangulation according to
    ! the 2-level ordering algorithm. The refined mesh is saved in 
    ! rdestTriangulation. There will be no correction of boundary
    ! vertices. Boundary parameter values are not handled here!

    ! The source triangulation to be refined
    TYPE(t_triangulation), INTENT(INOUT) :: rsourceTriangulation

    ! Destination triangulation structure that receives the refined mesg. 
    TYPE(t_triangulation), INTENT(OUT) :: rdestTriangulation
    
      ! local variables
      REAL(DP), DIMENSION(:,:), POINTER :: p_DcoordSource
      REAL(DP), DIMENSION(:,:), POINTER :: p_DcoordDest
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertAtElementSource
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertAtElementDest
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementSource
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertAtEdgeSource
      INTEGER(PREC_ELEMENTIDX) :: nquads,iel,iel1,iel2,iel3
      INTEGER(PREC_EDGEIDX) :: imt
      INTEGER(PREC_VERTEXIDX) :: ivt1,ivt2, ivtoffset, ivt
      INTEGER(I32), DIMENSION(2) :: Isize
      INTEGER :: nnve, ive
      REAL(DP) :: x,y
      
      ! Get the arrays with information of the source mesh.
      CALL storage_getbase_double2D (rsourceTriangulation%h_DvertexCoords,&
          p_DcoordSource)
      CALL storage_getbase_int2D (rsourceTriangulation%h_IverticesAtElement,&
          p_IvertAtElementSource)
      CALL storage_getbase_int2D (rsourceTriangulation%h_IedgesAtElement,&
          p_IedgesAtElementSource)
      CALL storage_getbase_int2D (rsourceTriangulation%h_IverticesAtEdge,&
          p_IvertAtEdgeSource)
      
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
      nnve = tria_getNNVE(rsourceTriangulation)
      
      IF ((nnve .LT. TRIA_NVETRI2D) .OR. (nnve .GT. TRIA_NVEQUAD2D)) THEN
      
        CALL output_line (&
            '2-level refinement supports only triangular and quad meshes!', &
            OU_CLASS_ERROR,OU_MODE_STD,'tria_refineMesh2lv2D')
        CALL sys_halt()
        
      END IF
      
      ! Initialise the basic mesh data in rdestTriangulation:
      
      ! 2D mesh
      rdestTriangulation%ndim = NDIM2D
      
      ! Every element is divided into 4 subelements
      rdestTriangulation%NEL = 4 * rsourceTriangulation%NEL
      rdestTriangulation%InelOfType(:) = 4 * rsourceTriangulation%InelOfType(:)

      ! We expect NVT+NMT+nquads new points.
      rdestTriangulation%NVT = &
          rsourceTriangulation%NVT + rsourceTriangulation%NMT + nquads
      
      ! Allocate memory for the new vertex coordinates and
      ! get the pointers to the coordinate array
      Isize = (/NDIM2D,INT(rdestTriangulation%NVT,I32)/)
      CALL storage_new2D ('tria_refineMesh2lv2D', 'DCORVG', Isize, ST_DOUBLE, &
          rdestTriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double2D(&
          rdestTriangulation%h_DvertexCoords,p_DcoordDest)
      
      ! Ok, let's start the refinement. In the first step, we copy the
      ! corner coordinates of the coarse mesh to the fine mesh; they
      ! don't change during the refinement.
      DO ivt=1,rsourceTriangulation%NVT
        p_DcoordDest(1,ivt) = p_DcoordSource(1,ivt)
        p_DcoordDest(2,ivt) = p_DcoordSource(2,ivt)
      END DO
      
      ! Each edge produces an edge midpoint which is stored as new
      ! point in the fine mesh. To calculate the coordinates, take
      ! the mean of the coordinates in the coarse mesh.
      ivtoffset = rsourceTriangulation%NVT
      DO imt=1,rsourceTriangulation%NMT
        ivt1 = p_IvertAtEdgeSource (1,imt)
        ivt2 = p_IvertAtEdgeSource (2,imt)
        p_DcoordDest(1,ivtoffset+imt) = &
            0.5_DP * ( p_DcoordSource (1,ivt1) + p_DcoordSource (1,ivt2) )
        p_DcoordDest(2,ivtoffset+imt) = &
            0.5_DP * ( p_DcoordSource (2,ivt1) + p_DcoordSource (2,ivt2) )
      END DO
      
      ! Allocate memory for IverticesAtElement and get a pointer to it.
      ! Fill the array with zero, so we won't have problems when mixing
      ! triangles into a quad mesh.
      Isize = (/nnve,INT(rdestTriangulation%NEL,I32)/)
      CALL storage_new2D ('tria_refineMesh2lv2D', 'KVERT', Isize, ST_INT, &
          rdestTriangulation%h_IverticesAtElement, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int2D(&
          rdestTriangulation%h_IverticesAtElement,p_IvertAtElementDest)
    
      ! Are there quads in the mesh? They produce additional element midpoints
      ! which get numbers NVT+NMT+1..*
      IF (nnve .GT. TRIA_NVETRI2D) THEN
      
        ! Loop over the elements to find the quads.
        ivtoffset = rsourceTriangulation%NVT + rsourceTriangulation%NMT
        DO iel = 1,rsourceTriangulation%NEL
        
          IF (p_IvertAtElementSource (TRIA_NVEQUAD2D,iel) .NE. 0) THEN
          
            ! New element midpoint
            ivtoffset = ivtoffset+1
            
            ! Sum up the coordinates of the corners to get the midpoint.
            x = 0.0_DP
            y = 0.0_DP
            DO ive = 1,TRIA_NVEQUAD2D
              x = x + p_DcoordSource (1,p_IvertAtElementSource(ive,iel))
              y = y + p_DcoordSource (2,p_IvertAtElementSource(ive,iel))
            END DO

            ! Store the midpoint
            p_DcoordDest(1,ivtoffset) = 0.25_DP*x
            p_DcoordDest(2,ivtoffset) = 0.25_DP*y
          
          END IF
        
        END DO
      
      END IF
      
      ! Ok, that was easy up to here. Now the interesting part: 
      ! the two-level ordeting.
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
      !  - The vertex numbers in the coars grid are transferred to the fine grid
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
      ! Let's loop over the vertices on the coarse grid. Everyone produces
      ! four elements on the fine grid:
      !
      ! In nquads we count the number of quads +NVT+NMT we reach.
      ! That's the number of the midpoint of that element!
      nquads = rsourceTriangulation%NVT+rsourceTriangulation%NMT
      
      IF (nnve .EQ. TRIA_NVETRI2D) THEN
        ! Pure triangle mesh
        DO iel = 1,rsourceTriangulation%NEL
      
          ! Determine number of subelements.
          iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
          iel2 = iel1+1
          iel3 = iel1+2
          
          ! Step 1: Initialise IverticesOnElement for element IEL
          p_IvertAtElementDest(1,iel) = p_IedgesAtElementSource (1,iel)
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (2,iel)
          p_IvertAtElementDest(3,iel) = p_IedgesAtElementSource (3,iel)
          
          ! Step 2: Initialise IverticesOnElement for element IEL1
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (1,iel)
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (1,iel)
          p_IvertAtElementDest(3,iel1) = p_IedgesAtElementSource (3,iel)

          ! Step 3: Initialise IverticesOnElement for element IEL2
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (2,iel)
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (2,iel)
          p_IvertAtElementDest(3,iel2) = p_IedgesAtElementSource (1,iel)

          ! Step 4: Initialise IverticesOnElement for element IEL3
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (3,iel)
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (3,iel)
          p_IvertAtElementDest(3,iel3) = p_IedgesAtElementSource (2,iel)
        
        END DO
      
      ELSE IF (nquads .EQ. rsourceTriangulation%NEL) THEN
        
        ! Pure QUAD mesh
        DO iel = 1,rsourceTriangulation%NEL
      
          ! Determine number of subelements.
          iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
          iel2 = iel1+1
          iel3 = iel1+2
          
          ! In nquads we count the number of the quad +NVT+NMT we process.
          ! That's the number of the midpoint of that element!
          ! As we reached a new quad, we increase nquads
          nquads = nquads+1
          
          ! Step 1: Initialise IverticesOnElement for element IEL
          p_IvertAtElementDest(1,iel) = p_IvertAtElementSource (1,iel)
          p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (1,iel)
          p_IvertAtElementDest(3,iel) = nquads
          p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource (4,iel)
          
          ! Step 2: Initialise IverticesOnElement for element IEL1
          p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (2,iel)
          p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (2,iel)
          p_IvertAtElementDest(3,iel1) = nquads
          p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource (1,iel)
        
          ! Step 3: Initialise IverticesOnElement for element IEL2
          p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (3,iel)
          p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (3,iel)
          p_IvertAtElementDest(3,iel2) = nquads
          p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource (2,iel)

          ! Step 4: Initialise IverticesOnElement for element IEL3
          p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (4,iel)
          p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (4,iel)
          p_IvertAtElementDest(3,iel3) = nquads
          p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource (3,iel)
        
        END DO ! iel
        
      ELSE
      
        ! Triangles and quads mixed.
      
        DO iel = 1,rsourceTriangulation%NEL
      
          ! Is that a triangle or a quad?  
          IF (p_IvertAtElementSource(TRIA_NVEQUAD2D,iel) .EQ. 0) THEN
          
            ! Triangular element.
            !
            ! Determine number of subelements.
            iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
            iel2 = iel1+1
            iel3 = iel1+2
            
            ! Step 1: Initialise IverticesOnElement for element IEL
            p_IvertAtElementDest(1,iel) = p_IedgesAtElementSource (1,iel)
            p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (2,iel)
            p_IvertAtElementDest(3,iel) = p_IedgesAtElementSource (3,iel)
            
            ! Step 2: Initialise IverticesOnElement for element IEL1
            p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (1,iel)
            p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (1,iel)
            p_IvertAtElementDest(3,iel1) = p_IedgesAtElementSource (3,iel)

            ! Step 3: Initialise IverticesOnElement for element IEL2
            p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (2,iel)
            p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (2,iel)
            p_IvertAtElementDest(3,iel2) = p_IedgesAtElementSource (1,iel)

            ! Step 4: Initialise IverticesOnElement for element IEL3
            p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (3,iel)
            p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (3,iel)
            p_IvertAtElementDest(3,iel3) = p_IedgesAtElementSource (2,iel)
          
          ELSE
          
            ! Quadrilateral element
            !
            ! Determine number of subelements.
            iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
            iel2 = iel1+1
            iel3 = iel1+2
            
            ! In nquads we count the number of the quad +NVT+NMT we process.
            ! That's the number of the midpoint of that element!
            ! As we reached a new quad, we increase nquads
            nquads = nquads+1
            
            ! Step 1: Initialise IverticesOnElement for element IEL
            p_IvertAtElementDest(1,iel) = p_IvertAtElementSource (1,iel)
            p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (1,iel)
            p_IvertAtElementDest(3,iel) = nquads
            p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource (4,iel)
            
            ! Step 2: Initialise IverticesOnElement for element IEL1
            p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (2,iel)
            p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (2,iel)
            p_IvertAtElementDest(3,iel1) = nquads
            p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource (1,iel)
          
            ! Step 3: Initialise IverticesOnElement for element IEL2
            p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (3,iel)
            p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (3,iel)
            p_IvertAtElementDest(3,iel2) = nquads
            p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource (2,iel)

            ! Step 4: Initialise IverticesOnElement for element IEL3
            p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (4,iel)
            p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (4,iel)
            p_IvertAtElementDest(3,iel3) = nquads
            p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource (3,iel)
          
          END IF
        
        END DO ! iel
        
      END IF
      
      ! The last step of setting up the raw mesh on the finer level:
      ! Set up InodalProperty. But that's the most easiest thing: Simply
      ! copy the nodal property array from the coarse mesh to the fine mesh.
      ! The nodal information of the edges (number NVT+1..NVT+NMT)
      ! that way converts into the nodal information about the new vertices
      ! on the fine mesh!
      CALL storage_copy (rsourceTriangulation%h_InodalProperty,&
          rdestTriangulation%h_InodalProperty)
    
    END SUBROUTINE

    ! ---------------------------------------------------------------
  
    SUBROUTINE tria_refineBdry2lv2D(rsourceTriangulation,rdestTriangulation,rboundary)

    ! This routine refines the boundary definition of rsourceTriangulation
    ! according to the 2-level ordering algorithm to generate a new 
    ! IverticesAtBoundary. 
    ! If rboundary is specified, the parameter values of the boundary vertices are 
    ! updated and the coordinates of the boundary points are corrected according 
    ! to the analytic boundarty. the boundary vertices are not sorted for their
    ! parameter value!

    ! The source triangulation to be refined
    TYPE(t_triangulation), INTENT(IN) :: rsourceTriangulation

    ! Destination triangulation structure that receives the refined mesg. 
    TYPE(t_triangulation), INTENT(INOUT) :: rdestTriangulation
    
    ! OPTIONAL: Defintion of analytic boundary.
    ! If specified, the coordinates of the new boundary vertices are
    ! recomputed according to the analytic boundary.
    TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary

      ! local variables
      REAL(DP), DIMENSION(:), POINTER :: p_DvertParamsSource
      REAL(DP), DIMENSION(:), POINTER :: p_DedgeParamsSource
      REAL(DP), DIMENSION(:), POINTER :: p_DvertParamsDest
      REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordDest
      INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IvertAtBoundartySource
      INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IedgesAtBoundartySource
      INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IvertAtBoundartyDest
      INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdxSource
      INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdxDest
      INTEGER :: ivbd,ibct
      
      ! Get the definition of the boundary vertices and -edges.
      CALL storage_getbase_int (rsourceTriangulation%h_IverticesAtBoundary,&
          p_IvertAtBoundartySource)
      CALL storage_getbase_int (rsourceTriangulation%h_IedgesAtBoundary,&
          p_IedgesAtBoundartySource)
      CALL storage_getbase_int (rsourceTriangulation%h_IboundaryCpIdx,&
          p_IboundaryCpIdxSource)

      ! The number of boundary vertices is doubled.
      rdestTriangulation%NVBD = 2*rsourceTriangulation%NVBD
      rdestTriangulation%NBCT = rsourceTriangulation%NBCT 
          
      ! Create new arrays in the fine grid for the vertices and indices.
      CALL storage_new ('tria_refineBdry2lv2D', &
          'KVBD', INT(rdestTriangulation%NVBD,I32), &
          ST_INT, rdestTriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)

      CALL storage_new ('tria_generateBasicBoundary', &
          'KBCT', INT(rdestTriangulation%NBCT+1,I32), &
          ST_INT, rdestTriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT)

      CALL storage_getbase_int (rdestTriangulation%h_IverticesAtBoundary,&
          p_IvertAtBoundartyDest)
      CALL storage_getbase_int (rdestTriangulation%h_IboundaryCpIdx,&
          p_IboundaryCpIdxDest)

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
      
      DO ibct = 1,SIZE(p_IboundaryCpIdxDest)
        p_IboundaryCpIdxDest(ibct) = 2*(p_IboundaryCpIdxSource(ibct)-1) + 1
      END DO
      
      DO ivbd = 0,SIZE(p_IvertAtBoundartySource)-1
        p_IvertAtBoundartyDest(1+2*ivbd) = p_IvertAtBoundartySource(1+ivbd)
        p_IvertAtBoundartyDest(2+2*ivbd) = p_IedgesAtBoundartySource(1+ivbd)
      END DO
      
      ! Let's see if parameter values of boundary vertices are available.
      IF ((rsourceTriangulation%h_DvertexParameterValue .NE. ST_NOHANDLE) .AND. &
          (rsourceTriangulation%h_DedgeParameterValue .NE. ST_NOHANDLE)) THEN

        ! Also interleave the parameter values of the vertices and edge midpoints.
        CALL storage_getbase_double (rsourceTriangulation%h_DvertexParameterValue,&
            p_DvertParamsSource)
        CALL storage_getbase_double (rsourceTriangulation%h_DedgeParameterValue,&
            p_DedgeParamsSource)
        
        ! Create a new array for the boundary vertex parameter values
        ! and fill it with the data from the coarse grid.
        CALL storage_new ('tria_refineBdry2lv2D', &
            'DVBDP', INT(rdestTriangulation%NVBD,I32), &
            ST_DOUBLE, rdestTriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double (rdestTriangulation%h_DvertexParameterValue,&
            p_DvertParamsDest)
            
        DO ivbd = 0,SIZE(p_IvertAtBoundartySource)-1
          p_DvertParamsDest(1+2*ivbd) = p_DvertParamsSource(1+ivbd)
          p_DvertParamsDest(2+2*ivbd) = p_DedgeParamsSource(1+ivbd)
        END DO
        
        ! If the analytic boundary is given, compute the coordinates of the
        ! boundary vertices from that.
        IF (PRESENT(rboundary)) THEN
          
          ! Get the array with the vertex coordinates.
          ! We want to correct the coordinates of the boundary points
          ! according to the analytic boundary.
          CALL storage_getbase_double2d (rdestTriangulation%h_DvertexCoords,&
              p_DcornerCoordDest)
            
          ! Loop through the boundary points and canculate the correct
          ! coordinates. 
          DO ibct = 1,rdestTriangulation%NBCT
            DO ivbd = p_IboundaryCpIdxDest(ibct),p_IboundaryCpIdxDest(ibct+1)-1
              CALL boundary_getCoords(rboundary,ibct,p_DvertParamsDest(ivbd),&
                  p_DcornerCoordDest(1,p_IvertAtBoundartyDest(ivbd)),&
                  p_DcornerCoordDest(2,p_IvertAtBoundartyDest(ivbd)))
            END DO
          END DO
          
        END IF
        
      END IF
    
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<function>

  INTEGER FUNCTION tria_searchBoundaryNode(inode,rtriangulation)

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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: inode

  ! The triangulation structure where to search the boundary node.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</input>

!<result>
  ! If inode is a boundary vertex: The index of the inode in IboundaryVertexPos.
  ! If inode is a boundary edge: The index of the inode in IboundaryEdgePos.
  ! =0 if inode was not found (e.g. because inode is not on the boundary e.g.).
!</result>

!</function>

    ! local variables
    INTEGER :: ibct
    INTEGER(I32), DIMENSION(:), POINTER :: p_InodalProperty
    INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_InodePos
    INTEGER :: ipos,ileft,iright
    
    ! Get the boundary component of the vertex
    CALL storage_getbase_int (rtriangulation%h_InodalProperty,p_InodalProperty)
    ibct = p_InodalProperty(inode)

    ! We can quit if ibct=0: The node is not on the boundary.
    IF (ibct .EQ. 0) THEN
      tria_searchBoundaryNode = 0
      RETURN
    END IF
    
    ! Do we have a vertex or an edge number?
    IF (inode .LE. rtriangulation%NVT) THEN
      ! Search in the IboundaryVertexPos array
      CALL storage_getbase_int2d (rtriangulation%h_IboundaryVertexPos,p_InodePos)
    ELSE
      ! Search in the IboundaryEdgePos array
      CALL storage_getbase_int2d (rtriangulation%h_IboundaryEdgePos,p_InodePos)
    END IF

    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Use bisection search to find the node in the array.
    ileft = p_IboundaryCpIdx(ibct)
    iright = p_IboundaryCpIdx(ibct+1)-1
    
    IF (p_InodePos(1,ileft) .EQ. inode) THEN
      ! Return the index in the node array.
      tria_searchBoundaryNode = p_InodePos(2,ileft)
    ELSE IF (p_InodePos(1,iright) .EQ. inode) THEN
      ! Return the index in the node array.
      tria_searchBoundaryNode = p_InodePos(2,iright)
    ELSE
      DO WHILE (ileft .LT. iright)
        ipos = (ileft+iright)/2
        IF (p_InodePos(1,ipos) .GT. inode) THEN
          iright = ipos
        ELSE IF (p_InodePos(1,ipos) .LT. inode) THEN
          ileft = ipos
        ELSE
          ! We found the node. Return the index in the node array.
          tria_searchBoundaryNode = p_InodePos(2,ipos)
          EXIT
        END IF
      END DO
    END IF
    
  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_compress2LevelOrdHierarchy (rtriangulationFine,rtriangulationCoarse)
  
!<description>
  ! This routine can be used to save memory after a refinement with the
  ! 2-level-ordering algorithm was applied. It's typically applied 'backwards'
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
  !   DO I=NLMAX-1,NLMIN,-1             
  !     CALL tria_compress2LevelOrdHierarchy (rtria(i+1),rtria(i))
  !   END DO
  !
  ! Afterwards, the finest mesh contains all coordinates and the coarse grid
  ! coordinates are part of the fine grid coordinates.
  !
  ! WARNING: Use this routine with care. It does not check whether
  ! rtriangulationCoarse,rtriangulationFine are 'compatible' to each other.
!</description>

!<input>
  ! Fine grid triangulation.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulationFine
!</input>

!<inputoutput>
  ! Coarse grid triangulation where redundant data should be removed from.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulationCoarse
!</inputoutput>
  
!</subroutine>
  
    ! Release the vertex coordinates array in the coarse triangulation -- as long
    ! as it's not a copy of someone else...
    IF (IAND(rtriangulationCoarse%iduplicationFlag,TR_SHARE_DVERTEXCOORDS) .EQ. 0) THEN
      CALL storage_free (rtriangulationCoarse%h_DvertexCoords)
    END IF 
    
    ! Share the same handle.
    rtriangulationCoarse%h_DvertexCoords = rtriangulationFine%h_DvertexCoords
    
    ! Mark the DvertexCoords array in the source mesh as 'being a copy
    ! of someone else', so the array is not released when the coarse
    ! mesh is released!
    rtriangulationCoarse%iduplicationFlag = &
      IOR(rtriangulationCoarse%iduplicationFlag,TR_SHARE_DVERTEXCOORDS)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_getPointsOnEdge (rtriangulation,Dcoords,npointsPerEdge, &
      DparValue,rboundary)
  
!<description>
  ! This routine can be used to calculate the coordinates of a fixed number
  ! of regularly distributed (inner-edge) points per edge.
  !
  ! npointsPerEdge specifies how many points per edge should be calculated.
  ! E.g. if npointsPerEdge=3, the routine will place 3 points
  ! on each edge and calculate their coordinates. By default, the points will
  ! be distributed regularly, so in this example at the 'relative position'
  ! {0.25, 0.5, 0.75) of the edge. If the caller wants to specify a different
  ! position of the points manually, the DparValue array can be specified.
  !
  ! Another example: By calling this routine with npointsPerEdge=1, the 
  ! routine will calculate the midpoints of all edges.
!</description>

!<input>
  ! Triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! Number of points per edge to generate
  INTEGER, INTENT(IN) :: npointsPerEdge
  
  ! OPTIONAL: Array with parameter values of the points on the edge.
  ! DIMENSION(npointsPerEdge). DparValue is a value in the range [0,1]
  ! and specifies the 'relative position' or 'parameter value' of each 
  ! point on the edge. If not specified, tria_getPointsOnEdge assumes a regular
  ! distribution of points, defined by the number of points on the edge.
  ! E.g. if npointsPerEdge=2 and DparValue is not specified, 
  ! the routine assumes 3 inner points on the edge corresponding to 
  ! DparValue=/(0.333333,0.666666)/.
  ! If specified, the caller can specify the exact parameter values of
  ! the three points on the edge, e.g. DparValue=/(0.25,0.75)/.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DparValue
  
  ! OPTIONAL: Definition of the domain.
  ! If specified, the routine will calculate the points on the boundary
  ! edges using the definition on the boundary. If not specified, the routine
  ! calculates the points on the edges based only on the triangulation.
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary
!</input>

!<output>
  ! Coordinates of the regularly distributed points on all the edges
  ! of the triangulation.
  !   DIMENSION(space-dimension, npointsPerEdge * #edges)
  ! Here, Dcoords(:,1..npointsPerEdge) receives the coordinates of the 
  ! points on the first edge, Dcoords(:,npointsPerEdge+1:2*npointsPerEdge)
  ! the coordinates on edge number 2 etc.
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dcoords
!</output>

!</subroutine>

    ! local variables      
    REAL(DP), DIMENSION(npointsPerEdge) :: Dparameters
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dcoords
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
    INTEGER(PREC_EDGEIDX) :: iedge
    INTEGER(PREC_VERTEXIDX) :: ipointpos,ipoint1,ipoint2
    INTEGER :: idim,ipoint
    INTEGER :: ibdc,ibdedge
    REAL(DP) :: dpar1,dpar2
    
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_EDGEIDX), DIMENSION(:), POINTER :: p_IedgesAtBoundary
    INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    
    ! If DparValue is specified, take that. Otherwise, create the parameter
    ! values of the (inner-edge) points manually.
    IF (PRESENT(DparValue)) THEN
      IF (SIZE(DparValue) .LT. npointsPerEdge) THEN
        CALL output_line ('DparValue not large enough!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'tria_getPointsOnEdge')
        CALL sys_halt()
      END IF
      Dparameters(:) = DparValue(1:npointsPerEdge)
    END IF
    
    ! Get the triangulation stuff
    IF (rtriangulation%ndim .EQ. 0) THEN
      CALL output_line ('Triangulation not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_getPointsOnEdge')
      CALL sys_halt()
    END IF

    IF (rtriangulation%h_IverticesAtEdge .EQ. 0) THEN
      CALL output_line ('IverticesAtEdge not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_getPointsOnEdge')
      CALL sys_halt()
    END IF
    
    CALL storage_getbase_double2d(rtriangulation%h_DvertexCoords,p_Dcoords)
    CALL storage_getbase_int2d(rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
  
    ! Loop through all edges
    DO iedge = 0,rtriangulation%NMT-1
    
      ! Where do the set of points start in Dcoords?
      ipointpos = iedge*npointsPerEdge
      
      ! Endpoints of that edge?
      ipoint1 = p_IverticesAtEdge(1,iedge)
      ipoint2 = p_IverticesAtEdge(2,iedge)
    
      ! Calculate the points on the edge
      DO ipoint = 1,npointsPerEdge
    
        ! Calculate the point coordinates
        DO idim = 1,UBOUND(Dcoords,1)
          
          Dcoords(idim,ipoint+ipointpos) = &
            p_Dcoords(idim,ipoint1) * Dparameters(ipoint) + &
            p_Dcoords(idim,ipoint2) * (1.0_DP-Dparameters(ipoint))
        
        END DO
      
      END DO
      
    END DO
    
    ! Is the boundary structure given?
    IF (PRESENT(rboundary)) THEN
    
      ! 2D?
      IF (rtriangulation%ndim .EQ. NDIM2D) THEN
      
        ! Ok, we can calculate the correct position of the points on the boundary!
        ! Get the array with the boundary vertices and their parameter values.
        CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
            p_IboundaryCpIdx)
        CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
            p_IverticesAtBoundary)
        CALL storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
            p_IedgesAtBoundary)
        CALL storage_getbase_double (rtriangulation%h_DvertexParameterValue,&
            p_DvertexParameterValue)
            
        ! Loop through the boundary components and the points in each component
        DO ibdc = 1,rtriangulation%NBCT
        
          DO ibdedge = p_IboundaryCpIdx(ibdc),p_IboundaryCpIdx(ibdc+1)-1
          
            ! Get the boundary edge
            iedge = p_IedgesAtBoundary(ibdedge)
            
            ! Get the parameter value of the points adjacent to that edge.
            dpar1 = p_DvertexParameterValue(ibdedge)
            IF (ibdedge .NE. p_IboundaryCpIdx(ibdc+1)-1) THEN
              dpar2 = p_DvertexParameterValue(ibdedge+1)
            ELSE
              ! Last edge ends with maximum parameter value on that boundary component
              dpar2 = boundary_dgetMaxParVal(rboundary,ibdc)
            END IF
            
            ! Where do the set of points start in Dcoords?
            ipointpos = iedge*npointsPerEdge
            
            ! Calculate the points on the edge
            DO ipoint = 1,npointsPerEdge
          
              ! Calculate the point coordinates
              CALL boundary_getCoords(rboundary, ibdc, &
                  dpar1 * Dparameters(ipoint) + dpar2 * (1.0_DP-Dparameters(ipoint)),  &
                  Dcoords(1,ipoint+ipointpos), Dcoords(2,ipoint+ipointpos))
                  
            END DO
          
          END DO
        
        END DO
      
      END IF
    
    END IF
  
  END SUBROUTINE

END MODULE

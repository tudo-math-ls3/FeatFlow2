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
!# 4.) tria_done
!#     -> Cleans up a triangulation structure, releases memory from the heap.
!#
!# 5.) tria_rawQuadToTri
!#     -> Converts a raw mesh into a triangular mesh.
!#
!# 6.) tria_duplicate
!#     -> Creates a duplicate / backup of a triangulation.
!#        Some information may be shared between two triangulation structures.
!#
!# 7.) tria_restore
!#     -> Restores a triangulation previously backed up with tria_backup.
!#
!# 8.) tria_recoverHandles
!#     -> Recovers destroyed handles of a triangulation structure from a 
!#        backup.
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
  
  ! Maximum number of corner-vertices in each element.
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
  INTEGER, PARAMETER :: TRIA_NVETRIA2D = 3

  ! Number of vertices per element for quadrilateral element shapes in 2D.
  INTEGER, PARAMETER :: TRIA_NVEQUAD2D = 4
  
!</constantblock>

!<constantblock description="KIND values for triangulation data">
  
  ! kind value for indexing the points in a triangulation
  INTEGER, PARAMETER :: PREC_POINTIDX   = I32

  ! kind value for indexing the edges in a triangulation
  INTEGER, PARAMETER :: PREC_EDGEIDX    = I32

  ! kind value for indexing the elements in a triangulation
  INTEGER, PARAMETER :: PREC_ELEMENTIDX = I32

!</constantblock>
  
!<constantblock description="Duplication flags. Specifies which information is shared \
!                            between triangulation structures">

  INTEGER(I32), PARAMETER :: TR_SHARE_DCORNERCOORDINATES     = 2** 0  ! DCORVG 
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
    INTEGER(PREC_POINTIDX), DIMENSION(TRIA_MAXNVE2D) :: Icorners
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
    INTEGER(PREC_POINTIDX)   :: NVT = 0
    
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
  
    ! Number of vertices per edge; normally = 0.
    ! If a regular distribution of vertices on edges is given,
    ! NVPED saves the number of vertices on each edge;
    ! e.g. 1 if midpoints on edges exist in p_DfreeVertexCoordinates.
    INTEGER             :: nverticesPerEdge = 0
    
    ! Total number of vertices on edges; normally = 0.
    ! Total number of vertices on all edges, realized in p_DfreeVertexCoordinates. 
    ! E.g. if midpoints on edges exist, there is NVEDT=NMT.
    INTEGER(PREC_POINTIDX)   :: nVerticesOnAllEdges = 0
    
    ! Number of inner-element vertices; normally = 0.
    ! If a regular distribution of vertices in the inner of 
    ! each element is given, NIELV saves the number of vertices 
    ! in the inner of each element; e.g. 1 if element-midpoints 
    ! exist in p_DfreeVertexCoordinates.
    INTEGER             :: nverticesInEachElement = 0

    ! Total number of vertices in elements; normally = 0.
    ! Total number of vertices on all elements, realized in p_DfreeVertexCoordinates. 
    ! E.g. if element-midpoints exist in DCORMG, there is NIEVT=NEL.
    INTEGER(PREC_POINTIDX)   :: nverticesInAllElements = 0
    
    ! Number of additional vertices; normally = 0.
    ! Can be set <> 0 if there are any additional vertices 
    ! realized in p_DfreeVertexCoordinates, that don't belong to a regular 
    ! distribution of vertices in corners, on edges or on elements.
    INTEGER(PREC_POINTIDX)    :: nadditionalVertices = 0
  
    ! A list of all corner(!)-vertices of the elements in the triangulation.
    ! Handle to 
    !       p_RcornerCoordinates = array [1..NDIM2D,1..NVT] of double
    ! with
    !   p_DcornerCoordinates(1,.) = X-coordinate.
    !   p_DcornerCoordinates(2,.) = Y-coordinate.
    ! This is a handle to the old DCORVG-array.
    INTEGER        :: h_DcornerCoordinates = ST_NOHANDLE
    
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
    ! Number on the  is sorted for the vertex number, 
    ! thus allowing quick access to the index of a vertex
    ! in p_IverticesAtBoundary.
    INTEGER          :: h_IboundaryVertexPos = ST_NOHANDLE

    ! Inverse index array to p_IedgesAtBoundary. 
    ! Handle to 
    !       p_IboundaryEdgePos = array [1..2,1..NMBD] of integer.
    ! p_IboundaryEdgePos(1,.) contains a vertex number on the 
    ! boundary and p_IboundaryEdgePos(2,.) the appropriate index
    ! of this vertex inside of the p_IedgesAtBoundary-array. 
    ! Number on the  is sorted for the vertex number, 
    ! thus allowing quick access to the index of a vertex
    ! in p_IedgesAtBoundary.
    INTEGER          :: p_IboundaryEdgePos = ST_NOHANDLE
    
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
                                rtriangulation%h_DcornerCoordinates)
  
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
    CALL checkAndCopy(idupflag, TR_SHARE_DCORNERCOORDINATES,&
          rtriangulation%h_DcornerCoordinates, &
          rbackupTriangulation%h_DcornerCoordinates)

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
          rtriangulation%p_IboundaryEdgePos, &
          rbackupTriangulation%p_IboundaryEdgePos)
    
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
    CALL checkAndCopy(idupflag, TR_SHARE_DCORNERCOORDINATES,&
          rtriangulation%h_DcornerCoordinates, &
          rbackupTriangulation%h_DcornerCoordinates)

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
          rtriangulation%p_IboundaryEdgePos, &
          rbackupTriangulation%p_IboundaryEdgePos)
    
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

  SUBROUTINE tria_recoverHandles (rtriangulation,rbackupTriangulation)
  
!<description>
  ! With this routine, destroyed handles of shared information between
  ! two triangulation structures can be recovered. TRIA is assumed to
  ! be a valid triangulation structure. The handles of all arrays
  ! of information that TRIA shares with TRIADS are copied to TRIADS.
  !
  ! This routine is typically used in refinement processes when e.g.
  ! coordinates of grid points are shared between two levels of
  ! refinement. Upon refining the coarser grid, the handle in the
  ! source grid gets invalid, while the handle of the fine grid is
  ! the new valid one. TRICSH can now be used to copy the handle
  ! back from the fine grid triangulation structure to the coarse
  ! grid triangulation structure.
  !
  ! Example: DcornerCoordinates is shared between all refinement levels.
  ! Then the caller can create the grid by:
  !
  !   CALL GENTRI (IMETH=2, filename='anything.tri')  -> read coarse gr.
  !   DO I=NLMIN+1,NLMAX                              -> refine to NLMAX
  !     CALL GENTRI (TRIA(I-1),TRIA(I),IMETH=1,IDPFLG=1)
  !     --> generates TRIA(I), whereas destroys LCORVG on level I-1!
  !   END DO
  !   DO I=NLMAX-1,NLMIN,-1             -> write LCORVG of lv. n+1 into
  !     CALL tria_recoverHandles (TRIA(I),TRIA(I+1))  -> h_DcornerCoordinates of level n
  !   END DO
  !
  ! Afterwards, the information is again shared correctly between all
  ! levels.
  !
  ! WARNING: Use this routine with care. It does not check whether
  ! rtriangulation and rbackupTriangulation are 'compatible' to each other
  ! and may throw away old handles which may lead to memory leaks!
  ! Only apply this routine to handles of memory block which are really
  ! invalid!
!</description>

!<input>
  ! Valid triangulation structure, which is a modified copy of 
  ! rtriangulation.
  TYPE(t_triangulation), INTENT(IN) :: rbackupTriangulation
!</input>

!<inputoutput>
  ! Destination triangulation.
  ! The handles of all arrays that are shared with rbackupTriangulation
  ! are copied to rtriangulation. Old handles in rtriangulation are assumed
  ! to be invalid and will be overwritten by those of rbackupTriangulation.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>
  
    ! local variables
    INTEGER(I32) :: idupFlag

    idupFlag = rtriangulation%iduplicationFlag
      
    ! Call checkAndCopy for all the arrays. this will either copy the handle
    ! or copy the content of the array.
    
    ! Bit  0: DCORVG 
    CALL checkAndCopy(idupflag, TR_SHARE_DCORNERCOORDINATES,&
          rtriangulation%h_DcornerCoordinates, &
          rbackupTriangulation%h_DcornerCoordinates)

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
          rtriangulation%p_IboundaryEdgePos, &
          rbackupTriangulation%p_IboundaryEdgePos)
    
    ! Bit  1: DCORMG 
    CALL checkAndCopy(idupflag, TR_SHARE_DFREEVERTEXCOORDINATES,&
          rtriangulation%h_DfreeVertexCoordinates, &
          rbackupTriangulation%h_DfreeVertexCoordinates)
    
  CONTAINS
  
    SUBROUTINE checkAndCopy (idupFlag,ibitfield,idesthandle,isourcehandle)
    
    ! Checks if idupFlag has all bits ibitfield set.
    ! If yes, the handle isourcehandle is copied to idesthandle.
    
    INTEGER(I32), INTENT(IN) :: ibitfield
    INTEGER(I32), INTENT(IN) :: idupFlag
    INTEGER, INTENT(IN) :: isourcehandle
    INTEGER, INTENT(INOUT) :: idesthandle
    
      IF (IAND(idupFlag,ibitfield) .EQ. ibitfield) THEN
        idesthandle = isourcehandle
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
    CALL checkAndRelease(idupflag, TR_SHARE_DCORNERCOORDINATES,&
          rtriangulation%h_DcornerCoordinates)

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
          rtriangulation%p_IboundaryEdgePos)
    
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

  SUBROUTINE tria_rawQuadToTri (rtriangulation)
  
!<description>
  ! This routine converts a 2D mesh into a triangular 2D mesh. All quads are 
  ! converted to triangles.
  ! Warning: This routine can only applied to a 'raw' mesh, e.g. a mesh which 
  !  comes directly from a .TRI file. It's not advisable to apply this routine to
  !  a 'standard' mesh that contains already further mesh information 
  !  (adjacencies,...), as any information about possible two-level ordering
  !  is lost!
!</description>

!<inputoutput>
  ! The triangulation structure to be converted. Is replaced by a triangular
  ! mesh.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>

    INTEGER(PREC_ELEMENTIDX) :: i
    INTEGER :: icount
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER :: h_IverticesAtElementTri
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementTri
    INTEGER(I32), DIMENSION(2) :: Isize
   
    ! There are some things missing... i.e. the calculation of adjacencies
    PRINT *,'Conversion to triangular mesh not yet fully implemented!'
    STOP
   
    ! For this routine we currently assume that there are only triangles
    ! and quads in the triangulation. Might be a matter of change in 
    ! the future...
   
    ! Get the points-at-element array
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
    ! Count the quads; they are recogniseable by having the 4th corner vertex
    ! number nonzero.
    icount = 0
    DO i=1,rtriangulation%NEL
      IF (p_IverticesAtElement(4,i) .NE. 0) icount = icount+1
    END DO
    
    ! Create a new p_IverticesAtElement array for the triangular mesh.
    Isize = (/TRIA_MAXNVE2D,icount+rtriangulation%NEL/)
    CALL storage_new2D ('tria_quadToTri', 'KVERTTRI', Isize, ST_INT, &
                        h_IverticesAtElementTri,ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int2d (h_IverticesAtElementTri,p_IverticesAtElementTri)
    
    ! Convert the array
    CALL quadToTriang_aux1 (rtriangulation%NEL, icount,&
                            p_IverticesAtElement, p_IverticesAtElementTri)
    
    ! Replace the old array by the new, release the old one.
    CALL storage_free (rtriangulation%h_IverticesAtElement)
    rtriangulation%h_IverticesAtElement = h_IverticesAtElementTri
    
    ! That's it.

  CONTAINS

    SUBROUTINE quadToTriang_aux1 (nel, nelquad, IverticesAtElement, Kvert_triang)

  !<description>
    ! Purpose: Convert quad mesh to triangular mesh
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
        
      ! Copy the old IverticesAtElement two times, once to the first half and
      ! once the quads in it to the second half of Kvert_triang:
      DO i=1,nel
        Kvert_triang(:,i) = IverticesAtElement(:,i)
      END DO
      
      j = 0
      DO i=1,nel
        IF (IverticesAtElement(4,i) .NE. 0) THEN
          j=j+1
          Kvert_triang(:,nel+j) = IverticesAtElement(:,i)
        END IF
      END DO

      ! Correct Kvert_triang:
      DO i=1,nel
        ! Set the 4th entry in the first half of Kvert_triang to 0.
        ! So the first triangle in each QUAD element consists of the
        ! first three vertices 1..3.
        Kvert_triang(4,i) = 0
      END DO

      DO i=1,nelquad
        ! The second triangle in each quad consists of vertices 1,3,4.
        ! So in the 2nd half, shift entry 3->2 and 4->3 and set the 4th
        ! entry to 0.
        Kvert_triang(2,nel+i) = Kvert_triang(3,nel+i)
        Kvert_triang(3,nel+i) = Kvert_triang(4,nel+i)
        Kvert_triang(4,nel+i) = 0
      END DO

      ! That's it.
      
    END SUBROUTINE  
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_readTriFile2D(rboundary, rtriangulation, sfilename)

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
  ! NEL,NVT,NMT,NBCT,NVBD,
  ! DcornerCoordinates, IverticesAtElement, InodalProperty, 
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
  ! An rboundary object specifying the underlying domain.
  TYPE(t_boundary), INTENT(IN) :: rboundary

  ! The name of the .prm file to read.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
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
    CALL tria_genRawBoundary2D (rboundary,rtriangulation)

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
  ! NEL,NVT,NMT,NBCT,
  ! DcornerCoordinates, IverticesAtElement and InodalProperty.
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
        rtriangulation%h_DcornerCoordinates, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    CALL storage_getbase_double2D(&
        rtriangulation%h_DcornerCoordinates,p_Ddata2D)
        
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

  SUBROUTINE tria_genRawBoundary2D (rboundary,rtriangulation)

!<description>  
  ! Auxiliary routine of tria_readTriFile2D.
  ! This routine initialises basic boundary arrays and cleans up 
  ! a basic triangulation. That means:
  ! -> NVBD is calculated
  ! -> IboundaryCpIdx is created and generated
  ! -> IverticesAtBoundary is created and generated.
  !    The vertices are ordered for the boundary component according
  !    to IboundaryCpIdx but not ordered for their parameter value.
  ! -> DvertexParameterValue is created and generated (unsorted)
  ! -> The parameter values of the boundary vertices are extracted
  !    from the first coordinate in DcornerCoordinates and put into
  !    the DvertexParameterValue array.
  ! -> Based on the parameter value of each boundary vertex, the 
  !    DcornerCoordinates array receives the actual coordinates of 
  !    the boundary vertices.
!</description>
  
!<input>
  ! The parametrisation that specifies the coordinates of the boundary points
  TYPE(t_boundary), INTENT(IN) :: rboundary
!</input>
  
!<inputoutput>
  ! Triangulation to be initialised with basic data.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
    REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_POINTIDX) :: ivbd,ivt
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

    ! Allocate memory for IverticesAtBoundary and DvertexParameterValue.
    CALL storage_new ('tria_generateBasicBoundary', &
        'KVBD', INT(rtriangulation%NVBD,I32), &
        ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    CALL storage_new ('tria_generateBasicBoundary', &
        'DVBD', INT(rtriangulation%NVBD,I32), &
        ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
    
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    CALL storage_new ('tria_generateBasicBoundary', &
        'KBCT', INT(rtriangulation%NBCT+1,I32), &
        ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    CALL storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    CALL storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
        
    CALL storage_getbase_double2D (&
        rtriangulation%h_DcornerCoordinates,p_DcornerCoordinates)
        
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
        
        ! Store the parameter value; it's saved in DcornerCoordinates(1,.)
        p_DvertexParameterValue (ivbd) = p_DcornerCoordinates(1,ivt)
        
        ! Replace the coordinates in DcornerCoordinates by those
        ! given by the parametrisation.
        CALL boundary_getCoords(rboundary, ibct, p_DvertexParameterValue (ivbd), &
            p_DcornerCoordinates(1,ivt), p_DcornerCoordinates(2,ivt))
        
      END IF
    END DO
    
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
!</description>

!<input>
  ! OPTIONAL: Boundary structure that defines the domain.
  ! Must be specified for 2D domains.
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary
!</input>

!<inputoutput>
  ! Triangulation structure to be initialised
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    SELECT CASE (rtriangulation%ndim)
    CASE (NDIM1D)
    CASE (NDIM2D)
      
      IF (.NOT. PRESENT(rboundary)) THEN
        CALL output_line ('rboundary not available!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'tria_initStandardMeshFromRaw')
        STOP
      END IF
      
      ! Generate all standard arrays for 2D meshes.
      CALL tria_sortBoundaryVertices2D   (rtriangulation)
      CALL tria_genElementsAtVertex2D    (rtriangulation)
      CALL tria_genNeighboursAtElement2D (rtriangulation)
      CALL tria_genEdgesAtElement2D      (rtriangulation)
      CALL tria_genElementsAtEdge2D      (rtriangulation)
      CALL tria_genVerticesAtEdge2D      (rtriangulation)
      CALL tria_genEdgeNodalProperty2D   (rtriangulation)
      CALL tria_genElementVolume2D       (rtriangulation)
      CALL tria_genElementsAtBoundary2D  (rtriangulation)
      CALL tria_genEdgesAtBoundary2D     (rtriangulation)
      CALL tria_genEdgeParameterValue2D  (rtriangulation,rboundary)
    CASE (NDIM3D)
    CASE DEFAULT
      CALL output_line ('Triangulation structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_generateStandardMeshFromRaw')
      STOP
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
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_POINTIDX) :: istart,iend
    INTEGER :: ibct,ivbd
    INTEGER :: hresort
    INTEGER(I32), DIMENSION(:), POINTER :: p_Iresort

    ! Get pointers to the arrays
    CALL storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    CALL storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
        
    CALL storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! Sort the p_DcornerCoordinates (sub)-arrays for increasing
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
      CALL sort_dp(p_DvertexParameterValue(istart:iend),Imapping=p_Iresort(istart:iend))
      
    END DO

    ! Resort the vertices according to the calculated mapping.
    ! Afterwards, the vertices are ordered according to the increasing
    ! parameter value.
    p_Iresort(:) = p_IverticesAtBoundary(p_Iresort(:))
    p_IverticesAtBoundary(:) = p_Iresort(:)
    
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
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>
  
!<result>
  ! Maximum number of vertices per element.
!</result>
  
!</function>

    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_Idata2D
    
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
    INTEGER(PREC_POINTIDX) :: ivt,isize
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    
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
    
    ! Set the first element in p_IverticesAtElement to 1. Then, add
    ! all the length information together to get the index array.
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_POINTIDX) :: ivt,ivtneighbour
    
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_POINTIDX) :: NVT
    INTEGER(PREC_EDGEIDX) :: iedge
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      STOP
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      STOP
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      STOP
    END IF

    IF (rtriangulation%NMT .EQ. 0) THEN
      CALL output_line ('Edge information (NMT) not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge2D')
      STOP
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

          iedge = p_IedgesAtElement (ive,iel)-NVT
          p_IelementsAtEdge(1,iedge) = iel
          p_IelementsAtEdge(2,iedge) = 0
          
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive,nnve
    INTEGER(PREC_POINTIDX) :: ivtneighbour
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_EDGEIDX) :: iedge
    INTEGER(PREC_POINTIDX) :: NVT
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      STOP
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      STOP
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genVerticesAtEdge2D')
      STOP
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_EDGEIDX) :: iedge
    INTEGER(I32), DIMENSION(2) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
      STOP
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
      STOP
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
        
        ! Check the neighbour element. If the neighbour element has a number less
        ! than our current element number (including 0 which is the case for
        ! boundary edges), we 'see' the edge the first time and so
        ! we give it a number.
        IF (p_IneighboursAtElement(ive,iel) .LT. iel) THEN
        
          iedge = iedge + 1
          
          ! Add NVT to iedge to get the edge number
          p_IedgesAtElement(ive,iel) = iedge
        
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER :: ive
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_POINTIDX) :: Isize

    ! Is everything here we need?
    IF (rtriangulation%h_InodalProperty .EQ. ST_NOHANDLE) THEN
      CALL output_line ('InodalPropertys not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      STOP
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      STOP
    END IF
    
    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgeNodalProperty2D')
      STOP
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
    REAL(DP), DIMENSION(:), POINTER :: p_DelementVolume
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_ELEMENTIDX) :: isize
    REAL(DP) :: dtotalVolume
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE2D) :: Dpoints

    ! Is everything here we need?
    IF (rtriangulation%h_DcornerCoordinates .EQ. ST_NOHANDLE) THEN
      CALL output_line ('h_DcornerCoordinates not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      STOP
    END IF

    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement  not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementVolume2D')
      STOP
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
    CALL storage_getbase_double2D (rtriangulation%h_DcornerCoordinates,&
        p_DcornerCoordinates)
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    CALL storage_getbase_double (rtriangulation%h_DelementVolume,&
        p_DelementVolume)
        
    dtotalVolume = 0.0_DP
        
    ! Currently, we support triangules and quads.
    IF (UBOUND(p_IverticesAtElement,1) .EQ. TRIA_NVETRIA2D) THEN

      ! Calculate the element volume for all elements
      DO iel=1,rtriangulation%NEL
        ! triangular element
        Dpoints(1:NDIM2D,1:TRIA_NVETRIA2D) = &
            p_DcornerCoordinates(1:NDIM2D,p_IverticesAtElement(1:TRIA_NVETRIA2D,iel))
        p_DelementVolume(iel) = gaux_getArea_tria2D(Dpoints)
        
        dtotalVolume = dtotalVolume+p_DelementVolume(iel)
      END DO
    
    ELSE

      ! Calculate the element volume for all elements
      DO iel=1,rtriangulation%NEL
      
        IF (p_IverticesAtElement(4,iel) .EQ. 0) THEN
          ! triangular element
          Dpoints(1:NDIM2D,1:TRIA_NVETRIA2D) = &
              p_DcornerCoordinates(1:NDIM2D,p_IverticesAtElement(1:TRIA_NVETRIA2D,iel))
          p_DelementVolume(iel) = gaux_getArea_tria2D(Dpoints)
        ELSE
          ! quad element
          Dpoints(1:NDIM2D,1:TRIA_NVEQUAD2D) = &
              p_DcornerCoordinates(1:NDIM2D,p_IverticesAtElement(1:TRIA_NVEQUAD2D,iel))
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER :: ive, ibct, ivbd, iadjElement, nnve
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_POINTIDX) :: isize, ivt

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IelementsAtVertexIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IelementsAtVertexIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IelementsAtVertex .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IelementsAtVertex not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtBoundary2D')
      STOP
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
      CALL storage_new ('tria_genElementsAtBoundary2D', 'KMID', &
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
  ! IedgesAtBoundary (KEBD). 
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
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IedgesAtBoundary
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    INTEGER :: ive, ibct, ivbd, nnve
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_POINTIDX) :: isize, ivt

    ! Is everything here we need?
    IF (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IedgesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IelementsAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IelementsAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
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
      CALL storage_new ('tria_genEdgesAtBoundary2D', 'KMID', &
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
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genEdgeParameterValue2D(rtriangulation,rboundary)

!<description>
  ! This routine generates the parameter value of edge midpoints
  ! on the boundary DedgeParameterValue (DMBDP). Initialises NMBD.
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
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
    INTEGER(PREC_POINTIDX), DIMENSION(:), POINTER :: p_IboundaryCpIdx
    
    INTEGER :: ibct, ivbd
    INTEGER(PREC_POINTIDX) :: isize
    REAL(DP) :: dpar1,dpar2

    ! Is everything here we need?
    IF (rtriangulation%h_DvertexParameterValue .EQ. ST_NOHANDLE) THEN
      CALL output_line ('DvertexParameterValue not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IverticesAtBoundary .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IverticesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
    END IF

    IF (rtriangulation%h_IboundaryCpIdx .EQ. ST_NOHANDLE) THEN
      CALL output_line ('IboundaryCpIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary2D')
      STOP
    END IF

    ! Get the arrays.
    CALL storage_getbase_double (rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)
    CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    
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
    
      ! On each boundary component, loop through all vertices
      DO ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-2
      
        ! Get the parameter value of the current vertex and its neighbour.
        dpar1 = p_DvertexParameterValue(ivbd)
        dpar2 = p_DvertexParameterValue(ivbd)
        
        ! The edge parameter value is the mean.
        p_DedgeParameterValue(ivbd) = 0.5_DP*(dpar1+dpar2)
      
      END DO
      
      ! The 'last' vertex is a special case as there is no real neighbour.
      ! The parameter value of the 'neighbour' is the maximum parameter
      ! value of the boundary component + the parameter value of the
      ! first vertex of the boundary component.
      !
      ! Note that ivbd points to p_IboundaryCpIdx(ibct+1)-1 by
      ! Fortran standard!
      dpar1 = p_DvertexParameterValue(ivbd)
      dpar2 = p_DvertexParameterValue(p_IboundaryCpIdx(ibct)) + &
          boundary_dgetMaxParVal(rboundary, ibct)
      p_DedgeParameterValue(ivbd) = 0.5_DP*(dpar1+dpar2)
    
    END DO
    
    rtriangulation%NMBD = rtriangulation%NVBD
    
  END SUBROUTINE

END MODULE

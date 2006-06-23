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
!# 2.) tria_done
!#     -> Cleans up a triangulation structure, releases memory from the heap.
!#
!# 3.) tria_quadToTri
!#     -> Converts an arbitrary mesh into a triangular mesh.
!#
!# </purpose>
!##############################################################################

MODULE triangulation

  USE fsystem
  USE storage
  USE basicgeometry
  USE linearalgebra

  IMPLICIT NONE
  
  INCLUDE 'stria.inc'

!<constants>

!<constantblock description="Triangulation constants">
  
  ! Maximum number of corner-vertices in each element.
  ! We set this to 4 to allow triangle and quadrilateral polygons.
  ! May be set to higher vales in the future for the support of
  ! isoparametric elements!
  ! This is the old NNVE.
  INTEGER, PARAMETER :: TRIA_MAXNVE2D = 4

  ! Maximum number of edges in each element.
  ! We set this to 4 to allow triangle and quadrilateral polygons.
  ! This is the old NNVE, too.
  INTEGER, PARAMETER :: TRIA_MAXNME2D = TRIA_MAXNVE2D
  
!</constantblock>

!<constantblock description="KIND values for triangulation data">
  
  ! kind value for indexing the points in a triangulation
  INTEGER, PARAMETER :: PREC_POINTIDX   = I32

  ! kind value for indexing the edges in a triangulation
  INTEGER, PARAMETER :: PREC_EDGEIDX    = I32

  ! kind value for indexing the elements in a triangulation
  INTEGER, PARAMETER :: PREC_ELEMENTIDX = I32

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
    ! e.g. 1 if midpoints on edges exist in p_RfreeVertexCoordinates.
    INTEGER             :: nverticesPerEdge = 0
    
    ! Total number of vertices on edges; normally = 0.
    ! Total number of vertices on all edges, realized in p_RfreeVertexCoordinates. 
    ! E.g. if midpoints on edges exist, there is NVEDT=NMT.
    INTEGER(PREC_POINTIDX)   :: nVerticesOnAllEdges = 0
    
    ! Number of inner-element vertices; normally = 0.
    ! If a regular distribution of vertices in the inner of 
    ! each element is given, NIELV saves the number of vertices 
    ! in the inner of each element; e.g. 1 if element-midpoints 
    ! exist in p_RfreeVertexCoordinates.
    INTEGER             :: nverticesInEachElement = 0

    ! Total number of vertices in elements; normally = 0.
    ! Total number of vertices on all elements, realized in p_RfreeVertexCoordinates. 
    ! E.g. if element-midpoints exist in DCORMG, there is NIEVT=NEL.
    INTEGER(PREC_POINTIDX)   :: nverticesInAllElements = 0
    
    ! Number of additional vertices; normally = 0.
    ! Can be set <> 0 if there are any additional vertices 
    ! realized in p_RfreeVertexCoordinates, that don't belong to a regular 
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
    ! p_RneighbourElement(IEL)%Ineighbours(.) describes the elements adjacent 
    ! to IEL along the edges (p_RedgesOnElement(IEL)%Iedges(.)-NVT).
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
    
    ! Array with area of each element. 
    ! Handle to 
    !       p_DelementArea = array [1..NEL+1] of double.
    ! p_DelementArea [NEL+1] gives the total area of the domain.
    INTEGER         :: h_DelementArea = ST_NOHANDLE
    
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
    !       p_RfreeVertexCoordinates = array [1..NDIM2D,1..NVT] of double
    ! Array containing the coordinates of all vertices on edges,
    ! inner element vertices and additional nodes in the geometry.
    !
    ! p_RfreeVertexCoordinates(1..nVerticesOnAllEdges) contains the coordinates 
    ! of the regular distributed vertices on edges. 
    !
    ! p_RfreeVertexCoordinates(nVerticesOnAllEdges + 1 ..
    !                          nVerticesOnAllEdges + nVerticesOnAllElements) 
    ! contains the coordinates of the regular distributed 
    ! inner-element vertices. 
    !
    ! p_RfreeVertexCoordinates(nVerticesOnAllEdges + nVerticesOnAllElements + 1 ..
    !                          nVerticesOnAllEdges + nVerticesOnAllElements + nadditionalVertices)
    ! contains the coordinates of any additional vertices that do not 
    ! belong to regular distributed vertices on edges or on elements.
    !
    ! This is an extended version of the old DCORMG array.
    ! The DCORMG-array is originally designed to collect the midpoints 
    ! that are associated to the edges. In the new style, this behaviour 
    ! is only a special case that happens if NVEDT=NMT is set.
    ! The new style allowS p_RfreeVertexCoordinates to save all vertex 
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
    INTEGER           :: h_RfreeVertexCoordinates = ST_NOHANDLE
    
  END TYPE

!<typeblock>
  
!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_wrp_tria2Structure (TRIA, p_rtriangulation)
  
  USE afc_util
  
!<description>
  ! Wrapper routine. Accepts an 'old' triangulation structure array of CC2D
  ! and converts it completely to a triangulation structure. All 'old'
  ! information in the triangulation structure is overwritten.
  ! If p_rtriangulation is NULL(), a new structure will be created. 
  ! Otherwise, the existing structure is recreated/updated.
!</description>
  
!<input>
  ! The old triangulation structure array that should be converted to
  ! the new triangulation structure.
  INTEGER, DIMENSION(SZTRIA), INTENT(IN) :: TRIA
!</input>
  
!<inputoutput>
  ! The triangulation structure which will be overwritten by the information
  ! in TRIA.
  TYPE(t_triangulation), POINTER      :: p_rtriangulation
!</inputoutput>
  
!</subroutine>

  ! Do we have a structure?
  IF (.NOT. ASSOCIATED(p_rtriangulation)) THEN
    ALLOCATE(p_rtriangulation)
  ELSE
    ! Release the old structure without removing it from the heap.
    CALL tria_done (p_rtriangulation,.TRUE.)
  END IF

  ! We take the ownership of the triangulation structure - for compatibility
  ! to old FEAT 1.x:
  p_rtriangulation%Itria = TRIA
  
  ! Copy static entries
  p_rtriangulation%NVT                      = TRIA(ONVT  )
  p_rtriangulation%NMT                      = TRIA(ONMT  )
  p_rtriangulation%NEL                      = TRIA(ONEL  )
  p_rtriangulation%NBCT                     = TRIA(ONBCT )
  p_rtriangulation%NVBD                     = TRIA(ONVBD )
  p_rtriangulation%NMBD                     = TRIA(ONVBD ) ! NMBD=NVBD !!!
  p_rtriangulation%nverticesPerEdge         = TRIA(ONVPED)
  p_rtriangulation%nVerticesOnAllEdges      = TRIA(ONVEDT)
  p_rtriangulation%nverticesInEachElement   = TRIA(ONVEDT)
  p_rtriangulation%nverticesInAllElements   = TRIA(ONIEVT)
  p_rtriangulation%nadditionalVertices      = TRIA(ONANT )
  
  ! The duplication flag stays at 0 - as for now, all
  ! arrays are created from the feat arrays as new arrays, and
  ! so they are not a copy of another array.
  
  p_rtriangulation%iduplicationFlag = 0
  
  ! *******************************************************
  ! Copy DCORVG, create p_RcornerCoordinates.
  
  CALL copy_featarray_double2d ('DCORVG',2,INT(p_rtriangulation%NVT),TRIA(OLCORVG),&
                                p_rtriangulation%h_DcornerCoordinates)
  
  ! *******************************************************
  ! Copy KVERT, create p_RverticesAtElement.
  
  CALL copy_featarray_int2d ('KVERT',4,INT(p_rtriangulation%NEL),TRIA(OLVERT),&
                             p_rtriangulation%h_IverticesAtElement)

  ! *******************************************************
  ! Copy KMID, create p_RedgesAtElement.
  
  CALL copy_featarray_int2d ('KMID',4,INT(p_rtriangulation%NEL),TRIA(OLMID),&
                             p_rtriangulation%h_IedgesAtElement)

  ! *******************************************************
  ! Copy KADJ, create p_RneighboursAtElement.
  
  CALL copy_featarray_int2d ('KADJ',4,INT(p_rtriangulation%NEL),TRIA(OLADJ),&
                             p_rtriangulation%h_IneighboursAtElement)

  ! *******************************************************
  ! Copy KMEL, create p_IelementsAtEdge.
  
  CALL copy_featarray_int2d ('KMEL',2,INT(p_rtriangulation%NMT),TRIA(OLMEL),&
                             p_rtriangulation%h_IelementsAtEdge)

  ! *******************************************************
  ! Copy KEAN, create p_IverticesAtEdge.
  
  CALL copy_featarray_int2d ('KEAN',2,INT(p_rtriangulation%NMT),TRIA(OLEAN),&
                             p_rtriangulation%h_IverticesAtEdge)

  ! *******************************************************
  ! Copy KNPR, create p_InodalProperty.
  
  CALL copy_featarray_int1d ('KNPR',INT(p_rtriangulation%NVT+p_rtriangulation%NMT),&
                             TRIA(OLNPR),p_rtriangulation%h_InodalProperty)
                             
  ! *******************************************************
  ! Copy KAREA, create hpDelementArea.
  
  CALL copy_featarray_double1d ('KAREA',INT(p_rtriangulation%NEL+1),&
                                TRIA(OLAREA),p_rtriangulation%h_DelementArea)
                             
  ! *******************************************************
  ! Initialise the new KVEL, create p_IelementsAtVertexIdx/p_IelementsAtVertex.
                             
  CALL translate_KVEL (TRIA(ONVEL),INT(p_rtriangulation%NVT),TRIA(OLVEL), &
                       p_rtriangulation%h_IelementsAtVertex,&
                       p_rtriangulation%h_IelementsAtVertexIdx)

  ! *******************************************************
  ! Copy KBCT, create p_IboundaryCpIdx.
  
  CALL copy_featarray_int1d ('KBCT',(p_rtriangulation%NBCT+1),&
                             TRIA(OLBCT),p_rtriangulation%h_IboundaryCpIdx)
  
  ! *******************************************************
  ! Copy KVBD, create p_IverticesAtBoundary.
  
  CALL copy_featarray_int1d ('KVBD',INT(p_rtriangulation%NVBD),&
                             TRIA(OLVBD),p_rtriangulation%h_IverticesAtBoundary)

  ! *******************************************************
  ! Copy KMBD, create p_IedgesAtBoundary.
  
  CALL copy_featarray_int1d ('KMBD',(p_rtriangulation%NMBD),&
                             TRIA(OLMBD),p_rtriangulation%h_IedgesAtBoundary)

  ! *******************************************************
  ! Copy KEBD, create p_IelementsAtBoundary.
  
  CALL copy_featarray_int1d ('KEBD',INT(p_rtriangulation%NMBD),&
                             TRIA(OLEBD),p_rtriangulation%h_IelementsAtBoundary)

  ! *******************************************************
  ! Copy DVBDP, create p_DvertexParameterValue.
  
  CALL copy_featarray_double1d ('DVBDP',INT(p_rtriangulation%NVBD),&
                                TRIA(OLVBDP),p_rtriangulation%h_DvertexParameterValue)

  ! *******************************************************
  ! Copy DMBDP, create p_DedgeParameterValue.
  
  CALL copy_featarray_double1d ('DMBDP',INT(p_rtriangulation%NMBD),&
                                TRIA(OLMBDP),p_rtriangulation%h_DedgeParameterValue)


  ! *******************************************************
  ! Copy DCORMG, create p_RfreeVertexCoordinates.
  
  CALL copy_featarray_double2d ('DCORMG',2,INT(p_rtriangulation%NMT),TRIA(OLCORMG),&
                                p_rtriangulation%h_RfreeVertexCoordinates)
  

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
                          idim1, ST_INT, ihandle, &
                          ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the array
      CALL storage_new ('tria_wrp_tria2Structure', name, &
                        idim1,ST_INT, ihandle, &
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
                          idim1, ST_DOUBLE, ihandle, &
                          ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the array
      CALL storage_new ('tria_wrp_tria2Structure', name, &
                        idim1,ST_DOUBLE, ihandle, &
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
                          nentries, ST_INT, ihandle, &
                          ST_NEWBLOCK_NOINIT)
        CALL storage_new ('tria_wrp_tria2Structure', 'KADJIDX', &
                          NVT+1, ST_INT, ihandleidx, &
                          ST_NEWBLOCK_NOINIT)
      END IF
    ELSE
      ! Allocate the new array
      IF (ihandleidx .NE. 0) CALL storage_free (ihandleidx)
      CALL storage_new ('tria_wrp_tria2Structure', 'KADJ', &
                        nentries, ST_INT, ihandle, &
                        ST_NEWBLOCK_NOINIT)
      CALL storage_new ('tria_wrp_tria2Structure', 'KADJIDX', &
                        NVT+1, ST_INT, ihandleidx, &
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

  SUBROUTINE tria_done (p_rtriangulation,bkeepStructure)
  
!<description>
  ! This routine cleans up a triangulation structure.
  ! All memory allocated by handles in the structure is released from the heap.
  ! The routine does not clean up old FEAT 1.x legacy arrays in the
  ! rtriangulation%Itria substructure!
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the structure rtriangulation itself is not 
  ! released from memory. If set to FALSE or not existent (the usual setting), 
  ! the structure p_rboundary will also be removed from the heap after 
  ! cleaning up.
  LOGICAL, INTENT(IN), OPTIONAL :: bkeepStructure
!</input>

!<inputoutput>
  ! The triangulation structure to be cleaned up.
  TYPE(t_triangulation), POINTER :: p_rtriangulation
!</inputoutput>
  
!</subroutine>

    INTEGER :: idupflag
    
    IF (.NOT. ASSOCIATED(p_rtriangulation)) RETURN
    
    idupflag = p_rtriangulation%iduplicationFlag
    
    ! Bit  8: KMM    is a copy of another structure
    ! ... does not exist!?!

    ! Just release all allocated handles.
    ! Take care of which handles are duplicates from other structures - 
    ! these must not be released, as we are not the owner of them!
    
    ! Bit  0: DCORVG is a copy of another structure
    CALL checkAndRelease(idupflag, 0,p_rtriangulation%h_DcornerCoordinates)

    ! Bit  2: KVERT  is a copy of another structure
    CALL checkAndRelease(idupflag, 2,p_rtriangulation%h_IverticesAtElement)

    ! Bit  3: KMID   is a copy of another structure
    CALL checkAndRelease(idupflag, 3,p_rtriangulation%h_IedgesAtElement)
    
    ! Bit  4: KADJ   is a copy of another structure
    CALL checkAndRelease(idupflag, 4,p_rtriangulation%h_IneighboursAtElement)

    ! Bit  6: KMEL   is a copy of another structure
    CALL checkAndRelease(idupflag, 6,p_rtriangulation%h_IelementsAtEdge)

    ! Bit 16: KEAN   is a copy of another structure
    CALL checkAndRelease(idupflag,16,p_rtriangulation%h_IverticesAtEdge)

    ! Bit  7: KNPR   is a copy of another structure
    CALL checkAndRelease(idupflag, 7,p_rtriangulation%h_InodalProperty)

    ! Bit 19: DAREA  is a copy of another structure
    CALL checkAndRelease(idupflag,19,p_rtriangulation%h_DelementArea)

    ! Bit  5: KVEL   is a copy of another structure
    CALL checkAndRelease(idupflag, 5,p_rtriangulation%h_IelementsAtVertexIdx)
    CALL checkAndRelease(idupflag, 5,p_rtriangulation%h_IelementsAtVertex)

    ! Bit 11: KBCT   is a copy of another structure
    CALL checkAndRelease(idupflag,11,p_rtriangulation%h_IboundaryCpIdx)

    ! Bit  9: KVBD   is a copy of another structure
    CALL checkAndRelease(idupflag, 9,p_rtriangulation%h_IverticesAtBoundary)

    ! Bit 14: KMBD   is a copy of another structure
    CALL checkAndRelease(idupflag,14,p_rtriangulation%h_IedgesAtBoundary)

    ! Bit 10: KEBD   is a copy of another structure
    CALL checkAndRelease(idupflag,10,p_rtriangulation%h_IelementsAtBoundary)

    ! Bit 12: DVBDP  is a copy of another structure
    CALL checkAndRelease(idupflag,12,p_rtriangulation%h_DvertexParameterValue)

    ! Bit 13: DMBDP  is a copy of another structure
    CALL checkAndRelease(idupflag,13,p_rtriangulation%h_DedgeParameterValue)

    ! Bit 17: KVBDI  is a copy of another structure
    CALL checkAndRelease(idupflag,17,p_rtriangulation%h_IboundaryVertexPos)

    ! Bit 18: KMBDI  is a copy of another structure
    CALL checkAndRelease(idupflag,18,p_rtriangulation%p_IboundaryEdgePos)
    
    ! Bit  1: DCORMG is a copy of another structure
    CALL checkAndRelease(idupflag, 1,p_rtriangulation%h_RfreeVertexCoordinates)
    
    ! Clean up the rest of the structure

    p_rtriangulation%iduplicationFlag = 0
    p_rtriangulation%NVT = 0
    p_rtriangulation%NMT = 0
    p_rtriangulation%NEL = 0
    p_rtriangulation%NBCT = 0
    p_rtriangulation%NVBD = 0
    p_rtriangulation%NMBD = 0
    p_rtriangulation%nverticesPerEdge = 0
    p_rtriangulation%nVerticesOnAllEdges = 0
    p_rtriangulation%nverticesInEachElement = 0
    p_rtriangulation%nverticesInAllElements = 0
    p_rtriangulation%nadditionalVertices = 0

    ! Deallocate the structure (if we are allowed to), finish.
    IF (.NOT. PRESENT(bkeepStructure)) THEN
      DEALLOCATE(p_rtriangulation)
    ELSE
      IF (.NOT. bkeepStructure) DEALLOCATE(p_rtriangulation)
    END IF

    ! That's it...

  CONTAINS
  
    ! **********************************************************
    ! Release handle ihandle if bit ibit in idubFlag is not set.
    ! Otherwise, ihandle is set to ST_NOHANDLE.    
    SUBROUTINE checkAndRelease (idupFlag,ibit,ihandle)
    
    INTEGER, INTENT(IN) :: ibit
    INTEGER(I32), INTENT(IN) :: idupFlag
    INTEGER, INTENT(INOUT) :: ihandle
    
      IF (IAND(idupFlag,2**ibit) .EQ. 0) THEN
        IF (ihandle .NE. ST_NOHANDLE) CALL storage_free(ihandle)
      ELSE
        ihandle = ST_NOHANDLE
      END IF
      
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>    

  SUBROUTINE tria_quadToTri (rtriangulation)
  
!<description>
  ! This routine converts a 2D mesh into a triangular 2D mesh. All quads are 
  ! converted to triangles.
  ! Warning: This should be applied only to the coarse grid; any information
  !  about possible two-level ordering is lost!
!</description>

!<inputoutput>
  ! The triangulation structure to be converted. Is replaced by a triangular
  ! mesh.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>

    INTEGER(PREC_ELEMENTIDX) :: i,icount
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER :: h_IverticesAtElementTri
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementTri
    INTEGER, DIMENSION(2) :: Isize
   
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
    INTEGER,INTENT(IN) :: nel

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
  
!</subroutine>
  
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

END MODULE

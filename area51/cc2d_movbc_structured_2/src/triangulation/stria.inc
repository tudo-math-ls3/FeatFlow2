!***********************************************************************
! This file contains the indices for the STRIA-structure that describes
! a triangulation. The structure has the following form:
!
! STRIA = record
!   // Status flags
!
!   TRIFL   : integer
!   IDPFLG  : integer
!   IPRMFL  : integer
!
!   // Dimension flags
! 
!   NDIM    : integer
!   NEL     : integer
!   NVT     : integer
!   ...
!   NANT    : integer
! 
!   // Handles
!  
!   LCORVG  : integer
!   LCORMG  : integer
!   ...
!   ILMBDP  : integer
!
!   // User-defined variables
!
!   TRIUD   : array [1..32] of integer
! end
!
! It's realised as a simple integer array:
!
!   TRIA = array [1..SZTRIA] of integer
!
! The "index/offset" constants defined in this file correspond to the
! position of the desired information in the array, i.e. for example:
!
!   NVE = TRIA [ONVE]
!
! P.S.: Some elements in the array are not used, they are reserved
!       for future use. Furthermore we don't care about the fact that
!       the basic F77 compiler only uses 6 letters. Some "names"
!       have more than 6 letters but those compilers using only the
!       first 6 should ignore the remaining without danger of
!       confound two index names - as all index names differ in the
!       first 6 letters.
!
! A grid consists of a number of (corner-)vertices, a number of edges
! connecting these vertices and a number of elements, each surrounded
! by vertices and corners. Optionally the grid can have a number
! of edge midpoints, each associated with an edge.
!
! Remark: The grid is totally independent of any discretization!
!  It's a pure coincidence if grid vertices/edges correspond to
!  any quantities of Finite Element spaces!
!
! The following arrays and quantities characterize a grid:
!
! Constants:
! ----------
! NNVE   - Maximum number of vertices for each element
!   Defined in cbasictria.inc. Contains the maximum
!   number of corner-vertices per element - 3 or 4, depending on whether
!   the element is triangle or quadrilateral. Gives at the same
!   time the maximum number of edges on an element, as each node
!   corresponds to one edge adjacent to it.
!
! Structure-specific variables
! ----------------------------
! TRIFL   - triangulation type flag.
!            =0: triangulation is in old FEAT style
!           <>0: triangulation is in new, extended style
!                In this case TRIFL is a user-defined flag
!                that can identify which additional information is saved
!                in the structure. TRIFL should them be used as a
!                bitfield with Bit0=1 to indicate extended style.
! IDPFLG  - integer
!           Duplication flag. Bitfield. Used by TRIDUP/TRIRST/TRIDEL to 
!           mark which information of a triangulation structure 
!           coincides with those of another triangulation structure and
!           must not be released from memory when deleting a copy of a
!           triangulation structure.
!           When a bit is set to 1, the corresponding array is
!           maintained by another triangulation structure and must not
!           be deleted by TRIDEL. When the bit is 0, the array is a real
!           copy of another array and must be deleted in TRIDEL.
!           Bit  0: Handle of DCORVG is a copy of another structure
!           Bit  1: Handle of DCORMG is a copy of another structure
!           Bit  2: Handle of KVERT  is a copy of another structure
!           Bit  3: Handle of KMID   is a copy of another structure
!           Bit  4: Handle of KADJ   is a copy of another structure
!           Bit  5: Handle of KVEL   is a copy of another structure
!           Bit  6: Handle of KMEL   is a copy of another structure
!           Bit  7: Handle of KNPR   is a copy of another structure
!           Bit  8: Handle of KMM    is a copy of another structure
!           Bit  9: Handle of KVBD   is a copy of another structure
!           Bit 10: Handle of KEBD   is a copy of another structure
!           Bit 11: Handle of KBCT   is a copy of another structure
!           Bit 12: Handle of DVBDP  is a copy of another structure
!           Bit 13: Handle of DMBDP  is a copy of another structure
!           Bit 14: Handle of KMBD   is a copy of another structure
!           Bit 15: Handle of KXNPR  is a copy of another structure
!           Bit 16: Handle of KEAN   is a copy of another structure
!           Bit 17: Handle of KVBDI  is a copy of another structure
!           Bit 18: Handle of KMBDI  is a copy of another structure
!           Bit 19: Handle of DAREA  is a copy of another structure
!           Bit 20: Handle of KFBVT  is a copy of another structure
!           Bit 21: Handle of KFBMT  is a copy of another structure
! IPRMFL  - Parametrization flag.
!           This is a (user defined) identifier of a parametrization that
!           allowes the user to choose the correct parametrization, this
!           structure refers to.
!           (At the moment this is not used and should be kept = 0!
!            This corresponds to the standard parametrization).
!
! Mandatory Dimensions:
! ---------------------
! NDIM    - Number of dimensions in the triangulation structure;
!           usually =2 or =3
! NEL     - Number of elements
! NVT     - Number of vertices, i.e. corners of elements
! NMT     - Number of edges
! NVE     - Number of (corner-)vertices per element
! NVEL    - maximum number of elements meeting in a (corner-)vertex
! NBCT    - Number of (real) boundary components the triangulation
!           structure corresponds to
! NVBD    - Total number of corner vertices on the (real) boundary
!
! Mandatory Arrays:
! -----------------
! The following arrays are mandatory as they explicitly describe the
! structure of the grid. The arrays marked with "*" are extended arrays
! not included in the original FEAT structure.
!
! DCORVG  - array [1..2,1..NVT] of double precision
!   Array with X/Y-coordinates of all (corner-)vertices in the
!   geometry. For a vertex IVT
!    - DCORVG(1,IVT) is the X-coordinate
!    - DCORVG(2,IVT) is the Y-coordinate
! KVERT   - array [1..NNVE,1..NEL] of integer
!   Gives for each element the node numbers of the corner-vertices
!   in mathematically positive sense.
! KMID    - array [1..NNVE,1..NEL] of integer
!   Gives for each element the numbers of the edges in mathematically
!   positive sense. KMID(I,IEL) is exactly the number of the edge
!   adjacent to the node KVERT(I,IEL) on element IEL. To be able to
!   distinguish a number of an edge from a vertex number, edges are
!   numbered in the range NVT+1..NVT+NMT. Some arrays below use
!   the "real" edge number; this can be obtained by subtracting NVT
!   from the value of KMID.
! KADJ    - array [1..NNVE,1..NEL] of integer
!   Contains for each element the numbers of adjacent elements
!   in mathematically positive sense. KADJ(I,IEL) is the element
!   adjacent to IEL on the edge (KMID(I,IEL)-NVT).
! KVEL    - array [1..NVEL,1..NVT] of integer
!   KVEL(.,IVT) contains the numbers of the elements meeting in a
!   vertex. If there are less than NVEL elements meeting in a vertex,
!   KVEL(x..NVEL,IVT) is filled up with 0.
! KMEL    - array [1..2,1..NMT] of integer
!   KMEL(.,IMT) contains the numbers of the two elements adjacent
!   to an edge IMT. For boundary edges, KMEL(2,IMT) is set to 0.
! KNPR    - array [1..NVT+NMT] of integer
!   Standard nodal property. KNPR(1,.) defines for each vertex 
!   (1..NVT) and each edge (NVT+1..NVT+NMT) its function inside of the
!   geometry. Generally said, the range of the KNPR-array 
!   characterizes the type of the node:
!   = 0    : The vertex/edge is an inner node
!   > 0    : For vertices: The vertex is a boundary node on the real
!              boundary. KNPR(.) defines the number of the boundary
!              component.
!            For edges: The edge is an edge on the real boundary.
!              KNPR(.) defines the number of the vertex preceding
!              the edge. Therefore, KNPR(KNPR(.)) defines the number
!              of the boundary component of that edge.
!   < 0,
!   >= -NEL: The vertex/edge is an invalid node to element -KNPR(.)
!            (-> hanging node, not implemented for now)
!   < -NEL : The vertex/edge is a fictitious boundary node.
!            -KNPR(.)-NEL defines the number of the fictitious
!            boundary component.
! KXNPR  *- array [1..2,1..NVT+NMT] of integer
!   Extended Nodal property. KXNPR(1,.) defines for each vertex 
!   (1..NVT) and each edge (NVT+1..NVT+NMT) its function inside of the
!   geometry. KXNPR(1,INODE) is a flag-bitfield + identifier for the
!   type of the node INODE, while KXNPR(2,INODE) contains a numerical
!   information tag depending on the bitfield.
!   KXNPR(1,INODE) is build as follows:
!   Bit0..Bit7: The first 8 bit give a value of 0..255 to characterize
!               the type of the node numerically. This is normally
!               not used (and thus set to 0), but can be defined
!               if a special node characterization is necessary.
!   Bit  8: =0: The node is an inner node according to the
!     (256)     parametrization/geometry.
!           =1: The node is a real boundary node according to the
!               parametrization/geometry.
!               KXNPR(2,I) = number of the boundary component that
!                            contains node I
!   Bit  9: =0: The node is a standard node in the mesh
!     (512) =1: The node is an invalid node to a neighbour element.
!               KXNPR(2,INODE) = number of the element, this node is
!                                invalid to.
!               Remarks: 1.) "Invalid node" normally means "hanging 
!                node". As hanging nodes are not implemented in 
!                FEAT, this bit is not used for now and only reserved
!                for future use.
!               2.) If the node is inside of a fictitious boundary
!                component, use the neighbour element to obtain
!                information about the fict. bdry. component.
!   Bit 10+11 : Reserved
!   The following bits 12+ are not directly specifying geometry
!   information these bits are used to indicate the type of the node
!   in the triangulation, so that the application knows e.g. which type
!   of boundary conditions to impose. When used at all, the bits are
!   typically used in the following manner:
!   Bit 12: =0: Usual setting. Treat the node as inner node or
!    (4096)     Neumann-type boundary node.
!           =1: Treat the node as Dirichlet node (inner or boundary).
!               The solution in this node is prescribed by a given
!               (analytic) function.
!   Bit 13: =0: The node is not lying in a fictitious boundary
!    (8192)     component. 
!           =1: The node belongs currently to a fictitious boundary
!               component.
!               If the node is not a real boundary node:
!                 KXNPR(2,I) = number of the fictitious boundary 
!                              component that contains node I
!               For boundary nodes inside of a fictitious boundary
!               component, KXNPR(2,I) still shows the number of the
!               real boundary component.
!   Bit 14..20,32: reserved (for future types of boundary conditions)
!   Bit 21..31: user defineable
! KEAN   *- array [1..2,1..NMT] of integer
!   Number of adjacent node to an edge. For every edge I, KEAN(1,I)
!   defines the number of the node preceding this edge and
!   KEAN(2,I) the number of the node following the edge.
!   Remark: I must be the "real" edge number, i.e. when taking the
!   value from KMID, NVT must be subtracted.
! 
! The following structures define the mandatory information
! describing the boundary:
!
! KVBD    - array [1..NVBD] of integer
!   Array containing the numbers of all corner vertices on the 
!   boundary. Takes values in the range 1..NVT.
! KVBDI  *- array [1..2,1..NVBD] of integer
!   Inverse index array to KVBD. KVBDI(1,.) contains a vertex 
!   number on the boundary and KVBDI(2,.) the appropriate index
!   of this vertex inside of the KVBD-array. KVBDI is sorted
!   for the vertex number, thus allowing quick access to the
!   index of a vertex in KVBD.
! KMBD    - array [1..NVBD] of integer
!   Array containing the numbers of all edges on the boundary.
!   For compatibility reasons this takes values in the range
!   NVT+1..NVT+NMT. There are as many edges on the boundary
!   as vertices, so NVBD = number of vertices on the boundary
!   = nummber of edges on the boundary.
! KMBDI   - array [1..2,1..NVBD] of integer
!   Inverse index array to KMBD. KMBDI(1,.) contains an edge
!   number on the boundary and KMBDI(2,.) the appropriate index
!   of this edge inside of the KMBD-array.
!   For compatibility reasons KMBDI(1,.) has a value in the range
!   NVT+1..NVT+NMT). 
!   KMBDI is sorted for the edge number, thus allowing quick
!   access to the index of an edge in KMBD.
! KEBD    - array [1..NVBD] of integer
!   Array containing the numbers of elements adjacent to the
!   boundary. The element KEBD(IVBD) contains the boundary
!   edge KMBD(IVBD) and the node KVBD(IVBD).
!   Remark that only those elements are listed here that have
!   at least two adjacent vertices on the boundary; elements with
!   only one vertex on the boundary (e.g. an apex) are not
!   maintained in this array.
! KBCT    - array [1..NBCT+1] of integer 
!   Index vector for KVBD/KMBD/DVBDP. For a (real) boundary 
!   component IBCT all corner nodes for that boundary component are 
!   saved in KVBD ( KBCT(IBCT)..KBCT(IBCT+1)-1 ).
!   All boundary edges of this boundary component are saved in 
!   KMBD ( KBCT(IBCT)..KBCT(IBCT+1)-1 ).
!   KBCT(NBCT+1) points to NVBD+1 for easier access to the last
!   boundary component.
! DVBDP   - array [1..NVBD] of double precision
!   DVBDP(I) contains the parameter value of boundary node I,
!   which corresponds to the corner vertex KVBD(I).
!
! Optional information
! --------------------
! The following information is somehow "optional" as they are not really
! necessary for the definition of the grid, but for easier access to
! different quantities. If they are not available, the corresponding
! number/handle should be set to 0.
!
! Optional Dimensions:
! --------------------
!
!   NVPED   - Number of vertices per edge; normally = 0.
!             If a regular distribution of vertices on edges is given,
!             NVPED saves the number of vertices on each edge;
!             e.g. 1 if midpoints on edges exist in DCORMG.
!   NVEDT   - Total number of vertices on edges; normally = 0.
!             Total number of vertices on all edges, realized in DCORMG. 
!             E.g. if midpoints on edges exist, there is NVEDT=NMT.
!   NVEDT   - Number of inner-element vertices; normally = 0.
!             If a regular distribution of vertices in the inner of 
!             each element is given, NIELV saves the number of vertices 
!             in the inner of each element; e.g. 1 if element-midpoints 
!             exist in DCORMG.
!   NIEVT   - Total number of vertices on edges; normally = 0.
!             Total number of vertices on all edges, realized in DCORMG. 
!             E.g. if element-midpoints exist in DCORMG, there is 
!             NIEVT=NEL.
!   NANT    - Number of additional nodes; normally = 0.
!             Can be set <> 0 if there are any additional vertices 
!             realized in DCORMG, that don't belong to a regular 
!             distribution of vertices in corners, on edges or on 
!             elements.
!   NFBVT   - Number of fictitious boundary vertices; normally = 0
!   NFBMT   - Number of fictitious boundary edges; normally = 0
!
! Optional Arrays:
! ----------------
!
! DMBDP   - array [1..NVBD] of double precision
!   DMBDP(.,IVBD) contains all parameter values of the boundary
!   midpoints.
! DCORMG  - array [1..2,1..NVEDT+NIEVT+NANT] of double precision
!   Array containing the X/Y-coordinates of all vertices on edges,
!   inner element vertices and additional nodes in the geometry.
!   DCORMG(1,.) contains the X-, DCORMG(2,.) the Y-coordinates.
!   DCORMG(.,1..NVEDT) contains the coordinates of the regular 
!   distributed vertices on edges. 
!   DCORMG(.,NVEDT+1..NVEDT+NIEVT) contains the coordinates of the 
!   regular distributed inner-element vertices. 
!   DCORMG(.,NVEDT+NIEVT+1..NVEDT+NIEVT+NANT) contains the coordinates 
!   of any additional vertices that do not belong to regular distributed
!   vertices on edges or on elements.
! DAREA   - array [1..NEL+1] of double
!   Array with area of each element. DAREA [NEL+1] gives the total area
!   of the domain.
!
! KFBVT   - array [1..NFBVT] of integer
!   Array with all fictitious boundary vertices.
!   Normally not used, can be build with GENFBA.
! KFBMT   - array [1..NFBMT] of integer
!   Array with all fictitious boundary edges. To be able to
!   distinguish a number of an edge from a vertex number, edges are
!   numbered in the range NVT+1..NVT+NMT.
!   Normally not used, can be build with GENFBA.
!
! Remark: User defined information
!  At index ITRIUD the structure STRIA provides a set of 32 integer
!  values that can be defined by the user. The content of these
!  variables are not maintained by the standard triangulation routines.
!  The user is responsible for handling its content.
!
! Remark: Fictitious boundary nodes
!  If the grid contains fictitious boundary vertices/edges, the
!  corresponding value in KXNPR must be set appropriately. Furthermore,
!  the vertex/edge can be added to the list KFBVT/KFBMT, which contains
!  a list of all fictitious boundary vertices/edges. These arrays can be
!  build from KXNPR with GENFBA if desired.
!
! Remark: Additional nodes
!  The DCORMG-array is originally designed to collect the midpoints 
!  that are associated to the edges. In the new style, this behaviour 
!  is only a special case that happens if NVEDT=NMT is set.
!  The new style allowes DCORMG to save all vertex coordinates that 
!  are not corner vertices. This includes regular distributed vertices
!  on edges (either midpoints or points at 1/3 and 2/3 of the edge, 
!  or...), on elements (either midpoints or the 4 Gauss-points, or...)
!  as well as additional vertices the user wants to be realized, that
!  don't belong to any of the previous two groups.
!
!  There is a numbering introduced to number all the nodes in DCORMG
!  as well as DCORVG when accessing them in auxiliary routines:
!  - 1..NVT 
!    represent all corner vertices.
!  - NVT+1..NVT+NVEDT 
!    represent all regular distributed vertices on edges.
!  - NVT+NVEDT+1..NVT+NVEDT+NIEVT 
!    represent all regular distributed vertices in the inner of elements
!  - NVT+NVEDT+NIEVT..NVT+NVEDT+NIEVT+NANT 
!    represent all additional nodes that don't belong to regular 
!    distributed vertices on edges or inside of the elements.
!  NVT+NVEDT+NIEVT+NANT is also the total number of nodes returned
!  by TNDCNT. There is a special handling of the node-numbering
!  if DCORMG does not exist, to maintain compatibility. In this case
!  TNDCNT returns NVT - as only NVT vertices exist in the triangulation
!  structure. Furthermore only in this case the following alternative
!  node numbering is used:
!  - 1..NVT 
!    represent all corner vertices.
!  - NVT+1..NVT+NMT
!    represent all edge-midpoits
!  - NVT+NMT+1..NVT+NMT+NEL
!    represent all element midpoints
!  Note that in this case the coordinates are not realized in DCORMG!
!  The vertices are only virtual, i.e. it's possible to access their
!  coordinates with NDE2XY, but that's all!
!
! The numbers of regularly distributed vertices on edges/elements can
! be directly calculated with an appropriate formula. E.g. let's assume,
! we have n regularly distributed vertices on each edge (excluding the
! starting/ending point). Then the corresponding vertices on 
! egde E (=1..NMT) have the numbers:
!       (NVT-1)+(NMT-1) + (E-1)*n + 1
!    .. (NVT-1)+(NMT-1) + (E-1)*n + n
! the formula for regular distributed vertices on elements is similar.
!
! A more detailed description about the node numbering can be found in
! structria.f.
!***********************************************************************

      INTEGER ONEL,ONVT,ONMT,ONVE,ONVEL,ONBCT,ONVBD,ONANT,OTRIFL,OLVBDI
      INTEGER ONVPED,ONVEDT,ONIELV,ONIEVT,OIPRMFL,ONDIM
      INTEGER OLCORVG,OLCORMG,OLVERT,OLMID,OLADJ,OLVEL,OLMEL,OLNPR,OLMM
      INTEGER OLVBD,OLEBD,OLBCT,OLVBDP,OLMBDP,OLMBD,OLXNPR,OLEAN
      INTEGER OLMBDI,OLAREA,OIDPFLG,OLFBVT,OLFBMT,ONFBVT,ONFBMT,OTRIUD
     
      INTEGER SZTRIA

      PARAMETER (OTRIFL  =  1)
      PARAMETER (OIDPFLG =  2)
      PARAMETER (OIPRMFL =  3)

      PARAMETER (ONDIM   =  5)
      PARAMETER (ONEL    =  6)
      PARAMETER (ONVT    =  7)
      PARAMETER (ONMT    =  8)
      PARAMETER (ONVE    =  9)
      PARAMETER (ONVEL   = 10)
      PARAMETER (ONBCT   = 11)
      PARAMETER (ONVBD   = 12)
      PARAMETER (ONVPED  = 13)
      PARAMETER (ONVEDT  = 14)
      PARAMETER (ONIELV  = 15)
      PARAMETER (ONIEVT  = 16)
      PARAMETER (ONANT   = 17)
      
      PARAMETER (ONFBVT  = 18)
      PARAMETER (ONFBMT  = 19)
    
! Indices 22..29 reserved for future use
    
      PARAMETER (OLCORVG =  30)
      PARAMETER (OLCORMG =  31)
      PARAMETER (OLVERT  =  32)
      PARAMETER (OLMID   =  33)
      PARAMETER (OLADJ   =  34)
      PARAMETER (OLVEL   =  35)
      PARAMETER (OLMEL   =  36)
      PARAMETER (OLNPR   =  37)
      PARAMETER (OLMM    =  38)
      PARAMETER (OLVBD   =  39)
      PARAMETER (OLEBD   =  40)
      PARAMETER (OLBCT   =  41)
      PARAMETER (OLVBDP  =  42)
      PARAMETER (OLMBDP  =  43)
      PARAMETER (OLMBD   =  44)
      PARAMETER (OLXNPR  =  45)
      PARAMETER (OLEAN   =  46)
      PARAMETER (OLVBDI  =  47)
      PARAMETER (OLMBDI  =  48)
      PARAMETER (OLAREA  =  49)
      
      PARAMETER (OLFBVT  =  50)
      PARAMETER (OLFBMT  =  51)
      
! Indices 54..64 reserved for future use
! Indices 65..96 user-defineable.

      PARAMETER (OTRIUD  =  65)
      
      PARAMETER (SZTRIA  =  96)


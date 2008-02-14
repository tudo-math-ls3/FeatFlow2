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
!#  1.) tria_wrp_tria2Structure
!#      -> Wrapper. Create a FEAT 2.0 triangulation structure from a
!#         FEAT 1.x triangulation structure.
!#
!#  2.) tria_readTriFile2D
!#      -> Reads a .TRI file and creates a 'raw' 2D mesh with only basic 
!#         information.
!#
!#  3.) tria_readTriFile3D
!#      -> Reads a .TRI file and creates a 'raw' 3D mesh with only basic 
!#         information.
!#
!#  4.) tria_initStandardMeshFromRaw
!#      -> Generates all standard arrays for a mesh, i.e. converts a 'raw' mesh
!#         (as set up by tria_readTriFile2D e.g.) to a standard mesh.
!#
!#  5.) tria_refine2LevelOrdering
!#      -> Refines a mesh according to the 2-level ordering algorithm.
!#         Creates a 'raw' fine mesh from a 'standard' coarse mesh.
!#
!#  6.) tria_compress2LevelOrdHierarchy
!#      -> Can be used to compress a mesh hierarchy created with the 2-level
!#         ordering. Saves memory. Shares vertex coordinates between fine
!#         and coarse mesh.
!#
!#  7.) tria_quickRefine2LevelOrdering
!#      -> Refines a mesh multiple times according to the 2-level ordering 
!#         algorithm. Creates a 'raw' fine mesh from a 'raw' or 'standard' 
!#         coarse mesh.
!#
!#  8.) tria_done
!#      -> Cleans up a triangulation structure, releases memory from the heap.
!#
!#  9.) tria_rawGridToTri
!#      -> Converts a raw mesh into a triangular mesh.
!#
!# 10.) tria_duplicate
!#      -> Creates a duplicate / backup of a triangulation.
!#         Some information may be shared between two triangulation structures.
!#
!# 11.) tria_restore
!#      -> Restores a triangulation previously backed up with tria_backup.
!#
!# 12.) tria_searchBoundaryNode
!#      -> Search for the position of a boundary vertex / edge on the boundary.
!#
!# 13.) tria_getPointsOnEdge
!#      -> For all edges in a triangulation, calculate the coordinates 
!#         of a number of points on each edge
!#
!# 14.) tria_getNVE = tria_getNVE_direct / tria_getNVE_indirect
!#      -> Get the number of vertices/edges on an element
!#
!# 15.) tria_createRawTria1D
!#      -> Creates a 'raw' 1D triangulation $[a,b]$ with $n$ sub-intervals 
!#         of the same length
!#
!# 16.) tria_infoStatistics
!#      -> Prints out statistics about a mesh.
!#
!# 17.) tria_exportTriFile
!#      -> Exports a triangulation structure to a .TRI file.
!#
!# 18.) tria_getNeighbourVertex
!#      -> Calculates the vertex number of the neighbour vertex of a 
!#         vertex on an edge.
!#
!# Auxiliary routines:
!#
!#  1.) tria_readRawTriangulation1D / tria_readRawTriangulation2D
!#      -> Reads basic data from a TRI file.
!#
!#  2.) tria_genRawBoundary2D
!#     -> Generates basic information about boundary vertices.
!#
!#  3.) tria_sortBoundaryVertices2D
!#      -> Sorts the boundary vertices for increasing parameter values.
!#
!#  4.) tria_genElementsAtVertex2D
!#      -> Generates the IelementsAtVertex array for a 2D triangulation
!#
!#  5.) tria_genNeighboursAtElement1D
!#      -> Generates the IneighboursAtElement array for a 1D triangulation
!#
!#  6.) tria_genNeighboursAtElement2D
!#      -> Generates the IneighboursAtElement array for a 2D triangulation
!#
!#  7.) tria_genEdgesAtElement2D
!#      -> Generates the IedgesAtElement array for a 2D triangulation
!#
!#  8.) tria_genElementsAtEdge2D
!#      -> Generates the IelementsAtEdge array for a 2D triangulation
!#
!#  9.) tria_genVerticesAtEdge2D
!#      -> Generates the IverticesAtEdge array for a 2D triangulation
!#
!# 10.) tria_genEdgeNodalProperty2D
!#      -> Generates the InodalProperty array part for all edges
!#
!# 11.) tria_genElementVolume1D
!#      -> Generates the DelementVolume array for a 1D triangulation
!#
!# 12.) tria_genElementVolume2D
!#      -> Generates the DelementVolume array for a 2D triangulation
!#
!# 13.) tria_genElementsAtBoundary2D
!#      -> Generates the IelementsAtBoundary array for a 2D triangulation
!#
!# 14.) tria_genEdgesAtBoundary2D
!#      -> Generates the IedgesAtBoundary array for a 2D triangulation
!#
!# 15.) tria_genEdgeParameterValue2D
!#      -> Generates the DedgeParameterValue array for a 2D triangulatioon
!# 
!# 16.) tria_genBoundaryVertexPos2D
!#      -> Generates the IboundaryVertexPos2D array for a 2D triangulatioon
!#
!# 17.) tria_genBoundaryEdgePos2D
!#      -> Generates the IboundaryEdgePos2D array for a 2D triangulatioon
!#
!# 18.) tria_readRawTriangulation3d
!#      -> Read a 3D .TRI file
!#
!# 19.) tria_genRawBoundary3d
!#      -> Generate basic boundary information in 3D
!#
!# 20.) tria_refineMesh2lv3D
!#      -> Refine a 3D mesh with 2-level ordering
!#
!# 21.) tria_refineBdry2lv3D
!#      -> Refine the 3D boundary
!#
!# 22.) tria_genFacesAtBoundary
!#      -> Generate the array with face numbers on the boundary
!#
!# 23.) tria_genEdgeNodalProperty3d
!#      -> Generate the edge nodal property in 3D
!#
!# 24.) tria_genFaceNodalProperty3d
!#      -> Generate the face nodal property in 3D
!#
!# 25.) tria_genFacesAtVertex
!#      -> Generate the arrays with the faces on the boundary
!#
!# 26.) tria_genFacesAtEdge
!#      -> Generate the arrays with the faces at an edge
!#
!# 27.) tria_genEdgesAtFace
!#      -> Generate the arrays with the edges at a face
!#
!# 28.) tria_genElementsAtFace
!#      -> Generate the arrays with the elements at a face
!#
!# 29.) tria_genVerticesAtFace
!#      -> Generate the arrays with the vertices at a face
!#
!# 30.) tria_genFacesAtElement
!#      -> Generate the arrays with the faces at an element
!#
!# 31.) tria_genVerticesAtEdge3D
!#      -> Generate the arrays describing the vertices at an edge in 3D
!#
!# 32.) tria_genElementsAtEdge3D
!#      -> Generate the arrays describing the elements at an edge in 3D
!#
!# 33.) tria_genElementsAtVertex3D
!#      -> Generate the 'elements at vertex' arrays
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
  
  ! Maximum number of corner-vertices in each 1D element.
  INTEGER, PARAMETER :: TRIA_MAXNVE1D = 2

  ! Maximum number of corner-vertices in each 2D element.
  ! We set this to 4 to allow triangle and quadrilateral element shapes.
  ! May be set to higher vales in the future for the support of
  ! isoparametric elements!
  ! This is the old NNVE.
  INTEGER, PARAMETER :: TRIA_MAXNVE2D = 4
  
  ! Maximum number of edges in each 2D element.
  ! We set this to 4 to allow triangle and quadrilateral element shapes.
  ! This is the old NNVE, too.
  INTEGER, PARAMETER :: TRIA_MAXNME2D = TRIA_MAXNVE2D
  
  ! Number of vertices per line element in 1D.
  INTEGER, PARAMETER :: TRIA_NVELINE1D = 2
  
  ! Number of vertices per element for triangular element shapes in 2D.
  INTEGER, PARAMETER :: TRIA_NVETRI2D  = 3

  ! Number of vertices per element for quadrilateral element shapes in 2D.
  INTEGER, PARAMETER :: TRIA_NVEQUAD2D = 4
  
  ! number of elements of a 3d connector
  integer, parameter :: TRIA_NCONNECT3D = 5
  
  
!</constantblock>

!<constantblock description="KIND values for triangulation data">
  
  ! kind value for indexing the vertices in a triangulation
  INTEGER, PARAMETER :: PREC_VERTEXIDX  = I32

  ! Alternative name for PREC_VERTEXIDX
  INTEGER, PARAMETER :: PREC_POINTIDX   = PREC_VERTEXIDX

  ! kind value for indexing the edges in a triangulation
  INTEGER, PARAMETER :: PREC_EDGEIDX    = I32
  
  ! kind value for indexing the faces in a triangulation
  INTEGER, PARAMETER :: PREC_FACEIDX    = I32

  ! kind value for indexing the elements in a triangulation
  INTEGER, PARAMETER :: PREC_ELEMENTIDX = I32

!</constantblock>
  
!<constantblock description="Flags to be specified as cflags in the refinement routines.">
  
  ! After refinement, those points of quad elements on the fine mesh which 
  ! were formally element midpoints on the coarse mesh were recalculated
  ! by taking the mean of the corners.
  ! This helps avoiding tangled elements when there is a 'hole' in the domain.
  INTEGER(I32), PARAMETER :: TRIA_R2LV_AVERAGEMIDPOINTS  = 2**0

  ! After refinement, the coordinates of points on the boundary are
  ! recalculated using their parameter values (only 2D)
  INTEGER(I32), PARAMETER :: TRIA_R2LV_RECALCCOORDSONBD  = 2**1

  ! Standard parameter settings for 2-level refinement
  INTEGER(I32), PARAMETER :: TRIA_R2LV_STANDARD = TRIA_R2LV_AVERAGEMIDPOINTS + &
                                                  TRIA_R2LV_RECALCCOORDSONBD

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
  INTEGER(I32), PARAMETER :: TR_SHARE_IREFINEMENTPATCH       = 2**20   
  INTEGER(I32), PARAMETER :: TR_SHARE_ICOARSEGRIDELEMENT     = 2**21  
  
  INTEGER(I32), PARAMETER :: TR_SHARE_IVERTICESATFACE        = 2**22  ! KVAR 
  INTEGER(I32), PARAMETER :: TR_SHARE_IFACESATELEMENT        = 2**23  ! KAREA
  INTEGER(I32), PARAMETER :: TR_SHARE_IELEMENTSATFACE        = 2**24  ! K???
  INTEGER(I32), PARAMETER :: TR_SHARE_IEDGESATFACE           = 2**23 
  INTEGER(I32), PARAMETER :: TR_SHARE_IFACESATEDGEIDX        = 2**24
  INTEGER(I32), PARAMETER :: TR_SHARE_IFACESATEDGE           = 2**25
  INTEGER(I32), PARAMETER :: TR_SHARE_IFACESATVERTEXIDX      = 2**26
  INTEGER(I32), PARAMETER :: TR_SHARE_IFACESATVERTEX         = 2**27
  INTEGER(I32), PARAMETER :: TR_SHARE_IFACESATBOUNDARY       = 2**28
  
  INTEGER(I32), PARAMETER :: TR_SHARE_IBOUNDARYCPEDGESIDX    = 2**29
  INTEGER(I32), PARAMETER :: TR_SHARE_IBOUNDARYCPFACESIDX    = 2**30  
  
  
  
  ! Share everything
  INTEGER(I32), PARAMETER :: TR_SHARE_ALL = NOT(0_I32)

!</constantblock>
  
!<constantblock description="Format tags for TRI file formats.">

  ! Standard TRI file format, compatible to FEAT1.
  ! For 2D triangulations, Vertex coordinates of boundary vertices are
  ! replaced by parameter values.
  INTEGER(I32), PARAMETER :: TRI_FMT_STANDARD          = 0
  
  ! Standard TRI file format, but the vertex coordinates are 
  ! exported 'as they are', not as parameter values.
  INTEGER(I32), PARAMETER :: TRI_FMT_NOPARAMETRISATION = 2**0
  
!</constantblock>
  
!</constants>


!<types>

!<typeblock>

  ! Each 1D element consitst of at most TRIA_MAXNVE1D points.
  ! Each point has a number, which is usually an integer value.
  TYPE t_elementCorners1D
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE1D) :: Icorners
  END TYPE
  
  ! Each 2D element consists of at most TRIA_MAXNVE2D points.
  ! Each point has a number, which is usually an integer value.
  TYPE t_elementCorners2D
    INTEGER(PREC_VERTEXIDX), DIMENSION(TRIA_MAXNVE2D) :: Icorners
  END TYPE

  ! Each 2D element contains at most TRIA_MAXNME2D edges.
  ! Each edge has a number, which is usually an integer value.
  TYPE t_elementEdges2D
    INTEGER(PREC_EDGEIDX), DIMENSION(TRIA_MAXNME2D) :: Iedges
  END TYPE
  
  ! Each 1D element contains at most TRIA_MAXNVE1D neighbour elements,
  ! each meeting the element in a vertice.
  ! Each neighbour element has a number, which is usually an integer value.
  TYPE t_elementNeighbours1D
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNVE1D) :: Ineighbours
  END TYPE
  
  
  ! Each 2D element contains at most NMAXEDGES neighbour elements,
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
    ! Bit 20: IrefinementPatch/IrefinementPatchIndex is a copy of another structure
    ! Bit 21: IcoarseGridElement is a copy of another structure
    INTEGER(I32)             :: iduplicationFlag
  
    ! Dimension of the triangulation.
    ! NDIM1D=1D triangulation, NDIM2D=2D triangulation
    ! NDIM2D=3D triangulation, 0=not initialised
    INTEGER                  :: ndim = 0
  
    ! Number of points in the domain corresponding to corners of elements;
    ! coincides with SIZE(RcornerCoordinates)
    INTEGER(PREC_VERTEXIDX)   :: NVT = 0
    
    ! Number of edges in the domain belonging to elements
    INTEGER(PREC_EDGEIDX)    :: NMT = 0
    
    ! total number of faces in the domain
    INTEGER(PREC_EDGEIDX)    :: NAT = 0
    
    
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
    
    ! the number of faces on the boundary
    INTEGER             :: NABD = 0
    
    ! Number of edges on the boundary; coincides with NVBD
    ! for 2D domains.
    INTEGER             :: NMBD = 0

    ! Maximum number of vertices per element.
    ! 3 for 2D triangles, 4 for 2D quads or 3D tetrahedrons, 12 for 3D hexas
    INTEGER             :: NNVE = 0
    
    ! Maximum number of edges per element.
    ! 3 for 2D triangles, 4 for 2D quads or 3D tetrahedrons, 12 for 3D hexas.
    INTEGER             :: NNEE = 0
    
    ! Maximum number of areas per element. One hexa has e.g. 6 areas.
    INTEGER             :: NNAE = 0
    
    ! Maximum number of vertices per face. 3 for 3D tetraheral meshes,
    ! 4 for 3D hexahedral meshes. Unused in 2D.
    INTEGER             :: NNVA = 0
    
    ! Number of elements with a defined number of vertices per element.
    ! InelOfType(TRIA_NVELINE1D) = number of lines in the mesh (1D).
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
    INTEGER        :: h_DvertexCoords = ST_NOHANDLE
    
    ! Vertices Adjacent to an Element.
    ! Handle to h_IverticesAtElement=array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the corner-vertices
    ! in mathematically positive sense.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IverticesAtElement(4,.)=0 for a triangle in a quad mesh.
    ! This is a handle to the old KVERT array.
    INTEGER        :: h_IverticesAtElement = ST_NOHANDLE

    ! Edges Adjacent to an Element.
    ! Handle to 
    !       p_IedgesAtElement = array [1..NVE,1..NEL] of integer
    ! For each element the node numbers of the edges following the
    ! corner vertices in mathematically positive sense.
    ! This is the old KMID array.
    ! On pure triangular meshes, there is NVE=3. On mixed or pure quad
    ! meshes, there is NVE=4. In this case, there is 
    ! IedgesAtElement(4,.)=0 for a triangle in a quad mesh.
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
    
    ! Elements Adjacent to an Edge. Only 2D.
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
    
    ! This array defines a patch index and is created during the
    ! refinement. For a mesh that does not stem from a refinement,
    ! this array is undefined. For a mesh that stems from a refinement,
    ! this array defines a pointer into p_IrefinementPatch.
    ! Using p_IrefinementPatchIndex in combination with p_IrefinementPatch
    ! allows to get the element numbers of the fine grid elements
    ! that were created from a coarse grid element.
    ! Handle to 
    !       p_IrefinementPatchIndex=array [1..NELcoarse+1] of integer.
    ! So if IEL is the number of a coarse grid element, the fine grid 
    ! elements are to be found in
    ! p_IrefinementPatch ( p_IrefinementPatchIndex(IEL)..p_IrefinementPatchIndex(IEL+1)-1 )
    ! By subtracting
    !     p_IrefinementPatchIndex(IVT+1)-p_IrefinementPatchIndex(IVT)
    ! one can get the number of elements that a coarse grid element
    ! was refined to.
    INTEGER        :: h_IrefinementPatchIndex = ST_NOHANDLE
    
    ! This array defines patches that were created during the
    ! refinement. For a mesh that does not step from a refinement,
    ! this array is undefined. For a mesh that stems from a refinement,
    ! this array defines for every coarse grid element the number
    ! of the fine grid elements that were created from that coarse
    ! grid element.
    ! Handle to
    !       p_IrefinementPatches = array(1..*) of integer
    ! p_IrefinementPatch ( p_IrefinementPatchIndex(IEL)..p_IrefinementPatchIndex(IEL+1)-1 )
    ! contains the numbers of the elements, a coarse grid element
    ! was divided into.
    INTEGER        :: h_IrefinementPatch = ST_NOHANDLE
    
    ! If a mesh stems from a refinement of a coarse mesh, this array
    ! defines for every element on the fine mesh the element
    ! number on the coarse mesh where the fine grid element comes from.
    ! Handle to
    !       p_IcoarseGridElement = array(1..NEL) of integer
    ! For a mesh that does not come from a refinement, this handle is 
    ! undefined.
    INTEGER        :: h_IcoarseGridElement = ST_NOHANDLE
    
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

    ! In 2d this is a reference to the IboundaryCpIdx array.
    ! (Note: By du0plication flag, both arrays are identical and have
    ! the same handles).
    ! In 3d, this array is in index into the p_IedgesAtBoundary array that
    ! works like the p_IboundaryCpIdx for the vertices... see above.
    INTEGER          :: h_IboundaryCpEdgesIdx = ST_NOHANDLE
    
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
    INTEGER          :: h_IedgesAtBoundary = ST_NOHANDLE
    
    ! Ihis array is an index into the p_IfacesAtBoundary array it
    ! works like the p_IboundaryCpIdx for the vertices... see above.
    INTEGER          :: h_IboundaryCpFacesIdx = ST_NOHANDLE

    ! Faces adjacent to the boundary. Only 3D, undefined in 2D.
    ! Handle to 
    !       p_IfacesAtBoundary = array [1..NMBD] of integer.
    ! This array contains a list of all edges on the (real) boundary
    ! with increasing number.
    ! The boundary edges of boundary component i are saved at
    !        p_IboundaryCpFacesIdx(i)..p_IboundaryCpFacesIdx(i+1)-1.
    INTEGER          :: h_IfacesAtBoundary = ST_NOHANDLE

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
    
    ! Handle to h_IelementsAtEdgeIdx3d=array [1..NMT+1] of integer.
    ! Index array for h_IelementsAtEdge3d of length NMT+1 for describing the
    ! elements attached to an edge. for edge IVE, the array
    ! p_IelementsAtEdge3d contains the indices of the elements around this
    ! edge at the array positions
    !     p_IelementsAtEdgeIdx3d(IVE)..p_IelementsAtEdgeIdx3d(IVE+1)-1.
    ! By subtracting
    !     p_IelementsAtEdgeIdx3d(IVE+1)-p_IelementsAtEdgeIdx3d(IVE)
    ! One can get the number of elements attached to edge IVE. Only 3D.
    INTEGER        :: h_IelementsAtEdgeIdx3d = ST_NOHANDLE
    
    ! Elements Adjacent to an Edge. Only 3D.
    ! Array containing the Elements Adjacent to an edge.
    ! Handle to 
    !       p_IelementsAtEdge3d = array(1..*) of integer
    ! p_IelementsAtEdge3D ( p_IelementsAtEdgeIdx3d(IVT)..p_IelementsAtEdgeIdx3d(IVT+1)-1 )
    ! contains the number of the adjacent element in an edge.
    INTEGER        :: h_IelementsAtEdge3d = ST_NOHANDLE
    
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
    INTEGER        :: h_IfacesAtEdgeIdx = ST_NOHANDLE

    ! Array containing the Faces Adjacent to an Edge.
    ! Handle to 
    !       p_IfacesAtEdge = array(1..*) of integer
    ! p_IfacesAtEdge ( p_IfacesAtEdgeIdx(IVT)..p_IfacesAtEdgeIdx(IVT+1)-1 )
    ! contains the number of the adjacent faces in an edge.
    INTEGER        :: h_IfacesAtEdge    = ST_NOHANDLE
    
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
    INTEGER        :: h_IfacesAtVertexIdx = ST_NOHANDLE

    ! Array containing the Faces Adjacent to a Vertex. Only 3D.
    ! Handle to 
    !       p_IfacesAtVertex = array(1..*) of integer
    ! p_IfacesAtVertex ( p_IfacesAtVertexIdx(IVT)..p_IfacesAtVertexIdx(IVT+1)-1 )
    ! contains the number of the adjacent faces in a vertex.
    INTEGER        :: h_IfacesAtVertex    = ST_NOHANDLE
     
    ! Faces adjacent to an element. Only 3D.
    ! Handle to 
    !       p_IfacesAtElement = array [1..NNAE,1..NEL] of integer
    ! For each element the node numbers of the edges following the
    ! corner vertices in mathematically positive sense.
    ! This is the old KMID array.
    ! On pure tethrahedral meshes, there is NNAE=4. On mixed or pure hexahedral
    ! meshes, there is NNAE=8. In this case, there is 
    ! IedgesAtElement(5:8,.)=0 for a tethrahedral in a hexahedral mesh.
    ! To be able to distinguish a number of an edge from a vertex number, 
    ! edges are numbered in the range NVT+NMT+1..NVT+NMT+NAT.
    INTEGER        :: h_IfacesAtElement = ST_NOHANDLE
    
    ! Vertices Adjacent to a face.
    ! Handle to 
    !       p_IverticesAtFace = array [1..NNVA,1..NEL] of integer
    ! For each element the node numbers of the edges following the
    ! corner vertices in mathematically positive sense.
    ! On pure tetrahedral meshes, there is NVA=3. On mixed or pure 
    ! hexahedral meshes, there is NVA=4. In this case, there is 
    ! IverticesAtFace(4,.)=0 for a tetrahedral in a hexahedral mesh.
    INTEGER        :: h_IverticesAtFace = ST_NOHANDLE 
    
    ! Elements Adjacent to a Face. Only 3D.
    ! Handle to 
    !       p_IelementsAtEdge = array [1..2,1..NAT] of integer.
    ! The numbers of the two elements adjacent to an edge IMT in 2D. 
    ! For boundary edges, p_IelementsOnEdge(2,IMT) is set to 0.
    ! This is the old KMEL array.
    INTEGER        :: h_IelementsAtFace = ST_NOHANDLE

    ! Edges Adjacent to a Face. Only 3D.
    ! Handle to 
    !       p_IedgesAtFace = array [1..NNVA,1..NAT] of integer.
    ! The numbers of the edges adjacent to a face in 3D. 
    INTEGER        :: h_IedgesAtFace

  END TYPE

!</typeblock>

!<typeblock>
  ! a connector connects to adjacent cells (i.e. a face in 3d)
  ! structure used to generate 3d connectivity
  type t_connector3d
      ! the array stores at 1-4 the vertices of a face
      ! at 5 the element the face belongs to
      ! at 6 it stores the local face number
      integer, dimension(6) :: I_conData
  end type
!</typeblock>
  
!</types>

  INTERFACE tria_getNVE
    MODULE PROCEDURE tria_getNVE_direct
    MODULE PROCEDURE tria_getNVE_indirect
  END INTERFACE

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

    !$OMP PARALLEL DO PRIVATE(i)
    DO j=0,idim2-1
      DO i=0,idim1-1
        p_array(i+1,j+1) = KWORK(kpos+idim1*j+i)
      END DO
    END DO
    !$OMP END PARALLEL DO
    
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
    
    !$OMP PARALLEL DO
    DO i=0,idim1-1
      p_array(i+1) = KWORK(kpos+i) !p_array2(i)
    END DO
    !$OMP END PARALLEL DO
    
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
    
    !$OMP PARALLEL DO PRIVATE(i)
    DO j=0,idim2-1
      DO i=0,idim1-1
        p_array(i+1,j+1) = DWORK(kpos+idim1*j+i)
      END DO
    END DO
    !$OMP END PARALLEL DO
    
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
    !$OMP PARALLEL DO
    DO i=0,idim1-1
      p_array(i+1) = DWORK(kpos+i)
    END DO
    !$OMP END PARALLEL DO
    
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
      rbackupTriangulation%NAT                    = rtriangulation%NAT
      rbackupTriangulation%NBCT                   = rtriangulation%NBCT
      rbackupTriangulation%NVBD                   = rtriangulation%NVBD
      rbackupTriangulation%NMBD                   = rtriangulation%NMBD
      rbackupTriangulation%NNVE                   = rtriangulation%NNVE
      rbackupTriangulation%NNEE                   = rtriangulation%NNEE
      rbackupTriangulation%NNAE                   = rtriangulation%NNAE
      rbackupTriangulation%NNVA                   = rtriangulation%NNVA
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

    ! Bit 20: IrefinementPatch
    CALL checkAndCopy(idupflag, TR_SHARE_IREFINEMENTPATCH,&
          rtriangulation%h_IrefinementPatch, &
          rbackupTriangulation%h_IrefinementPatch)
    CALL checkAndCopy(idupflag, TR_SHARE_IREFINEMENTPATCH,&
          rtriangulation%h_IrefinementPatchIndex, &
          rbackupTriangulation%h_IrefinementPatchIndex)

    ! Bit 21: IcoarseGridElement 
    CALL checkAndCopy(idupflag, TR_SHARE_ICOARSEGRIDELEMENT,&
          rtriangulation%h_IcoarseGridElement, &
          rbackupTriangulation%h_IcoarseGridElement)
    
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
    ! Bit  5: KMEL  3d Version          
    CALL checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdgeIdx3d)
          
    ! Bit  5: KMEL  3d Version      
    CALL checkAndRelease(idupflag, TR_SHARE_IELEMENTSATEDGE,&
          rtriangulation%h_IelementsAtEdge3d)

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

    ! Bit 20: IrefinementPatch          
    CALL checkAndRelease(idupflag, TR_SHARE_IREFINEMENTPATCH, &
           rtriangulation%h_IrefinementPatch)
    CALL checkAndRelease(idupflag, TR_SHARE_IREFINEMENTPATCH, &
           rtriangulation%h_IrefinementPatchIndex)

    ! Bit 21: IcoarseGridElement
    CALL checkAndRelease(idupflag, TR_SHARE_ICOARSEGRIDELEMENT, &
           rtriangulation%h_IcoarseGridElement)

    ! Bit 22: KVAR          
    CALL checkAndRelease(idupflag, TR_SHARE_IVERTICESATFACE, &
           rtriangulation%h_IverticesAtFace)
           
    ! Bit 23: KAREA           
    CALL checkAndRelease(idupflag, TR_SHARE_IFACESATELEMENT, &
           rtriangulation%h_IFacesAtElement)
    
    ! Bit 24: K????       
    CALL checkAndRelease(idupflag, TR_SHARE_IELEMENTSATFACE, &
           rtriangulation%h_IelementsAtFace)
           
    ! Bit 23: K????       
    call checkAndRelease(idupflag, TR_SHARE_IELEMENTSATFACE, &
           rtriangulation%h_IelementsAtFace)
           
    ! Bit 24: K????       
    call checkAndRelease(idupflag, TR_SHARE_IEDGESATFACE, &
           rtriangulation%h_IedgesAtFace)
           
    ! Bit 25: K????       
    call checkAndRelease(idupflag, TR_SHARE_IFACESATEDGEIDX, &
           rtriangulation%h_IfacesAtEdgeIdx)
       
    call checkAndRelease(idupflag, TR_SHARE_IFACESATVERTEXIDX, &
           rtriangulation%h_IfacesAtVertexIdx)

    call checkAndRelease(idupflag, TR_SHARE_IFACESATVERTEX, &
           rtriangulation%h_IfacesAtVertex)
    
    call checkAndRelease(idupflag, TR_SHARE_IFACESATEDGE, &
           rtriangulation%h_IfacesAtEdge)

    call checkAndRelease(idupflag, TR_SHARE_IFACESATBOUNDARY, &
           rtriangulation%h_IfacesAtBoundary)
           
    CALL checkAndRelease(idupflag,TR_SHARE_IBOUNDARYCPEDGESIDX,&
          rtriangulation%h_IboundaryCpEdgesIdx)

    CALL checkAndRelease(idupflag,TR_SHARE_IBOUNDARYCPFACESIDX,&
          rtriangulation%h_IboundaryCpFacesIdx)
           
           
    
    ! Clean up the rest of the structure

    rtriangulation%iduplicationFlag = 0
    rtriangulation%ndim = 0
    rtriangulation%NVT = 0
    rtriangulation%NNVE = 0
    rtriangulation%NNAE = 0
    rtriangulation%NNVA = 0
    rtriangulation%NNEE = 0
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

    ! Check if quadrilaterals exists at all
    IF (icount .EQ. 0) RETURN
    
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
    rtriangulation%nnve = 3
    
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
  ! The name of the .tri file to read.
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
        rtriangulation%NNVE,rtriangulation%NBCT
        
    nve = rtriangulation%NNVE

    ! Comment: 'DCORVG'
    READ (iunit,*)

    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM2D, NVT)
    Isize = (/NDIM2D,INT(rtriangulation%NVT,I32)/)
    CALL storage_new2D ('tria_read_tri2D', 'DCORVG', Isize, ST_DOUBLE, &
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    CALL storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
        
    ! Read the data from the file, store it in the array.
    ! read data into p_Ddata:
    ! first read nvt x-coordinates into p_Ddata(1,ivt)
    ! then read nvt  y-coordinates into p_Ddata(2,ivt)
    READ (iunit,*) ((p_Ddata2D(idim,ivt),idim=1,NDIM2D), ivt=1,rtriangulation%NVT)

    ! Comment: 'KVERT'
    READ (iunit,*)

    

    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/nve,INT(rtriangulation%NEL,I32)/)
    CALL storage_new2D ('tria_read_tri2D', 'KVERT', Isize, ST_INT, &
        rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    CALL storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)

    ! read ive=1 indices to nve into p_Idata2D(ive,iel) where iel=1 to NEL
    READ (iunit,*) ((p_Idata2D(ive,iel),ive=1,nve), iel=1,rtriangulation%NEL)

    ! Loop through the elements and determine how many elements
    ! of each element type we have.
    rtriangulation%InelOfType(:) = 0
    DO iel=1,rtriangulation%nel
      ! start at the last index of element iel down to the first
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
    ivbd = 0
    !$OMP PARALLEL DO REDUCTION(+:ivbd)
    DO ivt=1,rtriangulation%NVT
      !IF (p_InodalProperty(ivt) .NE. 0) rtriangulation%NVBD = rtriangulation%NVBD+1
      IF (p_InodalProperty(ivt) .NE. 0) ivbd = ivbd+1
    END DO
    !$OMP END PARALLEL DO
    rtriangulation%NVBD = ivbd

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
      !$OMP PARALLEL DO PRIVATE(ibct,ivbd)
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
      !$OMP END PARALLEL DO
      
    ELSE
    
      ! No parametrisation available, the array with boundary parameter values 
      ! is not generaterd.
      !
      ! Check all vertices to find out, which vertices are on the boundary.
      !$OMP PARALLEL DO PRIVATE(ibct,ivbd)
      DO ivt=1,rtriangulation%NVT
        IF (p_InodalProperty(ivt) .NE. 0) THEN
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
        END IF
      END DO
      !$OMP END PARALLEL DO
      
    END IF
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_createRawTria1D(rtriangulation,dleft,dright,nintervals)
  
!<description>
  ! This routine creates a 'raw' 1D triangulation with ninvervals
  ! sub-intervals of same length.
!</description>

!<input>
  ! The left end of the interval. Must be < dright.
  REAL(DP), INTENT(IN) :: dleft
  
  ! The right end of the interval. Must be > dleft.
  REAL(DP), INTENT(IN) :: dright
  
  ! OPTIONAL: The number of sub-intervals to create. If given, ninvervals
  ! must be > 0. If not given, one interval is created.
  INTEGER, OPTIONAL, INTENT(IN) :: nintervals
!</input>

!<output>
  ! The triangulation.
  TYPE(t_triangulation), INTENT(OUT) :: rtriangulation
!</output>

!</subroutine>

  ! Some local variables
  INTEGER :: i,nintv
  REAL(DP) :: t, s
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dcoords
  INTEGER, DIMENSION(:,:), POINTER :: p_Iverts
  INTEGER, DIMENSION(:), POINTER :: p_Idata
  INTEGER, DIMENSION(2) :: Isize
    
    ! Check parameters
    IF (dleft .GE. dright) THEN
      PRINT *, "tria_createRawTria1D: dleft must be less than dright!"
      STOP
    END IF
    
    ! Set number of intervals
    nintv = 1
    IF (PRESENT(nintervals)) THEN
      IF (nintervals .GT. 0) nintv = nintervals
    END IF
  
    ! Set the triangulation's dimension
    rtriangulation%ndim = NDIM1D

    ! We have nintv+1 vertices
    rtriangulation%NVT = nintv+1
    
    rtriangulation%NNVE = 2

    ! Allocate vertices
    Isize = (/1, nintv+1/)
    CALL storage_new2D('tria_createRawTria1D', 'DCORVG', Isize, &
        ST_DOUBLE, rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_Dcoords)
    
    ! Initialise vertices
    p_Dcoords(1,1) = dleft
    s = 1.0_DP / REAL(nintv, DP)
    DO i=2, nintv
      t = REAL(i-1, DP) * s
      p_Dcoords(1,i) = (1.0_DP - t) * dleft + t * dright
    END DO
    p_Dcoords(1,nintv+1) = dright
    
    ! And we have nintv elements
    rtriangulation%NEL = nintv
    rtriangulation%InelOfType(TRIA_NVELINE1D) = nintv

    ! Allocate elements
    Isize = (/2, nintv/)
    CALL storage_new2D('tria_createRawTria1D', 'KVERT', Isize, &
        ST_INT, rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int2d(rtriangulation%h_IverticesAtElement, p_Iverts)
    
    ! Initialise elements
    DO i=1, nintv
      p_Iverts(1,i) = i
      p_Iverts(2,i) = i+1
    END DO
    
    ! There is one boundary component - the interval ends
    rtriangulation%NBCT = 1
    
    ! Allocate memory for boundary components
    CALL storage_new ('tria_createRawTria1D', 'KBCT', 3, ST_INT, &
        rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the boundary components
    CALL storage_getbase_int(rtriangulation%h_IboundaryCpIdx, p_Idata)
    p_Idata = (/ 1, 2, 3 /)

    ! There is one vertice per boundary component
    rtriangulation%NVBD = 2
    
    ! Allocate memory for boundary components
    CALL storage_new ('tria_createRawTria1D', 'KVBD', 2, ST_INT, &
        rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the boundary components
    CALL storage_getbase_int(rtriangulation%h_IverticesAtBoundary, p_Idata)
    p_Idata = (/ 1, nintv+1 /)

    ! Allocate memory for nodal property
    CALL storage_new ('tria_createRawTria1D', 'KNPR', nintv+1, ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    CALL storage_getbase_int(rtriangulation%h_InodalProperty,p_Idata)
    
    ! Set up nodal property
    p_Idata(1) = 1
    p_Idata(nintv+1) = 2

    ! That's it
      
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
      ! Generate all standard arrays for 1D meshes.
      CALL tria_genElementsAtVertex2D    (rtriangulation)
      CALL tria_genNeighboursAtElement1D (rtriangulation)
      CALL tria_genElementVolume1D       (rtriangulation)

      CALL tria_sortBoundaryVertices2D   (rtriangulation)
      CALL tria_genElementsAtBoundary2D  (rtriangulation)
      CALL tria_genBoundaryVertexPos2D   (rtriangulation)
       
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
    ! vertices at element info provided by tri-File
    call tria_genElementsAtVertex3D   (rtriangulation)
    call tria_genNeighboursAtElement3D(rtriangulation)
    call tria_genEdgesAtElement3D     (rtriangulation)
    call tria_genElementsAtEdge3D     (rtriangulation)
    call tria_genVerticesAtEdge3D     (rtriangulation)
    
    ! faces have global numbers
    ! nvt+nmt+1 = first face
    call tria_genFacesAtElement       (rtriangulation)
    call tria_genVerticesAtFace       (rtriangulation)
    call tria_genElementsAtFace       (rtriangulation) 
    call tria_genEdgesAtFace          (rtriangulation)
    call tria_genFacesAtEdge          (rtriangulation)
    call tria_genFacesAtVertex        (rtriangulation)
    
    !----BOUNDARY------
    call tria_genFacesAtBoundary      (rtriangulation)
    call tria_genEdgesAtBoundary3D    (rtriangulation)

    !----Properties----!
    call tria_genEdgeNodalProperty3d  (rtriangulation)
    call tria_genFaceNodalProperty3d  (rtriangulation)

    
    ! CALL tria_genElementVolume3d       (rtriangulation)
    
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
    !$OMP PARALLEL DO
    DO ivbd=1,rtriangulation%NVBD
      p_Iresort(ivbd) = ivbd
    END DO
    !$OMP END PARALLEL DO
    
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
    !$OMP PARALLEL DO
    DO ivbd=1,SIZE(p_IverticesAtBoundary)
      p_Iresort(ivbd) = p_IverticesAtBoundary(p_Iresort(ivbd))
    END DO
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO
    DO ivbd=1,SIZE(p_IverticesAtBoundary)
      p_IverticesAtBoundary(ivbd) = p_Iresort(ivbd)
    END DO
    !$OMP END PARALLEL DO
    
    CALL storage_free (hresort)

  END SUBROUTINE

!************************************************************************

!<function>

  INTEGER FUNCTION tria_getNNVE(rtriangulation)

!<description>
  ! Determines the maximum number of vertices per element in
  ! rtriangulation.
  !
  ! DEPRECATED: Alternatively, the rtriangulation%NNVE-variable can be 
  ! used. tria_getNNVE determines NNVE without rtriangulation%NNVE by
  ! checking the array size.
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

!<function>

  INTEGER FUNCTION tria_getNNVA(rtriangulation)

!<description>
  ! Determines the maximum number of vertices per face in
  ! rtriangulation.
!</description>

!<inputoutput>
  ! The triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
!</inputoutput>
  
!<result>
  ! Maximum number of vertices per face.
!</result>
  
!</function>

    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_Idata2D
    
    ! NNVE is given by the first dimension of KVERT
    CALL storage_getbase_int2D(&
        rtriangulation%h_IverticesAtFace,p_Idata2D)
        
    tria_getNNVA = UBOUND(p_Idata2D,1)

  END FUNCTION

  ! ***************************************************************************

!<function>

  PURE INTEGER FUNCTION tria_getNVE_indirect (IvertEdgAtElement,iel)
  
!<description>
  ! This routine calculates the number of vertices/edges on element iel.
!</description>

!<input>
  ! This may be either the IverticesAtElement or the IedgesAtElement
  ! array of a triangulation
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IvertEdgAtElement
  
  ! Number of the element whose NVE should be calculated
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
!</input>

!<result>
  ! Number of vertices on element iel.
!</result>

!</subroutine>

    INTEGER :: i

    DO i=UBOUND(IvertEdgAtElement,1),1,-1
      IF (IvertEdgAtElement(i,iel) .NE. 0) THEN
        tria_getNVE_indirect = i
        RETURN
      END IF
    END DO
  
  END FUNCTION

  ! ***************************************************************************

!<function>

  INTEGER FUNCTION tria_getNVE_direct (rtriangulation,iel)
  
!<description>
  ! This routine calculates the number of vertices/edges on element iel.
!</description>

!<input>
  ! Triangulation structure
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! Number of the element whose NVE should be calculated
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
!</input>

!<result>
  ! Number of vertices on element iel.
!</result>

!</subroutine>

    INTEGER :: i
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    
    CALL storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)

    DO i=UBOUND(p_IverticesAtElement,1),1,-1
      IF (p_IverticesAtElement(i,iel) .NE. 0) THEN
        tria_getNVE_direct = i
        RETURN
      END IF
    END DO
  
  END FUNCTION

  !************************************************************************

!<subroutine>

  ELEMENTAL SUBROUTINE tria_getNeighbourVertex(ivertex,ivt1,ivt2,ineighbour)

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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: ivertex
  
  ! Vertex number of one vertex adjacent to an edge.
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: ivt1

  ! Vertex number of the other vertex adjacent to that edge.
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: ivt2
!</input>

!<output>
  ! Vertex number of the neighbour of ivertex.
  INTEGER(PREC_VERTEXIDX), INTENT(OUT) :: ineighbour
!</output>

!</subroutine>

    ! Note: Directly implementing this formula into the program
    ! code brings more speed :-)
    ! But to have a reference not to forget the formula, we have
    ! this routine...
    
    ineighbour = ivt1 + ivt2 - ivertex
    
  END SUBROUTINE 

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
    INTEGER(PREC_VERTEXIDX) :: ivt,isize,isize2
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
    
    nnve = rtriangulation%NNVE

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
      CALL storage_getsize (rtriangulation%h_IelementsAtVertex, isize2)
      IF (isize .NE. isize2) THEN
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
        
        ! Remember the element number and increase the pointer in the
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

  SUBROUTINE tria_genNeighboursAtElement1D(rtriangulation)

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
    INTEGER :: nnve
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    
    INTEGER(PREC_ELEMENTIDX) :: iel1, iel2, ivi1, ivi2
    INTEGER(PREC_VERTEXIDX) :: ivt, ive
    
    nnve = rtriangulation%NNVE

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
   
    ! Get the index array that tells us how many elements are adjacent to
    ! each vertex.
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertex,&
        p_IelementsAtVertex)

    ! We have to look at all vertices to find neighbouring information.
    !
    ! Loop through all vertices.
    DO ive = 1, rtriangulation%NVT
    
      ! Get both elements at this vertice
      ivi1 = p_IelementsAtVertexIdx(ive)
      ivi2 = p_IelementsAtVertexIdx(ive+1)
      
      ! If (ivi2 - ivi1) < 2, then we have no neighbour elements on
      ! this vertice
      IF ((ivi2 - ivi1) .LT. 2) CYCLE
      
      ! Get indices of elements on this vertice
      iel1 = p_IelementsAtVertex(ivi1)
      iel2 = p_IelementsAtVertex(ivi1+1)
      
      ! Add iel2 as a neighbour to iel1
      DO ivt = 1, nnve
        IF (p_IverticesAtElement(ivt, iel1) .EQ. ive) THEN
          p_IneighboursAtElement(ivt, iel1) = iel2
          EXIT
        END IF
      END DO
    
      ! Add iel1 as a neighbour to iel2
      DO ivt = 1, nnve
        IF (p_IverticesAtElement(ivt, iel2) .EQ. ive) THEN
          p_IneighboursAtElement(ivt, iel2) = iel1
          EXIT
        END IF
      END DO

    END DO ! ive
    
    ! That's it

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

    nnve = rtriangulation%NNVE

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
      
      END DO ! ive
      
    END DO ! iel
    
    ! Next, we have to look at all edges to find neighbouring information.
    !
    ! Loop through all edges adjacent to all vertices.
    !$OMP PARALLEL DO PRIVATE(ivt,ivtneighbour,iel,iedgeneighbour)
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
    !$OMP END PARALLEL DO
    
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
    
    ! in 3d this is a list because there no fixed number
    ! of elements at an edge
    ! Do we have (enough) memory for that array?
    IF (rtriangulation%h_IelementsAtEdge .EQ. ST_NOHANDLE) THEN
      Isize = (/2,rtriangulation%NMT/)
      CALL storage_new2D ('tria_genElementsAtEdge2D', 'KMID', &
          Isize, ST_INT, &
          rtriangulation%h_IelementsAtEdge, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize2D (rtriangulation%h_IelementsAtEdge, Isize)
      IF (Isize(2) .NE. rtriangulation%NMT) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genElementsAtEdge2D', &
            rtriangulation%NMT, rtriangulation%h_IelementsAtEdge, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int2D (rtriangulation%h_IelementsAtEdge,p_IelementsAtEdge)
    
    NVT = rtriangulation%NVT
    
    ! Loop through all elements and all edges on the elements
    !$OMP PARALLEL DO PRIVATE(ive,iedge)
    ! all elements
    DO iel = 1,UBOUND(p_IedgesAtElement,2)
      ! loop 1-4
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
    !$OMP END PARALLEL DO
    
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
      IF (Isize(2) .NE. rtriangulation%NMT) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genVerticesAtEdge2D', &
            rtriangulation%NMT, rtriangulation%h_IverticesAtEdge, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    
    nnve = UBOUND(p_IedgesAtElement,1)
    NVT = rtriangulation%NVT
    
    ! Loop through all elements and all edges on the elements
    !$OMP PARALLEL DO PRIVATE(ive,iedge,ivtneighbour)
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
    !$OMP END PARALLEL DO
    
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
        Isize(2) = rtriangulation%NEL
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
    INTEGER(PREC_VERTEXIDX) :: isize

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
    IF (isize .LT. rtriangulation%NVT+rtriangulation%NMT) THEN
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
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    CALL storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
        p_IneighboursAtElement)
        
    ! Initialise the nodal property with 0 by default.
    CALL lalg_clearVectorInt (&
        p_InodalProperty(rtriangulation%NVT+1:rtriangulation%NVT+rtriangulation%NMT))
    
    ! Loop through all elements and all edges on the elements
    !$OMP PARALLEL DO PRIVATE(ive)
    DO iel = 1,UBOUND(p_IedgesAtElement,2)
    
      DO ive = 1,UBOUND(p_IedgesAtElement,1)
      
        ! Stop if we handled all edges; this is important if there are triangles
        ! in a quad mesh e.g.
        IF (p_IedgesAtElement(ive,iel) .EQ. 0) EXIT
        
        ! The edge nodal property is initialised with 0 by default -- inner edge.
        ! Is there a neighbour? If yes, we have an inner edge. If not, this is 
        ! a boundary edge. 
        IF (p_IneighboursAtElement(ive,iel) .EQ. 0) THEN
        
          ! Get the number of the boundary component from the vertex preceeding
          ! the edge and store it as information for the edge. 
          p_InodalProperty(p_IedgesAtElement(ive,iel)) = &
              p_InodalProperty(p_IverticesAtElement(ive,iel))
        
        END IF
        
      END DO
      
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_genElementVolume1D(rtriangulation)

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
    REAL(DP), DIMENSION(TRIA_MAXNVE1D) :: Dpoints
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
      CALL storage_new ('tria_genElementVolume1D', 'DAREA', &
          INT(rtriangulation%NEL+1,I32), ST_DOUBLE, &
          rtriangulation%h_DelementVolume, ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize (rtriangulation%h_DelementVolume, isize)
      IF (isize .NE. rtriangulation%NEL+1) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genElementVolume1D', &
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
        
    ! Calculate the element volume for all elements
    !$OMP PARALLEL DO PRIVATE(ive,Dpoints) REDUCTION(+:dtotalVolume)
    DO iel=1,rtriangulation%NEL
    
      ! line element
      DO ive=1,TRIA_NVELINE1D
        Dpoints(ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
      END DO
      p_DelementVolume(iel) = ABS(Dpoints(1) - Dpoints(2))
      
      dtotalVolume = dtotalVolume+p_DelementVolume(iel)
    END DO
    !$OMP END PARALLEL DO
    
    ! Store the total volume in the last element of DelementVolume
    p_DelementVolume(rtriangulation%NEL+1) = dtotalVolume
    
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
      !$OMP PARALLEL DO PRIVATE(ive, Dpoints) REDUCTION(+:dtotalVolume)
      DO iel=1,rtriangulation%NEL
        ! triangular element
        DO ive=1,TRIA_NVETRI2D
          Dpoints(1,ive) = p_DvertexCoords(1,p_IverticesAtElement(ive,iel))
          Dpoints(2,ive) = p_DvertexCoords(2,p_IverticesAtElement(ive,iel))
        END DO
        p_DelementVolume(iel) = gaux_getArea_tria2D(Dpoints)
        
        dtotalVolume = dtotalVolume+p_DelementVolume(iel)
      END DO
      !$OMP END PARALLEL DO
    
    ELSE

      ! Calculate the element volume for all elements
      !$OMP PARALLEL DO PRIVATE(ive,Dpoints) REDUCTION(+:dtotalVolume)
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
      !$OMP END PARALLEL DO
      
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
    ! get the vertices at the element
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    ! get the neighbours at the element
    CALL storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
    ! get the vertices on the boundary
    CALL storage_getbase_int (rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
    ! get the index array for the boundary component
    CALL storage_getbase_int (rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    ! get the elements at a vertex
    CALL storage_getbase_int (rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
    ! get the elements at vertex index array
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

    nnve = rtriangulation%NNVE

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
      CALL storage_getsize (rtriangulation%h_IedgesAtBoundary, isize)
      IF (isize .NE. rtriangulation%NVBD) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_genEdgesAtBoundary2D', &
            rtriangulation%NVBD, rtriangulation%h_IedgesAtBoundary, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF
    END IF
    
    CALL storage_getbase_int (rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)

    nnve = rtriangulation%NNVE

    ! Loop through all boundary components
    !$OMP PARALLEL DO PRIVATE(ivbd,ivt,iel,ive)
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
    !$OMP END PARALLEL DO
    
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
    !$OMP PARALLEL DO
    DO ivbd = 1,rtriangulation%NVBD
      p_IboundaryVertexPos(1,ivbd) = p_IverticesAtBoundary(ivbd)
      p_IboundaryVertexPos(2,ivbd) = ivbd
    END DO
    !$OMP END PARALLEL DO
    
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
    !$OMP PARALLEL DO
    DO ivbd = 1,rtriangulation%NMBD
      p_IboundaryEdgePos(1,ivbd) = p_IedgesAtBoundary(ivbd)
      p_IboundaryEdgePos(2,ivbd) = ivbd
    END DO
    !$OMP END PARALLEL DO
    
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
      !$OMP PARALLEL DO PRIVATE(dpar1,dpar2)
      DO ivbd = p_IboundaryCpIdx(ibct),p_IboundaryCpIdx(ibct+1)-2
      
        ! Get the parameter value of the current vertex and its neighbour.
        dpar1 = p_DvertexParameterValue(ivbd)
        dpar2 = p_DvertexParameterValue(ivbd+1)
        
        ! The edge parameter value is the mean.
        p_DedgeParameterValue(ivbd) = 0.5_DP*(dpar1+dpar2)
      
      END DO
      !$OMP END PARALLEL DO
      
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

  SUBROUTINE tria_quickRefine2LevelOrdering(nfine,rtriangulation,rboundary,cflags)

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
  ! For a value <= 0, nothing will happen.
  INTEGER, INTENT(IN) :: nfine

  ! OPTIONAL: Boundary structure that defines the parametrisation of the boundary.
  ! If specified, the coordinates of the new boundary vertices are
  ! recomputed according to the analytic boundary.
  TYPE(t_boundary), INTENT(IN), OPTIONAL :: rboundary

  ! OPTIONAL: Bitfield of TRIA_R2LV_xxxx constants that allow to specify
  ! options for the refinement. If not specified, TRIA_R2LV_STANDARD
  ! is used as default.
  INTEGER(I32), INTENT(IN), OPTIONAL :: cflags
!</input>

!<inputoutput>
  ! The triangulation to be refined; can be a 'raw' or a 'standard' mesh.
  ! Is overwritten by the refined mesh.
  TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
!</inputoutput>

!</subroutine>
 
    INTEGER :: ifine

    SELECT CASE(rtriangulation%ndim)
    CASE (NDIM1D)
      ! Refine nfine times:
      DO ifine = 1,nfine
      
        ! Create missing arrays in the source mesh.
        ! Create only those arrays we need to make the refinement process
        ! as fast as possible.
        IF (rtriangulation%h_IelementsAtVertex .EQ. ST_NOHANDLE) &
          CALL tria_genElementsAtVertex2D    (rtriangulation)
        IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) &
          CALL tria_genNeighboursAtElement2D (rtriangulation)

        CALL tria_sortBoundaryVertices2D   (rtriangulation)

        IF (rtriangulation%h_IelementsAtBoundary .EQ. ST_NOHANDLE) &
          CALL tria_genElementsAtBoundary2D  (rtriangulation)
       
        ! Refine the mesh, replace the source mesh.
        CALL tria_refine2LevelOrdering(rtriangulation,rboundary=rboundary,&
            cflags=cflags)
            
      END DO

    CASE (NDIM2D)
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
        IF (reallocRefinedInodalProperty2D (rtriangulation)) &
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
        CALL tria_refine2LevelOrdering(rtriangulation,rboundary=rboundary,&
            cflags=cflags)
            
      END DO
      
      CASE (NDIM3D)
       DO ifine = 1,nfine     
       
        IF (rtriangulation%h_IelementsAtVertex .EQ. ST_NOHANDLE) &
          CALL tria_genElementsAtVertex3D    (rtriangulation)
          
   
        IF (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) &
          CALL tria_genNeighboursAtElement3D (rtriangulation)
          
        IF (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) &
          CALL tria_genEdgesAtElement3D      (rtriangulation)

        IF (rtriangulation%h_IelementsAtEdge3d .EQ. ST_NOHANDLE) &
          CALL tria_genElementsAtEdge3D      (rtriangulation)
          
        IF (rtriangulation%h_IverticesAtEdge .EQ. ST_NOHANDLE) &
          CALL tria_genVerticesAtEdge3D      (rtriangulation)
 
        IF (rtriangulation%h_IfacesAtElement .EQ. ST_NOHANDLE) &
          CALL tria_genFacesAtElement        (rtriangulation)
           
        IF (rtriangulation%h_IverticesAtFace .EQ. ST_NOHANDLE) &
          CALL tria_genVerticesAtFace        (rtriangulation)          
          
        IF (rtriangulation%h_IelementsAtFace .EQ. ST_NOHANDLE) &
          CALL tria_genElementsAtFace        (rtriangulation)             
          
        IF (rtriangulation%h_IfacesAtBoundary .EQ. ST_NOHANDLE) &
          CALL tria_genFacesAtBoundary       (rtriangulation)             

        IF (rtriangulation%h_IedgesAtFace .EQ. ST_NOHANDLE) &
          CALL tria_genEdgesAtFace           (rtriangulation)             

          
        IF (rtriangulation%h_IfacesAtEdge .EQ. ST_NOHANDLE) &
          CALL tria_genFacesAtEdge           (rtriangulation)             

        IF (rtriangulation%h_IedgesAtBoundary .EQ. ST_NOHANDLE) &
          CALL tria_genEdgesAtBoundary3d     (rtriangulation)
          
        call tria_genEdgeNodalProperty3d  (rtriangulation)
        
        call tria_genFaceNodalProperty3d  (rtriangulation)
          
          
                       

!    call tria_genFacesAtEdge          (rtriangulation)          

!    call tria_genNeighboursAtElement3D(rtriangulation)
!    call tria_genEdgesAtElement3D     (rtriangulation)
!    call tria_genElementsAtEdge3D     (rtriangulation)
!    call tria_genVerticesAtEdge3D     (rtriangulation)
        
        
!    call tria_genElementsAtVertex3D   (rtriangulation)
!    call tria_genFacesAtElement3D     (rtriangulation)
!    call tria_genNeighboursAtElement3D(rtriangulation)
!    call tria_genEdgesAtElement3D     (rtriangulation)
!    call tria_genElementsAtEdge3D     (rtriangulation)
!    call tria_genVerticesAtEdge3D     (rtriangulation)
!    
!    ! faces have global numbers
!    ! nvt+nmt+1 = first face
!    call tria_genFacesAtElement       (rtriangulation)
!    call tria_genVerticesAtFace       (rtriangulation)
!    call tria_genElementsAtFace       (rtriangulation) 
!    call tria_genEdgesAtFace          (rtriangulation)
!    call tria_genFacesAtEdge          (rtriangulation)
!    call tria_genFacesAtVertex        (rtriangulation)
   
       
        CALL tria_refine2LevelOrdering(rtriangulation,rboundary=rboundary,&
            cflags=cflags)
       END DO
    END SELECT
    
  CONTAINS
  
    LOGICAL FUNCTION reallocRefinedInodalProperty2D (rtriangulation)
    
    ! Reallocates the InodalProperty-array in rtriangulation such that
    ! is provides enough space for the refined mesh.
    
    TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
    
    ! Return value: whether the memory was reallocated or not.
    
      ! local variables
      INTEGER(PREC_VERTEXIDX) :: nnodes,isize
      INTEGER(I32), DIMENSION(:), POINTER :: p_InodalProperty
      
      ! Calculate the number of nodes in the mesh.
      ! This is: #vertices 
      !         +#edges (as every edge generates a new vertex)
      !         +#quads (as every quad generates a midpoint)
      nnodes = rtriangulation%NVT + rtriangulation%NMT + &
          rtriangulation%InelOfType(TRIA_NVEQUAD2D)
          
      ! Reallocate the memory if necessary.
      ! Copy the old content as we mustn't destroy the old nodal 
      ! property tags of the vertices.
      ! New elements are filled with zero = specify inner vertices.
      CALL storage_getsize (rtriangulation%h_InodalProperty, isize)
      IF (isize .LT. nnodes) THEN
        CALL storage_realloc ('tria_genEdgeNodalProperty2D', &
            nnodes, rtriangulation%h_InodalProperty, &
            ST_NEWBLOCK_NOINIT, .TRUE.)
        
        IF (rtriangulation%InelOfType(TRIA_NVEQUAD2D) .NE. 0) THEN
          ! Fill the last InelOfType(TRIA_NVEQUAD2D) entries of
          ! the array by 0. These elements will generate the
          ! nodal property for the element midpoints.
          ! Edge midpoints are recomputed anyway, so we don't
          ! have to initialise that part of the array!
          CALL storage_getbase_int (&
            rtriangulation%h_InodalProperty,p_InodalProperty)
          CALL lalg_clearVectorInt (&
            p_InodalProperty(nnodes-rtriangulation%InelOfType(TRIA_NVEQUAD2D)+1:))
        END IF
        
        reallocRefinedInodalProperty2D = .TRUE.
      ELSE
        reallocRefinedInodalProperty2D = .FALSE.
      END IF
      
    END FUNCTION
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_refine2LevelOrdering(&
      rsourceTriangulation,rdestTriangulation,rboundary,cflags)

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
  
  ! OPTIONAL: Bitfield of TRIA_R2LV_xxxx constants that allow to specify
  ! options for the refinement. If not specified, TRIA_R2LV_STANDARD
  ! is used as default.
  INTEGER(I32), INTENT(IN), OPTIONAL :: cflags
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
    INTEGER(I32) :: cflagsAct
    
    cflagsAct = TRIA_R2LV_STANDARD
    IF (PRESENT(cflags)) cflagsAct = cflags
    
    ! Call the correct submethod depending on the dimension.
    SELECT CASE (rsourceTriangulation%ndim)
    CASE (NDIM1D)
      ! Refine the basic mesh
      CALL tria_refineMesh2lv1D(rsourceTriangulation,rdestTria)
    
    CASE (NDIM2D)
      ! Refine the basic mesh
      CALL tria_refineMesh2lv2D(rsourceTriangulation,rdestTria)
      
      ! Refine the boundary
      CALL tria_refineBdry2lv2D(rsourceTriangulation,rdestTria,&
          IAND(cflagsAct,TRIA_R2LV_RECALCCOORDSONBD) .NE. 0,rboundary)
      
      IF (IAND(cflagsAct,TRIA_R2LV_AVERAGEMIDPOINTS) .NE. 0) THEN
        ! Recalculate corner points of quads that were element midpoints
        ! on the coarse mesh.
        CALL tria_averageMidpoints2D(rsourceTriangulation,rdestTria)
      END IF
      
    CASE (NDIM3D)
      ! refine the basic mesh
      call tria_refineMesh2lv3D(rsourceTriangulation, rdestTria)
      
      ! Refine the boundary
      call tria_refineBdry2lv3D(rsourceTriangulation,rdestTria,rboundary)
      
      
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
  
    SUBROUTINE tria_refineMesh2lv1D(rsourceTriangulation,rdestTriangulation)

    ! This routine refines the given 1D mesh rsourceTriangulation according to
    ! the 2-level ordering algorithm. The refined mesh is saved in 
    ! rdestTriangulation. There will be no correction of boundary
    ! vertices. Boundary parameter values are not handled here!

    ! The source triangulation to be refined
    TYPE(t_triangulation), INTENT(INOUT) :: rsourceTriangulation

    ! Destination triangulation structure that receives the refined mesh. 
    TYPE(t_triangulation), INTENT(OUT) :: rdestTriangulation
    
      ! local variables
      REAL(DP), DIMENSION(:,:), POINTER :: p_DcoordSource
      REAL(DP), DIMENSION(:,:), POINTER :: p_DcoordDest
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertAtElementSource
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertAtElementDest
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IrefinementPatchIndex
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IrefinementPatch
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IcoarseGridElement
      INTEGER(I32), DIMENSION(:), POINTER :: p_InodalPropSource, p_InodalPropDest
      INTEGER(PREC_ELEMENTIDX) :: iel,iel2
      INTEGER(PREC_VERTEXIDX) :: ivt1,ivt2, ivtoffset, ivt
      INTEGER(I32), DIMENSION(2) :: Isize
      INTEGER :: nnve, ive
      REAL(DP) :: x,y
      
      ! Get the arrays with information of the source mesh.
      CALL storage_getbase_double2D (rsourceTriangulation%h_DvertexCoords,&
          p_DcoordSource)
      CALL storage_getbase_int2D (rsourceTriangulation%h_IverticesAtElement,&
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
      Isize = (/NDIM1D,INT(rdestTriangulation%NVT,I32)/)
      CALL storage_new2D ('tria_refineMesh2lv1D', 'DCORVG', Isize, ST_DOUBLE, &
          rdestTriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double2D(&
          rdestTriangulation%h_DvertexCoords,p_DcoordDest)
      
      ! Allocate memory for the refinement information arrays.
      ! These arrays define for every coarse grid element the fine
      ! grid elements and for every fine grid element the coarse
      ! grid element where it comes from.
      CALL storage_new ('tria_refineMesh2lv1D', 'h_IrefinementPatchIndex', &
          rsourceTriangulation%NEL+1, ST_INT, &
          rdestTriangulation%h_IrefinementPatchIndex, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatchIndex,p_IrefinementPatchIndex)

      CALL storage_new ('tria_refineMesh2lv1D', 'h_IrefinementPatch', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IrefinementPatch, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatch,p_IrefinementPatch)
    
      CALL storage_new ('tria_refineMesh2lv1D', 'h_IcoarseGridElement', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IcoarseGridElement, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(&
          rdestTriangulation%h_IcoarseGridElement,p_IcoarseGridElement)

      ! The p_IrefinementPatchIndex array can directly be initialised
      ! as every coarse grid element gets 4 fine grid elements.
      DO iel = 0,rsourceTriangulation%NEL
        p_IrefinementPatchIndex(1+iel) = 1 + iel*2
      END DO
      
      ! Also the p_IcoarseGridElement array can directly be initialised.
      ! The coarse grid element number of element 1..NEL(coarse) stays
      ! the same. The number of the other coarse grid elements can
      ! be calculated by a formula.
      DO iel = 1,rsourceTriangulation%NEL
        p_IcoarseGridElement(iel) = iel
      END DO

      DO iel = rsourceTriangulation%NEL+1,rdestTriangulation%NEL
        p_IcoarseGridElement(iel) = iel-rsourceTriangulation%NEL
      END DO
    
      ! Ok, let's start the refinement. In the first step, we copy the
      ! corner coordinates of the coarse mesh to the fine mesh; they
      ! don't change during the refinement.
      !$OMP PARALLEL DO
      DO ivt=1,rsourceTriangulation%NVT
        p_DcoordDest(1,ivt) = p_DcoordSource(1,ivt)
      END DO
      !$OMP END PARALLEL DO
      
      ! Each line produces an line midpoint which is stored as new
      ! point in the fine mesh. To calculate the coordinates, take
      ! the mean of the coordinates in the coarse mesh.
      ivtoffset = rsourceTriangulation%NVT
      !$OMP PARALLEL DO PRIVATE(ivt1,ivt2)
      DO iel=1, rsourceTriangulation%NEL
        ivt1 = p_IvertAtElementSource (1,iel)
        ivt2 = p_IvertAtElementSource (2,iel)
        p_DcoordDest(1,ivtoffset+iel) = &
            0.5_DP * ( p_DcoordSource (1,ivt1) + p_DcoordSource (1,ivt2) )
      END DO
      !$OMP END PARALLEL DO
      
      ! Allocate memory for IverticesAtElement and get a pointer to it.
      ! Fill the array with zero, so we won't have problems when mixing
      ! triangles into a quad mesh.
      Isize = (/nnve,INT(rdestTriangulation%NEL,I32)/)
      CALL storage_new2D ('tria_refineMesh2lv1D', 'KVERT', Isize, ST_INT, &
          rdestTriangulation%h_IverticesAtElement, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int2D(&
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
      
      !$OMP PARALLEL DO PRIVATE(iel2)
      DO iel = 1,rsourceTriangulation%NEL
    
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
      
      END DO
      !$OMP END PARALLEL DO
      
      ! Create a new nodal property array
      CALL storage_new('tria_refineMesh2lv1D', 'h_InodalProperty',&
          rdestTriangulation%NVT, ST_INT,&
          rdestTriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
      
      ! Get the nodal property arrays
      CALL storage_getbase_int(rsourceTriangulation%h_InodalProperty, &
                               p_InodalPropSource)
      CALL storage_getbase_int(rdestTriangulation%h_InodalProperty, &
                               p_InodalPropDest)
      
      ! Copy the nodal property of the coarse mesg
      DO ivt=1, rsourceTriangulation%NVT
        p_InodalPropDest(ivt) = p_InodalPropSource(ivt)
      END DO
      
      ! We also need to copy the boundary vertices array to the finer level.
      rdestTriangulation%NVBD = rsourceTriangulation%NVBD
      CALL storage_copy (rsourceTriangulation%h_IverticesAtBoundary,&
          rdestTriangulation%h_IverticesAtBoundary)
      
      rdestTriangulation%NBCT = rsourceTriangulation%NBCT
      CALL storage_copy (rsourceTriangulation%h_IboundaryCpIdx,&
          rdestTriangulation%h_IboundaryCpIdx)
    
    END SUBROUTINE
    
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
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IrefinementPatchIndex
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IrefinementPatch
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IcoarseGridElement
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
      nnve = rsourceTriangulation%NNVE
      
      IF ((nnve .LT. TRIA_NVETRI2D) .OR. (nnve .GT. TRIA_NVEQUAD2D)) THEN
      
        CALL output_line (&
            '2-level refinement supports only triangular and quad meshes!', &
            OU_CLASS_ERROR,OU_MODE_STD,'tria_refineMesh2lv2D')
        CALL sys_halt()
        
      END IF
      
      ! Initialise the basic mesh data in rdestTriangulation:
      
      ! 2D mesh
      rdestTriangulation%ndim = NDIM2D
      
      rdestTriangulation%nnve = nnve
      
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
      !$OMP PARALLEL DO
      DO ivt=1,rsourceTriangulation%NVT
        p_DcoordDest(1,ivt) = p_DcoordSource(1,ivt)
        p_DcoordDest(2,ivt) = p_DcoordSource(2,ivt)
      END DO
      !$OMP END PARALLEL DO
      
      ! Each edge produces an edge midpoint which is stored as new
      ! point in the fine mesh. To calculate the coordinates, take
      ! the mean of the coordinates in the coarse mesh.
      ivtoffset = rsourceTriangulation%NVT
      !$OMP PARALLEL DO PRIVATE(ivt1,ivt2)
      DO imt=1,rsourceTriangulation%NMT
        ivt1 = p_IvertAtEdgeSource (1,imt)
        ivt2 = p_IvertAtEdgeSource (2,imt)
        p_DcoordDest(1,ivtoffset+imt) = &
            0.5_DP * ( p_DcoordSource (1,ivt1) + p_DcoordSource (1,ivt2) )
        p_DcoordDest(2,ivtoffset+imt) = &
            0.5_DP * ( p_DcoordSource (2,ivt1) + p_DcoordSource (2,ivt2) )
      END DO
      !$OMP END PARALLEL DO
      
      ! Allocate memory for IverticesAtElement and get a pointer to it.
      ! Fill the array with zero, so we won't have problems when mixing
      ! triangles into a quad mesh.
      Isize = (/nnve,INT(rdestTriangulation%NEL,I32)/)
      CALL storage_new2D ('tria_refineMesh2lv2D', 'KVERT', Isize, ST_INT, &
          rdestTriangulation%h_IverticesAtElement, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int2D(&
          rdestTriangulation%h_IverticesAtElement,p_IvertAtElementDest)

      ! Allocate memory for the refinement information arrays.
      ! These arrays define for every coarse grid element the fine
      ! grid elements and for every fine grid element the coarse
      ! grid element where it comes from.
      CALL storage_new ('tria_refineMesh2lv2D', 'h_IrefinementPatchIndex', &
          rsourceTriangulation%NEL+1, ST_INT, &
          rdestTriangulation%h_IrefinementPatchIndex, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatchIndex,p_IrefinementPatchIndex)

      CALL storage_new ('tria_refineMesh2lv2D', 'h_IrefinementPatch', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IrefinementPatch, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(&
          rdestTriangulation%h_IrefinementPatch,p_IrefinementPatch)
    
      CALL storage_new ('tria_refineMesh2lv2D', 'h_IcoarseGridElement', &
          rdestTriangulation%NEL, ST_INT, &
          rdestTriangulation%h_IcoarseGridElement, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(&
          rdestTriangulation%h_IcoarseGridElement,p_IcoarseGridElement)
    
      ! The p_IrefinementPatchIndex array can directly be initialised
      ! as every coarse grid element gets 4 fine grid elements.
      DO iel = 0,rsourceTriangulation%NEL
        p_IrefinementPatchIndex(1+iel) = 1 + iel*4
      END DO
      
      ! Also the p_IcoarseGridElement array can directly be initialised.
      ! The coarse grid element number of element 1..NEL(coarse) stays
      ! the same. The number of the other coarse grid elements can
      ! be calculated by a formula.
      DO iel = 1,rsourceTriangulation%NEL
        p_IcoarseGridElement(iel) = iel
      END DO

      DO iel = rsourceTriangulation%NEL+1,rdestTriangulation%NEL
        p_IcoarseGridElement(iel) = (iel-rsourceTriangulation%NEL-1)/3 + 1
      END DO
    
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
      ! Let's loop over the vertices on the coarse grid. Everyone produces
      ! four elements on the fine grid:
      !
      ! In nquads we count the number of quads +NVT+NMT we reach.
      ! That's the number of the midpoint of that element!
      nquads = rsourceTriangulation%NVT+rsourceTriangulation%NMT
      
      IF (nnve .EQ. TRIA_NVETRI2D) THEN
      
        ! Pure triangle mesh
        !$OMP PARALLEL DO PRIVATE(iel1,iel2,iel3)
        DO iel = 1,rsourceTriangulation%NEL
      
          ! Determine number of subelements.
          iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
          iel2 = iel1+1
          iel3 = iel1+2
          
          ! Save them
          p_IrefinementPatch(1+(iel-1)*4) = iel
          p_IrefinementPatch(1+(iel-1)*4+1) = iel1
          p_IrefinementPatch(1+(iel-1)*4+2) = iel2
          p_IrefinementPatch(1+(iel-1)*4+3) = iel3
          
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
        !$OMP END PARALLEL DO
      
      ELSE IF (nquads .EQ. rsourceTriangulation%NEL) THEN
        
        ! Pure QUAD mesh
        !$OMP PARALLEL DO PRIVATE(iel1,iel2,iel3,nquads)
        DO iel = 1,rsourceTriangulation%NEL
      
          ! Determine number of subelements.
          iel1 = rsourceTriangulation%NEL+3*(iel-1)+1
          iel2 = iel1+1
          iel3 = iel1+2
          
          ! Save them
          p_IrefinementPatch(1+(iel-1)*4) = iel
          p_IrefinementPatch(1+(iel-1)*4+1) = iel1
          p_IrefinementPatch(1+(iel-1)*4+2) = iel2
          p_IrefinementPatch(1+(iel-1)*4+3) = iel3

          ! In nquads we count the number of the quad +NVT+NMT we process.
          ! That's the number of the midpoint of that element!
          ! As we reached a new quad, we increase nquads
          
          !nquads = nquads+1 !replaced because of OpenMP directive
          nquads = rsourceTriangulation%NVT + rsourceTriangulation%NMT + iel
          
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
        !$OMP END PARALLEL DO
        
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
            
            ! Save them
            p_IrefinementPatch(1+(iel-1)*4) = iel
            p_IrefinementPatch(1+(iel-1)*4+1) = iel1
            p_IrefinementPatch(1+(iel-1)*4+2) = iel2
            p_IrefinementPatch(1+(iel-1)*4+3) = iel3

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
            
            ! Save them
            p_IrefinementPatch(1+(iel-1)*4) = iel
            p_IrefinementPatch(1+(iel-1)*4+1) = iel1
            p_IrefinementPatch(1+(iel-1)*4+2) = iel2
            p_IrefinementPatch(1+(iel-1)*4+3) = iel3

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
  
    SUBROUTINE tria_refineBdry2lv2D(rsourceTriangulation,rdestTriangulation,&
        brecalcBoundaryCoords,rboundary)

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
    
    ! Recalculate the coordinates of all boundary vertices according to
    ! their parameter value.
    LOGICAL, INTENT(IN) :: brecalcBoundaryCoords
    
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
      
      !$OMP PARALLEL DO
      DO ibct = 1,SIZE(p_IboundaryCpIdxDest)
        p_IboundaryCpIdxDest(ibct) = 2*(p_IboundaryCpIdxSource(ibct)-1) + 1
      END DO
      !$OMP END PARALLEL DO
      !$OMP PARALLEL DO
      DO ivbd = 0,SIZE(p_IvertAtBoundartySource)-1
        p_IvertAtBoundartyDest(1+2*ivbd) = p_IvertAtBoundartySource(1+ivbd)
        p_IvertAtBoundartyDest(2+2*ivbd) = p_IedgesAtBoundartySource(1+ivbd)
      END DO
      !$OMP END PARALLEL DO
      
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
            
        !$OMP PARALLEL DO
        DO ivbd = 0,SIZE(p_IvertAtBoundartySource)-1
          p_DvertParamsDest(1+2*ivbd) = p_DvertParamsSource(1+ivbd)
          p_DvertParamsDest(2+2*ivbd) = p_DedgeParamsSource(1+ivbd)
        END DO
        !$OMP END PARALLEL DO
        
        ! If the analytic boundary is given, compute the coordinates of the
        ! boundary vertices from that.
        IF (PRESENT(rboundary) .AND. brecalcBoundaryCoords) THEN
          
          ! Get the array with the vertex coordinates.
          ! We want to correct the coordinates of the boundary points
          ! according to the analytic boundary.
          CALL storage_getbase_double2d (rdestTriangulation%h_DvertexCoords,&
              p_DcornerCoordDest)
            
          ! Loop through the boundary points and calculate the correct
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

    ! ---------------------------------------------------------------
  
    SUBROUTINE tria_averageMidpoints2D(rsourceTriangulation,rdestTriangulation)

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
    TYPE(t_triangulation), INTENT(IN) :: rsourceTriangulation

    ! Destination triangulation that specifies the fine mesh after
    ! 2-level refinement
    TYPE(t_triangulation), INTENT(INOUT) :: rdestTriangulation

      ! local variables
      INTEGER :: i
      REAL(DP) :: dx,dy
      INTEGER(PREC_ELEMENTIDX) :: iel,ivtoffset !,ielbd
      !INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtBoundaryCoarse
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementCoarse
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElementFine
      INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElementCoarse
      REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoordsFine
      INTEGER(PREC_VERTEXIDX) :: ipt

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
      ! twice that way, but we don't care...
      !
      ! Unfortunately, due to the incompleteness of the triangulation structure,
      ! this approach does not work for mixed mesh. We have to loop through
      ! all elements on the coarse mesh to find the numbers of the quads :-(
      
      CALL storage_getbase_int2d (rsourceTriangulation%h_IverticesAtElement,&
          p_IverticesAtElementCoarse)
          
      ! If this is a pure triangle mesh, there is nothing to do.
      IF (UBOUND(p_IverticesAtElementCoarse,1) .LE. TRIA_NVETRI2D) RETURN
      
      ! CALL storage_getbase_int (rsourceTriangulation%h_IelementsAtBoundary,&
      !     p_IelementsAtBoundaryCoarse)

      CALL storage_getbase_int2d (rsourceTriangulation%h_IedgesAtElement,&
          p_IedgesAtElementCoarse)
      CALL storage_getbase_int2d (rdestTriangulation%h_IverticesAtElement,&
          p_IverticesAtElementFine)
      CALL storage_getbase_double2d (rdestTriangulation%h_DvertexCoords,&
          p_DvertexCoordsFine)

      ivtoffset = rsourceTriangulation%NVT + rsourceTriangulation%NMT
      
      DO iel = 1,rsourceTriangulation%NEL

        ! Get the element number of the coarse mesh element
        ! iel = p_IelementsAtBoundaryCoarse(ielbd)
      
        ! Triangle or quad?
        IF (p_IverticesAtElementCoarse(TRIA_NVEQUAD2D,iel) .NE. 0) THEN
        
          ! New quad. Increase ivtoffset, this is then the number of the
          ! point that was the midpoint in the coarse mesh.
          ivtoffset = ivtoffset + 1
        
          ! Ok, that's a quad. The midpoint must be recalculated.
          ! We do this based on the 'edge midpoints' which are now
          ! vertices in the fine mesh. The coordinates of these points
          ! may differ from the coordinates calculated by linear 
          ! interpolation due to boundary adjustment!
          dx = 0.0_DP
          dy = 0.0_DP
          DO i=1,TRIA_NVEQUAD2D
            ipt = p_IedgesAtElementCoarse(i,iel)
            dx = dx + p_DvertexCoordsFine(1,ipt)
            dy = dy + p_DvertexCoordsFine(2,ipt)
          END DO
          
          ! Save the vertex coordinates. 
          p_DvertexCoordsFine(1,ivtoffset) = 0.25_DP*dx
          p_DvertexCoordsFine(2,ivtoffset) = 0.25_DP*dy
        
        END IF
      
      END DO
    
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE tria_searchBoundaryNode(inode,rtriangulation,iindex)

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

!<output>
  ! If inode is a boundary vertex: The index of the inode in IboundaryVertexPos.
  ! If inode is a boundary edge: The index of the inode in IboundaryEdgePos.
  ! =0 if inode was not found (e.g. because inode is not on the boundary e.g.).
  INTEGER, INTENT(OUT) :: iindex
!</output>

!</subroutine>

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
      iindex = 0
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
      iindex = p_InodePos(2,ileft)
    ELSE IF (p_InodePos(1,iright) .EQ. inode) THEN
      ! Return the index in the node array.
      iindex = p_InodePos(2,iright)
    ELSE
      DO WHILE (ileft .LT. iright)
        ipos = (ileft+iright)/2
        IF (p_InodePos(1,ipos) .GT. inode) THEN
          iright = ipos
        ELSE IF (p_InodePos(1,ipos) .LT. inode) THEN
          ileft = ipos
        ELSE
          ! We found the node. Return the index in the node array.
          iindex = p_InodePos(2,ipos)
          EXIT
        END IF
      END DO
    END IF
    
  END SUBROUTINE

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
    !$OMP PARALLEL DO PRIVATE(ipointpos,ipoint1,ipoint2,ipoint,idim)
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
    !$OMP END PARALLEL DO
    
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

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tria_infoStatistics (rtriangulation,bheadline,ilevel)
  
!<description>
  ! Prints out statistical information of the given triangulation to the
  ! terminal. The output is formatted as a table with an optional headline
  ! and an optional level identifier.
!</description>

!<input>
  ! Triangulation structure.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation
  
  ! OPTIONAL: Print out a headline above the statistical data.
  ! =FALSE: don't print = standard.
  ! =TRUE: print a headline.
  LOGICAL, INTENT(IN), OPTIONAL :: bheadline
  
  ! OPTIONAL: Level identifier.
  ! If specified, an additional column 'Level' is added to the front of
  ! the statistics table. ilevel is printed to this column.
  INTEGER, INTENT(IN), OPTIONAL :: ilevel
!</input>

!</subroutine>

    SELECT CASE (rtriangulation%NDIM)
    CASE (NDIM2D)
      IF (PRESENT(bheadline)) THEN
        IF (bheadline) THEN
          ! Print a headline
          IF (PRESENT(ilevel)) THEN
            CALL output_line(&
              'Lv. dim.       NVT        NMT        NEL    NBCT' &
            //'    NVBD     #trias      #quads')
            CALL output_line(&
              '------------------------------------------------' &
            //'-------------------------------')
          ELSE
            CALL output_line(&
              'dim.       NVT        NMT        NEL    NBCT' &
            //'    NVBD     #trias      #quads')
            CALL output_line(&
              '--------------------------------------------' &
            //'-------------------------------')
          END IF
        END IF
      END IF

      ! Print out the statistics
      IF (PRESENT(ilevel)) CALL output_line (TRIM(sys_si(ilevel,3))//' ',&
          bnolinebreak=.TRUE.)
      
      CALL output_line (&
          TRIM(sys_si(rtriangulation%NDIM,4)) &
        //TRIM(sys_si(rtriangulation%NVT,11)) &
        //TRIM(sys_si(rtriangulation%NMT,11)) &
        //TRIM(sys_si(rtriangulation%NEL,11)) &
        //TRIM(sys_si(rtriangulation%NBCT,8)) &
        //TRIM(sys_si(rtriangulation%NVBD,8)) &
        //TRIM(sys_si(rtriangulation%InelOfType(TRIA_NVETRI2D),11)) &
        //TRIM(sys_si(rtriangulation%InelOfType(TRIA_NVEQUAD2D),11)) )
    END SELECT

  END SUBROUTINE

  !************************************************************************

!<subroutine>

  SUBROUTINE tria_exportTriFile(rtriangulation, sfilename, ctriFormat)

!<description>
    ! This routine exports a triangulation into a .TRI file.
!</description>

!<input>
    ! Triangulation structure, to be exported
    TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation

    ! The name of the .tri file to write.
    CHARACTER(LEN=*), INTENT(IN) :: sfilename
    
    ! OPTIONAL: Format tag of the TRI file to export.
    ! TRI_FMT_STANDARD: Standard TRI file format, compatible to FEAT1.
    !    Vertex coordinates of boundary vertices are in 2D replaced
    !    by parameter values.
    ! TRI_FMT_NOPARAMETRISATION: Standard TRI file format, but the 
    !    vertex coordinates are exported 'as they are', not as parameter
    !    values.
    ! If not specified, TRI_FMT_STANDARD is assumed.
    INTEGER(I32), INTENT(IN), OPTIONAL :: ctriFormat
!</input>
!</subroutine>

    ! Depending on the dimension of the triangulation, call the corresponding
    ! export routine!
    SELECT CASE(rtriangulation%ndim)
    CASE (NDIM1D)
      CALL output_line ('1D TRI file export not implemented!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_exportTriFile')
      CALL sys_halt()
      
    CASE (NDIM2D)
      CALL tria_exportTriFile2D(rtriangulation, sfilename, ctriFormat)
      
    CASE (NDIM3D)
      CALL output_line ('3D TRI file export not implemented!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_exportTriFile')
      CALL sys_halt()
      
    CASE DEFAULT
      CALL output_line ('Triangulation structure not properly initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_exportTriFile')
      CALL sys_halt()
    END SELECT

  END SUBROUTINE tria_exportTriFile

  !************************************************************************

!<subroutine>

  SUBROUTINE tria_exportTriFile2D(rtriangulation, sfilename, ctriFormat)

!<description>
    ! Auxiliary routine. This routine exports a 2D triangulation into a .TRI file.
!</description>

!<input>
    ! Triangulation structure, to be exported
    TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation

    ! The name of the .tri file to write.
    CHARACTER(LEN=*), INTENT(IN) :: sfilename

    ! OPTIONAL: Format tag of the TRI file to export.
    ! TRI_FMT_STANDARD: Standard TRI file format, compatible to FEAT1.
    !    Vertex coordinates of boundary vertices are in 2D replaced
    !    by parameter values.
    ! TRI_FMT_NOPARAMETRISATION: Standard TRI file format, but the 
    !    vertex coordinates are exported 'as they are', not as parameter
    !    values.
    ! If not specified, TRI_FMT_STANDARD is assumed.
    INTEGER(I32), INTENT(IN), OPTIONAL :: ctriFormat
!</input>
!</subroutine>

    ! Local variables
    REAL(DP), DIMENSION(:,:), POINTER     :: p_Ddata2D
    REAL(DP), DIMENSION(:), POINTER     :: p_Ddata
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_Idata2D
    INTEGER(I32), DIMENSION(:), POINTER   :: p_Idata
    INTEGER(I32), DIMENSION(:), POINTER   :: p_IverticesAtBoundary
    INTEGER(I32), DIMENSION(:), POINTER   :: p_InodalProperty
    INTEGER(I32), DIMENSION(:), POINTER   :: p_IboundaryCpIdx
    CHARACTER(SYS_STRLEN) :: ckmmstr
    INTEGER(I32) :: ivt, iel, ivbd
    INTEGER      :: idim, ive, ibct
    INTEGER      :: iunit
    LOGICAL      :: bnoParameters
    
    bnoParameters = .FALSE.
    IF (PRESENT(ctriFormat)) &
      bnoParameters = ctriFormat .EQ. TRI_FMT_NOPARAMETRISATION
    
    ! Open the file
    CALL io_openFileForWriting(sfilename, iunit, SYS_REPLACE)

    ! Comment: Header
    WRITE (iunit,*) 'Coarse mesh exported by FeatFlow2 exporter'
    IF (.NOT. bnoParameters) THEN
      WRITE (iunit,*) 'Parametrisation PARXC, PARYC, TMAXC'
      !WRITE (iunit,*) 'Parametrisierung PARXC, PARYC, TMAXC'
    ELSE
      WRITE (iunit,*) 'No parametrisation'
    END IF

    ! Write NEL,NVT,NMT,NVE,NBCT to the file
    WRITE (iunit,*) rtriangulation%NEL,rtriangulation%NVT,rtriangulation%NMT,&
        rtriangulation%NNVE,rtriangulation%NBCT,'NEL NVT NMT NVE NBCT'

    ! Write: 'DCORVG'
    WRITE (iunit,*) 'DCORVG'

    ! Get the pointers to the coordinate array
    CALL storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
    
    ! Write the data to the file
    !
    IF (.NOT. bnoParameters) THEN
      ! Standard FEAT1 format. For boundary vertices, the parameter value is
      ! exported instead of the vertex coordinate!
      CALL storage_getbase_int(&
          rtriangulation%h_InodalProperty,p_InodalProperty)
      CALL storage_getbase_double(&
          rtriangulation%h_DvertexParameterValue,p_Ddata)
      
      DO ivt = 1, rtriangulation%NVT
        IF (p_InodalProperty(ivt) .GT. 0) THEN
          ! Get the index of the vertex in the IverticesAtBoundary array.
          CALL tria_searchBoundaryNode(ivt,rtriangulation,ivbd)
          
          ! Write the parameter value. 2nd entry is =0.
          WRITE (iunit,*) p_Ddata(ivbd),0.0_DP
        ELSE
          WRITE (iunit,*) (p_Ddata2D(idim,ivt),idim=1,NDIM2D)
        END IF
      END DO
    ELSE
      ! Parameterless format.
      ! The coordinates of the vertices are exportes 'as they are' independent
      ! of whether they are on the boundary or not.
      DO ivt = 1, rtriangulation%NVT
        WRITE (iunit,*) (p_Ddata2D(idim,ivt),idim=1,NDIM2D)
      END DO
    END IF

    ! Write: 'KVERT'
    WRITE (iunit,*) 'KVERT'

    ! Get the pointer to the IverticesAtElement array and read the array
    CALL storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)

    ! Write the data to the file
    DO iel = 1, rtriangulation%NEL
      WRITE (iunit,*) (p_Idata2D(ive,iel),ive=1,SIZE(p_Idata2D,1))
    END DO

    ! Write: 'KNPR'
    WRITE (iunit,*) 'KNPR'

    ! Get the pointer to the InodalProperty array
    CALL storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_Idata)
    
    ! Write the data
    DO ivt = 1, rtriangulation%NVT
      WRITE (iunit,*) p_Idata(ivt)
    END DO

    ! Write: 'KMM'
    WRITE (iunit,*) 'KMM'
    
    ! Get the pointer to the IboundaryCpIdx and IverticesAtBoundary arrays
    CALL storage_getbase_int(&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    CALL storage_getbase_int(&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)

    ! Write the data
    ckmmstr = ''
    DO ibct = 1, rtriangulation%NBCT
      ckmmstr = ADJUSTL(TRIM(ckmmstr))//' '//&
          TRIM(sys_siL(p_IverticesAtBoundary(p_IboundaryCpIdx(ibct)),10))//' '//&
          TRIM(sys_siL(p_IverticesAtBoundary(p_IboundaryCpIdx(ibct+1)-1),10))
    END DO
    WRITE (iunit,*) TRIM(ckmmstr)

    ! Close the file, finish
    CLOSE(iunit)
    
  END SUBROUTINE tria_exportTriFile2D
  
!====================================================================
!
!       ++++      ++++
!           +     +   +  
!           +     +    +
!        +++      +    +
!           +     +    + 
!           +     +   +  
!       ++++      ++++
!           
!====================================================================

!<subroutine>
  subroutine tria_readTriFile3D(rtriangulation, sfilename, rboundary)

!<description>
  ! This routine reads a .TRI file of a 3D triangulation into memory
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
  ! The name of the .tri file to read.
  character(len=*), intent(in) :: sfilename

  ! OPTIONAL: An rboundary object specifying the underlying domain.
  ! If not specified, the routine assumes that the TRI file does not specify
  ! boundary parameter values, i.e. the point coordinates in the TRI file
  ! are all real coordinates. The array DvertexParameterValue is not
  ! generated in this case.
  type(t_boundary), intent(in), optional :: rboundary
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

    ! we create a 3d triangulation here.
    rtriangulation%ndim = NDIM3D

    ! read the basic mesh
    call tria_readRawTriangulation3d(iunit,rtriangulation)

    ! create the basic boundary information
    call tria_genRawBoundary3d (rtriangulation,rboundary)

    ! close the file, finish
    close(iunit)
  
  end subroutine
  
  !************************************************************************  
!<subroutine>  
  subroutine tria_readRawTriangulation3d(iunit,rtriangulation)
  
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
  real(dp), dimension(:,:), pointer :: p_Ddata3d
  integer(i32), dimension(:,:), pointer :: p_Idata3d
  integer(i32), dimension(:), pointer :: p_Idata    
  
  ! integer mesh parameters
  integer(i32) :: NVT, NEL, NBCT, NVE, NEE, NAE
  
  integer(i32), dimension(2) :: Isize
  
  ! some counter variables
  integer(i32) :: idim, ivt,ive,iel
  
  
  ! start to read the tri file the first two lines are comments
  read(iunit,*)
  read(iunit,*)

  ! read variables
  read(iunit,*) NEL, NVT, NBCT, NVE, NEE, NAE

  ! assign the variables in the structure    
  rtriangulation%NVT  = NVT
  rtriangulation%NEL  = NEL
  rtriangulation%NBCT = NBCT
  rtriangulation%NNEE  = NEE
  rtriangulation%NNAE  = NAE
  rtriangulation%NNVE  = NVE
  
  ! Vertices per face
  IF (NAE .EQ. 4) THEN
    rtriangulation%NNVA = 3
  ELSE
    rtriangulation%NNVA = 4
  END IF
  
  ! skip Comment: 'DCORVG'
  read (iunit,*) 
  ! allocate array of coordinates
  ! allocate(p_ddata3d(3,NVT))
  Isize = (/NDIM3D, INT(rtriangulation%NVT, i32)/)
  
  call storage_new2D('tria_readRawTriangulation3d', 'DCORVG', &
       Isize, ST_DOUBLE, rtriangulation%h_DvertexCoords, &
       ST_NEWBLOCK_NOINIT)
       
  ! get a pointer to the allocated memory, store it in p_Ddata3d       
  call storage_getbase_double2D(&
       rtriangulation%h_DvertexCoords, p_Ddata3d)
       
  ! Read the data from the file, store it in the array.
  ! read data into p_Ddata3d :
  ! first read nvt x-coordinates into p_Ddata3d(1,ivt)
  ! then read nvt  y-coordinates into p_Ddata3d(2,ivt)
  ! then read nvt  z-coordinates into p_Ddata3d(3,ivt)
  read (iunit,*) ((p_Ddata3d(idim,ivt),idim=1,NDIM3D), ivt=1,NVT)
    
  ! skip Comment: 'KVERT'
  read (iunit,*)
  
  ! Allocate memory for IverticesAtElement
  ! build the old KVERT...
  ! 2d array of size(NVE, NEL)
  Isize = (/nve,int(nel,i32)/)
  call storage_new2D('tria_readRawTriangulation3d', 'KVERT',&
       Isize, ST_INT, rtriangulation%h_IverticesAtElement, &
       ST_NEWBLOCK_NOINIT)
       
  ! get a pointer to the memory we just allocated
  call storage_getbase_int2d(&
       rtriangulation%h_IverticesAtElement, p_Idata3d)
       
  ! read ive=1 indices to nve into p_Idata3D(ive,iel) where iel=1 to NEL     
  read (iunit,*) ((p_idata3d(ive,iel),ive=1,nve), iel=1,NEL)
    
  ! skip Comment: 'KNPR'
  read (iunit,*)            
  
  ! Allocate memory for InodalProperty
  call storage_new('tria_readRawTriangulation3d', 'KNPR', &
       INT(NVT,i32), ST_INT, &
       rtriangulation%h_InodalProperty, &
       ST_NEWBLOCK_ZERO)  
       
  ! get a pointer to the memory
  call storage_getbase_int(&
       rtriangulation%h_InodalProperty, p_Idata)
  
  ! Read the data   
  read (iunit,*) (p_idata(ivt),ivt=1,NVT)   
  
  ! done reading raw mesh data ...
  
  end subroutine
  
  !************************************************************************  
  
!<subroutine>  
  subroutine tria_genRawBoundary3d(rtriangulation,rboundary)
  
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
  real(dp), dimension(:,:), pointer :: p_DvertexCoords
  integer(prec_vertexidx), dimension(:), POINTER :: p_IboundaryCpIdx
  integer(prec_vertexidx), dimension(:), POINTER :: p_IverticesAtBoundary
  integer(prec_vertexidx) :: ivbd,ivt
  integer :: ibct
  integer(i32), dimension(:), pointer :: p_InodalProperty
  
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
    if(p_InodalProperty(ivt) /= 0) ivbd = ivbd + 1
  end do
  
  ! assign number of vertices on the boundary
  rtriangulation%NVBD = ivbd
  
  ! Allocate memory for IverticesAtBoundary.
  call storage_new('tri_genRawBoundary3d', &
       'KVBD', int(rtriangulation%NVBD, i32), &
       ST_INT, rtriangulation%h_IverticesAtBoundary, &
       ST_NEWBLOCK_NOINIT)
  
  ! allocate memory for the boundary compnent index vector and
  ! init with zeros
  call storage_new('tri_genRawBoundary3d', &
       'KBCT', int(rtriangulation%NBCT+1,i32), &
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
    if(p_InodalProperty(ivt) /= 0) then
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
    if(p_InodalProperty(ivt) /= 0) then
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
  
  
  end subroutine
  
!************************************************************************

!<subroutine>
    subroutine tria_genElementsAtVertex3D(rtriangulation)    
    
!<description>
  ! This routine generates the array IelementsAtVertex.
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement.
  ! If necessary, new memory is allocated.
!</description>
    
!<inputoutput>
  ! The triangulation structure to be updated.
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
    
!</subroutine>    
    ! local variables
    integer(PREC_ELEMENTIDX), dimension(:), pointer :: p_IelementsAtVertexIdx
    integer(PREC_ELEMENTIDX), dimension(:), pointer :: p_IelementsAtVertex
    integer(PREC_VERTEXIDX) , dimension(:,:), pointer :: p_idata3d    
    
    integer(PREC_ELEMENTIDX) :: iel
    
    integer :: ive, nnve, haux1
    
    integer(PREC_VERTEXIDX)  :: ivt, Isize, Isize2
    
    integer(i32), dimension(:), pointer :: p_Iaux1
    
    ! allocate memory for the p_IelementsAtVertexIdx array
    if(rtriangulation%h_IelementsAtVertexIdx == ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtVertex3D', 'IelementsAtVertexIdx', &
          int(rtriangulation%NVT+1,I32), ST_INT, &
          rtriangulation%h_IelementsAtVertexIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IelementsAtVertexIdx, isize)
      if (isize .NE. rtriangulation%NVT+1) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtVertex3D', &
            int(rtriangulation%NVT+1,I32), rtriangulation%h_IelementsAtVertexIdx, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      end if
    
    
    end if
    
   ! Get the index array.
    call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
        p_IelementsAtVertexIdx)
        
    ! Get some data arrays about the vertices.
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_idata3d)

    ! Fill the index array with zero.
    call storage_clear (rtriangulation%h_IelementsAtVertexIdx)
    
    ! shorthand to number of elements at vertex
    nnve = rtriangulation%NNVE
    
    ! first we calculate the number of elements at each vertex simply by counting
    
    ! loop over all elements
    do iel=1, rtriangulation%NEL
      ! loop over all vertices at the element
      do ive = 1, nnve
        ! ivt is the ive-th vertex at element iel
        ivt = p_idata3d(ive,iel)
        
        ! handle triangle case
        if(ivt .eq. 0) exit
        
        ! increase the number of elements by one
        p_IelementsAtVertexIdx(ivt+1) = p_IelementsAtVertexIdx(ivt+1) + 1
        
      end do ! end ive
      
    end do ! end iel
    
    ! set the first index to 1
    p_IelementsAtVertexIdx(1) = 1
    
    ! in the next step we sum up the number of elements at two
    ! successive vertices to create the index array
    ! thus at the penultimate position of p_IelementsAtVertexIdx
    ! we find the length of p_IelementsAtVertex
    do ivt = 2, rtriangulation%nvt+1
        p_IelementsAtVertexIdx(ivt) = p_IelementsAtVertexIdx(ivt) + p_IelementsAtVertexIdx(ivt-1)
    end do
    
    ! set the size
    Isize = p_IelementsAtVertexIdx(rtriangulation%NVT+1)-1
   
    
    ! Isize contains now the length of the array where we store the adjacency
    ! information (IelementsAtVertex).
    ! Do we have (enough) memory for that array?
    if (rtriangulation%h_IelementsAtVertex == ST_NOHANDLE) then
      call storage_new ('tria_genElementsAtVertex3D', 'IelementsAtVertex', &
          int(Isize,I32), ST_INT, &
          rtriangulation%h_IelementsAtVertex, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize (rtriangulation%h_IelementsAtVertex, Isize2)
      if (Isize /= Isize2) then
        ! If the size is wrong, reallocate memory.
        call storage_realloc ('tria_genElementsAtVertex3D', &
            int(isize,I32), rtriangulation%h_IelementsAtVertex, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      end if
    end if
    
    ! get the pointer to the array
    call storage_getbase_int(rtriangulation%h_IelementsAtVertex, p_IelementsAtVertex)
    
    
    ! Duplicate the p_IelementsAtVertexIdx array. We use that as pointer and index
    ! when new elements at a vertex are found.
    haux1 = ST_NOHANDLE
    call storage_copy (rtriangulation%h_IelementsAtVertexIdx,haux1)
    call storage_getbase_int (haux1,p_Iaux1)
    
    ! loop over all elements
    do iel = 1, rtriangulation%nel
        ! loop over all vertices of the element
        do ive = 1,nnve
        
          ! ivt is the ive-th vertex at element iel
          ivt = p_idata3d(ive,iel)
          
          ! handle triangle case
          if( ivt .eq. 0) exit

          ! store the adjacency information at position p_Iaux1(ivt)        
          p_IelementsAtVertex( p_Iaux1(ivt) ) = iel
          ! increase the position of the next element in p_Iaux1(ivt)
          p_Iaux1(ivt) = p_Iaux1(ivt) + 1
        
        end do ! end iel
        
    end do ! end ive
    
    call storage_free(haux1)
    
    end subroutine ! end tria_genElementsAtVertex3D

!************************************************************************    
!<subroutine>      
  recursive subroutine tria_mergesort(p_ConnectorList, l, r, pos)
   
!<description>
  ! this routine sorts a connector list
  ! it is used as an auxilliary routine
  ! during the Neighbours at elements routine
!</description>
    
!<input>
  ! the array positions l...r will be sorted
  ! the sorting key is element 'pos' of the connector
  integer :: l,r,pos
!</input>
  
!<inputoutput>
  ! the list of connectors    
  type(t_connector3d), dimension(:), pointer :: p_ConnectorList
!</inputoutput>

!</subroutine>      
    
  ! local variables
  integer :: m
    
      if(l < r) then
        
          m = l + (r-l)/2
            
          call tria_mergesort(p_ConnectorList,l,m, pos)
          call tria_mergesort(p_ConnectorList,m+1,r, pos)
          call tria_merge(p_ConnectorList,l,m,r,pos)
        
      end if
    
  end subroutine ! tria_mergesort
        
!************************************************************************      

!<subroutine> 
  subroutine tria_merge(p_ConnectorList, l, m, r, pos)
    
!<description>
  ! 
  ! standard auxilliary routine in the mergesort algorithm
  ! 
!</description>
    
!<input> 
  ! the array positions l...r will be sorted
  ! the sorting key is element 'pos' of the connector
  integer :: l,r,m,pos
!</input>

!<inputoutput>
  ! the list of connectors     
  type(t_connector3d), dimension(:), pointer :: p_ConnectorList
!</inputoutput>

!</subroutine>     

  ! local variables
  integer :: i,j,n1,n2,k
    
  type(t_connector3d), dimension(:), pointer :: p_L, p_R
    
  ! function body
  
    
  ! init counters
  n1 = m - l + 1
    
  n2 = r - m 
    
  k = l    
    
  ! allocate memory for merging
  allocate(p_L(n1))
  allocate(p_R(n2))
    
  ! fill left array
  do i=1,n1
      p_L(i) = p_ConnectorList(l+i-1)
  end do
    
  ! fill right array
  do j=1,n2
      p_R(j) = p_ConnectorList(m+j)
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
      if(p_L(i)%I_conData(pos) <= p_R(j)%I_conData(pos)) then
          p_ConnectorList(k) = p_L(i)
          i = i + 1
          k = k + 1
      else
          p_ConnectorList(k) = p_R(j)                 
        j = j + 1
        k = k + 1
    end if
    
  end do
    
  ! copy the remaining entries of p_L (if present)
  do
  if(i > n1) exit
    
      p_ConnectorList(k) = p_L(i)
      ! increment counters
      k = k + 1
      i = i + 1
        
  end do

  ! copy the remaining entries of p_R (if present)
  do
  if(j > n2) exit
    
      p_ConnectorList(k) = p_R(j)
      ! increment counters
      k = k + 1
      j = j + 1
        
  end do

  ! done merging
    
  ! free p_L and p_R
  deallocate(p_L)
  deallocate(p_R)
    
  end subroutine ! end tria_merge 
  
!************************************************************************   
      
!<subroutine>      
  subroutine tria_buildConnectorList(p_IConnectList, rtriangulation)

!<description>
  ! this routine builds the connector list used 
  ! in the neighbours at elements 
  ! routine
!</description>

!<input>
  type(t_triangulation) :: rtriangulation
!</input>
    
!<inputoutput>    
  ! the list of connectors this routine is supposed to build
  type(t_connector3d), dimension(:), pointer :: p_IConnectList    
!</inputoutput>  
  
!</subroutine>    

  ! local variables
  integer :: i,j,k, NVFACE
  
  integer, dimension(:,:), pointer :: p_idata3d
    
  ! function body
  
  ! Get some data arrays about the vertices.
  call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
      p_idata3d)
  
    
  ! allocate memory for NEL connectors
  allocate(p_IConnectList(rtriangulation%NEL*rtriangulation%NNAE))
    
  ! number of face
  NVFACE = rtriangulation%NNAE
    
  ! loop through all hexahedrons
  do i=1, rtriangulation%NEL
      ! build connectors for each hexahedron
      
      !=========================================================  
      ! first face
      j=1
      do k=1,4
          p_IConnectList( (i-1) * NVFACE + j)%I_conData(k) = &
          p_idata3d(k,i)
      end do
      ! save the number of the element this face was found from
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(k)=i
      ! assign the local face number
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(6)=1
      j=j+1
      !=========================================================
      ! sixth face
      do k=5,8
          p_IConnectList( (i-1) * NVFACE + j)%I_conData(k-4) = &
          p_idata3d(k,i)
      end do
        
      ! save the number of the element this face was found from
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(5)=i
      
      ! assign the local face number
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(6)=6
      j=j+1
      !=========================================================
      ! second face
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(1) = &
      p_idata3d(1,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(2) = &
      p_idata3d(2,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(3) = &
      p_idata3d(5,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(4) = &
      p_idata3d(6,i)
        
      ! save the number of the element this face was found from
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(5)=i
      
      ! assign the local face number
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(6)=2
        
      ! increment counter
      j=j+1
      
      !=========================================================  
      ! fourth face
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(1) = &
      p_idata3d(4,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(2) = &
      p_idata3d(3,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(3) = &
      p_idata3d(7,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(4) = &
      p_idata3d(8,i)
        
      ! save the number of the element this face was found from
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(5)=i

      ! assign the local face number
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(6)=4
        
      ! increment counter
      j=j+1
        
      !=========================================================  
      ! third face
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(1) = &
      p_idata3d(2,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(2) = &
      p_idata3d(3,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(3) = &
      p_idata3d(6,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(4) = &
      p_idata3d(7,i)
        
      ! save the number of the element this face was found from
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(5)=i
        
      ! assign the local face number
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(6)=3
        
      ! increment counter
      j=j+1
        
      !=========================================================  
      ! fifth face
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(1) = &
      p_idata3d(1,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(2) = &
      p_idata3d(4,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(3) = &
      p_idata3d(8,i)
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(4) = &
      p_idata3d(5,i)
        
      ! save the number of the element this face was found from
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(5)=i
      
      ! assign the local face number
      p_IConnectList( (i-1) * NVFACE + j)%I_conData(6)=5
        
      !=========================================================  
    
  end do
    
  ! done...
  end subroutine ! end tria_buildConnectorList
  
!************************************************************************   
  
!<subroutine>  
  subroutine tria_genNeighboursAtElement3D(rtriangulation)

!<inputoutput>
  ! The triangulation structure to be updated.
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>

  ! the array this routine is supposed to build
  integer, dimension(:,:), pointer :: p_IneighboursAtElement3
     
!</subroutine>
     
  ! local variables
  integer :: iElements, nae, iel, j
  
  integer(i32), dimension(2) :: Isize

  ! a pointer to a list of connectors
  type(t_connector3d), dimension(:), pointer :: p_IConnectList
    
  ! function body ...
  
  nae = rtriangulation%NNAE
  
  ! Do we have (enough) memory for that array?
  if (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) then
    Isize = (/nae,int(rtriangulation%NEL,I32)/)
    call storage_new2D ('tria_genNeighboursAtElement3D', 'KADJ', &
        Isize, ST_INT, &
        rtriangulation%h_IneighboursAtElement, ST_NEWBLOCK_NOINIT)
  else
    call storage_getsize2D (rtriangulation%h_IneighboursAtElement, Isize)
    if (Isize(2) .NE. rtriangulation%NEL) then
      ! If the size is wrong, reallocate memory.
      call storage_realloc ('tria_genNeighboursAtElement3D', &
          rtriangulation%NEL, rtriangulation%h_IneighboursAtElement, &
          ST_NEWBLOCK_NOINIT, .FALSE.)
    end if
  end if
  
  ! get a pointer to the memory just allocated
  call storage_getbase_int2d(rtriangulation%h_IneighboursAtElement, &
  p_IneighboursAtElement3)
  
  p_IneighboursAtElement3(:,:) = 0  
    
  ! first build the connector list
  call tria_buildConnectorList(p_IConnectList,rtriangulation)
    
  iElements = rtriangulation%NEL*rtriangulation%NNAE
    
  ! ConnectorList is build, now sort it
  call tria_sortElements3dInt(p_IConnectList,iElements)
  call tria_sortElements3d(p_IConnectList,iElements)
  
  ! assign the neighbours at elements
  ! traverse the connector list
  do iel = 2,iElements
  
    ! check for equivalent connectors... that means:
    ! check if all 4 vertices that define the face are equal
    j = 0
    do while(p_IConnectList(iel-1)%I_conData(j+1) == &
         p_IConnectList(iel)%I_conData(j+1) )
         ! increment counter
         j=j+1 
    end do
    
    ! assign information
    if(j==4) then
      ! iel is a neighbour of iel-1 at the p_IConnectList(iel-1)%I_conData(6) face
      p_IneighboursAtElement3(p_IConnectList(iel-1)%I_conData(6), &
      p_IConnectList(iel-1)%I_conData(5)) = &
      p_IConnectList(iel)%I_conData(5)
      
      ! iel-1 is a neighbour of iel at the p_IConnectList(iel)%I_conData(6) face
      p_IneighboursAtElement3(p_IConnectList(iel)%I_conData(6), &
      p_IConnectList(iel)%I_conData(5)) = &
      p_IConnectList(iel-1)%I_conData(5)
    end if
  
  end do 
  
  
  ! free list of connectors
  deallocate(p_IConnectList)
  
  ! done... easy...
    
  end subroutine ! tria_genNeighboursAtElement3D
  
!************************************************************************   
  
!<subroutine>    
  subroutine tria_sortElements3dInt(p_IConnectList, iElements)

!<description>
  ! this routine builds the connector list used 
  ! in the neighbours at elements 
  ! routine
!</description>
    
  ! parameter values

!<inputoutput>        
  type(t_connector3d), dimension(:), pointer :: p_IConnectList
!</inputoutput>          
    
!<input>    
  integer, intent(in) :: iElements
!</input>
    
!</subroutine>
    
  ! local variables
  integer, dimension(:), pointer :: p_IEntries
  
  integer :: i

  ! create a sorted numbering in all connectors    
  do i=1,iElements
    
      p_IEntries =>p_IConnectList(i)%I_conData(1:4)
        
      call tria_quickSortInt(p_IEntries, 1, 4)
    
  end do
    
  end subroutine ! end tria_sortElements
    
!************************************************************************     
  
!<subroutine>   
  recursive subroutine tria_quickSortInt(p_IConnector, l, r)
    
      
!<description>
  ! this routine builds the connector list used 
  ! in the neighbours at elements 
  ! routine
!</description>
    
  integer, dimension(:), pointer :: p_IConnector
    
!<input>    
  integer, intent(in) :: l,r
!</input>  

!</subroutine>   

  integer :: m, lpos, rpos, pivot
  
  integer :: temp  
    
  if( l < r) then
    
      ! counter from left
      lpos = l+1
        
      ! counter from right
      rpos = r
        
      ! assign the pivot element
      pivot = p_IConnector(l)
        
      do while(lpos <= rpos)
            
          ! we find an element less than pivot => increment
          if(p_IConnector(lpos) <= pivot) then
            lpos=lpos+1
          ! we find an element greater than => decrement
          else if(p_IConnector(rpos) > pivot) then
            rpos=rpos-1
          ! ok swap the elements  
          else
              ! swap connectors
              temp = p_IConnector(lpos)
                
              p_IConnector(lpos) = p_IConnector(rpos)
                
              p_IConnector(rpos) = temp
          end if
        
      end do
        
      ! swap pivot element            
      temp = p_IConnector(l)
        
      p_IConnector(l) = p_IConnector(rpos)
        
      p_IConnector(rpos) = temp
        
      ! assign new pivot position
      m=rpos
      ! recursively sort the two subarrays
      call tria_quickSortInt(p_IConnector, l, m-1)
      call tria_quickSortInt(p_IConnector, m+1,r)
    
  end if !
    
  end subroutine !end tria_quickSortInt
  
!====================================================================
!
!   this subroutine establishes the lexicographic ordering on the
!   list of connectors in 3d
!
!====================================================================
!<subroutine>
  subroutine tria_sortElements3d(p_IConnectList,iElements)
  
!<description>
  ! this routine builds the connector list used 
  ! in the neighbours at elements 
  ! routine
!</description>

!<input>    
  integer, intent(in) :: iElements
!</input>  

!<inputoutput>        
  type(t_connector3d), dimension(:), pointer :: p_IConnectList
!</inputoutput>

!</subroutine>

  ! local
  integer :: j
    
  do j=TRIA_NCONNECT3D,1,-1
      call tria_mergesort(p_IConnectList, 1, iElements, j)
  end do
    
  end subroutine ! end tria_sortElements
  
!====================================================================  
  
!<subroutine>  
  subroutine tria_genEdgesAtElement3D(rtriangulation)    
!<description>
  ! this routine creates the information
  ! which edges belong to an element
  ! and stores them in IedgesAtElement
  ! Furthermore the edge numbering is calculated. That means,
  ! it generates information about the edges adjacent to
  ! each element IedgesAtElement (KMID) and calculates the correct NMT.
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IneighboursAtElement.
  ! If necessary, new memory is allocated.  
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
    type(t_triangulation), intent(inout) :: rtriangulation
!</inputoutput>
  
!</subroutine>
  ! local variables
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IneighboursAtElement
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx
  integer, dimension(:), pointer :: p_IelementsAtVertex
  integer(prec_edgeidx), dimension(:,:), pointer :: p_IedgesAtElement
  integer(prec_elementidx) :: iel,ied, iSCElement
  integer(prec_edgeidx) :: iedge
  integer(i32), dimension(2) :: Isize
  
  integer(prec_vertexidx) :: iVertex1, iVertex2, iVGlobal1, iVGlobal2
  integer(prec_vertexidx) :: iVertexAtEdge1,iVertexAtEdge2
  
  ! list of local edges
  integer, dimension(2,12) :: Iedges
  
  integer :: i
  
  ! initialize the local edges
  Iedges= reshape((/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/),&
  (/2,12/))
  
  ! check if the arrays are present      
  if (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) then
    call output_line ('IverticesAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
    call sys_halt()
  end if

  if (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) then
    call output_line ('IneighboursAtElement not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtElement2D')
    call sys_halt()
  end if

  ! Get the arrays.
  call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
  
  call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
  ! get the elements at vertex array
  call storage_getbase_int(rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
  ! Get the index array.
  call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
       p_IelementsAtVertexIdx)
  
    
  ! size of I_edgesAtElement
  Isize=(/rtriangulation%NNEE, rtriangulation%NEL/)
 
  ! allocate memory 
  call storage_new2D ('tria_genEdgesAtElement2D', 'KMID', &
        Isize, ST_INT, &
        rtriangulation%h_IedgesAtElement, ST_NEWBLOCK_NOINIT)
  
  ! Fill IedgesAtElement with 0. That's important in case some
  ! elements in the array are not tackled when searching for edges
  ! (e.g. in meshes where triangles and quads are mixed).
  call storage_clear (rtriangulation%h_IedgesAtElement)

  ! get the pointer to the memory
  call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    
  ! iedge counts the edges and specifies the last given edge number.
  ! The numbering starts with NVT, the first edge gets the number NVT+1.
  iedge = rtriangulation%NVT
  
  ! loop over all elements
  do iel=1,rtriangulation%NEL
  
    ! loop over all edges at the element
    do ied=1, rtriangulation%NNEE
    
      ! get the local vertex indices of the current edge
      iVertex1 = Iedges(1,ied)
      iVertex2 = Iedges(2,ied)
    
      ! get the corresponding global vertex numbers
      iVGlobal1 = p_IverticesAtElement(iVertex1,iel)
      iVGlobal2 = p_IverticesAtElement(iVertex2,iel)
      
      ! get the element with the lowest element number
      ! that contains the edge consisting of the
      ! vertices iVGlobal1 and iVGlobal2
      call tria_smallestCommonElement(rtriangulation, &
                   iVGlobal1, iVGlobal2,iel,iSCElement) 
                   
      
      ! if the smallest common element index greater or equal to
      ! the current iel the edge does not yet have an index
      ! so assign an index to the edge                   
      if(iSCElement >= iel) then
      
      ! increment edge number
      iedge = iedge + 1
      
      ! assign the edge number
      p_IedgesAtElement(ied,iel) = iedge
      
      else
      ! the smallest common element index is less than the current iel
      ! the edge already has a number, so                   
      ! search for edge (iVertex1,iVertex2) in the smallest common element
      
      do i=1,rtriangulation%NNEE
      
          ! get the indices of the current edge of iSCElement
          iVertexAtEdge1=p_IVerticesAtElement(Iedges(1,i),iSCElement)
          iVertexAtEdge2=p_IVerticesAtElement(Iedges(2,i),iSCElement)
          
          ! for better readability ... ;)
          ! if iVGlobal1 == iVertexAtEdge1 and iVGlobal2 == iVertexAtEdge2
          if( (iVGlobal1 == iVertexAtEdge1) & 
              .and. (iVGlobal2 == iVertexAtEdge2) .or. &
              ! or iVGlobal1==iVertexAtEdge2 and iVGlobal2==iVertexAtEdge1
              (iVGlobal1 == iVertexAtEdge2) .and. &
              (iVGlobal2 == iVertexAtEdge1)) then
          
          ! assign the edge number
          p_IedgesAtElement(ied,iel) = p_IedgesAtElement(i,iSCElement)
          
          ! the edge was found we can break out of the loop
          exit
              
          end if
      
      end do ! end i
      
      end if
    
    end do ! end ied
  
  end do ! end iel
  
  ! Save the correct NMT.
  rtriangulation%NMT = iedge-rtriangulation%NVT
  
  end subroutine ! end tria_genEdgesAtElement3D
  
!====================================================================    
!<subroutine>
  subroutine tria_smallestCommonElement(rtriangulation, iVGlobal1, iVGlobal2, iel, iSCElement)       
  
!<description>
! this routine computes the element with the lowest element number
! that contains the edge consisting of vertices iVGlobal1 and iVGlobal2
! and returns it in iSCElement
!</description>

!<input>  
  integer(prec_vertexidx), intent(in) :: iVGlobal1, iVGlobal2, iel
  type(t_triangulation), intent(in) :: rtriangulation  
!</input>
  
!<output>  
  integer(prec_elementidx), intent(out) :: iSCElement  
!</output>
  
!</subroutine>

! local variables
  integer, dimension(:), pointer :: p_IelementsAtVertexIdx
  integer, dimension(:), pointer :: p_IelementsAtVertex
  
  integer(prec_elementidx) :: ive11, ive12, ive21, ive22
  
  ! used for looping
  integer(prec_elementidx) :: i,j
  
  integer(prec_elementidx) :: iel1, iel2

  ! get the elements at vertex array
  call storage_getbase_int(rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
  ! Get the index array.
  call storage_getbase_int (rtriangulation%h_IelementsAtVertexIdx,&
       p_IelementsAtVertexIdx)
       
  iSCElement = iel
       
  ive11 = p_IelementsAtVertexIdx(iVGlobal1)
  ive12 = p_IelementsAtVertexIdx(iVGlobal1+1)-1
  
  ive21 = p_IelementsAtVertexIdx(iVGlobal2)       
  ive22 = p_IelementsAtVertexIdx(iVGlobal2+1)-1
  
  ! loop over all the elements attached to vertex iVGloval1
  do i = ive11, ive12
    
    ! current element at vertex iVGlobal1
    iel1 = p_IelementsAtVertex(i)
    
    ! loop over all elements attached to vertex iVGlobal2
    do j = ive21, ive22
    
      ! get the current element at vertex iVGlobal2
      iel2 = p_IelementsAtVertex(j) 
      
      ! check the vertices share element iel1
      if( ( iel2 == iel1) .and. (iel2 < iSCElement  ) ) then
        iSCElement = iel2      
      end if
      
    end do ! end j
  end do ! end i       

  ! done...

  end subroutine ! end tria_smallestCommonElement

!====================================================================    
!<subroutine>  
  subroutine tria_genElementsAtEdge3D(rtriangulation)
  
!<description>
  ! This routine generates information about the elements adjacent to each 
  ! edge IelementsAtEdge (KMID). 
  ! For this purpose, the following arrays are used:
  ! IverticesAtElement, IneighboursAtElement.
  ! NMT must be set up correctly.
  ! If necessary, new memory is allocated.
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
    
  
!</subroutine>  
  ! local variables
  integer(prec_elementidx), dimension(:,:), pointer :: p_IneighboursAtElement
  integer(prec_elementidx), dimension(:,:), pointer :: p_IedgesAtElement
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IverticesAtElement
  integer(prec_edgeidx), dimension(:), pointer :: p_IelementsAtEdgeIdx3d
  integer(prec_edgeidx), dimension(:), pointer :: p_IelementsAtEdge3d
  integer(prec_edgeidx), dimension(:), pointer :: p_Iaux1
  integer(prec_elementidx) :: iel
  integer(prec_edgeidx) :: NMT
  integer(prec_edgeidx) :: iedge,iglobalEdge
  integer(i32) :: Isize,haux1

  ! Is everything here we need?
  if (rtriangulation%h_IverticesAtElement .EQ. ST_NOHANDLE) then
    call output_line ('IverticesAtElement not available!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
    call sys_halt()
  end if

  if (rtriangulation%h_IedgesAtElement .EQ. ST_NOHANDLE) then
    call output_line ('IedgesAtElement not available!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
    call sys_halt()
  end if

  if (rtriangulation%h_IneighboursAtElement .EQ. ST_NOHANDLE) then
    call output_line ('IneighboursAtElement not available!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
    call sys_halt()
  end if

  if (rtriangulation%NMT .EQ. 0) then
    call output_line ('Edge information (NMT) not initialised!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'tria_genElementsAtEdge3D')
    call sys_halt()
  end if

  ! Get the arrays.
  call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
  call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
  call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
      p_IneighboursAtElement)
      
  if(rtriangulation%h_IelementsAtEdgeIdx3d == ST_NOHANDLE) then
    call storage_new ('tria_genElementsAtEdge3D', 'IelementsAtVertexIdx', &
        int(rtriangulation%NMT+1,I32), ST_INT, &
        rtriangulation%h_IelementsAtEdgeIdx3d, ST_NEWBLOCK_NOINIT)
  end if      
    
  ! Fill the index array with zero.
  call storage_clear (rtriangulation%h_IelementsAtEdgeIdx3d)
    
    
  call storage_getbase_int(rtriangulation%h_IelementsAtEdgeIdx3d,p_IelementsAtEdgeIdx3d)
    
  ! nvt shorthand  
  NMT = rtriangulation%NMT
  
  ! first we calculate the number of elements at each edge simply by counting
    
  ! loop over all elements
  do iel = 1, rtriangulation%NEL
    ! loop over all local edges
    do iedge = 1, rtriangulation%NNEE

      ! iglobalEdge is the iedge-th edge at element iel
      iglobalEdge = p_IedgesAtElement(iedge,iel)-rtriangulation%NVT
            
      ! increase the number of elements by one
      p_IelementsAtEdgeIdx3d(iglobalEdge+1) = p_IelementsAtEdgeIdx3d(iglobalEdge+1) + 1
        
    end do ! end iedge
        
  end do ! end ive
    
  ! set the first index to 1
  p_IelementsAtEdgeIdx3d(1) = 1
    
  ! in the next step we sum up the number of elements at two
  ! successive edges to create the index array,
  ! thus at the penultimate position of p_IelementsAtEdgeIdx3d
  ! we find the length of p_IelementsAtEdge3d
  do iedge = 2, NMT+1
      p_IelementsAtEdgeIdx3d(iedge) = p_IelementsAtEdgeIdx3d(iedge) + p_IelementsAtEdgeIdx3d(iedge-1)
  end do
    
  ! set the size
  Isize = p_IelementsAtEdgeIdx3d(rtriangulation%NMT+1)-1
  
  ! allocate memory
  if(rtriangulation%h_IelementsAtEdge3d == ST_NOHANDLE) then
    call storage_new ('tria_genElementsAtEdge3D', 'IelementsAtEdge3d', &
        int(Isize), ST_INT, &
        rtriangulation%h_IelementsAtEdge3d, ST_NEWBLOCK_NOINIT)
  end if      
    
  ! Fill the index array with zero.
  call storage_clear (rtriangulation%h_IelementsAtEdge3d)
    
  ! get the pointer  
  call storage_getbase_int(rtriangulation%h_IelementsAtEdge3d,p_IelementsAtEdge3d)

  ! Duplicate the p_IelementsAtVertexIdx array. We use that as pointer and index
  ! when new elements at a vertex are found.
  haux1 = ST_NOHANDLE
  call storage_copy (rtriangulation%h_IelementsAtEdgeIdx3d,haux1)
  call storage_getbase_int (haux1,p_Iaux1)  

  do iel = 1, rtriangulation%NEL
    ! loop over all local edges
    do iedge = 1, rtriangulation%NNEE

      ! iglobalEdge is the iedge-th edge at element iel
      iglobalEdge = p_IedgesAtElement(iedge,iel)-rtriangulation%NVT
       
      ! store the adjacency information at position p_Iaux1(ivt)        
      p_IelementsAtEdge3d( p_Iaux1(iglobalEdge) ) = iel
      ! increase the position of the next element in p_Iaux1(ivt)
      p_Iaux1(iglobalEdge) = p_Iaux1(iglobalEdge) + 1
        
    end do ! end iedge
        
  end do ! end ive
    
  call storage_free(haux1)  
  
  end subroutine ! end tria_genElementsAtEdge3D
  
!====================================================================    
!<subroutine>  
  subroutine tria_genVerticesAtEdge3D(rtriangulation)  
!<description>
! this routine assigns the Vertices at Edge array
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
    
  
!</subroutine>  
  
  ! local variables
  integer(prec_elementidx), dimension(:,:), pointer :: p_IedgesAtElement
  
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IverticesAtElement
  
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IverticesAtEdge
  
  integer(prec_edgeidx), dimension(:), pointer :: p_IelementsAtEdgeIdx3d
  integer(prec_edgeidx), dimension(:), pointer :: p_IelementsAtEdge3d
  integer(prec_elementidx) :: iel 
  integer(prec_edgeidx) :: NMT,ilocEdge
  integer(prec_edgeidx) :: iedge,iglobalEdge
  integer(prec_vertexidx) :: iVertexGlobal1, iVertexGlobal2
  
  ! list of local edges
  integer, dimension(2,12) :: Iedges
  
  integer, dimension(2)    :: Isize
  
  ! initialize the local edges
  Iedges= reshape((/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/),&
  (/2,12/))
  ! nmt shorthand
  NMT = rtriangulation%NMT

  ! get the arrays
  call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
  call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
  call storage_getbase_int(rtriangulation%h_IelementsAtEdgeIdx3d,p_IelementsAtEdgeIdx3d)
  call storage_getbase_int(rtriangulation%h_IelementsAtEdge3d,p_IelementsAtEdge3d)
  
  Isize = (/2,rtriangulation%NMT/)
  
  ! allocate memory
  if(rtriangulation%h_IverticesAtEdge == ST_NOHANDLE) then
    call storage_new2d('tria_genElementsAtEdge3D', 'IverticesAtEdge', &
        Isize, ST_INT, &
        rtriangulation%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
  end if      
  
  call storage_getbase_int2D (rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)

  
  ! loop over all edges (with global edge numbers)
  do iedge=1,rtriangulation%NMT
      ! get an element that is attached to that edge
      iel = p_IelementsAtEdge3d(p_IelementsAtEdgeIdx3d(iedge))
      ! loop over all local edges
      do ilocEdge = 1,rtriangulation%NNEE
        ! get the global edge number
        iglobalEdge = p_IedgesAtElement(ilocEdge,iel)-rtriangulation%NVT
        ! check if this edge's global number equal to the current edge
        if(iedge == iglobalEdge) then
          ! get the global indices of the vertices attached to that edge
          iVertexGlobal1 = p_IverticesAtElement(Iedges(1,ilocEdge),iel)
          iVertexGlobal2 = p_IverticesAtElement(Iedges(2,ilocEdge),iel)
          
          ! assign the vertex numbers            
          p_IverticesAtEdge(1,iedge)=iVertexGlobal1
          p_IverticesAtEdge(2,iedge)=iVertexGlobal2
          
          exit
          
        end if ! iedge == iglobalEdge
      
      end do ! end iel
      
  end do ! end iedge  
  
  end subroutine ! end tria_genVerticesAtEdge3D

!==================================================================== 
   
!<subroutine>  
  subroutine tria_genFacesAtElement(rtriangulation)
!<description>
  ! this routine builds the FacesAtElement array and
  ! assigns global face numbers
!</description>
  
  ! it is super easy just use IneighboursAtElement(iface,iel)

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
  
!</subroutine>  

  ! local variables
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IverticesAtElement
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IfacesAtElement
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IneighboursAtElement
  
  
  integer(prec_vertexidx)  :: ifaceNeighbour
  integer(prec_elementidx) :: iel,iface, ifaceGlobal,ineighbour, ifaceNumber

  ! list of local face numbers
  integer, dimension(4,6) :: Ifaces
  
  integer, dimension(2) :: Isize
  
  ! initialize the local face numbers
  Ifaces= reshape((/1,2,3,4, 1,2,5,6, 2,3,6,7, 3,4,7,8, 1,4,5,8, 5,6,7,8/),&
  (/4,6/))
  
  ! face numbering starts at NVT+NMT+1
  ifaceGlobal = rtriangulation%NVT+rtriangulation%NMT
  
  ! allocate memory and get pointers
  call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
  
  call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
      p_IneighboursAtElement)
      
  Isize = (/rtriangulation%NNAE,rtriangulation%NEL/)
  
  ! allocate memory
  if(rtriangulation%h_IfacesAtElement == ST_NOHANDLE) then
    call storage_new2d('tria_genElementsAtEdge3D', 'IfacesAtElement', &
        Isize, ST_INT, &
        rtriangulation%h_IfacesAtElement, ST_NEWBLOCK_NOINIT)
  end if      
  
  ! get the pointer
  call storage_getbase_int2D(rtriangulation%h_IfacesAtElement,p_IfacesAtElement)
  
  ! loop over all elements
  do iel=1,rtriangulation%NEL
    ! loop over all local faces
    do iface = 1,rtriangulation%NNAE
      
      ! check if a face number was already assigned
      if(p_IneighboursAtElement(iface,iel) == 0 .or. &
         p_IneighboursAtElement(iface,iel) >  iel) then

        ! a face number was not yet assigned
        ! increment face number
        ifaceGlobal = ifaceGlobal + 1
      
        ! assign the global face number
        p_IfacesAtElement(iface,iel) = ifaceGlobal
      
      else
      
        ! a face number was already assigned
        
        ! get the element number of the neighbour
        ineighbour = p_IneighboursAtElement(iface,iel)
        
        ifaceNeighbour = 1
        do 
         if(iel == p_IneighboursAtElement(ifaceNeighbour,ineighbour)) exit
         ifaceNeighbour = ifaceNeighbour + 1
        end do
        
        ifaceNumber = p_IfacesAtElement(ifaceNeighbour,ineighbour)
        
        ! assign the global face number
        p_IfacesAtElement(iface,iel) = ifaceNumber
      
      end if
    
    end do ! end iface
  
  end do ! end iel
  
  ! number of faces in total
  rtriangulation%NAT = ifaceGlobal-rtriangulation%NMT - &
                       rtriangulation%NVT 
  
  end subroutine ! end tria_genFacesAtElement
  
  
!====================================================================      
!<subroutine>  
  subroutine tria_genVerticesAtFace(rtriangulation)
!<description>
  ! this routine builds the VerticesAtFace array and
!</description>
  
!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
  
!</subroutine>    

  ! local variables
  integer(prec_vertexidx), dimension(:,:), pointer  :: p_IverticesAtElement
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IfacesAtElement
  integer(prec_elementidx), dimension(:,:), pointer :: p_IverticesAtFace
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IneighboursAtElement
  
  
  integer(prec_vertexidx)  :: ifaceNeighbour
  integer(prec_elementidx) :: iel,iface, ifaceGlobal,ineighbour, ifaceNumber

  ! list of local face numbers
  integer, dimension(4,6) :: Ifaces
  
  integer, dimension(2) :: Isize
  
  ! initialize the local face numbers
  Ifaces= reshape((/1,2,3,4, 1,2,5,6, 2,3,6,7, 3,4,7,8, 1,4,5,8, 5,6,7,8/),&
  (/4,6/))
  
  ! allocate memory and get pointers
  call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
  
  call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
      p_IneighboursAtElement)
      
  call storage_getbase_int2D (rtriangulation%h_IfacesAtElement,&
      p_IfacesAtElement)
      
  ! init isize  
  Isize = (/4,rtriangulation%NAT/)
  
  ! allocate memory
  if(rtriangulation%h_IverticesAtFace == ST_NOHANDLE) then
    call storage_new2d('tria_genElementsAtEdge3D', 'IverticesAtFace', &
        Isize, ST_INT, &
        rtriangulation%h_IverticesAtFace, ST_NEWBLOCK_NOINIT)
  end if      
  
  ! get the pointer
  call storage_getbase_int2D(rtriangulation%h_IverticesAtFace,p_IverticesAtFace)
    
  ifaceGlobal = 0
  
  ! loop over all elements
  do iel=1,rtriangulation%NEL
    ! loop over all local faces
    do iface = 1,rtriangulation%NNAE
      
      ! check if a face number was already assigned
      if(p_IneighboursAtElement(iface,iel) == 0 .or. &
         p_IneighboursAtElement(iface,iel) >  iel) then

        ! a face number was not yet assigned
        ! increment face number
        ifaceGlobal = ifaceGlobal + 1
      
        ! assign the vertices at this face
        p_IverticesAtFace(:,ifaceGlobal) = p_IverticesAtElement(ifaces(:,iface),iel)
      
         
      else
      
        ! a face number was already assigned
        
        ! get the element number of the neighbour
        ineighbour = p_IneighboursAtElement(iface,iel)
        
        ifaceNeighbour = 1
        
        ! number of the face
        do 
         if(iel == p_IneighboursAtElement(ifaceNeighbour,ineighbour)) exit
         ifaceNeighbour = ifaceNeighbour + 1
        end do
        
        ! index of the face in the array is obtained by nmt and nvt
        ifaceNumber = p_IfacesAtElement(ifaceNeighbour,ineighbour) - &
                      rtriangulation%NMT - rtriangulation%NVT
        
        ! assign the vertices at this face
        p_IverticesAtFace(:,ifaceNumber) = p_IverticesAtElement(ifaces(:,iface),iel)
      
      end if
    
    end do ! end iface
  
  end do ! end iel

  end subroutine ! end tria_genVerticesAtFace
  
!====================================================================      
!<subroutine>  
  subroutine tria_genElementsAtFace(rtriangulation) 
!<description>
  ! this routine builds the ElementsAtFace array 
!</description>
  
!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>

!</subroutine>
  
  ! local variables
  integer(prec_elementidx), dimension(:,:), pointer :: p_IfacesAtElement
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IneighboursAtElement
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IelementsAtFace
  
  
  integer(prec_faceidx) :: iface
  integer(prec_elementidx) :: iel, ifaceNumber
  
  integer, dimension(2) :: Isize

  ! size of target array  
  Isize = (/2,rtriangulation%NAT/)
  
  ! get pointers to some needed connectivity information
  call storage_getbase_int2D (rtriangulation%h_IneighboursAtElement,&
      p_IneighboursAtElement)
      
  call storage_getbase_int2D (rtriangulation%h_IfacesAtElement,&
      p_IfacesAtElement)

  ! allocate memory
  if(rtriangulation%h_IelementsAtFace == ST_NOHANDLE) then
    call storage_new2d('tria_genElementsAtFace', 'IelementsAtFace', &
        Isize, ST_INT, &
        rtriangulation%h_IelementsAtFace, ST_NEWBLOCK_NOINIT)
  end if      
  
  call storage_getbase_int2D (rtriangulation%h_IelementsAtFace,&
      p_IelementsAtFace)
  

  ! loop over all elements
  do iel=1,rtriangulation%NEL
    ! loop over all faces of this element
    do iface=1,rtriangulation%NNAE
    
      ! if there is no neighbour at this element
      if(p_IneighboursAtElement(iface,iel) == 0) then
        
        ifaceNumber = p_IfacesAtElement(iface,iel) - &
                      rtriangulation%NVT - &
                      rtriangulation%NMT
                      
        p_IelementsAtFace(1,ifaceNumber) = iel
        p_IelementsAtFace(2,ifaceNumber) = 0
      
      else if (p_IneighboursAtElement(iface,iel) < iel) then
        
        ! There is a neighbour and it has a smaller number than the current element --
        ! so we haven't had that face! Store the two adjacent elements.
        
        ifaceNumber = p_IfacesAtElement(iface,iel) - &
                      rtriangulation%NVT - &
                      rtriangulation%NMT
        
        p_IelementsAtFace(1,ifaceNumber) = iel
        
        p_IelementsAtFace(2,ifaceNumber) = &
        p_IneighboursAtElement(iface,iel)      
      
      end if
    
    end do ! end iface
  end do ! end iel  
  
  end subroutine ! end tria_genElementsAtFace
  
  
!====================================================================          
!<subroutine>    
  subroutine tria_genEdgesAtFace(rtriangulation)  
!<description>
  ! this routine builds the EdgesAtFace array 
!</description>
  
!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>
  
!</subroutine>
  ! local variables
  
  integer(prec_edgeidx), dimension(:,:), pointer    :: p_IedgesAtElement
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IelementsAtFace
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IedgesAtFace
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IfacesAtElement
  
  integer :: iface, iel, ilocalFace,iglobalFace
  
  ! list of local edges
  integer, dimension(2) :: Isize
  
  Isize = (/4,rtriangulation%NAT/)
  
  ! allocate memory
  if(rtriangulation%h_IedgesAtFace == ST_NOHANDLE) then
    call storage_new2d('tria_genEdgesAtFace', 'IedgesAtFace', &
        Isize, ST_INT, &
        rtriangulation%h_IedgesAtFace, ST_NEWBLOCK_NOINIT)
  end if      
  
  
  ! get pointers to some needed connectivity information
  call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,&
      p_IedgesAtElement)
      
  call storage_getbase_int2D (rtriangulation%h_IelementsAtFace,&
      p_IelementsAtFace)
      
  call storage_getbase_int2D (rtriangulation%h_IfacesAtElement,&
      p_IfacesAtElement)

  call storage_getbase_int2D (rtriangulation%h_IedgesAtFace,&
      p_IedgesAtFace)

      
  ! loop over all Faces
  do iface=1,rtriangulation%NAT
  
    ! get a face this element is conncected to
    iel = p_IelementsAtFace(1,iface)
    
    ! determine which local face is iface 
    do ilocalFace=1,rtriangulation%NNAE
     iglobalFace = p_IfacesAtElement(ilocalFace,iel) - &
                   rtriangulation%NVT - &
                   rtriangulation%NMT
     if(iglobalFace == iface) exit
    end do
    
    ! now I need the edge numbers for this face
    ! get them by edgesAtElement
    
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
    
    end select
    
  end do ! end iface
     
  end subroutine ! end tria_genEdgesAtFace
  
!====================================================================        
  
!<subroutine>    
  subroutine tria_genFacesAtEdge(rtriangulation) 
  
!<description>
  ! this routine builds the FacesAtEdge array 
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>

!</subroutine>
    
  ! local parameters
  integer(prec_elementidx), dimension(:,:), pointer :: p_IedgesAtFace
  
  integer(prec_edgeidx), dimension(:), pointer :: p_IfacesAtEdgeIdx
  
  integer(prec_edgeidx), dimension(:), pointer :: p_IfacesAtEdge
  
  integer :: iface, iglobalEdge, iedge
  
  integer(i32) :: haux1
  
  integer(prec_edgeidx), dimension(:), pointer :: p_Iaux1
  
  ! list of local edges
  integer, dimension(2) :: Isize
  
  Isize = (/4,rtriangulation%NAT/)
  
  ! allocate memory
  if(rtriangulation%h_IfacesAtEdgeIdx == ST_NOHANDLE) then
    call storage_new ('tria_genElementsAtEdge3D', 'IfacesAtEdgeIdx', &
        int(rtriangulation%NMT+1,I32), ST_INT, &
        rtriangulation%h_IfacesAtEdgeIdx, ST_NEWBLOCK_NOINIT)
  end if      
    
  ! Fill the index array with zero.
  call storage_clear (rtriangulation%h_IfacesAtEdgeIdx)
    
  ! get the pointer  
  call storage_getbase_int(rtriangulation%h_IfacesAtEdgeIdx,p_IfacesAtEdgeIdx)
    
  ! Fill the index array with zero.
  call storage_clear(rtriangulation%h_IfacesAtEdgeIdx)
    
  ! get the pointer  
  call storage_getbase_int(rtriangulation%h_IfacesAtEdgeIdx,p_IfacesAtEdgeIdx)
  
  ! get the pointer  
  call storage_getbase_int2d(rtriangulation%h_IedgesAtFace,p_IedgesAtFace)
  
  
  p_IfacesAtEdgeIdx(1) = 1;
  
  ! create the index array
  do iface=1, rtriangulation%NAT
  
    ! increase the facecount at these edges by one
    do iedge=1,4
    iglobalEdge = p_IedgesAtFace(iedge,iface) - &
                  rtriangulation%NVT
                
    p_IfacesAtEdgeIdx(iglobalEdge+1) = &
    p_IfacesAtEdgeIdx(iglobalEdge+1) + 1 
    end do ! end iedge
  
  end do ! end iface

  ! create the actual index array  
  do iedge = 2,rtriangulation%NMT+1
    p_IfacesAtEdgeIdx(iedge) = p_IfacesAtEdgeIdx(iedge) + p_IfacesAtEdgeIdx(iedge-1);
  end do ! end iedge
  
  Isize(1) = p_IfacesAtEdgeIdx(rtriangulation%NMT+1)-1
  
  ! allocate memory
  if(rtriangulation%h_IfacesAtEdge == ST_NOHANDLE) then
    call storage_new ('tria_genFacesAtEdge3D', 'IfacesAtEdge', &
        int(Isize(1)), ST_INT, &
        rtriangulation%h_IfacesAtEdge, ST_NEWBLOCK_NOINIT)
  end if      
  
  ! get the pointer
  call storage_getbase_int(rtriangulation%h_IfacesAtEdge, p_IfacesAtEdge)
  
  haux1 = ST_NOHANDLE
  call storage_copy (rtriangulation%h_IfacesAtEdgeIdx,haux1)
  call storage_getbase_int (haux1,p_Iaux1) 
  
  ! assign the connectivity info
  do iface=1, rtriangulation%NAT
  
    ! increase the facecount at these edges by one
    do iedge=1,4
    
      ! iglobalFace is the iedge-th edge at face iface
      iglobalEdge = p_IedgesAtFace(iedge,iface) - &
             rtriangulation%NVT

       
      ! store the adjacency information at position p_Iaux1(ivt)        
      p_IfacesAtEdge( p_Iaux1(iglobalEdge) ) = iface
      ! increase the position of the next element in p_Iaux1(ivt)
      p_Iaux1(iglobalEdge) = p_Iaux1(iglobalEdge) + 1    
     
    end do ! end iedge
  
  end do ! end iface
  
  call storage_free(haux1)   
  
  end subroutine ! end tria_genFacesAtEdge

!====================================================================        
  
!<subroutine>    
  subroutine tria_genFacesAtVertex(rtriangulation)
!<description>
  ! this routine builds the FacesAtVertex list 
!</description>

!<inputoutput>  
  type(t_triangulation), intent(inout) :: rtriangulation  
!</inputoutput>


!</subroutine>
  ! local variables
  
  integer(prec_elementidx), dimension(:,:), pointer :: p_IverticesAtFace
  
  integer(prec_elementidx), dimension(:), pointer :: p_IfacesAtVertexIdx
  integer(prec_elementidx), dimension(:), pointer :: p_IfacesAtVertex
  
  integer :: iface, ivt, iGlobalVertex
  
  integer(i32) :: haux1
  
  integer(prec_edgeidx), dimension(:), pointer :: p_Iaux1
  
  ! list of local edges
  integer :: Isize
  
  Isize = rtriangulation%NVT+1
  
  ! allocate memory
  if(rtriangulation%h_IfacesAtVertexIdx == ST_NOHANDLE) then
    call storage_new ('tria_genElementsAtEdge3D', 'IfacesAtVertexIdx', &
        int(Isize,I32), ST_INT, &
        rtriangulation%h_IfacesAtVertexIdx, ST_NEWBLOCK_NOINIT)
  end if      
    
  ! Fill the index array with zero.
  call storage_clear (rtriangulation%h_IfacesAtVertexIdx)
    
  ! get the pointer  
  call storage_getbase_int(rtriangulation%h_IfacesAtVertexIdx,p_IfacesAtVertexIdx)
  
  call storage_getbase_int2D(rtriangulation%h_IverticesAtFace,p_IverticesAtFace)
    
  ! Fill the index array with zero.
  call storage_clear(rtriangulation%h_IfacesAtVertexIdx)
    
  ! set the first value to 1
  p_IfacesAtVertexIdx(1) = 1;
  
  ! create the index array
  do iface=1,rtriangulation%NAT
    do ivt=1,4
      ! get the global vertex number
      iGlobalVertex = p_IverticesAtFace(ivt,iface)
      ! increment the facecount for this vertex
      p_IfacesAtVertexIdx(iGlobalVertex+1) = &
      p_IfacesAtVertexIdx(iGlobalVertex+1) + 1
    end do ! end ivt
  end do ! end iface

  
  ! create the actual index array  
  do iface = 2,rtriangulation%NVT+1
    p_IfacesAtVertexIdx(iface) = p_IfacesAtVertexIdx(iface) + p_IfacesAtVertexIdx(iface-1);
  end do ! end iface
  
  Isize = p_IfacesAtVertexIdx(rtriangulation%NVT+1)-1
  
  ! allocate memory
  if(rtriangulation%h_IfacesAtVertex == ST_NOHANDLE) then
    call storage_new ('tria_genFacesAtVertex', 'IfacesAtVertex', &
        int(Isize), ST_INT, &
        rtriangulation%h_IfacesAtVertex, ST_NEWBLOCK_NOINIT)
  end if      

  ! get the pointer  
  call storage_getbase_int(rtriangulation%h_IfacesAtVertex,p_IfacesAtVertex)

  ! build the auxilliary array  
  haux1 = ST_NOHANDLE
  call storage_copy (rtriangulation%h_IfacesAtVertexIdx,haux1)
  call storage_getbase_int (haux1,p_Iaux1) 
  
  ! assign the connectivity info
  do iface=1, rtriangulation%NAT
  
    ! get the global vertex number and assign the faces
    do ivt=1,4
    
      ! iglobalFace is the iedge-th edge at face iface
      iGlobalVertex = p_IverticesAtFace(ivt,iface)
      
      ! store the adjacency information at position p_Iaux1(ivt)        
      p_IfacesAtVertex( p_Iaux1(iGlobalVertex) ) = iface
      ! increase the position of the next element in p_Iaux1(ivt)
      p_Iaux1(iGlobalVertex) = p_Iaux1(iGlobalVertex) + 1    
     
    end do ! end iedge
  
  end do ! end iface
  
  call storage_free(haux1)   
  
  end subroutine ! end tria_genFacesAtVertex

!====================================================================        
  
!<subroutine>
  subroutine tria_refineMesh2lv3D(rsourceTriangulation,rdestTriangulation)
!<description>
  ! This routine refines the given 2D mesh rsourceTriangulation according to
  ! the 2-level ordering algorithm. The refined mesh is saved in 
  ! rdestTriangulation. There will be no correction of boundary
  ! vertices. Boundary parameter values are not handled here!
!</description>

!<inputoutput>
  ! The source triangulation to be refined
  type(t_triangulation), intent(INOUT) :: rsourceTriangulation
!</inputoutput>

!<output>
  ! Destination triangulation structure that receives the refined mesg. 
  type(t_triangulation), intent(OUT) :: rdestTriangulation
!</output>
  
!</subroutine>      
  ! local variables
  
  real(dp), dimension(:,:), pointer :: p_DcoordSource
  real(dp), dimension(:,:), pointer :: p_DcoordDest
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IvertAtElementSource
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IvertAtElementDest
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IedgesAtElementSource
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IvertAtEdgeSource
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IverticesAtFace
  integer(prec_vertexidx), dimension(:,:), pointer :: p_IfacesAtElement
  integer(prec_vertexidx), dimension(:), pointer :: p_IrefinementPatch
  integer(prec_vertexidx), dimension(:), pointer :: p_IrefinementPatchIndex
  
  integer(prec_vertexidx), dimension(:), pointer :: p_IcoarseGridElement
  
  integer(prec_elementidx) :: nquads,iel,iel1,iel2,iel3
  integer(prec_edgeidx) :: imt
  integer(prec_vertexidx) :: ivt1, ivt2, ivtoffset, ivt, iae
  integer(prec_vertexidx) :: iel4, iel5, iel6, iel7, midPointOfIel
  integer(prec_vertexidx) :: midpointFaceA, midpointFaceB, midpointFaceC
  integer(prec_vertexidx) :: midpointFaceD, midpointFaceE, midpointFaceF
  
  integer(i32), dimension(2) :: Isize
  integer :: nnve, ive, nve
  real(dp) :: x,y,z
  
  ! Get the arrays with information of the source mesh.
  call storage_getbase_double2D (rsourceTriangulation%h_DvertexCoords,&
      p_DcoordSource)
  call storage_getbase_int2D (rsourceTriangulation%h_IverticesAtElement,&
      p_IvertAtElementSource)
  call storage_getbase_int2D (rsourceTriangulation%h_IedgesAtElement,&
      p_IedgesAtElementSource)
  call storage_getbase_int2D (rsourceTriangulation%h_IverticesAtEdge,&
      p_IvertAtEdgeSource)

  call storage_getbase_int2D (rsourceTriangulation%h_IverticesAtFace,&
      p_IverticesAtFace)
      
  call storage_getbase_int2d(rsourceTriangulation%h_IfacesAtElement, &
      p_IfacesAtElement)
      
  
      
      
      
      
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
      
  nquads = rsourceTriangulation%NEL
  nnve = rsourceTriangulation%NNVE
      
  ! Initialise the basic mesh data in rdestTriangulation:
      
  ! 3D mesh
  rdestTriangulation%ndim = NDIM3D
  
  rdestTriangulation%nnve = nnve
     
  ! Every element is divided into 8 subelements
  rdestTriangulation%NEL = 8 * rsourceTriangulation%NEL
  rdestTriangulation%InelOfType(:) = 8 * rsourceTriangulation%InelOfType(:)

  ! We have now found rdestTriangulation%NEL, so we can:
  ! Allocate memory for the refinement information arrays.
  ! These arrays define for every coarse grid element the fine
  ! grid elements and for every fine grid element the coarse
  ! grid element where it comes from.
  call storage_new ('tria_refineMesh2lv3D', 'h_IrefinementPatchIndex', &
      rsourceTriangulation%NEL+1, ST_INT, &
      rdestTriangulation%h_IrefinementPatchIndex, ST_NEWBLOCK_ZERO)
  call storage_getbase_int(&
      rdestTriangulation%h_IrefinementPatchIndex,p_IrefinementPatchIndex)

  call storage_new ('tria_refineMesh2lv3D', 'h_IrefinementPatch', &
      rdestTriangulation%NEL, ST_INT, &
      rdestTriangulation%h_IrefinementPatch, ST_NEWBLOCK_ZERO)
  call storage_getbase_int(&
      rdestTriangulation%h_IrefinementPatch,p_IrefinementPatch)
      
  call storage_new ('tria_refineMesh2lv3D', 'h_IcoarseGridElement', &
      rdestTriangulation%NEL, ST_INT, &
      rdestTriangulation%h_IcoarseGridElement, ST_NEWBLOCK_ZERO)
  call storage_getbase_int(&
      rdestTriangulation%h_IcoarseGridElement,p_IcoarseGridElement)
      

  ! The p_IrefinementPatchIndex array can directly be initialised
  ! as every coarse grid element gets 4 fine grid elements.
  do iel = 0,rsourceTriangulation%NEL
    p_IrefinementPatchIndex(1+iel) = 1 + iel*8
  end do
  
  ! Also the p_IcoarseGridElement array can directly be initialised.
  ! The coarse grid element number of element 1..NEL(coarse) stays
  ! the same. The number of the other coarse grid elements can
  ! be calculated by a formula.
  do iel = 1,rsourceTriangulation%NEL
    p_IcoarseGridElement(iel) = iel
  end do
  
  ! assign the other element numbers
  ! the number for an element is iel1 = rsourceTriangulation%NEL+7*(iel-1)+1  
  do iel = rsourceTriangulation%NEL+1,rdestTriangulation%NEL
    p_IcoarseGridElement(iel) = (iel-rsourceTriangulation%NEL-1)/7 + 1
  end do

  ! New value of NVT
  ! We expect NVT+NMT+NAT+nquads new points.
  rdestTriangulation%NVT = &
      rsourceTriangulation%NVT + &
      rsourceTriangulation%NMT + &
      rsourceTriangulation%NAT + &
      nquads
      
  ! Allocate memory for the new vertex coordinates and
  ! get the pointers to the coordinate array
  Isize = (/NDIM3D,int(rdestTriangulation%NVT,I32)/)
  call storage_new2D ('tria_refineMesh2lv3D', 'DCORVG', Isize, ST_DOUBLE, &
       rdestTriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
  call storage_getbase_double2D(&
       rdestTriangulation%h_DvertexCoords,p_DcoordDest)
      
  ! Ok, let's start the refinement. In the first step, we copy the
  ! corner coordinates of the coarse mesh to the fine mesh; they
  ! don't change during the refinement.
  do ivt=1,rsourceTriangulation%NVT
     p_DcoordDest(1,ivt) = p_DcoordSource(1,ivt)
     p_DcoordDest(2,ivt) = p_DcoordSource(2,ivt)
     p_DcoordDest(3,ivt) = p_DcoordSource(3,ivt)
  end do
      
  ! Each edge produces an edge midpoint which is stored as new
  ! point in the fine mesh. To calculate the coordinates, take
  ! the mean of the coordinates in the coarse mesh.
  ivtoffset = rsourceTriangulation%NVT
  do imt=1,rsourceTriangulation%NMT
     ivt1 = p_IvertAtEdgeSource (1,imt)
     ivt2 = p_IvertAtEdgeSource (2,imt)
     p_DcoordDest(1,ivtoffset+imt) = &
     0.5_DP * ( p_DcoordSource (1,ivt1) + p_DcoordSource (1,ivt2) )
     p_DcoordDest(2,ivtoffset+imt) = &
     0.5_DP * ( p_DcoordSource (2,ivt1) + p_DcoordSource (2,ivt2) )
     p_DcoordDest(3,ivtoffset+imt) = &
     0.5_DP * ( p_DcoordSource (3,ivt1) + p_DcoordSource (3,ivt2) )     
  end do
  
  ivtoffset = rsourceTriangulation%NVT+rsourceTriangulation%NMT

  ! each midpoint of a face produces a vertex on the next level
  do iae=1,rsourceTriangulation%NAT
    p_DcoordDest(:,ivtoffset+iae) = 0
    do ive=1,4
      ! sum up all x y z coordinates
        ivt1 = p_IverticesAtFace(ive,iae)
        p_DcoordDest(1,ivtoffset+iae) = &
        p_DcoordDest(1,ivtoffset+iae) + &
        p_DcoordSource(1,ivt1)
        
        p_DcoordDest(2,ivtoffset+iae) = &
        p_DcoordDest(2,ivtoffset+iae) + &
        p_DcoordSource(2,ivt1)
        
        p_DcoordDest(3,ivtoffset+iae) = &
        p_DcoordDest(3,ivtoffset+iae) + &
        p_DcoordSource(3,ivt1)
    end do ! end ive
    
    p_DcoordDest(1,ivtoffset+iae) = &
    p_DcoordDest(1,ivtoffset+iae) * 0.25_DP
    p_DcoordDest(2,ivtoffset+iae) = &
    p_DcoordDest(2,ivtoffset+iae) * 0.25_DP
    p_DcoordDest(3,ivtoffset+iae) = &
    p_DcoordDest(3,ivtoffset+iae) * 0.25_DP
  end do ! end iae
     
      
  ! Allocate memory for IverticesAtElement and get a pointer to it.
  ! Fill the array with zero, so we won't have problems when mixing
  ! triangles into a quad mesh.
  nve = rsourceTriangulation%NNVE
  Isize = (/nve,int(rdestTriangulation%NEL,I32)/)
  call storage_new2D ('tria_refineMesh2lv3D', 'KVERT', Isize, ST_INT, &
       rdestTriangulation%h_IverticesAtElement, ST_NEWBLOCK_ZERO)
  call storage_getbase_int2D(&
       rdestTriangulation%h_IverticesAtElement,p_IvertAtElementDest)
    
  ! increase the offset.
  ivtoffset = rsourceTriangulation%NVT + rsourceTriangulation%NMT + &
              rsourceTriangulation%NAT
              
  do iel = 1,rsourceTriangulation%NEL
        
    ! New element midpoint
    ivtoffset = ivtoffset+1
           
    ! Sum up the coordinates of the corners to get the midpoint.
    x = 0.0_DP
    y = 0.0_DP
    z = 0.0_DP
    do ive = 1,nve
      x = x + p_DcoordSource (1,p_IvertAtElementSource(ive,iel))
      y = y + p_DcoordSource (2,p_IvertAtElementSource(ive,iel))
      z = z + p_DcoordSource (3,p_IvertAtElementSource(ive,iel))
    end do

    ! Store the midpoint
    p_DcoordDest(1,ivtoffset) = x * 0.125_DP
    p_DcoordDest(2,ivtoffset) = y * 0.125_DP
    p_DcoordDest(3,ivtoffset) = z * 0.125_DP
        
  end do
      
      
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
      ! Let's loop over the vertices on the coarse grid. Everyone produces
      ! four elements on the fine grid:
      !
      ! In nquads we count the number of quads +NVT+NMT we reach.
      ! That's the number of the midpoint of that element!
      
   do iel = 1,rsourceTriangulation%NEL
      
     ! Determine number of subelements.
     iel1 = rsourceTriangulation%NEL+7*(iel-1)+1
     iel2 = iel1+1
     iel3 = iel1+2
     iel4 = iel1+3
     iel5 = iel1+4
     iel6 = iel1+5
     iel7 = iel1+6
     
     ! Save them
     p_IrefinementPatch(1+(iel-1)*8) = iel
     p_IrefinementPatch(1+(iel-1)*8+1) = iel1
     p_IrefinementPatch(1+(iel-1)*8+2) = iel2
     p_IrefinementPatch(1+(iel-1)*8+3) = iel3
     p_IrefinementPatch(1+(iel-1)*8+4) = iel4
     p_IrefinementPatch(1+(iel-1)*8+5) = iel5
     p_IrefinementPatch(1+(iel-1)*8+6) = iel6
     p_IrefinementPatch(1+(iel-1)*8+7) = iel7

     
     ! the vertex number of the midpoint is
     ! the old face number
     midpointFaceA = p_IfacesAtElement(1,iel)
     midpointFaceB = p_IfacesAtElement(2,iel)
     midpointFaceC = p_IfacesAtElement(3,iel)
     midpointFaceD = p_IfacesAtElement(4,iel)
     midpointFaceE = p_IfacesAtElement(5,iel)
     midpointFaceF = p_IfacesAtElement(6,iel)
     
     ! In nquads we count the number of the quad +NVT+NMT+NAT we process.
     ! That's the number of the midpoint of that element!
     ! As we reached a new quad, we increase nquads
     
     ! nquads should be the index of the vertex that is the midpoint of
     ! the old hexahedron
     midPointOfIel = rsourceTriangulation%NVT + rsourceTriangulation%NMT &
            + rsourceTriangulation%NAT + iel
          
     ! Step 1: Initialise IverticesOnElement for element IEL
     ! the first vertex of the original hexahedron
     p_IvertAtElementDest(1,iel) = p_IvertAtElementSource (1,iel)
     ! the index of the first edge is now the number of the 2nd vertex
     p_IvertAtElementDest(2,iel) = p_IedgesAtElementSource (1,iel)
     ! the midpoint of face A is the 3rd vertex 
     p_IvertAtElementDest(3,iel) = midpointFaceA
     ! the index of the 4th edge is now the index of the 4th vertex
     p_IvertAtElementDest(4,iel) = p_IedgesAtElementSource (4,iel)
     
     ! the index of the 5th edge is now the index of the 5th vertex     
     p_IvertAtElementDest(5,iel) = p_IedgesAtElementSource (5,iel)
     ! the index of the 6th vertex is the index of the midpoint of face B
     p_IvertAtElementDest(6,iel) = midpointFaceB
     ! the index of the 7th vertex is the index of the midpoint of the old
     ! hexahedron
     p_IvertAtElementDest(7,iel) = midPointOfIel
     ! the index of the 8th vertex is the index of the midpoint of face E
     p_IvertAtElementDest(8,iel) = midpointFaceE

         
     ! Step 2: Initialise IverticesOnElement for element IEL1
     p_IvertAtElementDest(1,iel1) = p_IvertAtElementSource (2,iel)
     p_IvertAtElementDest(2,iel1) = p_IedgesAtElementSource (2,iel)
     p_IvertAtElementDest(3,iel1) = midPointFaceA
     p_IvertAtElementDest(4,iel1) = p_IedgesAtElementSource (1,iel)

     p_IvertAtElementDest(5,iel1) = p_IedgesAtElementSource (6,iel)
     p_IvertAtElementDest(6,iel1) = midpointFaceC
     p_IvertAtElementDest(7,iel1) = midPointOfIel
     p_IvertAtElementDest(8,iel1) = midpointFaceB
        
     ! Step 3: Initialise IverticesOnElement for element IEL2
     p_IvertAtElementDest(1,iel2) = p_IvertAtElementSource (3,iel)
     p_IvertAtElementDest(2,iel2) = p_IedgesAtElementSource (3,iel)
     p_IvertAtElementDest(3,iel2) = midpointFaceA
     p_IvertAtElementDest(4,iel2) = p_IedgesAtElementSource (2,iel)
     
     p_IvertAtElementDest(5,iel2) = p_IedgesAtElementSource (7,iel)
     p_IvertAtElementDest(6,iel2) = midpointFaceD
     p_IvertAtElementDest(7,iel2) = midPointOfIel
     p_IvertAtElementDest(8,iel2) = midpointFaceC
     
     ! Step 4: Initialise IverticesOnElement for element IEL3
     p_IvertAtElementDest(1,iel3) = p_IvertAtElementSource (4,iel)
     p_IvertAtElementDest(2,iel3) = p_IedgesAtElementSource (4,iel)
     p_IvertAtElementDest(3,iel3) = midpointFaceA
     p_IvertAtElementDest(4,iel3) = p_IedgesAtElementSource (3,iel)

     p_IvertAtElementDest(5,iel3) = p_IedgesAtElementSource (8,iel)
     p_IvertAtElementDest(6,iel3) = midpointFaceE
     p_IvertAtElementDest(7,iel3) = midPointOfIel
     p_IvertAtElementDest(8,iel3) = midpointFaceD

     ! Step 5: Initialise IverticesOnElement for element IEL
     p_IvertAtElementDest(1,iel4) = p_IvertAtElementSource (5,iel)
     p_IvertAtElementDest(2,iel4) = p_IedgesAtElementSource (9,iel)
     p_IvertAtElementDest(3,iel4) = midpointFaceF
     p_IvertAtElementDest(4,iel4) = p_IedgesAtElementSource (12,iel)

     p_IvertAtElementDest(5,iel4) = p_IedgesAtElementSource (5,iel)
     p_IvertAtElementDest(6,iel4) = midpointFaceB
     p_IvertAtElementDest(7,iel4) = midPointOfIel
     p_IvertAtElementDest(8,iel4) = midpointFaceE
         
     ! Step 6: Initialise IverticesOnElement for element IEL1
     p_IvertAtElementDest(1,iel5) = p_IvertAtElementSource (6,iel)
     p_IvertAtElementDest(2,iel5) = p_IedgesAtElementSource (10,iel)
     p_IvertAtElementDest(3,iel5) = midpointFaceF
     p_IvertAtElementDest(4,iel5) = p_IedgesAtElementSource (9,iel)
     
     p_IvertAtElementDest(5,iel5) = p_IedgesAtElementSource (6,iel)
     p_IvertAtElementDest(6,iel5) = midpointFaceC
     p_IvertAtElementDest(7,iel5) = midPointOfIel
     p_IvertAtElementDest(8,iel5) = midpointFaceB
     
        
     ! Step 7: Initialise IverticesOnElement for element IEL2
     p_IvertAtElementDest(1,iel6) = p_IvertAtElementSource (7,iel)
     p_IvertAtElementDest(2,iel6) = p_IedgesAtElementSource (11,iel)
     p_IvertAtElementDest(3,iel6) = midpointFaceF
     p_IvertAtElementDest(4,iel6) = p_IedgesAtElementSource (10,iel)
     
     p_IvertAtElementDest(5,iel6) = p_IedgesAtElementSource (7,iel)
     p_IvertAtElementDest(6,iel6) = midpointFaceD
     p_IvertAtElementDest(7,iel6) = midPointOfIel
     p_IvertAtElementDest(8,iel6) = midpointFaceC
     
     ! Step 8: Initialise IverticesOnElement for element IEL3
     p_IvertAtElementDest(1,iel7) = p_IvertAtElementSource(8,iel)
     p_IvertAtElementDest(2,iel7) = p_IedgesAtElementSource(12,iel)
     p_IvertAtElementDest(3,iel7) = midpointFaceF
     p_IvertAtElementDest(4,iel7) = p_IedgesAtElementSource (11,iel)
     
     p_IvertAtElementDest(5,iel7) = p_IedgesAtElementSource (8,iel)
     p_IvertAtElementDest(6,iel7) = midpointFaceE
     p_IvertAtElementDest(7,iel7) = midPointOfIel
     p_IvertAtElementDest(8,iel7) = midpointFaceD
     
     
        
   end do ! iel
      
      ! The last step of setting up the raw mesh on the finer level:
      ! Set up InodalProperty. But that's the most easiest thing: Simply
      ! copy the nodal property array from the coarse mesh to the fine mesh.
      ! The nodal information of the edges and faces (number NVT+1...NVT+NMT+NAT)
      ! that way converts into the nodal information about the new vertices
      ! on the fine mesh!
   call storage_copy (rsourceTriangulation%h_InodalProperty,&
          rdestTriangulation%h_InodalProperty)
  
    rdestTriangulation%NNEE = rsourceTriangulation%NNEE
    rdestTriangulation%NNAE = rsourceTriangulation%NNAE
    rdestTriangulation%NNVA = rsourceTriangulation%NNVA
  
  end subroutine
  
  
!==================================================================== 
 
!<subroutine> 
  subroutine tria_refineBdry2lv3D(rsourceTriangulation,rdestTriangulation,rboundary)
!<description>  
    ! This routine refines the boundary definition of rsourceTriangulation
    ! according to the 2-level ordering algorithm to generate a new 
    ! IverticesAtBoundary. 
    ! If rboundary is specified, the parameter values of the boundary vertices are 
    ! updated and the coordinates of the boundary points are corrected according 
    ! to the analytic boundary. the boundary vertices are not sorted for their
    ! parameter value!
!</description>  

    ! The source triangulation to be refined
    type(t_triangulation), intent(IN) :: rsourceTriangulation

    ! Destination triangulation structure that receives the refined mesg. 
    type(t_triangulation), intent(INOUT) :: rdestTriangulation
    
    ! OPTIONAL: Defintion of analytic boundary.
    ! If specified, the coordinates of the new boundary vertices are
    ! recomputed according to the analytic boundary.
    type(t_boundary), intent(IN), optional :: rboundary

!</subroutine> 

  ! local variables

      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IvertAtBoundartySource
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IedgesAtBoundartySource
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_InodalPropertySource
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IvertAtBoundartyDest
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_InodalPropertyDest
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IboundaryCpIdxSource
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IboundaryCpIdxDest
      
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IfacesAtBoundary
      integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IedgesAtBoundary
      integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtEdge
      integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtFace
      
      integer :: ivbd,ibct,isize,NMBD,NABD,NMT,NAT
      integer :: NVT,NEL,ivt,NVBD
      
      ! these names are just too long...
      NMBD = rsourceTriangulation%NMBD
      NABD = rsourceTriangulation%NABD
      NVT = rsourceTriangulation%NVT
      NMT = rsourceTriangulation%NMT
      NAT = rsourceTriangulation%NAT
      NEL = rsourceTriangulation%NEL
      
      ! Get the definition of the boundary vertices and -edges.
      call storage_getbase_int (rsourceTriangulation%h_IverticesAtBoundary,&
          p_IvertAtBoundartySource)
      call storage_getbase_int (rsourceTriangulation%h_IedgesAtBoundary,&
          p_IedgesAtBoundartySource)
          
      call storage_getbase_int (rsourceTriangulation%h_InodalProperty,&
          p_InodalPropertySource)

      call storage_getbase_int (rsourceTriangulation%h_IboundaryCpIdx,&
          p_IboundaryCpIdxSource)

      call storage_getbase_int (rsourceTriangulation%h_IfacesAtBoundary,&
          p_IfacesAtBoundary)

      call storage_getbase_int (rsourceTriangulation%h_IedgesAtBoundary,&
          p_IedgesAtBoundary)

      call storage_getbase_int2d (rsourceTriangulation%h_IverticesAtEdge,&
          p_IverticesAtEdge)

      call storage_getbase_int2d (rsourceTriangulation%h_IverticesAtFace,&
          p_IverticesAtFace)


      ! calculate the new value of NVBD
      rdestTriangulation%NVBD =  rsourceTriangulation%NVBD + &
                                 rsourceTriangulation%NMBD + &      
                                 rsourceTriangulation%NABD       
                                                
      ! shorter name                                          
      NVBD = rdestTriangulation%NVBD
                                                
      rdestTriangulation%NBCT = rsourceTriangulation%NBCT 
          
      ! Create new arrays in the fine grid for the vertices and indices.
      CALL storage_new ('tria_refineBdry2lv2D', &
          'KVBD', INT(rdestTriangulation%NVBD,I32), &
          ST_INT, rdestTriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)

      call storage_getbase_int (rdestTriangulation%h_IverticesAtBoundary,&
          p_IvertAtBoundartyDest)

      
      ! update the h_IverticesAtBoundary 
      
      ! assume we are at lvl "r" then:
      ! there are NVBD boundary vertices from lvl (r-1)
      ! there are NMBD new boundary vertices from the edges
      ! there are NABD new boundary vertices from the faces
      
      CALL storage_getsize (rdestTriangulation%h_InodalProperty,isize)
      IF (isize .NE. rdestTriangulation%NVT) THEN
        ! If the size is wrong, reallocate memory.
        CALL storage_realloc ('tria_refineBdry2lv3D', &
            INT(rdestTriangulation%NVT,I32),&
            rdestTriangulation%h_InodalProperty, &
            ST_NEWBLOCK_NOINIT, .FALSE.)
      END IF      
      
      call storage_getbase_int (rdestTriangulation%h_InodalProperty,&
          p_InodalPropertyDest)
      

      ! copy the nodal properties
      ! the vertex, edge and face properties of the source mesh
      ! are the vertex properties of the destination mesh
      call lalg_copyVectorint(p_InodalPropertySource(1:NVT+NMT+NAT),&
           p_InodalPropertyDest(1:NVT+NMT+NAT))
      
!      p_InodalPropertyDest(1:NVT+NMT+NAT) = &        
!      p_InodalPropertySource(1:NVT+NMT+NAT)

      ! there are still NEL vertices remaining but these
      ! are not boundary vertices so assign zeros there
      call lalg_clearVectorint( &
           p_InodalPropertyDest(NVT+NMT+NAT+1:NVT+NMT+NAT+NEL))
!      p_InodalPropertyDest(NVT+NMT+NAT+1:NVT+NMT+NAT+NEL)=0
      
      
      ! allocate memory
      if(rdestTriangulation%h_IboundaryCpIdx == ST_NOHANDLE) then
        call storage_new ('tria_genFacesAtBoundary', 'IboundaryCpIdx', &
            int(rdestTriangulation%NBCT+1,I32), ST_INT, &
            rdestTriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_NOINIT)
      end if      
                
      ivbd = 1
      do ivt=1,rdestTriangulation%NVT
        if(p_InodalPropertyDest(ivt) > 0) then
          p_IvertAtBoundartyDest(ivbd) = ivt
          ivbd = ivbd + 1
        end if
      end do
      
    call storage_getbase_int (rdestTriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdxDest)
      
    p_IboundaryCpIdxDest(:) = 0  
    ! the first element in p_IboundaryCpIdx is always 1
    p_IboundaryCpIdxDest(1) = 1
  
    ! assign the indices of the boundary vertices
    ! first save the number of vertices in each boundary component in
    ! p_IboundaryCpIdx(2:NBCT+1)
    do ivt = 1,rdestTriangulation%NVT
      if(p_InodalPropertyDest(ivt) /= 0) then
          ibct = p_InodalPropertyDest(ivt)
          p_IboundaryCpIdxDest(ibct+1) = p_IboundaryCpIdxDest(ibct+1) + 1
      end if
    end do
 
    ! now create the actual index array
    do ibct = 2, rdestTriangulation%NBCT+1
        p_IboundaryCpIdxDest(ibct) = p_IboundaryCpIdxDest(ibct)+p_IboundaryCpIdxDest(ibct-1)
    end do      
  
  end subroutine ! end tria_refineBdry2lv3D

!====================================================================

  function tria_BinSearch(p_Iarray, Ivalue, Ilbound, Iubound)
  
    integer(PREC_VERTEXIDX), dimension(:),intent(IN) :: p_Iarray
    integer(PREC_VERTEXIDX) :: Ivalue
    integer(PREC_VERTEXIDX), intent(IN) :: Ilbound, Iubound
    
    ! local variables
    integer :: Imid, tria_BinSearch,Ilboundloc, Iuboundloc
    
    Ilboundloc = Ilbound
    Iuboundloc = Iubound
    
    ! standard binary search scheme...
    do while( Ilboundloc <= Iuboundloc )
    
      Imid = (Ilboundloc + Iuboundloc) / 2
    
      if(p_Iarray(Imid) > Ivalue) then
        Iuboundloc = Imid-1
      else if(p_Iarray(Imid) < Ivalue) then
        Ilboundloc = Imid+1
      else ! found Ivalue
        tria_BinSearch = 1
        return
      end if
    
    end do ! end while
    
    ! Ivalue was not found...
    tria_BinSearch = 0
  
  end function ! end tria_BinSearch

!====================================================================

!<subroutine> 
  subroutine tria_genFacesAtBoundary      (rtriangulation)
!<description>
  ! This routine generates the facesAtBoundary information
  ! by using two conditions:
  !                          1) a boundary face consists solely of boundary vertices
  !                          2) it is not connected to more than one hexa(or whatever element)
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer(I32), dimension(:), pointer :: p_InodalProperty
    integer(PREC_ELEMENTIDX), dimension(:,:), pointer :: p_IelementsAtFace
    integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtFace
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IfacesAtBoundary
    integer(prec_vertexidx), dimension(:), pointer :: p_IboundaryCpFacesIdx    
    integer :: ive
    integer(PREC_ELEMENTIDX) :: iface, inumElements, iel, inumVertices
    integer(PREC_VERTEXIDX) :: iglobIndex, NFBD, findex, iBdyComp
    
    NFBD = 0
    
    ! Is everything here we need?
    if (rtriangulation%h_InodalProperty .EQ. ST_NOHANDLE) then
      call output_line ('InodalPropertys not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtBoundary')
      call sys_halt()
    end if

    if (rtriangulation%h_IelementsAtFace .EQ. ST_NOHANDLE) then
      call output_line ('IelementsAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtBoundary')
      call sys_halt()
    end if
    
    if (rtriangulation%h_IverticesAtFace .EQ. ST_NOHANDLE) then
      call output_line ('IverticesAtFace not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genFacesAtBoundary')
      call sys_halt()
    end if

    ! allocate memory
    if(rtriangulation%h_IboundaryCpFacesIdx == ST_NOHANDLE) then
      call storage_new ('tria_genFacesAtBoundary', 'IboundaryCpFacesIdx', &
          int(rtriangulation%NBCT+1,I32), ST_INT, &
          rtriangulation%h_IboundaryCpFacesIdx, ST_NEWBLOCK_NOINIT)
    end if      


    ! Get the arrays.
    call storage_getbase_int (rtriangulation%h_InodalProperty,p_InodalProperty)
    call storage_getbase_int2D (rtriangulation%h_IelementsAtFace,p_IelementsAtFace)
    call storage_getbase_int2D (rtriangulation%h_IverticesAtFace,p_IverticesAtFace)
    call storage_getbase_int (rtriangulation%h_IboundaryCpFacesIdx,p_IboundaryCpFacesIdx)
    
    ! the first index is one
    p_IboundaryCpFacesIdx(1)=1
    ! the remaining indices are initialized with zeros
    p_IboundaryCpFacesIdx(2:rtriangulation%NBCT+1)=0
    
    ! loop over all faces
    ! to count the number of faces on the boundary 
    ! we save this number in NFBD
    do iface=1,rtriangulation%NAT
    
      ! used to count the number of vertices
      inumVertices = 0
    
      ! count the number of boundary vertices at this face
      do ive = 1,4
        iglobIndex = p_IverticesAtFace(ive,iface)
        
        if(p_InodalProperty(iglobIndex) > 0) then
          inumVertices = inumVertices + 1
          ! save the number of the boundary component
          iBdyComp = p_InodalProperty(iglobIndex)
        end if
        
      end do ! end ive
    
      ! used to count the number of adjacent elements
      inumElements = 0
      
      ! count the number of adjacent elements
      do iel=1,2
        if(p_IelementsAtFace(iel,iface) > 0) then
          inumElements = inumElements + 1
        end if
      end do ! end iel

      if(inumVertices == 4 .and. inumElements < 2) then
      ! boundary face
      ! increase the number of entries for tthe corresponding
      ! boundary component
      p_IboundaryCpFacesIdx(iBdyComp+1)=&
      p_IboundaryCpFacesIdx(iBdyComp+1) + 1  
      else
      ! inner face
      end if
      
    end do ! end iface
    
    ! add up the entries to compute the actual index array
    do ive=2,rtriangulation%NBCT+1
    p_IboundaryCpFacesIdx(ive) = &
    p_IboundaryCpFacesIdx(ive) + p_IboundaryCpFacesIdx(ive-1)
    end do ! end ive
    
    ! number of faces on the boundary
    NFBD = p_IboundaryCpFacesIdx(rtriangulation%NBCT+1) - 1

    ! assign the number of faces on the boundary
    rtriangulation%NABD = NFBD
    
    ! shift the indices... we do the old trick
    p_IboundaryCpFacesIdx(2:rtriangulation%NBCT+1) = &
    p_IboundaryCpFacesIdx(1:rtriangulation%NBCT)

    ! allocate memory for the p_IfacesAtBoundary array
    if(rtriangulation%h_IfacesAtBoundary == ST_NOHANDLE) then
      call storage_new ('tria_genFacesAtBoundary', 'IfacesAtBoundary', &
          int(NFBD,I32), ST_INT, &
          rtriangulation%h_IfacesAtBoundary, ST_NEWBLOCK_NOINIT)
    end if      
    
    ! woohoo memory is allocates so get the pointer
    call storage_getbase_int (rtriangulation%h_IfacesAtBoundary,p_IfacesAtBoundary)
    
    ! auxilliary index that points to the current position in the
    ! p_IfacesAtBoundary array
    findex = 0
    
    ! loop over all faces
    do iface=1,rtriangulation%NAT
    
      ! used to count the number of vertices
      inumVertices = 0
    
      ! count the number of boundary vertices at this face
      do ive = 1,4
        iglobIndex = p_IverticesAtFace(ive,iface)
        
        if(p_InodalProperty(iglobIndex) > 0) then
          inumVertices = inumVertices + 1
        end if
        
      end do ! end ive
    
      ! used to count the number of adjacent elements
      inumElements = 0
      
      ! count the number of adjacent elements
      do iel=1,2
        if(p_IelementsAtFace(iel,iface) > 0) then
          inumElements = inumElements + 1
        end if
      end do ! end iel
      
      if(inumVertices == 4 .and. inumElements < 2) then
        ! we have identified a boundary face
        ! get the boundary component
        iBdyComp = p_InodalProperty(iglobIndex)
        ! so write its number in the array
        findex = p_IboundaryCpFacesIdx(iBdyComp+1)
        p_IfacesAtBoundary(findex) = iface
        ! increase the free position in the index array
        p_IboundaryCpFacesIdx(iBdyComp+1) = &
        p_IboundaryCpFacesIdx(iBdyComp+1) + 1
      else
      ! inner face
      end if
      
    end do ! end iface    
    
  end subroutine ! end tria_genFacesAtBoundary

!==================================================================== 

!<subroutine>  
  subroutine tria_genEdgesAtBoundary3d(rtriangulation)
!<description>
  ! This routine generates the edgesAtBoundary information
  ! by using a simple condition:
  ! -------a boundary edge has a boundary face attached-------
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! Local variables
    integer(I32), dimension(:), pointer :: p_InodalProperty
    
    integer(prec_vertexidx), dimension(:), pointer :: p_IfacesAtEdgeIdx
    integer(prec_vertexidx), dimension(:), pointer :: p_IfacesAtEdge
    integer(prec_vertexidx), dimension(:,:), pointer :: p_IverticesAtEdge
    
    integer(prec_vertexidx), dimension(:), pointer :: p_IedgesAtBoundary
    integer(prec_vertexidx), dimension(:), pointer :: p_Iaux
    integer(prec_vertexidx), dimension(:), pointer :: p_IfacesAtBoundary
    integer(prec_vertexidx), dimension(:), pointer :: p_IboundaryCpEdgesIdx
    integer :: ive, h_IbdyComponents, ibct, ibdyFace, inotFound
    integer(prec_elementidx) :: iface,iface1,iface2, ifaceIndex
    integer(prec_vertexidx) :: NMBD, eIndex, iBdyComp
    integer(prec_vertexidx), dimension(:), pointer :: p_IbdyComponents
    integer(prec_vertexidx), dimension(:), pointer :: p_IboundaryCpFacesIdx        
    
    NMBD = 0
    
    ! Is everything here we need?
    if (rtriangulation%h_InodalProperty .EQ. ST_NOHANDLE) then
      call output_line ('InodalPropertys not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3d')
      call sys_halt()
    end if

    if (rtriangulation%h_IverticesAtEdge .EQ. ST_NOHANDLE) then
      call output_line ('IverticesAtEdge not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3d')
      call sys_halt()
    end if


    if (rtriangulation%h_IfacesAtBoundary .EQ. ST_NOHANDLE) then
      call output_line ('IfacesAtBoundary not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3d')
      call sys_halt()
    end if


    if (rtriangulation%h_IfacesAtEdgeIdx .EQ. ST_NOHANDLE) then
      call output_line ('IfacesAtEdgeIdx not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3d')
      call sys_halt()
    end if

    if (rtriangulation%h_IfacesAtEdge .EQ. ST_NOHANDLE) then
      call output_line ('IfacesAtEdge not available!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'tria_genEdgesAtBoundary3d')
      call sys_halt()
    end if
    
    ! allocate memory
    if(rtriangulation%h_IboundaryCpEdgesIdx == ST_NOHANDLE) then
      call storage_new ('tria_genEdgesAtBoundary3d', 'IboundaryCpEdgesIdx', &
          int(rtriangulation%NBCT+1,I32), ST_INT, &
          rtriangulation%h_IboundaryCpEdgesIdx, ST_NEWBLOCK_NOINIT)
    end if      

    ! Get the arrays.
    call storage_getbase_int (rtriangulation%h_InodalProperty,p_InodalProperty)
    call storage_getbase_int (rtriangulation%h_IfacesAtBoundary,p_IfacesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IfacesAtEdgeIdx,p_IfacesAtEdgeIdx)
    call storage_getbase_int (rtriangulation%h_IfacesAtEdge,p_IfacesAtEdge)
    call storage_getbase_int (rtriangulation%h_IboundaryCpEdgesIdx, &
    p_IboundaryCpEdgesIdx)
    call storage_getbase_int2d(rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
    call storage_getbase_int (rtriangulation%h_IboundaryCpFacesIdx,p_IboundaryCpFacesIdx)    
    
    
    
    ! create an auxilliary array
    ! the auxilliary is (1:NAT) and
    ! and p_Iaux(i) = 0 if face i is an inner
    ! face and greater 0 if i is not a boundary face
    allocate(p_Iaux(rtriangulation%NAT))
    
    ! initialise it with zeros
    p_Iaux(:)=0

    ! another auxilliary array which is used as a pointer to
    ! the free positions for the different boundary components
    h_IbdyComponents = ST_NOHANDLE
    
    call storage_copy(rtriangulation%h_IboundaryCpFacesIdx, &
         h_IbdyComponents)
    
    call storage_getbase_int (h_IbdyComponents, &
         p_IbdyComponents)
    
    ! the first index is one
    p_IboundaryCpEdgesIdx(1)=1
    ! the remaining indices are initialized with zeros
    p_IboundaryCpEdgesIdx(2:rtriangulation%NBCT+1)=0
    
    ! iface1 points to the position in the
    ! facesAtBoundary array
    iface1=1
    
    ! in the first loop over the faces we
    ! count the number of edges on the boundary
    ! and save it as NMBD
    ! we also build the auxilliary array
    do iface=1,rtriangulation%NAT
    
      ! check if the maximum number of faces on the
      ! boundary is reached
      if(iface1 > rtriangulation%NABD) exit
      
      inotFound = 0
      ! check for this face in every boundary component
      ! until we have found it
      do ibct=1,rtriangulation%NBCT
        ! check if face iface a boundary face
        ibdyFace = p_IbdyComponents(ibct)
        ! if we are at or beyond the end of that bdy component
        ! then skip it
        if(ibdyFace >= p_IboundaryCpFacesIdx(ibct+1)) cycle
        if(p_IfacesAtBoundary(ibdyFace) == iface) then
          ! we have identified a boundary face
          p_Iaux(iface) = 1
          ! increase the number of boundary faces found
          iface1 = iface1 + 1
          ! increase the pointer to the boundary comp index
          p_IbdyComponents(ibct) = p_IbdyComponents(ibct) + 1
          exit ! break out of the ibct loop 
        else
          ! increase the not found counter
          inotFound = inotFound + 1
        end if
      end do ! end do ibct
      
      ! if not found
      if(inotFound == rtriangulation%NBCT) then
        p_Iaux(iface) = 0
      end if
    
    end do ! end iface
    
    ! check all edges
    do ive=1,rtriangulation%NMT
        
      ! we get the start and end indices into facesAtEdge array  
      iface1 = p_IfacesAtEdgeIdx(ive)
      iface2 = p_IfacesAtEdgeIdx(ive+1) - 1
      
      ! check all connected faces
      do iface=iface1,iface2
        ifaceIndex = p_IfacesAtEdge(iface)
        if(p_Iaux(ifaceIndex) == 1) then
          ! the edge is connected to a boundary face
          ! so it is a boundary edge
          ! increase the number of edges on the boundary
          NMBD = NMBD + 1
          ! get the boundary component number of 
          ! a vertex of the edge
          iBdyComp = p_IverticesAtEdge(1,ive)
          iBdyComp = p_InodalProperty(iBdyComp)
          p_IboundaryCpEdgesIdx(iBdyComp+1) = &
          p_IboundaryCpEdgesIdx(iBdyComp+1) + 1
          ! break out of this loop
          exit
        end if
      end do ! end iface
    
    end do ! end ive
    
    ! assign the number of faces on the boundary
    rtriangulation%NMBD = NMBD

    ! allocate memory
    if(rtriangulation%h_IedgesAtBoundary == ST_NOHANDLE) then
      call storage_new ('tria_genEdgesAtBoundary3d', 'IedgesAtBoundary', &
          int(NMBD,I32), ST_INT, &
          rtriangulation%h_IedgesAtBoundary, ST_NEWBLOCK_NOINIT)
    end if      
    
    ! get the pointer
    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,p_IedgesAtBoundary)
    
    eIndex = 0
    
    ! add up the entries to compute the actual index array
    do ive=2,rtriangulation%NBCT+1
    p_IboundaryCpEdgesIdx(ive) = &
    p_IboundaryCpEdgesIdx(ive) + p_IboundaryCpEdgesIdx(ive-1)
    end do ! end ive
    
    ! shift the indices...
    p_IboundaryCpEdgesIdx(2:rtriangulation%NBCT+1) = &
    p_IboundaryCpEdgesIdx(1:rtriangulation%NBCT)
    
    ! check all edges
    do ive=1,rtriangulation%NMT
    
      ! we need to get the indices of the faces
      ! connected to the edge
      iface1 = p_IfacesAtEdgeIdx(ive)
      iface2 = p_IfacesAtEdgeIdx(ive+1) - 1
      
      ! check all connected faces
      do iface=iface1,iface2
        ifaceIndex = p_IfacesAtEdge(iface)
        if(p_Iaux(ifaceIndex) == 1) then
          ! we have identified a boundary edge
          ! now get the corresponding boundary component          
          iBdyComp = p_IverticesAtEdge(1,ive)
          iBdyComp = p_InodalProperty(iBdyComp)
          
          eIndex = p_IboundaryCpEdgesIdx(iBdyComp+1)
          p_IedgesAtBoundary(eIndex) = ive
          ! increase the pointer into the edgesAtBoundary array
          p_IboundaryCpEdgesIdx(iBdyComp+1) = &
          p_IboundaryCpEdgesIdx(iBdyComp+1) + 1
          
          exit  
        end if
      end do ! end iface
    
    end do ! end iface
    
    
    ! free memory
    deallocate(p_Iaux)
    call storage_free(h_IbdyComponents)
    
  end subroutine ! end tria_genEdgesAtBoundary3d
  
!==================================================================== 

!<subroutine>
  subroutine tria_genEdgeNodalProperty3d(rtriangulation)
!<description>
  ! This routine generates the nodalproperty for the nodes that will be
  ! formed by the current set of edges.
  ! This is done by simply checking the pregenerated edgesAtBoundary
  ! information
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>

    ! local variables

    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_InodalProperty
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IboundaryCpEdgesIdx

    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IverticesAtBoundary      
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IfacesAtBoundary
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IedgesAtBoundary
    
    integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtEdge
      
    integer :: isize,IboundaryComponent,ive,NMBD,NABD,NMT,NAT
    integer :: NVT,NEL,NBCT, IcpIdx1, IcpIdx2
      
    ! these names are just too long...
    NMBD = rtriangulation%NMBD
    NABD = rtriangulation%NABD
    NVT = rtriangulation%NVT
    NMT = rtriangulation%NMT
    NAT = rtriangulation%NAT
    NEL = rtriangulation%NEL
      
    ! Get the definition of the boundary vertices and -edges.
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
        
    call storage_getbase_int (rtriangulation%h_IfacesAtBoundary,&
        p_IfacesAtBoundary)

    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
        
    call storage_getbase_int (rtriangulation%h_IboundaryCpEdgesIdx,&
        p_IboundaryCpEdgesIdx)
        
    call storage_getbase_int2d(rtriangulation%h_IverticesAtEdge,&
        p_IverticesAtEdge)
                                                
    NBCT = rtriangulation%NBCT 
    
    ! we check the size of the p_InodalProperty array
    ! its size should be NVT+NMT+NAT, because in the nodal
    ! property array we ant to store the properties of
    ! NVT vertices, NMT edges and NAT faces
    CALL storage_getsize (rtriangulation%h_InodalProperty,isize)
    IF (isize .NE. NVT+NMT+NAT) THEN
      ! If the size is wrong, reallocate memory.
      CALL storage_realloc ('tria_genEdgeNodalProperty3d', &
          INT(NVT+NMT+NAT,I32),&
          rtriangulation%h_InodalProperty,ST_NEWBLOCK_NOINIT)
    END IF      
    
    ! get the pointer
    call storage_getbase_int (rtriangulation%h_InodalProperty,&
        p_InodalProperty)


    ! Loop over all edges to see if the current edge ive is on 
    ! the boundary, if it is on the boundary we write the
    ! number of the boundary component on position NVT+ive 
    ! in the p_InodalProperty array
    do ive=1,rtriangulation%NMT
    
      ! get the boundary component index of the vertices
      ! that build the edge
      IboundaryComponent = p_IverticesAtEdge(1,ive)
      
      ! get the boundary component number  
      IboundaryComponent = p_InodalProperty( &
      IboundaryComponent)
      
      ! if the vertex is an inner vertex
      ! we can skip this edge
      if(IboundaryComponent == 0) then
        p_InodalProperty(rtriangulation%NVT+ive) = 0
        cycle
      end if
      
      IcpIdx1 = p_IboundaryCpEdgesIdx(IboundaryComponent)
      IcpIdx2 = p_IboundaryCpEdgesIdx(IboundaryComponent+1)-1
    
      ! check if edge is on border
      if(tria_BinSearch(p_IedgesAtBoundary,ive,IcpIdx1,IcpIdx2)&
          ==1) then
         ! the edge is on the boundary
         ! write that number to the p_InodalProperty array 
         p_InodalProperty(rtriangulation%NVT+ive) = &
         IboundaryComponent
       else
         ! the edge is an inner edge assign a zero
         p_InodalProperty(rtriangulation%NVT+ive) = 0
       end if
         
    end do ! end ive
           
  end subroutine ! end tria_genEdgeNodalProperty3d
  
!====================================================================
  
  subroutine tria_genFaceNodalProperty3d(rtriangulation)
!<description>
  ! This routine generates the faceNodal property information
  ! by simply checking the pregenerated facesAtBoundary array
!</description>

!<inputoutput>
  ! The triangulation structure to be updated.
  type(t_triangulation), intent(INOUT) :: rtriangulation
!</inputoutput>
  
!</subroutine>
      
    ! local variables
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_InodalProperty
    
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IboundaryCpFacesIdx    

    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IverticesAtBoundary      
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IfacesAtBoundary
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IedgesAtBoundary
    
    integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtFace
      
    integer :: isize,IboundaryComponent,NMBD,NABD,iface,NMT,NAT
    integer :: NVT,NEL,NBCT, IcpIdx1, IcpIdx2
      
    ! these names are just too long...
    NMBD = rtriangulation%NMBD
    NABD = rtriangulation%NABD
    NVT = rtriangulation%NVT
    NMT = rtriangulation%NMT
    NAT = rtriangulation%NAT
    NEL = rtriangulation%NEL
      
    ! Get the definition of the boundary vertices and -edges.
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary,&
        p_IverticesAtBoundary)
        
    call storage_getbase_int (rtriangulation%h_IfacesAtBoundary,&
        p_IfacesAtBoundary)

    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
        
    call storage_getbase_int (rtriangulation%h_IboundaryCpFacesIdx,&
        p_IboundaryCpFacesIdx)
        
    call storage_getbase_int2d(rtriangulation%h_IverticesAtFace,&
        p_IverticesAtFace)

    NBCT = rtriangulation%NBCT 
    
    ! we check the size of the p_InodalProperty array
    ! its size should be NVT+NMT+NAT, because in the nodal
    ! property array we ant to store the properties of
    ! NVT vertices, NMT edges and NAT faces
    CALL storage_getsize (rtriangulation%h_InodalProperty,isize)
    IF (isize .NE. NVT+NMT+NAT) THEN
      ! If the size is wrong, reallocate memory.
      CALL storage_realloc ('tria_genEdgeNodalProperty3d', &
          INT(NVT+NMT+NAT,I32),&
          rtriangulation%h_InodalProperty, &
          ST_NEWBLOCK_NOINIT)
    END IF      
    
    call storage_getbase_int (rtriangulation%h_InodalProperty,&
        p_InodalProperty)

    ! Loop over all faces and check if the current face iface
    ! is on the boundary, if it is on the boundary the
    ! entry NVT+NMT+iface in the p_InodalProperty array will
    ! contain the boundary component of the face iface
    do iface=1,rtriangulation%NAT
    
      ! get the boundary component index of the vertices
      ! that build the face
      IboundaryComponent = p_IverticesAtFace(1,iface)
         
      ! get the boundary component number from the
      ! p_InodalProperty array
      IboundaryComponent = p_InodalProperty( &
      IboundaryComponent)
      
      ! if the vertex is an inner vertex
      ! we can skip this face
      if(IboundaryComponent == 0) then
        p_InodalProperty(rtriangulation%NVT+NMT+iface) = 0
        cycle
      end if
      
      IcpIdx1 = p_IboundaryCpFacesIdx(IboundaryComponent)
      IcpIdx2 = p_IboundaryCpFacesIdx(IboundaryComponent+1)-1
      
      ! check if face is on border
      if(tria_BinSearch(p_IfacesAtBoundary,iface,IcpIdx1,IcpIdx2)&
         ==1) then
         
         ! write the boundary component number to the
         ! position NVT+NMT+iface
         p_InodalProperty(rtriangulation%NVT+NMT+iface) = &
         IboundaryComponent
       else
         ! the face is an inner face so assign zero
         p_InodalProperty(rtriangulation%NVT+NMT+iface) = 0
       end if
      
    end do ! end iface
  
  end subroutine ! end tria_genFaceNodalProperty3d
  
END MODULE



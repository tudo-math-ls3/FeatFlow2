!##############################################################################
!# ****************************************************************************
!# <name> ucd </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains various routines to export solution vector(s) to
!# visualisation file formats like GMV, AVS or MATLAB.
!# UCD stands for 'unstructured cell data'. The term is originally used for
!# describing the export format of AVS, but the approach how to provide
!# information for external postprocessing tools can also be generalised
!# into a much more abstract manner.
!#
!# The module contains the following subroutines:
!#
!#  1.) ucd_startGMV
!#      -> Starts a new GMV output file
!#
!#  2.) ucd_startAVS
!#      -> Starts a new AVS output file
!#
!#  3.) ucd_startVTK
!#      -> Starts a new VTK output file
!#
!#  4.) ucd_setAlternativeSource
!#      -> Allows to specify an alternative source file for parts of a
!#         triangulation or the whole triangulation. (Saves disc space
!#         if a sequence of output files are written.)
!#
!#  5.) ucd_addVariableVertexBased / ucd_addVariableElementBased
!#      -> Add a vertex- or element-based data vector to the output.
!#
!#  6.) ucd_addVarVertBasedVec / ucd_addVarElemBasedVec
!#      -> Add a vertex/element-based vector to the output.
!#
!#  7.) ucd_setVertexMaterial
!#      -> Allows to specify for every vertex/node (corner vertex, edge 
!#         midpoint,...) a material id.
!#
!#  8.) ucd_setCellMaterial
!#      -> Allows to specify for every cell a material id.
!#
!#  9.) ucd_setMaterials
!#      -> Allows to specify material names for the material id's that
!#         are assigned to cells/vertices in ucd_setVertexMaterial and
!#         ucd_setCellMaterial, respectively.
!#
!# 10.) ucd_setTracers
!#      -> Specify tracers by (X,Y) or (X,Y,Z) coordinates
!#
!# 11.) ucd_addTracerVariable
!#      -> Adds a tracer variable, i.e. a named value for every tracer.
!#
!# 12.) ucd_removeTracers
!#      -> Removes all tracer information.
!#
!# 13.) ucd_setSimulationTime
!#      -> Specify simulation time.
!#
!# 14.) ucd_addPolygon
!#      -> Adds a point-array as polygon of line-segments to the output.
!#
!# 15.) ucd_addCommentLine
!#      -> Add a comment line.
!#
!# 16.) ucd_addParameterList
!#      -> Add a configuration parameter list as comments to the output.
!#
!# 17.) ucd_setOutputNumberFormat
!#      -> Allows to specify a custom format of double precision numbers 
!#         in the output file.
!#
!# 18.) ucd_write
!#      -> Writes the output file to the hard disc, closes the file.
!#
!# 19.) ucd_release
!#      -> Releases a UCD export structure that was created by ucd_startXXXX.
!#
!# 20.) ucd_setFilename
!#      -> Change the filename of the output file.
!#
!# 21.) ucd_readGMV
!#      -> Read a GMV file.
!#
!# 22.) ucd_getVariable
!#      -> Get a data / information of a variable.
!#
!# To write a file (GMV, AVS, VTK, MATLAB or what else), one has to use the
!# following sequence of commands:
!#
!# a) ucd_startXXXX     - Create an output structure for a file
!#
!# b) ucd_addVariableVertexBased  
!#                      - Add a variable consisting of values in vertices
!#    ucd_addVariableElementBased 
!#                      - Add a variable consisting of values in elements
!#    ucd_addVarVertBasedVec
!#                      - Add a vector variable consisting of values
!#                        in vertices
!#
!# c) ucd_setAlternativeSource 
!#                      - Optional: Set an alternative input file for parts
!#                        of the mesh or the whole mesh. 
!#
!# d) ucd_write         - Write all pending data and close the file
!#
!# e) ucd_release       - Release the export structure
!#
!# ucd_startXXXX creates an export structure that identifies the output
!# file. This is the only point, where the caller must specify the type
!# of output. Afterwards, the ucd_setXXXX and ucd_addXXXX routines
!# can be used to specify all the output.\\
!#
!# Depending on the type of file (AVS, GMV or what else), the output may 
!# be written to the disc directly or first collected in memory. 
!# The final ucd_write command then finishes the output process 
!# and creates a valid GMV/AVS/...-file.\\
!#
!# The ucd_release command at the end must be used to release all memory 
!# that was used by the ucd_setXXXX and ucd_addXXXX routines.
!# </purpose>
!##############################################################################

MODULE ucd

  USE fsystem
  USE triangulation
  USE paramlist

  IMPLICIT NONE

!<constants>

!<constantblock description="Output format types">

  ! No output format
  INTEGER, PARAMETER :: UCD_FORMAT_NONE = 0
  
  ! GMV file format (ASCII)
  INTEGER, PARAMETER :: UCD_FORMAT_GMV  = 1
  
  ! AVS/Express file format (ASCII)
  INTEGER, PARAMETER :: UCD_FORMAT_AVS  = 2
  
  ! Visualization Toolkit (VTK) file format (ASCII)
  INTEGER, PARAMETER :: UCD_FORMAT_VTK  = 3

  ! GMV file format (binary)
  INTEGER, PARAMETER :: UCD_FORMAT_BGMV = 4
  
!</constantblock>

!<constantblock description="Data type flags">

  ! Element-based data
  INTEGER, PARAMETER :: UCD_BASE_ELEMENT = 0
  
  ! Vertex-based data
  INTEGER, PARAMETER :: UCD_BASE_VERTEX  = 1
  
!</constantblock>

!<constantblock description="Parameter constants for the VTK exporter">

  ! Export vector components as scalars
  INTEGER, PARAMETER :: UCD_PARAM_VTK_VECTOR_TO_SCALAR = 2**0
  
  ! Use quadratic cell elements when possible
  INTEGER, PARAMETER :: UCD_PARAM_VTK_USE_QUADRATIC    = 2**1

!</constantblock>

!<constantblock description="Constants for the cflags variable in ucd_startXXXX. Bitfield.">

  ! Standard flags. Write information in corner vertices only, linear interpolation.
  INTEGER(I32), PARAMETER :: UCD_FLAG_STANDARD            = 0
  
  ! Construct edge midpoints and write them as vertices to the output file.
  INTEGER(I32), PARAMETER :: UCD_FLAG_USEEDGEMIDPOINTS    = 2**0
  
  ! Construct element midpoints and write them as vertices to the output file.
  INTEGER(I32), PARAMETER :: UCD_FLAG_USEELEMENTMIDPOINTS = 2**1
  
  ! Quadratic buld interpolation. Information is given in corner vertices and
  ! edge midpoints. Implies UCD_FLAG_USEEDGEMIDPOINTS.
  INTEGER(I32), PARAMETER :: UCD_FLAG_BULBQUADRATIC       = 2**2
  
  ! Output of a linear interpolated solution on a once refined mesh. 
  ! Implies UCD_FLAG_USEEDGEMIDPOINTS and UCD_FLAG_USEELEMENTMIDPOINTS. 
  ! Cannot be used with UCD_FLAG_BULBQUADRATIC.
  INTEGER(I32), PARAMETER :: UCD_FLAG_ONCEREFINED         = 2**3

!</constantblock>

!<constantblock description="Specification flags for variables. Bitfield.">
  ! Standard 'scalar' variable.
  INTEGER(I32), PARAMETER :: UCD_VAR_STANDARD             = 0

  ! The variable specifies the X velocity of a velocity field.
  ! Implies UCD_VAR_VERTEXBASED.
  ! Cannot be used together with UCD_VAR_YVELOCITY or UCD_VAR_ZVELOCITY.
  INTEGER(I32), PARAMETER :: UCD_VAR_XVELOCITY            = 2**0

  ! The variable specifies the Y velocity of a velocity field.
  ! Implies UCD_VAR_VERTEXBASED.
  ! Cannot be used together with UCD_VAR_XVELOCITY or UCD_VAR_ZVELOCITY.
  INTEGER(I32), PARAMETER :: UCD_VAR_YVELOCITY            = 2**1

  ! The variable specifies the Z velocity of a velocity field.
  ! Implies UCD_VAR_VERTEXBASED.
  ! Cannot be used together with UCD_VAR_XVELOCITY or UCD_VAR_YVELOCITY.
  INTEGER(I32), PARAMETER :: UCD_VAR_ZVELOCITY            = 2**2
  
!</constantblock>

!<constantblock description="Constants for specifying alternative source files. Bitfield.">

  ! Whole triangulation is specified by a different source file.
  INTEGER(I32), PARAMETER :: UCD_ASRC_ALL                 = NOT(0)

  ! Coordinates of mesh points are specified in a different source file.
  INTEGER(I32), PARAMETER :: UCD_ASRC_POINTS              = 2**0

  ! Cell information is specified in a different source file.
  INTEGER(I32), PARAMETER :: UCD_ASRC_CELLS               = 2**1

  ! Material information is specified in a different source file.
  INTEGER(I32), PARAMETER :: UCD_ASRC_MATERIALS           = 2**2

  ! Polygons are specified in a difference source file.
  INTEGER(I32), PARAMETER :: UCD_ASRC_POLYGONS            = 2**3

!</constantblock>

!<constantblock description="Cell type constants for VTK exporter">

  ! A single vertex
  INTEGER, PARAMETER :: VTK_VERTEX                         =  1
  
  ! A set of vertices
  INTEGER, PARAMETER :: VTK_POLY_VERTEX                    =  2
  
  ! A line
  INTEGER, PARAMETER :: VTK_LINE                           =  3
  
  ! A line strip
  INTEGER, PARAMETER :: VTK_POLY_LINE                      =  4
  
  ! A triangle
  INTEGER, PARAMETER :: VTK_TRIANGLE                       =  5
  
  ! A triangle strip
  INTEGER, PARAMETER :: VTK_TRIANGLE_STRIP                 =  6
  
  ! A polygon
  INTEGER, PARAMETER :: VTK_POLYGON                        =  7
  
  ! A pixel
  INTEGER, PARAMETER :: VTK_PIXEL                          =  8
  
  ! A quadrilateral
  INTEGER, PARAMETER :: VTK_QUAD                           =  9
  
  ! A tetrahedron
  INTEGER, PARAMETER :: VTK_TETRA                          = 10
  
  ! A voxel (cube)
  INTEGER, PARAMETER :: VTK_VOXEL                          = 11
  
  ! A hexahedron
  INTEGER, PARAMETER :: VTK_HEXAHEDRON                     = 12
  
  ! A wedge
  INTEGER, PARAMETER :: VTK_WEDGE                          = 13
  
  ! A pyramid
  INTEGER, PARAMETER :: VTK_PYRAMID                        = 14
  
  ! A quadratic edge
  INTEGER, PARAMETER :: VTK_QUADRATIC_EDGE                 = 21
  
  ! A quadratic triangle
  INTEGER, PARAMETER :: VTK_QUADRATIC_TRIANGLE             = 22
  
  ! A quadratic quadrilateral
  INTEGER, PARAMETER :: VTK_QUADRATIC_QUAD                 = 23
  
  ! A quadratic tetrahedron
  INTEGER, PARAMETER :: VTK_QUADRATIC_TETRA                = 24
  
  ! A quadratic hexahedron
  INTEGER, PARAMETER :: VTK_QUADRATIC_HEXAHEDRON           = 25

!</constantblock>
  
!</constants>

  
!<types>

!<typeblock>
  
  ! UCD export structure. The structure is created by one of the ucd_startXXXX
  ! routines and configured to write output of type XXXX (e.g. GMV or AVS).
  ! After writing out the data, the structure can be released with ucd_release.
  
  TYPE t_ucdExport
  
    PRIVATE
    
    ! Output format. One of the UCD_FORMAT_XXXX constants.
    INTEGER :: coutputFormat = UCD_FORMAT_NONE
    
    ! Parameters for output routines. A combination of UCD_PARAM_XXX_YYYY,
    ! where XXX specifies the output format.
    INTEGER :: cparam = 0
    
    ! Name of the output file
    CHARACTER(LEN=SYS_STRLEN) :: sfilename = ""
    
    ! IO channel of the file
    INTEGER :: iunit = 0
    
    ! Export flags specifying the output
    INTEGER(I32) :: cflags            = 0
    
    ! Number of currently attached variables
    INTEGER :: nvariables        = 0
    
    ! Number of curently attached polygons
    INTEGER :: npolygons         = 0
    
    ! Number of currently attached tracer data fields
    INTEGER :: ntracerVariables  = 0
    
    ! Number of currently attached variable vectors
    INTEGER :: nvectors = 0

    ! The simulation time. =SYS_INFINITY if no simulation time is specified
    REAL(DP) :: dsimulationTime  = SYS_INFINITY
    
    ! Format of the simulation time. Fortran format string.
    CHARACTER(LEN=SYS_STRLEN) :: ssimTimeFormat = "(F20.5)"
    
    ! Format of the output of double-precision numbers. Fortran format string.
    CHARACTER(LEN=SYS_STRLEN) :: sdataFormat = "(E16.8E3)"
    
    ! An array containing the names of all the variables
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: p_SvariableNames => NULL()
    
    ! An array containing the names of all the variable vectors
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: p_SvarVecNames => NULL()
    
    ! Filename of file containing point coordinates.
    ! ""=no alternative source file.
    CHARACTER(LEN=SYS_STRLEN) :: saltFilePoints = ""

    ! Filename of file containing cell structure. 
    ! ""=no alternative source file.
    CHARACTER(LEN=SYS_STRLEN) :: saltFileCells = ""

    ! Filename of file containing material structure. 
    ! ""=no alternative source file.
    CHARACTER(LEN=SYS_STRLEN) :: saltFileMaterials = ""

    ! Filename of file containing polygons.
    ! ""=no alternative source file.
    CHARACTER(LEN=SYS_STRLEN) :: saltFilePolygons = ""
    
    ! A pointer to the underlying triangulation
    TYPE(t_triangulation), POINTER :: p_rtriangulation => NULL()
    
    ! A pointer to an array with specification flags. p_IvariableSpec(I)
    ! is a bitfield for variable I that specifies the type of the variable
    ! and how to handle it.
    INTEGER(I32), DIMENSION(:), POINTER :: p_IvariableSpec => NULL()
    
    ! A pointer to an array that specifies whether a variable is vertex based (1)
    ! or cell based (0).
    INTEGER, DIMENSION(:), POINTER :: p_IvariableBase => NULL()
    
    ! A pointer to an array that specifies the components of a vector variable.
    ! The first dimension of this vector is always 4.
    INTEGER, DIMENSION(:,:), POINTER :: p_Ivectors => NULL()
    
    ! A pointer to a list of handles of double precision pointers. 
    ! p_Hvariables(I) points to the data of variable I.
    INTEGER, DIMENSION(:), POINTER :: p_Hvariables => NULL()
    
    ! A pointer to a list of handles to polygon data.
    ! p_Hpolygon(I) points to a list of (X,Y) or (X,Y,Z) tags
    ! containing the points of a polygon of line segments.
    INTEGER, DIMENSION(:), POINTER :: p_Hpolygons => NULL()
    
    ! A pointer to a list of material identifier for polygonal data.
    ! Element i in this list specifies the material of the polygon.
    INTEGER :: hpolygonMaterial = ST_NOHANDLE
    
    ! A handle to a list of (X,Y) or (X,Y,Z) coordinates of tracers.
    INTEGER :: htracers = ST_NOHANDLE
    
    ! An array containing the names of all tracer variables
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: StracerVariable => NULL()
    
    ! A list of handles to tracer data. Each handle identifies an
    ! "array[1..#tracers] of double", which specifies data for each
    ! tracer. StracerVariable[i] os the name of the i'th array.
    INTEGER, DIMENSION(:), POINTER :: p_HtracerVariables => NULL()
    
    ! A handle to an array containing for every cell a cell material id.
    ! If not specified, every cell gets a default material id.
    INTEGER :: hIcellMaterial = ST_NOHANDLE

    ! A handle to an array containing for every node/vertex (corners, 
    ! midpoints,...) a material id.
    ! If not specified, every vertex/node gets a default material id.
    INTEGER :: hIvertexMaterial = ST_NOHANDLE
    
    ! Pointer to material names for vertex materials
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: SvertexMaterials => NULL()

    ! Pointer to material names for cell materials
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: ScellMaterials => NULL()
    
    ! Current length of the comment buffer
    INTEGER :: ncommentBufSize = 0
    
    ! Pointer to a buffer with comment lines.
    ! This is a character array. the lines are separated by the NEWLINE
    ! character.
    CHARACTER, DIMENSION(:), POINTER :: p_Scomments => NULL()
    
    ! Status variable: Number of vertices containing data in
    ! vertex based data arrays. All vertex-based variables in p_Hvariables
    ! have this length.
    INTEGER :: nvertices = 0
    
    ! Status variable: Number of cells containing data in
    ! cell based data arrays. All cell-based variables in p_Hvariables
    ! have this length.
    INTEGER :: ncells    = 0

    ! Status variable: Number of tracers.
    INTEGER :: ntracers  = 0
    
  END TYPE
  
!</typeblock>

!<typeblock>

  ! This structure is used for mesh refinement which is needed by some
  ! export routines. This structure is not meant to be used outside this
  ! module.

  TYPE t_ucdRefine
  
    ! Number of additional vertices
    INTEGER :: nvertices = 0
    
    ! Number of additional cells
    INTEGER :: ncells = 0
    
    ! Handle to additional vertice array
    INTEGER :: h_DvertexCoords = ST_NOHANDLE
    
    ! Handle to additional cell array
    INTEGER :: h_IverticesAtElement = ST_NOHANDLE

  END TYPE

!</typeblock>

!</types>

  PRIVATE :: ucd_moreVariables, ucd_moreVectors, ucd_morePolygons, &
      ucd_moreTracerVariables, ucd_refine

CONTAINS

!************************************************************************

!<subroutine>
  SUBROUTINE ucd_refine(rrefine, rtria, cflags)
  
!<description>
  ! This routine calculates the mesh refinement of a given triangulation
  ! and a combination of UCD_FLAG_XXXX flags.
  ! This routine is not meant to be called from outside this module.
!</description>

!<input>
  ! Specification of the underlying triangulation.
  TYPE(t_triangulation), INTENT(IN) :: rtria
  
  ! Bitfield that specifies the output.
  INTEGER(I32), INTENT(IN) :: cflags
!</input>

!<output>
  ! A mesh refinement
  TYPE(t_ucdRefine), INTENT(OUT) :: rrefine
!</output>

!</subroutine>

    ! Some local variables
    INTEGER :: i,j,k,ivt,off,numNewElem
    REAL(DP) :: dx
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords, p_DnewVerts
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertsAtEdge
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_InewVertsAtElement
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IvertsAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER, DIMENSION(2) :: I_dim
    LOGICAL :: bedgeMids, belemMids, brefined
    
    ! Initialize vertice / cell counts
    rrefine%nvertices = 0
    rrefine%ncells = 0
    
    ! Do we need to write edge / element midpoints?
    brefined  = (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)
    belemMids = (IAND(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. brefined
    bedgeMids = (IAND(cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
                (IAND(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. brefined
    
    ! If we use element midpoints or if the mesh is refines, we need to add
    ! the number of elements
    IF (belemMids) THEN
        
      rrefine%nvertices = rtria%NEL
      
    END IF
    
    ! Do we need to add edge midpoints?
    IF (bedgeMids) THEN
        
        rrefine%nvertices = rrefine%nvertices + rtria%NMT
    
    END IF
    
    ! If we do not allocate any vertices, we can leave this routine here
    IF (rrefine%nvertices .EQ. 0) RETURN
    
    ! Otherwise allocate the vertices
    IF (rtria%ndim .EQ. NDIM2D) THEN
      I_dim(1) = 2
    ELSE
      I_dim(1) = 3
    END IF
    I_dim(2) = rrefine%nvertices
    
    CALL storage_new2D("ucd_refine", "p_Dvertices", I_dim, ST_DOUBLE, &
        rrefine%h_DvertexCoords, ST_NEWBLOCK_ZERO)

    ! And get a pointer to them
    CALL storage_getbase_double2D(rrefine%h_DvertexCoords, p_DnewVerts)
    
    ! Get pointers to the triangulation's arrays
    CALL storage_getbase_double2D(rtria%h_DvertexCoords, p_DvertexCoords)
    CALL storage_getbase_int2D(rtria%h_IverticesAtEdge, p_IvertsAtEdge)
    CALL storage_getbase_int2D(rtria%h_IedgesAtElement, p_IedgesAtElement)
    CALL storage_getbase_int2D(rtria%h_IverticesAtElement, p_IvertsAtElement)
    
    ! Calculate the vertices
    off = 0
    IF (bedgeMids) THEN
      DO i = 1, rtria%NMT
        DO j = 1, UBOUND(p_DvertexCoords, 1)
          p_DnewVerts(j, i) = 0.5_DP * &
                             (p_DvertexCoords(j, p_IvertsAtEdge(1, i)) + &
                              p_DvertexCoords(j, p_IvertsAtEdge(2, i)))
        END DO
      END DO
      ! remember we already wrote some vertices
      off = rtria%NMT
    END IF
    
    ! Calculate element midpoints?
    IF (belemMids) THEN
    
      ! Go through all elements...
      DO i = 1, rtria%NEL
      
        ! ...and all coordinates...
        DO j = 1, UBOUND(p_DvertexCoords, 1)
        
          dx = 0.0_DP
          
          ! ...and through every vertex of the element
          DO k = 1, UBOUND(p_IvertsAtElement, 1)
          
            ivt = p_IvertsAtElement(k, i)
          
            ! Does this element only have k-1 vertices?
            IF (ivt .EQ. 0) EXIT
          
            dx = dx + p_DvertexCoords(j, ivt)
            
          END DO
          
          ! Store element midpoint
          p_DnewVerts(j, off+i) = dx / REAL(k-1, DP)
          
        END DO
        
      END DO
      
    END IF
    
    ! Do we also need to refine the elements?
    ! If not, then we can return here
    IF (.NOT. brefined) RETURN
    
    ! We first need to count the number of elements we will create for the
    ! refinement...
    ! Since we work in 2D here, every refinement takes 4 new elements...
    ! TODO: This code may need to be replaced for 3D grids...
    numNewElem = rtria%NEL * 4
    
    ! Allocate elements
    I_dim(1) = 4
    I_dim(2) = numNewElem
    CALL storage_new2D("ucd_refine", "p_DnewVertsAtElement", I_dim, &
        ST_INT, rrefine%h_IverticesAtElement, ST_NEWBLOCK_ZERO)

    ! And get a pointer to them
    CALL storage_getbase_int2D(rrefine%h_IverticesAtElement, p_InewVertsAtElement)

    ! Now go through all elements
    DO i = 1, rtria%NEL
    
      ! Count the number of vertices for this element
      DO k = 1, UBOUND(p_IvertsAtElement, 1)
        IF (p_IvertsAtElement(k, i) .EQ. 0) EXIT
      END DO
      
      ! Now this element has k-1 vertices
      k = k - 1
      
      SELECT CASE(k)
      CASE (3)
        ! Let a coarse grid triangle have 3 vertices (1, 2, 3) and 3 edges
        ! (a, b, c), then it will be refined into 4 new triangles (I, J, K, L)
        !
        !        1                        1
        !       / \                      / \
        !      /   \                    / J \
        !     a     c        ==>       a-----c
        !    /       \                / \ I / \
        !   /         \              / K \ / L \
        !  2-----b-----3            2-----b-----3
        !
        ! The triangles are indexed as follows:
        ! I => IEL
        ! J => NEL + 3*(IEL - 1) + 1
        ! K => NEL + 3*(IEL - 1) + 2
        ! L => NEL + 3*(IEL - 1) + 3
        ! Triangle I
        p_InewVertsAtElement(1, i) = rtria%NVT + p_IedgesAtElement(1, i)
        p_InewVertsAtElement(2, i) = rtria%NVT + p_IedgesAtElement(2, i)
        p_InewVertsAtElement(3, i) = rtria%NVT + p_IedgesAtElement(3, i)
        ! Triangle J
        j = rtria%NEL + 3*i - 2
        p_InewVertsAtElement(1, j) = p_IvertsAtElement(1, i)
        p_InewVertsAtElement(2, j) = rtria%NVT + p_IedgesAtElement(1, i)
        p_InewVertsAtElement(3, j) = rtria%NVT + p_IedgesAtElement(3, i)
        ! Triangle K
        j = j + 1
        p_InewVertsAtElement(1, j) = p_IvertsAtElement(2, i)
        p_InewVertsAtElement(2, j) = rtria%NVT + p_IedgesAtElement(2, i)
        p_InewVertsAtElement(3, j) = rtria%NVT + p_IedgesAtElement(1, i)
        ! Triangle L
        j = j + 1
        p_InewVertsAtElement(1, j) = p_IvertsAtElement(3, i)
        p_InewVertsAtElement(2, j) = rtria%NVT + p_IedgesAtElement(3, i)
        p_InewVertsAtElement(3, j) = rtria%NVT + p_IedgesAtElement(2, i)

      CASE (4)
        ! Let a coarse grid quadrilateral have 4 vertices (1, 2, 3, 4), 4 edges
        ! (a, b, c, d) and a midpoint (m), then it will be refined into 4 new
        ! quadrilaterals (I, J, K, L)
        !
        !  1---d---4           1---d---4
        !  |       |           | I | L |
        !  a   m   c    ==>    a---m---c
        !  |       |           | J | K |
        !  2---b---3           2---b---3
        !
        ! The quadrilaterals are indexed as follows:
        ! I => IEL
        ! J => NEL + 3*(IEL - 1) + 1
        ! K => NEL + 3*(IEL - 1) + 2
        ! L => NEL + 3*(IEL - 1) + 3
        
        ! Get the element's mid-point offset
        ivt = rtria%NVT + off + i
        ! Quad I
        p_InewVertsAtElement(1, i) = p_IvertsAtElement(1, i)
        p_InewVertsAtElement(2, i) = rtria%NVT + p_IedgesAtElement(1, i)
        p_InewVertsAtElement(3, i) = ivt
        p_InewVertsAtElement(4, i) = rtria%NVT + p_IedgesAtElement(4, i)
        ! Quad J
        j = rtria%NEL + 3*i - 2
        p_InewVertsAtElement(1, j) = p_IvertsAtElement(2, i)
        p_InewVertsAtElement(2, j) = rtria%NVT + p_IedgesAtElement(2, i)
        p_InewVertsAtElement(3, j) = ivt
        p_InewVertsAtElement(4, j) = rtria%NVT + p_IedgesAtElement(1, i)
        ! Quad K
        j = j+1
        p_InewVertsAtElement(1, j) = p_IvertsAtElement(3, i)
        p_InewVertsAtElement(2, j) = rtria%NVT + p_IedgesAtElement(3, i)
        p_InewVertsAtElement(3, j) = ivt
        p_InewVertsAtElement(4, j) = rtria%NVT + p_IedgesAtElement(2, i)
        ! Quad L
        j = j+1
        p_InewVertsAtElement(1, j) = p_IvertsAtElement(4, i)
        p_InewVertsAtElement(2, j) = rtria%NVT + p_IedgesAtElement(4, i)
        p_InewVertsAtElement(3, j) = ivt
        p_InewVertsAtElement(4, j) = rtria%NVT + p_IedgesAtElement(3, i)
        
      END SELECT
    
    END DO
    
    ! That's it
    
  END SUBROUTINE

!************************************************************************

!<subroutine>
  SUBROUTINE ucd_startGMV (rexport,cflags,rtriangulation,sfilename)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the GMV file
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  INTEGER(I32), INTENT(IN) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  TYPE(t_triangulation), INTENT(IN), TARGET :: rtriangulation
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  TYPE(t_ucdExport), INTENT(OUT) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(OUT) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_GMV
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    IF ((IAND(cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    END IF

    IF ((IAND(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      IF (rtriangulation%NDIM .EQ. NDIM2D) THEN
        rexport%ncells = rtriangulation%NEL*4
      ELSE
        rexport%ncells = rtriangulation%NEL*8
      END IF
    END IF
  
  END SUBROUTINE

!************************************************************************

!<subroutine>
  SUBROUTINE ucd_startBGMV (rexport,cflags,rtriangulation,sfilename)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the GMV file
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  INTEGER(I32), INTENT(IN) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  TYPE(t_triangulation), INTENT(IN), TARGET :: rtriangulation
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  TYPE(t_ucdExport), INTENT(OUT) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(OUT) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_BGMV
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    IF ((IAND(cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    END IF

    IF ((IAND(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      IF (rtriangulation%NDIM .EQ. NDIM2D) THEN
        rexport%ncells = rtriangulation%NEL*4
      ELSE
        rexport%ncells = rtriangulation%NEL*8
      END IF
    END IF
  
  END SUBROUTINE
  
!************************************************************************

!<subroutine>
  SUBROUTINE ucd_startAVS (rexport,cflags,rtriangulation,sfilename)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the AVS file
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  INTEGER(I32), INTENT(IN) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  TYPE(t_triangulation), INTENT(IN), TARGET :: rtriangulation
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  TYPE(t_ucdExport), INTENT(OUT) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(OUT) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_AVS
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    IF ((IAND(cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    END IF

    IF ((IAND(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      IF (rtriangulation%NDIM .EQ. NDIM2D) THEN
        rexport%ncells = rtriangulation%NEL*4
      ELSE
        rexport%ncells = rtriangulation%NEL*8
      END IF
    END IF
  
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE ucd_startVTK (rexport,cflags,rtriangulation,sfilename,cparam)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the VTK file
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  INTEGER(I32), INTENT(IN) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  TYPE(t_triangulation), INTENT(IN), TARGET :: rtriangulation
  
  ! OPTIONAL: Parameters for the VTK exporter
  INTEGER, OPTIONAL, INTENT(IN) :: cparam
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  TYPE(t_ucdExport), INTENT(OUT) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(OUT) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_VTK
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    
    ! Do we have any parameters?
    IF (PRESENT(cparam)) THEN
      rexport%cparam = cparam
    END IF
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    IF ((IAND(cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    END IF

    IF ((IAND(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
        (IAND(cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      IF (rtriangulation%NDIM .EQ. NDIM2D) THEN
        rexport%ncells = rtriangulation%NEL*4
      ELSE
        rexport%ncells = rtriangulation%NEL*8
      END IF
    END IF
  
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_release (rexport)

!<description>
  ! Releases all memory in the UCD export structure rexport.
!</description>

!<inputoutput>
  ! The ucd export structure that is to be released.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>
    INTEGER :: i
    
    ! Release all memory
    
    IF (ASSOCIATED(rexport%p_Scomments)) DEALLOCATE(rexport%p_Scomments)

    IF (ASSOCIATED(rexport%p_SvariableNames)) DEALLOCATE(rexport%p_SvariableNames)
    IF (ASSOCIATED(rexport%p_IvariableSpec)) DEALLOCATE(rexport%p_IvariableSpec)
    IF (ASSOCIATED(rexport%p_IvariableBase)) DEALLOCATE(rexport%p_IvariableBase)
    
    IF (ASSOCIATED(rexport%p_SvarVecNames)) DEALLOCATE(rexport%p_SvarVecNames)
    IF (ASSOCIATED(rexport%p_Ivectors)) DEALLOCATE(rexport%p_Ivectors)
    
    IF (ASSOCIATED(rexport%p_Hvariables)) THEN
      DO i=1,rexport%nvariables
        IF (rexport%p_Hvariables(i) .NE. ST_NOHANDLE) &
          CALL storage_free (rexport%p_Hvariables(i))
      END DO
      DEALLOCATE(rexport%p_Hvariables   )
    END IF
    
    IF (ASSOCIATED(rexport%p_Hpolygons    )) THEN
      DO i=1,rexport%npolygons
        IF (rexport%p_Hpolygons(i) .NE. ST_NOHANDLE) &
          CALL storage_free (rexport%p_Hpolygons(i))
      END DO
      DEALLOCATE(rexport%p_Hpolygons    )
    END IF
    
    IF (rexport%hpolygonMaterial .NE. ST_NOHANDLE) &
        CALL storage_free(rexport%hpolygonMaterial    )
    IF (rexport%hIcellMaterial    .NE. ST_NOHANDLE) &
        CALL storage_free(rexport%hIcellMaterial)
    IF (rexport%hIvertexMaterial    .NE. ST_NOHANDLE) &
        CALL storage_free(rexport%hIvertexMaterial)

    IF (ASSOCIATED(rexport%SvertexMaterials)) DEALLOCATE(rexport%SvertexMaterials)
    IF (ASSOCIATED(rexport%ScellMaterials)) DEALLOCATE(rexport%ScellMaterials)
    
    rexport%nvectors          = 0
    rexport%nvariables        = 0
    rexport%npolygons         = 0
    
    ! Release all tracer information
    CALL ucd_removeTracers (rexport)

    ! Don't deallocate the tringulation -- we are not the owner!!! :-)

  END SUBROUTINE

!************************************************************************

!<subroutine>
  SUBROUTINE ucd_setAlternativeSource (rexport,sfilename,caltFlags)
 
!<description>
  ! This command allows to specify an alternative source file for
  ! a part of the triangulation or the whole triangulation. 
  ! If no alternative source file is specified, the full mesh is 
  ! written to the output file. If an alternative source file is specified 
  ! for parts of a triangulation or for the full triangulation, a 
  ! reference to that file is written to the output file instead of
  ! the mesh. This allows to save disc space e.g. in a nonstationary
  ! simulation, where a mesh has only to be written into the first file
  ! while it usually stays fixed in all subsequent files.\\
  !
  ! sfilename is the filename of the source file containing alternative
  ! mesh information. caltFlags is a bitfield of type UCD_ASRC_xxxx
  ! flags that specify which parts of the triangulation are specified
  ! in that source file.
  ! The routine can be called multiple times with different filenames
  ! and a different flag-bitfield. This allows to specify different
  ! alternative source files for different parts of the triangulation.
  ! However, the most easiest way to use it is to specify 
  ! caltFlags=UCD_ASRC_ALL. In this case, the whole mesh is specified
  ! to be found in 'sfilename'.
  !
  ! Remark: Whether or not this technique is supported depends on 
  ! the output file format (GMV i.e. supports this). If not supported,
  ! the full mesh is written to the output file when ucd_write is called.
!</description>
 
!<input>
  ! Filename of the alternative input file
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
  
  ! Bitfield. Combination if UCD_ASRC_xxxx constants that specify which
  ! parts of the mesh are to be found in sfilename. UCD_ASRC_ALL specifies
  ! that the whole mesh is to be found in sfilename.
  INTEGER(I32), INTENT(IN) :: caltFlags
!</input>
  
!<inputoutput>
  ! The ucd export structure that specifies the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setAlternativeSource')
      CALL sys_halt()
    END IF
    
    ! Depending on the flags in caltFlags, put sfilename to the
    ! corresponding filename-strings in rexport that specify an alternative
    ! source for that part of the triangulation.
    !
    ! Point coordinates:
    IF (IAND(caltFlags,UCD_ASRC_POINTS) .NE. 0) THEN
      rexport%saltFilePoints = sfilename
    END IF

    ! Cell structure
    IF (IAND(caltFlags,UCD_ASRC_CELLS) .NE. 0) THEN
      rexport%saltFileCells = sfilename
    END IF

    ! Material structure
    IF (IAND(caltFlags,UCD_ASRC_MATERIALS) .NE. 0) THEN
      rexport%saltFileMaterials = sfilename
    END IF

    ! Polygons
    IF (IAND(caltFlags,UCD_ASRC_POLYGONS) .NE. 0) THEN
      rexport%saltFilePolygons = sfilename
    END IF
    
    ! The output routine must take care, that rexport%saltFileXXXX is
    ! correctly interpreted!
  
  END SUBROUTINE

!************************************************************************

!<subroutine>
  SUBROUTINE ucd_setMaterials (rexport,SmaterialsCells,SmaterialsVert)
  
!<description>
  ! This routine allows to specify names for material id's.
  ! SmaterialsCells is a list of material names for cells and
  ! SmaterialsVert a list of material names for vertices/nodes. 
  ! SmaterialsCell(i)/SmaterialsVert(i) is the 
  ! name that is to be assigned to material id i.
  !
  ! SmaterialVert is an optional parameter; if not specified, the
  ! vertex material names will coincide with the cell material names.
  !
  ! Whether or not named materials are supported, depends on the output
  ! format (GMV, AVS,...)
!</description>
 
!<input>
  ! Array with strings for the cell materials.
  ! The i'th string specifies a material id of material i.
  CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), INTENT(IN) :: SmaterialsCells
  
  ! OPTIONAL: Array with strings for the vertex/node materials.
  ! The i'th string specifies a material id of material i.
  ! If not specified, the same material names will be used for both,
  ! cells and vertices.
  CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), INTENT(IN), OPTIONAL :: SmaterialsVert
!</input>
  
!<inputoutput>
  ! The ucd export structure that specifies the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setMaterials')
      CALL sys_halt()
    END IF
    
    ! Make a copy of the strings. Use ALLOCATE/DEALLOCATE directly.
    IF (ASSOCIATED(rexport%ScellMaterials)) DEALLOCATE(rexport%ScellMaterials)
    ALLOCATE(rexport%ScellMaterials(SIZE(SmaterialsCells)))
    rexport%ScellMaterials = SmaterialsCells
    
    IF (ASSOCIATED(rexport%SvertexMaterials)) DEALLOCATE(rexport%SvertexMaterials)
    
    IF (PRESENT(SmaterialsVert)) THEN
      ALLOCATE(rexport%SvertexMaterials(SIZE(SmaterialsVert)))
      rexport%SvertexMaterials = SmaterialsVert
    ELSE
      ALLOCATE(rexport%SvertexMaterials(SIZE(SmaterialsCells)))
      rexport%SvertexMaterials = SmaterialsCells
    END IF
  
  END SUBROUTINE

!************************************************************************

!<subroutine>
  SUBROUTINE ucd_setCellMaterial (rexport,Imaterials)
  
!<description>
  ! This routine allows to specify for each cell a material id.
  ! Imaterials contains for each cell an integer with the corresponding
  ! id. How this id is visualised (if at all) depends on the 
  ! postprocessing tool (GMV, AVS,...).
!</description>
 
!<input>
  ! Array with as many elements as NEL in the triangulation. For every
  ! element i, Imaterials(i) specifies the element material id that
  ! should be assigned to that cell.
  INTEGER(I32), INTENT(IN), DIMENSION(:) :: Imaterials
!</input>
  
!<inputoutput>
  ! The ucd export structure that specifies the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_Idata
  INTEGER(PREC_ELEMENTIDX) :: NEL

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setCellMaterial')
      CALL sys_halt()
    END IF
    
    NEL = rexport%p_rtriangulation%NEL
    
    IF (SIZE(Imaterials) .LT. NEL) THEN
      CALL output_line ('Imaterials invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setCellMaterial')
      CALL sys_halt()
    END IF
    
    ! Copy that data and save it to the rexport structure.
    ! Create a new hImaterials handle if it does not exist.
    IF (rexport%hIcellMaterial .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('ucd_setCellMaterial','hIcellMaterial',&
          NEL,ST_INT,rexport%hIcellMaterial,ST_NEWBLOCK_NOINIT)
    END IF
    
    CALL storage_getbase_int (rexport%hIcellMaterial,p_Idata)
    CALL lalg_copyVectorInt (Imaterials(1:NEL),p_Idata(1:NEL))
  
  END SUBROUTINE
  
!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_setVertexMaterial (rexport,&
      ImaterialsVert, ImaterialsMid, ImaterialsElem)

!<description>
  ! This routine allows to specify for each vertex/node a material id.
  ! ImaterialsVert contains for each corner vertex an integer with the 
  ! corresponding id. How this id is visualised (if at all) depends on the 
  ! postprocessing tool (GMV, AVS,...).
  ! The optional ImaterialsMid and ImaterialsElem arrays allow to specify
  ! also for each edge midpoint and element midpoint a material id as well.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Array with as many elements as NVT in the triangulation. For every
  ! vertex i, ImaterialsVert(i) specifies the vertex material id that
  ! should be assigned to that vertex.
  INTEGER(I32), INTENT(IN), DIMENSION(:) :: ImaterialsVert

  ! OPTIONAL: Array with as many elements as NMT in the triangulation. 
  ! For every edge i, ImaterialsMid(i) specifies the material id
  ! that should be assigned to the corresponding edge midpoint.
  INTEGER(I32), INTENT(IN), DIMENSION(:), OPTIONAL :: ImaterialsMid

  ! OPTIONAL: Array with as many elements as NEL in the triangulation. 
  ! For every element i, ImaterialsElem(i) specifies the material id
  ! that should be assigned to the corresponding element midpoint.
  ! Note: The material of the element midpoint need not to coincide
  !  with the material of the element -- which is specified
  !  in ucd_setCellMaterial!
  INTEGER(I32), INTENT(IN), DIMENSION(:), OPTIONAL :: ImaterialsElem

!</input>

!</subroutine>

  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_Idata
  INTEGER(PREC_ELEMENTIDX) :: NEL
  INTEGER(PREC_EDGEIDX) :: NMT
  INTEGER(PREC_VERTEXIDX) :: NVT

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
      CALL sys_halt()
    END IF
    
    NVT = rexport%p_rtriangulation%NVT
    NMT = rexport%p_rtriangulation%NMT
    NEL = rexport%p_rtriangulation%NEL
    
    IF (SIZE(ImaterialsVert) .LT. NVT) THEN
      CALL output_line ('ImaterialsVert invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
      CALL sys_halt()
    END IF

    IF (PRESENT(ImaterialsMid)) THEN
      IF (SIZE(ImaterialsMid) .LT. NMT) THEN
        CALL output_line ('ImaterialsMid invalid!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
        CALL sys_halt()
      END IF
    END IF

    IF (PRESENT(ImaterialsElem)) THEN
      IF (SIZE(ImaterialsElem) .LT. NEL) THEN
        CALL output_line ('ImaterialsElem invalid!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
        CALL sys_halt()
      END IF
    END IF
    
    ! Create a new data array for the vertex materials if necessary.
    ! Fill it with 0, which is the default material id.
    
    ! Create a new hImaterials handle if it does not exist.
    IF (rexport%hIvertexMaterial .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('ucd_setVertexMaterial','hIvertexMaterial',&
          INT(rexport%nvertices,I32),ST_INT,rexport%hIvertexMaterial,ST_NEWBLOCK_ZERO)
    END IF

    CALL storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
    
    ! Copy that data and save it to the rexport structure.
    CALL lalg_copyVectorInt (ImaterialsVert(1:NVT),p_Idata(1:NVT))

    ! Copy edge midpoint data if available
    IF ((IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
        (IAND(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
        (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      IF (PRESENT(ImaterialsMid)) THEN
        CALL lalg_copyVectorInt( &
            ImaterialsMid(1:NMT),p_Idata(NVT+1:NVT+NMT))
      END IF
    END IF
    
    ! Copy element midpoint data if available
    IF ((IAND(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
        (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      IF (PRESENT(ImaterialsElem)) THEN
        CALL lalg_copyVectorInt( &
            ImaterialsElem(1:NEL),p_Idata(NVT+NMT+1:NVT+NMT+NEL))
      END IF
    END IF    

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_write (rexport)

!<description>
  ! Writes the output of UCD data into a postprocessig file.
  ! All pending data in rexport is written into the file. 
  ! If the file is identified by a filename, a new file is opened,
  ! data is written to it and the file is closed at the end.
!</description>

!<inputoutput>
  ! The ucd export structure that specifies the output.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    ! Check for special binary file output
    SELECT CASE(rexport%coutputFormat)
    CASE (UCD_FORMAT_BGMV)
      CALL ucd_writeBGMV (rexport)
      RETURN
    END SELECT

    ! If there is a new filename, open the output file.
    IF (rexport%sfilename .NE. '') THEN
      CALL io_openFileForWriting(rexport%sfilename, rexport%iunit, SYS_REPLACE)
      IF (rexport%iunit .LT. 0) THEN
        CALL output_line ('Cannot open file "'//TRIM(rexport%sfilename)&
                //'" for writing!', OU_CLASS_ERROR,OU_MODE_STD,'ucd_write')
        CALL sys_halt()
      END IF
    END IF
    
    IF (rexport%iunit .EQ. 0) THEN
      CALL output_line ('Cannot write UCD output: No output channel/filename!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_write')
      CALL sys_halt()
    END IF

    ! Start the output
    SELECT CASE(rexport%coutputFormat)
    CASE (UCD_FORMAT_GMV)
      CALL ucd_writeGMV (rexport)
    CASE (UCD_FORMAT_AVS)
      CALL ucd_writeAVS (rexport)
    CASE (UCD_FORMAT_VTK)
      CALL ucd_writeVTK (rexport)
    END SELECT
    
    IF (rexport%sfilename .NE. '') THEN
      ! Close the file if it was opened previously.
      CLOSE (rexport%iunit)
      rexport%iunit = 0
    END IF
    
  CONTAINS
    
    !****************************************************************
    
    SUBROUTINE ucd_writeGMV (rexport)

    ! Specific output routine for GMV output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startGMV
    ! subroutine before.
    
    ! The export structure with all information
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
    ! local variables
    INTEGER :: mfile,i,j,k,icoor
    INTEGER(PREC_VERTEXIDX) :: ivt,ivt1,ivt2,nnodes
    INTEGER(PREC_EDGEIDX) :: imt
    INTEGER(PREC_ELEMENTIDX) :: iel
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords,p_Ddata2D
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge 
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    REAL(DP) :: dx
    
      mfile = rexport%iunit

      CALL storage_getbase_double2d (rexport%p_Rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
      CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtElement,&
          p_IverticesAtElement)
      
      !----------------------------------------------------
      ! Write the GMV header
      WRITE(mfile,'(A)') 'gmvinput ascii'

      !----------------------------------------------------
      ! Simulation time
      IF (rexport%dsimulationTime .NE. SYS_INFINITY) THEN
        WRITE(mfile,'(A)',ADVANCE='NO')     'probtime '
        WRITE(mfile,rexport%ssimTimeFormat) rexport%dsimulationTime
      END IF

      !----------------------------------------------------
      ! Write all comments
      IF (rexport%ncommentBufSize .GT. 0) THEN
        WRITE(mfile,'(A)') 'comments'
        
        i = 1
        ! Find the end of the first line
        DO j=i,rexport%ncommentBufSize
          IF (rexport%p_Scomments(j) .EQ. NEWLINE) EXIT
        END DO
        
        ! Write out all lines, one after the other
        DO WHILE (j .LE. rexport%ncommentBufSize)
          ! Write the line (character i..j-1), continue with the next
          DO k=i,j-1
            WRITE(mfile,'(A)',ADVANCE='NO') rexport%p_Scomments(k)
          END DO
          WRITE(mfile,'(A)')
          
          ! Continue after the NEWLINE character, find the next NEWLINE.
          i = j+1
          DO j=i,rexport%ncommentBufSize
            IF (rexport%p_Scomments(j) .EQ. NEWLINE) EXIT
          END DO
        
        END DO
        
        WRITE(mfile,'(A)') 'endcomm'
      END IF
      
      !----------------------------------------------------
      ! Write the triangulation.
      !
      ! GMV output allows to specify an alternative source file for
      ! triangulation related data. Figure out if some of the data
      ! has to be taken from such an alternative file. If yes, we
      ! write the filename into the GMV file. If no, we write the
      ! triangulation data to the GMV.
      !
      ! Point coordinates:
      
      IF (rexport%saltFilePoints .NE. "") THEN
      
        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disc space!
        
        WRITE(mfile,'(A)')'nodes fromfile "'//TRIM(rexport%saltFilePoints)//'"'
      
      ELSE

        WRITE(mfile,'(A,I10)') 'nodes ',rexport%nvertices
        
        ! Loop through the X/Y/Z coordinates
        
        DO icoor = 1,MIN(UBOUND(p_DvertexCoords,1),3)
        
          DO ivt=1,rexport%p_Rtriangulation%NVT
            WRITE(mfile,rexport%sdataFormat) p_DvertexCoords(icoor,ivt)
          END DO

          ! Write coordinates of edge midpoints?
          IF ((IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
              (IAND(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
              (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
              
            CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                p_IverticesAtEdge)

            ! We construct them by hand.
            ! In a later implementation, one could take the coordinates
            ! of h_DfreecornerCoordinates...
            DO imt=1,rexport%p_Rtriangulation%NMT
              ivt1 = p_IverticesAtEdge(1,imt)
              ivt2 = p_IverticesAtEdge(2,imt)
              dx = 0.5_DP*(p_DvertexCoords(icoor,ivt1) + &
                           p_DvertexCoords(icoor,ivt2))
              WRITE(mfile,rexport%sdataFormat) dx
            END DO
              
          END IF

          ! Write coordinates of element midpoints?
          IF ((IAND(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
              (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
              
            CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                p_IverticesAtEdge)

            ! We construct them by hand.
            ! In a later implementation, one could take the coordinates
            ! of h_DfreecornerCoordinates...
            DO iel=1,rexport%p_Rtriangulation%NEL
              
              dx = 0.0_DP
              
              DO i=1,UBOUND(p_IverticesAtElement,1)
                ivt = p_IverticesAtElement(i,iel)
                IF (ivt .NE. 0) THEN
                  dx = dx + p_DvertexCoords(icoor,ivt)
                ELSE
                  ! We have only (i-1) vertices in that element; 
                  ! happens e.g. in triangles that are mixed into a quad mesh.
                  ! We stop here.
                  EXIT
                END IF
              END DO
              
              ! If all vertices of the element are touched, there is i=NVE+1.
              ! Divide by the number of vertices to get the coordinate of the 
              ! midpoint of the element.
              dx = dx / REAL(i-1,DP)
              
              WRITE(mfile,rexport%sdataFormat) dx
            END DO
              
          END IF
        
        END DO ! icoor
        
        ! If there are not enough coordinates, we must add 0's -- as
        ! GMV always expects 3D data.

        DO icoor = UBOUND(p_DvertexCoords,1)+1 , 3
        
          DO ivt=1,rexport%p_Rtriangulation%NVT
            WRITE(mfile,rexport%sdataFormat) 0.0_DP
          END DO

          ! Write coordinates of edge midpoints?
          IF ((IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
              (IAND(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
              (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
              
            DO imt=1,rexport%p_Rtriangulation%NMT
              WRITE(mfile,rexport%sdataFormat) 0.0_DP
            END DO
              
          END IF

          ! Write coordinates of element midpoints?
          IF ((IAND(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
              (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
              
            DO iel=1,rexport%p_Rtriangulation%NEL
              WRITE(mfile,rexport%sdataFormat) 0.0_DP
            END DO
              
          END IF
        
        END DO ! icoor
        
      END IF
      
      ! Mesh connectivity / Cells:
      
      IF (rexport%saltFileCells .NE. "") THEN
      
        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disc space!
        
        WRITE(mfile,'(A)')'cells fromfile "'//TRIM(rexport%saltFileCells)//'"'
      
      ELSE
      
        ! Write the connectivity to the mesh - i.e. the cells.

        WRITE(MFILE,'(A,I10)') 'cells ',rexport%ncells
        
        IF (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .EQ. 0) THEN
        
          SELECT CASE (rexport%p_rtriangulation%ndim)
          
          CASE (NDIM1D)
        
            ! Standard mesh.
            DO iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              DO i=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
              END DO
              
              ! We have i-1 vertices on that element -- so what is it?
              SELECT CASE (i-1)
              CASE (2)
                ! Line in 1D
                WRITE(mfile,'(A)') 'line 2'
                WRITE(mfile,'(3I8)') p_IverticesAtElement(1:2,iel)
              
              CASE DEFAULT
                CALL output_line ('Invalid element!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              END SELECT
                
            END DO
            
          CASE (NDIM2D)

            ! Standard mesh.
            DO iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              DO i=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
              END DO
              
              ! We have i-1 vertices on that element -- so what is it?
              SELECT CASE (i-1)
              
              CASE (3)
                ! Triangle
                WRITE(mfile,'(A)') 'tri 3'
                WRITE(mfile,'(3I8)') p_IverticesAtElement(1:3,iel)
                
              CASE (4)
                ! Quad
                WRITE(mfile,'(A)')'quad 4'
                WRITE(mfile,'(4I8)') p_IverticesAtElement(1:4,iel)
                
              CASE DEFAULT
                CALL output_line ('Invalid element!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              END SELECT
                
            END DO

          CASE (NDIM3D)
          
            ! Standard mesh.
            DO iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              DO i=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
              END DO
              
              ! We have i-1 vertices on that element -- so what is it?
              SELECT CASE (i-1)
              
              CASE (4)
                ! Tetrahedron
                WRITE(mfile,'(A)') 'ptet4 4'
                WRITE(mfile,'(4I8)') p_IverticesAtElement(1:4,iel)
                
              CASE (5)
                ! Pyramid
                WRITE(mfile,'(A)') 'ppyrmd5 5'
                WRITE(mfile,'(4I8)') p_IverticesAtElement(1:5,iel)

              CASE (6)
                ! Prism
                WRITE(mfile,'(A)') 'pprism6 6'
                WRITE(mfile,'(4I8)') p_IverticesAtElement(1:6,iel)

              CASE (8)
                ! Hexahedron
                WRITE(mfile,'(A)')'phex8 8'
                WRITE(mfile,'(8I8)') p_IverticesAtElement(1:8,iel)
                
              CASE DEFAULT
                CALL output_line ('Invalid element!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              END SELECT
                
            END DO

          END SELECT          
            
        ELSE

          ! 1x refined mesh
          
          ! The edge numbers give the numbers of the edge midpoints and thus
          ! the numbers of the new vertices on the once refined mesh.
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IedgesAtElement,&
              p_IedgesAtElement)
              
          SELECT CASE (rexport%p_rtriangulation%ndim)
          
          CASE (NDIM1D)
              
            DO iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              DO i=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
              END DO
              
              ! We have i-1 vertices on that element -- so what is it?
              SELECT CASE (i-1)
              CASE (2)
                ! Line in 1D.
                !
                ! The coarse grid element is
                !
                !   1 -----IEL----- 2
                !
                ! The once refined element is
                !
                !   1 -- IEL -- 1* -- NEL+IEL -- 2
                !
                ! Write the connectivity of element IEL
                WRITE(mfile,'(A)') 'line 2'
                WRITE(mfile,'(3I8)') p_IverticesAtElement(1,iel),p_IedgesAtElement(1,iel)

              END SELECT
              
            END DO

            DO iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              DO i=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
              END DO
              
              ! We have i-1 vertices on that element -- so what is it?
              SELECT CASE (i-1)
              CASE (2)
                ! Line in 1D.
                !
                ! Element "NEL+1"
                WRITE(mfile,'(A)') 'line 2'
                WRITE(mfile,'(2I8)') p_IedgesAtElement(1,iel),p_IverticesAtElement(2,iel)
                
              END SELECT
              
            END DO

          CASE (NDIM2D)
              
            DO iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              DO i=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
              END DO
              
              ! We have i-1 vertices on that element -- so what is it?
              SELECT CASE (i-1)

              CASE (3)
                ! Triangle.
                !
                ! Let a coarse grid triangle be locally numbered as:
                !
                !   2 
                !   |  \
                !   |    \
                !   | IEL  \
                !   |        \
                !   3----------1
                ! 
                ! Then the refinement process assigns the following numbers:
                !
                !   2_
                !   |  \_
                !   | NEL \_
                !   | +1     \_
                !   2*-------- 1* 
                !   | \_  IEL=1|  \_
                !   |NEL \_    | NEL \_
                !   |+3     \_ | +2     \_
                !   3----------3*---------1
                !
                ! So write the edge(-midpoint) numbers as corner numbers
                ! of the "main" element
                    
                WRITE(mfile,'(A)') 'tri 3'
                WRITE(mfile,'(3I8)') p_IedgesAtElement(1:3,iel)
                
              CASE (4)
                ! Quad
                !
                ! Let a coarse grid quad be locally numbered as
                !
                !  4-------3
                !  |       |
                !  |  IEL  |
                !  |       |
                !  1-------2
                !
                ! Then the refinement process assigns the following numbers:
                !
                !  4-----7-----3
                !  |NEL+3|NEL+2|
                !  |     |     |
                !  8-----9-----6
                !  |IEL=1|NEL+1|
                !  |     |     |
                !  1-----5-----2
                !
                ! So construct the corners of the 'smaller' elements from the
                ! corners of the coarse grid element, the numbers of the
                ! edge(-midpoint)s of the element and the number of the midpoint
                ! of the element -- which is defined as the element number itself.
                
                WRITE(mfile,'(A)')'quad 4'
                WRITE(mfile,'(4I8)') &
                    p_IverticesAtElement(1,iel),p_IedgesAtElement(1,iel),&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(4,iel)
                
              END SELECT
              
            END DO

            DO iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              DO i=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
              END DO
              
              ! We have i-1 vertices on that element -- so what is it?
              SELECT CASE (i-1)
                
              CASE (3)
                ! Triangle.
                    
                ! Element "NEL+1"
                WRITE(mfile,'(A)') 'tri 3'
                WRITE(mfile,'(3I8)') p_IverticesAtElement(2,iel), &
                    p_IedgesAtElement(2,iel),p_IedgesAtElement(1,iel)

                ! Element "NEL+2"
                WRITE(mfile,'(A)') 'tri 3'
                WRITE(mfile,'(3I8)') p_IverticesAtElement(1,iel), &
                    p_IedgesAtElement(1,iel),p_IedgesAtElement(3,iel)

                ! Element "NEL+3"
                WRITE(mfile,'(A)') 'tri 3'
                WRITE(mfile,'(3I8)') p_IverticesAtElement(3,iel), &
                    p_IedgesAtElement(3,iel),p_IedgesAtElement(2,iel)
                
              CASE (4)
                ! Quad
                !
                ! Element "NEL+1"
                WRITE(mfile,'(A)')'quad 4'
                WRITE(mfile,'(4I8)') &
                    p_IverticesAtElement(2,iel),p_IedgesAtElement(2,iel),&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(1,iel)

                ! Element "NEL+2"
                WRITE(mfile,'(A)')'quad 4'
                WRITE(mfile,'(4I8)') &
                    p_IverticesAtElement(3,iel),p_IedgesAtElement(3,iel),&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(2,iel)

                ! Element "NEL+3"
                WRITE(mfile,'(A)')'quad 4'
                WRITE(mfile,'(4I8)') &
                    p_IverticesAtElement(4,iel),p_IedgesAtElement(4,iel),&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(3,iel)
                
              END SELECT
              
            END DO
            
          CASE (NDIM3D)

            CALL output_line ('GMV export for 1x refined mesh in 3D'//&
                ' not implemented!', OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
            CALL sys_halt()
          
          END SELECT
          
        END IF
        
      END IF

      !----------------------------------------------------
      ! Write material names -- if specified
      !

      IF (rexport%saltFileMaterials .NE. "") THEN

        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disc space!
        
        WRITE(mfile,'(A)')'material fromfile "'//TRIM(rexport%saltFileMaterials)//'"'
        
      ELSE

        ! Cell material
        IF (ASSOCIATED(rexport%ScellMaterials) .and. &
            (rexport%hIcellMaterial .NE. ST_NOHANDLE)) THEN
          ! GMV only supports <= 1000 materials!
          WRITE(mfile,'(A,I10,I10)') 'material ',MIN(1000,SIZE(rexport%ScellMaterials)),0
          DO i=1,MIN(1000,SIZE(rexport%ScellMaterials))
            ! GMV supports only <= 8 characters and does not allow spaces
            ! in the material name. We replace all invalid spaces by "_".
            WRITE(mfile,'(A8)') &
                sys_charreplace(TRIM(rexport%ScellMaterials(i)),' ','_')
          END DO
          
          IF (rexport%hIcellMaterial .NE. ST_NOHANDLE) THEN
            ! Write a list of material id's. For every cell, we specify its
            ! material by the material number.
            CALL storage_getbase_int (rexport%hIcellMaterial,p_Idata)
            DO i=1,SIZE(p_Idata)
              WRITE(mfile,'(I4)') p_Idata(i)
            END DO
          END IF
        END IF
        
        ! Vertex materials; coincide with cell materials if not specified.
        IF (ASSOCIATED(rexport%SvertexMaterials) .and. &
            (rexport%hIvertexMaterial .NE. ST_NOHANDLE)) THEN
          ! GMV only supports <= 1000 materials!
          WRITE(mfile,'(A,I10,I10)') 'material ',MIN(1000,SIZE(rexport%SvertexMaterials)),1
          DO i=1,MIN(1000,SIZE(rexport%SvertexMaterials))
            ! GMV supports only <= 8 characters and does not allow spaces
            ! in the material name. We replace all invalid spaces by "_".
            WRITE(mfile,'(A8)') &
                sys_charreplace(TRIM(rexport%SvertexMaterials(i)),' ','_')
          END DO
          
          IF (rexport%hIvertexMaterial .NE. ST_NOHANDLE) THEN
            ! Write a list of material id's. For every vertex, we specify its
            ! material by the material number.
            CALL storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
            DO i=1,SIZE(p_Idata)
              WRITE(mfile,'(I4)') p_Idata(i)
            END DO
          END IF
        ELSE
          IF (ASSOCIATED(rexport%ScellMaterials) .and. &
              (rexport%hIvertexMaterial .NE. ST_NOHANDLE)) THEN
            ! GMV only supports <= 1000 materials!
            WRITE(mfile,'(A,I10,I10)') 'material ',MIN(1000,SIZE(rexport%ScellMaterials)),1
            DO i=1,MIN(1000,SIZE(rexport%ScellMaterials))
              ! GMV supports only <= 8 characters and does not allow spaces
              ! in the material name. We replace all invalid spaces by "_".
              WRITE(mfile,'(A8)') &
                  sys_charreplace(TRIM(rexport%ScellMaterials(i)),' ','_')
            END DO
          END IF
          
          IF (rexport%hIvertexMaterial .NE. ST_NOHANDLE) THEN
            ! Write a list of material id's. For every vertex, we specify its
            ! material by the material number.
            CALL storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
            DO i=1,SIZE(p_Idata)
              WRITE(mfile,'(I4)') p_Idata(i)
            END DO
          END IF
        END IF
      END IF
      
      IF (rexport%nvariables .GT. 0) THEN
        !----------------------------------------------------
        ! Write a velocity field -- if there is one
        !
        ! To write the velocity field, we have to look it up in the
        ! set of variables.
        ! Search for the X-velocity and write it out.
        ! If there is no X-velocity, there is no velocity at all.
        ! There may be a velocity field given in vertices or in cells,
        ! so we have to search twice!
        !
        ! Look for cell based velocity.
        DO i=1,rexport%nvariables
          IF ((IAND(rexport%p_IvariableSpec(i),UCD_VAR_XVELOCITY) .NE. 0) .AND. &
              (rexport%p_IvariableBase(i) .EQ. UCD_BASE_ELEMENT)) THEN
            
            ! Found it. Write it out.
            WRITE (mfile,'(A)') 'velocity 0'
            
            CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
            
            ! Don't be confused! ivt=number of cell, as we are in the 
            ! 'cell-oriented' case here!!!
            DO ivt=1,rexport%ncells
              WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
            END DO
            
            ! Find the Y-velocity
            DO j=1,rexport%nvariables
              IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_YVELOCITY) .NE. 0) .AND. &
                  (rexport%p_IvariableBase(j) .EQ. UCD_BASE_ELEMENT)) THEN
                
                ! Found it. Write it out.
                CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
                DO ivt=1,rexport%ncells
                  WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
                END DO
                
                EXIT
              END IF
            END DO
              IF (j .GT. rexport%nvariables) THEN
                ! Not found. Write out 0's instead.
                DO ivt=1,rexport%ncells
                  WRITE (MFILE,rexport%sdataFormat) 0.0_DP
                END DO
              END IF
              
              ! Find the Z-velocity
              DO j=1,rexport%nvariables
                IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_ZVELOCITY) .NE. 0) .AND. &
                    (rexport%p_IvariableBase(j) .EQ. UCD_BASE_ELEMENT)) THEN
                  
                  ! Found it. Write it out.
                  CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
                  DO ivt=1,rexport%ncells
                    WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
                  END DO
                  
                  EXIT
                END IF
              END DO
              IF (j .GT. rexport%nvariables) THEN
                ! Not found. Write out 0's instead.
                DO ivt=1,rexport%ncells
                  WRITE (MFILE,rexport%sdataFormat) 0.0_DP
                END DO
              END IF
              
            END IF
          END DO
          
          ! Look for vertex based velocity.
          DO i=1,rexport%nvariables
            IF ((IAND(rexport%p_IvariableSpec(i),UCD_VAR_XVELOCITY) .NE. 0) .AND. &
                (rexport%p_IvariableBase(i) .EQ. UCD_BASE_VERTEX)) THEN
              
              ! Found it. Write it out.
              WRITE (mfile,'(A)') 'velocity 1'
              
              CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
              DO ivt=1,rexport%nvertices
                WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
              END DO
              
              ! Find the Y-velocity
            DO j=1,rexport%nvariables
              IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_YVELOCITY) .NE. 0) .AND. &
                  (rexport%p_IvariableBase(j) .EQ. UCD_BASE_VERTEX)) THEN
                  
                ! Found it. Write it out.
                CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
                DO ivt=1,rexport%nvertices
                  WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
                END DO
                
                EXIT
              END IF
            END DO
            IF (j .GT. rexport%nvariables) THEN
              ! Not found. Write out 0's instead.
              DO ivt=1,rexport%nvertices
                WRITE (MFILE,rexport%sdataFormat) 0.0_DP
              END DO
            END IF
            
            ! Find the Z-velocity
            DO j=1,rexport%nvariables
              IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_ZVELOCITY) .NE. 0) .AND. &
                  (rexport%p_IvariableBase(j) .EQ. UCD_BASE_VERTEX)) THEN
                  
                ! Found it. Write it out.
                CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
                DO ivt=1,rexport%nvertices
                  WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
                END DO
                
                EXIT
              END IF
            END DO
            IF (j .GT. rexport%nvariables) THEN
              ! Not found. Write out 0's instead.
              DO ivt=1,rexport%nvertices
                WRITE (MFILE,rexport%sdataFormat) 0.0_DP
              END DO
            END IF
            
          END IF
        END DO
        
        !----------------------------------------------------
        ! Write all variables which are not velocities
        DO i=1,rexport%nvariables
          IF (IAND(rexport%p_IvariableSpec(i), &
              UCD_VAR_XVELOCITY+UCD_VAR_YVELOCITY+UCD_VAR_ZVELOCITY) .EQ. 0) THEN
          
            WRITE (MFILE,'(A)') 'variable'  
            
            IF (rexport%p_IvariableBase(i) .EQ. UCD_BASE_ELEMENT) THEN
              ! Cell based variable
              WRITE (MFILE,'(A,I5)') TRIM(rexport%p_SvariableNames(i)),0
              CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
              ivt1 = rexport%ncells
            ELSE
              ! Vertex based variable
              WRITE (MFILE,'(A,I5)') TRIM(rexport%p_SvariableNames(i)),1
              CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
              ivt1 = rexport%nvertices
            END IF

            DO ivt=1,ivt1
              WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
            END DO
            
            WRITE (MFILE,'(A)') 'endvars'      
          END IF
          
        END DO ! i

      END IF

      !----------------------------------------------------
      ! Write polygon data

      IF (rexport%saltFilePolygons .NE. "") THEN
      
        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disc space!
        
        WRITE(mfile,'(A)')'polygons fromfile "'//TRIM(rexport%saltFilePolygons)//'"'
      
      ELSE
        
        IF (ASSOCIATED(rexport%p_Hpolygons)) THEN
          
          ! At least one polygon
          WRITE (MFILE,'(A)') 'polygons'
          
          ! Materials
          CALL storage_getbase_int (rexport%hpolygonMaterial,p_Idata)
          
          ! Write all polygons.
          DO i=1,rexport%npolygons
            
            ! Coordinates        
            CALL storage_getbase_double2D (rexport%p_Hpolygons(i),p_Ddata2D)
            
            ! Write material, #points
            WRITE (MFILE,'(2I10)') p_Idata(i),UBOUND(p_Ddata2D,2)
            
            ! Either we have 2D or 3D coordinates. 
            ! Write coordinates of the points forming the line segments 
            ! of the polygon
            ! First all X-, then all Y- and at the end all Z-coordinates -- or 0.0.
            DO k=1,NDIM3D
              IF (UBOUND(p_Ddata2D,1) .GE. k) THEN
                DO j=1,UBOUND(p_Ddata2D,2)
                  WRITE (MFILE,'(E15.7)') p_Ddata2D(k,j)
                END DO
              ELSE
                DO j=1,UBOUND(p_Ddata2D,2)
                  WRITE (MFILE,'(E15.7)') 0.0_DP
                END DO
              END IF
            END DO
            
          END DO
          
          WRITE (mfile,'(A)') 'endpoly'
          
        END IF

      END IF

      !----------------------------------------------------
      ! Write tracer coordinates and data
      IF (rexport%ntracers .NE. 0) THEN
        WRITE(mfile,'(A,I10)') 'tracers ',rexport%ntracers
        
        CALL storage_getbase_double2d(rexport%htracers,p_Ddata2D)
        
        ! First write all X-coordinates, then Y-coordinates, 
        ! then Z-coordinates -- if specified.
        DO i=1,3
          IF (i .LE. UBOUND(p_Ddata,1)) THEN
            DO j=1,rexport%ntracers
              WRITE (mfile,rexport%sdataFormat) p_Ddata2D(j,i)
            END DO
          ELSE
            DO j=1,rexport%ntracers
              WRITE (mfile,rexport%sdataFormat) 0.0_DP
            END DO
          END IF
        END DO
        
        ! Write tracer variables if specified
        DO i=1,rexport%ntracerVariables
          WRITE (mfile,'(A32)') rexport%StracerVariable(i)
          
          CALL storage_getbase_double (rexport%p_HtracerVariables(i), p_Ddata)
          DO j=1,rexport%ntracers
            WRITE (mfile,rexport%sdataFormat) p_Ddata(j)
          END DO
        END DO
        
        WRITE (mfile,'(A)') 'endtrace'
        
      END IF

      !----------------------------------------------------
      ! Finally write the GMV footer, finish
      
      WRITE(mfile,'(A)') 'endgmv'
      
    END SUBROUTINE

    !****************************************************************
    
    SUBROUTINE ucd_writeBGMV (rexport)

    ! Specific output routine for GMV output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startBGMV
    ! subroutine before.
    
    ! The export structure with all information
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
    ! local variables
    INTEGER :: i,j,k,icoor
    INTEGER(PREC_VERTEXIDX) :: ivt,ivt1,ivt2,nnodes,nvt,nverts,matnum
    INTEGER(PREC_EDGEIDX) :: imt
    INTEGER(PREC_ELEMENTIDX) :: iel
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords,p_Ddata2D
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge 
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement

    INTEGER :: nod2ids(2), nod3ids(3), nod4ids(4), nod5ids(5), nod6ids(6), nod8ids(8)
    REAL, DIMENSION(:), ALLOCATABLE :: X,Y,Z,VAR
    REAL dx, dy,dz

    ! Open file for binary output
    CALL fgmvwrite_openfile(rexport%sfilename)
    
    CALL storage_getbase_double2d (rexport%p_Rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    
    !----------------------------------------------------
    ! Simulation time
    IF (rexport%dsimulationTime .NE. SYS_INFINITY) THEN
      CALL fgmvwrite_probtime(rexport%dsimulationTime)
    END IF

    !----------------------------------------------------
    ! Write the triangulation.
    !
    ! GMV output allows to specify an alternative source file for
    ! triangulation related data. Figure out if some of the data
    ! has to be taken from such an alternative file. If yes, we
    ! write the filename into the GMV file. If no, we write the
    ! triangulation data to the GMV.
    !
    ! Point coordinates:
    
    IF (rexport%saltFilePoints .NE. "") THEN
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!
      
      CALL fgmvwrite_nodes_fromfile(TRIM(rexport%saltFilePoints), INT(rexport%nvertices))
      
    ELSE

      ! Allocate temporal memory
      ALLOCATE(X(rexport%nvertices), Y(rexport%nvertices), Z(rexport%nvertices))
      
      SELECT CASE(rexport%p_Rtriangulation%ndim)
      
      CASE(NDIM1D)
        DO ivt=1,rexport%p_Rtriangulation%NVT
          X(ivt) =  REAL(p_DvertexCoords(1,ivt))
          Y(ivt) = 0.0E0
          Z(ivt) = 0.0E0
        END DO

        ! Store number of vertives already processed
        nvt = rexport%p_Rtriangulation%NVT

        ! Write coordinates of edge midpoints?
        IF ((IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
          
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          DO imt=1,rexport%p_Rtriangulation%NMT
            ivt1 = p_IverticesAtEdge(1,imt)
            ivt2 = p_IverticesAtEdge(2,imt)
            X(nvt+imt) = 0.5*(REAL(p_DvertexCoords(1,ivt1)) + &
                              REAL(p_DvertexCoords(1,ivt2)))
            Y(nvt+imt) = 0.0E0
            Z(nvt+imt) = 0.0E0
          END DO
          
          ! Store number of vertives already processed
          nvt = rexport%p_Rtriangulation%NVT + rexport%p_Rtriangulation%NMT
        END IF

        ! Write coordinates of element midpoints?
        IF ((IAND(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
          
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
              p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          DO iel=1,rexport%p_Rtriangulation%NEL
            
            dx = 0.0E0
            
            DO i=1,UBOUND(p_IverticesAtElement,1)
              ivt = p_IverticesAtElement(i,iel)
              IF (ivt .NE. 0) THEN
                dx = dx + REAL(p_DvertexCoords(1,ivt))
              ELSE
                ! We have only (i-1) vertices in that element; 
                ! happens e.g. in triangles that are mixed into a quad mesh.
                ! We stop here.
                EXIT
              END IF
            END DO
            
            ! If all vertices of the element are touched, there is i=NVE+1.
            ! Divide by the number of vertices to get the coordinate of the 
            ! midpoint of the element.
            dx = dx / REAL(i-1)
            
            X(nvt+iel) = dx
            Y(nvt+iel) = 0.0E0
            Z(nvt+iel) = 0.0E0
          END DO
          
        END IF


      CASE(NDIM2D)
        DO ivt=1,rexport%p_Rtriangulation%NVT
          X(ivt) =  REAL(p_DvertexCoords(1,ivt))
          Y(ivt) =  REAL(p_DvertexCoords(2,ivt))
          Z(ivt) = 0.0E0
        END DO

        ! Store number of vertives already processed
        nvt = rexport%p_Rtriangulation%NVT

        ! Write coordinates of edge midpoints?
        IF ((IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
          
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          DO imt=1,rexport%p_Rtriangulation%NMT
            ivt1 = p_IverticesAtEdge(1,imt)
            ivt2 = p_IverticesAtEdge(2,imt)
            X(nvt+imt) = 0.5*(REAL(p_DvertexCoords(1,ivt1)) + &
                              REAL(p_DvertexCoords(1,ivt2)))
            Y(nvt+imt) = 0.5*(REAL(p_DvertexCoords(2,ivt1)) + &
                              REAL(p_DvertexCoords(2,ivt2)))
            Z(nvt+imt) = 0.0E0
          END DO
          
          ! Store number of vertives already processed
          nvt = rexport%p_Rtriangulation%NVT + rexport%p_Rtriangulation%NMT
        END IF

        ! Write coordinates of element midpoints?
        IF ((IAND(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
          
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          DO iel=1,rexport%p_Rtriangulation%NEL
            
            dx = 0.0E0
            dy = 0.0E0
            
            DO i=1,UBOUND(p_IverticesAtElement,1)
              ivt = p_IverticesAtElement(i,iel)
              IF (ivt .NE. 0) THEN
                dx = dx + REAL(p_DvertexCoords(1,ivt))
                dy = dy + REAL(p_DvertexCoords(2,ivt))
              ELSE
                ! We have only (i-1) vertices in that element; 
                ! happens e.g. in triangles that are mixed into a quad mesh.
                ! We stop here.
                EXIT
              END IF
            END DO
            
            ! If all vertices of the element are touched, there is i=NVE+1.
            ! Divide by the number of vertices to get the coordinate of the 
            ! midpoint of the element.
            dx = dx / REAL(i-1)
            dy = dy / REAL(i-1)
            
            X(nvt+iel) = dx
            Y(nvt+iel) = dy
            Z(nvt+iel) = 0.0E0
          END DO
          
        END IF


      CASE(NDIM3D)
        DO ivt=1,rexport%p_Rtriangulation%NVT
          X(ivt) =  REAL(p_DvertexCoords(1,ivt))
          Y(ivt) =  REAL(p_DvertexCoords(2,ivt))
          Z(ivt) =  REAL(p_DvertexCoords(3,ivt))
        END DO

        ! Store number of vertives already processed
        nvt = rexport%p_Rtriangulation%NVT

        ! Write coordinates of edge midpoints?
        IF ((IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
          
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          DO imt=1,rexport%p_Rtriangulation%NMT
            ivt1 = p_IverticesAtEdge(1,imt)
            ivt2 = p_IverticesAtEdge(2,imt)
            X(nvt+imt) = 0.5*(REAL(p_DvertexCoords(1,ivt1)) + &
                              REAL(p_DvertexCoords(1,ivt2)))
            Y(nvt+imt) = 0.5*(REAL(p_DvertexCoords(2,ivt1)) + &
                              REAL(p_DvertexCoords(2,ivt2)))
            Z(nvt+imt) = 0.5*(REAL(p_DvertexCoords(3,ivt1)) + &
                              REAL(p_DvertexCoords(3,ivt2)))
          END DO

          ! Store number of vertives already processed
          nvt = rexport%p_Rtriangulation%NVT + rexport%p_Rtriangulation%NMT
        END IF

        ! Write coordinates of element midpoints?
        IF ((IAND(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
            (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
          
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          DO iel=1,rexport%p_Rtriangulation%NEL
            
            dx = 0.0E0
            dy = 0.0E0
            dz = 0.0E0
            
            DO i=1,UBOUND(p_IverticesAtElement,1)
              ivt = p_IverticesAtElement(i,iel)
              IF (ivt .NE. 0) THEN
                dx = dx + REAL(p_DvertexCoords(1,ivt))
                dy = dy + REAL(p_DvertexCoords(2,ivt))
                dz = dz + REAL(p_DvertexCoords(3,ivt))
              ELSE
                ! We have only (i-1) vertices in that element; 
                ! happens e.g. in triangles that are mixed into a quad mesh.
                ! We stop here.
                EXIT
              END IF
            END DO
            
            ! If all vertices of the element are touched, there is i=NVE+1.
            ! Divide by the number of vertices to get the coordinate of the 
            ! midpoint of the element.
            dx = dx / REAL(i-1)
            dy = dy / REAL(i-1)
            dz = dz / REAL(i-1)
            
            X(nvt+iel) = dx
            Y(nvt+iel) = dy
            Z(nvt+iel) = dz
          END DO
          
        END IF
        
      END SELECT
      
      ! Write nodes to file
      CALL fgmvwrite_node_data(rexport%nvertices, X, Y, Z)
      
      ! Deallocate temporal memory
      DEALLOCATE(X, Y, Z)

    END IF

    
    ! Mesh connectivity / Cells:
    
    IF (rexport%saltFileCells .NE. "") THEN
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!
      
      CALL fgmvwrite_cells_fromfile(TRIM(rexport%saltFileCells), INT(rexport%ncells))
      
    ELSE

      ! Write the connectivity to the mesh - i.e. the cells.
      CALL fgmvwrite_cell_header(INT(rexport%ncells))
      
      IF (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .EQ. 0) THEN
        
        SELECT CASE (rexport%p_rtriangulation%ndim)
          
        CASE (NDIM1D)
          
          ! Standard mesh.
          DO iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            DO i=1,UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
            END DO
            
            ! We have i-1 vertices on that element -- so what is it?
            SELECT CASE (i-1)
            CASE (2)
              ! Line in 1D
              nod2ids(1) = INT(p_IverticesAtElement(1,iel))
              nod2ids(2) = INT(p_IverticesAtElement(2,iel))
              CALL fgmvwrite_cell_type('line 2',2,nod2ids)
              
            CASE DEFAULT
              CALL output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            END SELECT
            
          END DO
          
        CASE (NDIM2D)
          
          ! Standard mesh.
          DO iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            DO i=1,UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
            END DO
            
            ! We have i-1 vertices on that element -- so what is it?
            SELECT CASE (i-1)
              
            CASE (3)
              ! Triangle
              nod3ids(1) = p_IverticesAtElement(1,iel)
              nod3ids(2) = p_IverticesAtElement(2,iel)
              nod3ids(3) = p_IverticesAtElement(3,iel)
              CALL fgmvwrite_cell_type('tri 3',3,nod3ids)
              
            CASE (4)
              ! Quad
              nod4ids(1) = p_IverticesAtElement(1,iel)
              nod4ids(2) = p_IverticesAtElement(2,iel)
              nod4ids(3) = p_IverticesAtElement(3,iel)
              nod4ids(4) = p_IverticesAtElement(4,iel)
              CALL fgmvwrite_cell_type('quad 4',4,nod4ids)
              
            CASE DEFAULT
              CALL output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            END SELECT
            
          END DO
          
        CASE (NDIM3D)
          
          ! Standard mesh.
          DO iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            DO i=1,UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
            END DO
            
            ! We have i-1 vertices on that element -- so what is it?
            SELECT CASE (i-1)
              
            CASE (4)
              ! Tetrahedron
              nod4ids(1) = p_IverticesAtElement(1,iel)
              nod4ids(2) = p_IverticesAtElement(2,iel)
              nod4ids(3) = p_IverticesAtElement(3,iel)
              nod4ids(4) = p_IverticesAtElement(4,iel)
              CALL fgmvwrite_cell_type('ptet4 4',4,nod4ids)

            CASE (5)
              ! Pyramid
              nod5ids(1) = p_IverticesAtElement(1,iel)
              nod5ids(2) = p_IverticesAtElement(2,iel)
              nod5ids(3) = p_IverticesAtElement(3,iel)
              nod5ids(4) = p_IverticesAtElement(4,iel)
              nod5ids(5) = p_IverticesAtElement(5,iel)
              CALL fgmvwrite_cell_type('ppyrmd5 5',5,nod5ids)
              
            CASE(6)
              ! Prism
              nod6ids(1) = p_IverticesAtElement(1,iel)
              nod6ids(2) = p_IverticesAtElement(2,iel)
              nod6ids(3) = p_IverticesAtElement(3,iel)
              nod6ids(4) = p_IverticesAtElement(4,iel)
              nod6ids(5) = p_IverticesAtElement(5,iel)
              nod6ids(6) = p_IverticesAtElement(6,iel)
              CALL fgmvwrite_cell_type('pprism6 6',6,nod6ids)

            CASE (8)
              ! Hexahedron
              nod8ids(1) = p_IverticesAtElement(1,iel)
              nod8ids(2) = p_IverticesAtElement(2,iel)
              nod8ids(3) = p_IverticesAtElement(3,iel)
              nod8ids(4) = p_IverticesAtElement(4,iel)
              nod8ids(5) = p_IverticesAtElement(5,iel)
              nod8ids(6) = p_IverticesAtElement(6,iel)
              nod8ids(7) = p_IverticesAtElement(7,iel)
              nod8ids(8) = p_IverticesAtElement(8,iel)
              CALL fgmvwrite_cell_type('phex8 8',8,nod8ids)
              
            CASE DEFAULT
              CALL output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            END SELECT
            
          END DO
          
        END SELECT
        
      ELSE

        ! 1x refined mesh
        
        ! The edge numbers give the numbers of the edge midpoints and thus
        ! the numbers of the new vertices on the once refined mesh.
        CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElement)
        
        SELECT CASE (rexport%p_rtriangulation%ndim)
          
        CASE (NDIM1D)
          
          DO iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            DO i=1,UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
            END DO
            
            ! We have i-1 vertices on that element -- so what is it?
            SELECT CASE (i-1)
            CASE (2)
              ! Line in 1D.
              !
              ! The coarse grid element is
              !
              !   1 -----IEL----- 2
              !
              ! The once refined element is
              !
              !   1 -- IEL -- 1* -- NEL+IEL -- 2
              !
              ! Write the connectivity of element IEL
              nod2ids(1) = p_IverticesAtElement(1,iel)
              nod2ids(2) = p_IedgesAtElement(1,iel)
              CALL fgmvwrite_cell_type('line 2',2,nod2ids)

            CASE DEFAULT
              CALL output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')             
            END SELECT
            
          END DO
          
          DO iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            DO i=1,UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
            END DO
            
            ! We have i-1 vertices on that element -- so what is it?
            SELECT CASE (i-1)
            CASE (2)
              ! Line in 1D.
              !
              ! Element "NEL+1"
              nod2ids(1) = p_IverticesAtElement(1,iel)
              nod2ids(2) = p_IedgesAtElement(1,iel)
              CALL fgmvwrite_cell_type('line 2',2,nod2ids)

            CASE DEFAULT
              CALL output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')   
            END SELECT
            
          END DO
          
        CASE (NDIM2D)
          
          DO iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            DO i=1,UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
            END DO
            
            ! We have i-1 vertices on that element -- so what is it?
            SELECT CASE (i-1)
              
            CASE (3)
              ! Triangle.
              !
              ! Let a coarse grid triangle be locally numbered as:
              !
              !   2 
              !   |  \
              !   |    \
              !   | IEL  \
              !   |        \
              !   3----------1
              ! 
              ! Then the refinement process assigns the following numbers:
              !
              !   2_
              !   |  \_
              !   | NEL \_
              !   | +1     \_
              !   2*-------- 1* 
              !   | \_  IEL=1|  \_
              !   |NEL \_    | NEL \_
              !   |+3     \_ | +2     \_
              !   3----------3*---------1
              !
              ! So write the edge(-midpoint) numbers as corner numbers
              ! of the "main" element
              
              nod3ids(1) = p_IedgesAtElement(1,iel)
              nod3ids(2) = p_IedgesAtElement(2,iel)
              nod3ids(3) = p_IedgesAtElement(3,iel)
              CALL fgmvwrite_cell_type('tri 3',3,nod3ids)
              
            CASE (4)
              ! Quad
              !
              ! Let a coarse grid quad be locally numbered as
              !
              !  4-------3
              !  |       |
              !  |  IEL  |
              !  |       |
              !  1-------2
              !
              ! Then the refinement process assigns the following numbers:
              !
              !  4-----7-----3
              !  |NEL+3|NEL+2|
              !  |     |     |
              !  8-----9-----6
              !  |IEL=1|NEL+1|
              !  |     |     |
              !  1-----5-----2
              !
              ! So construct the corners of the 'smaller' elements from the
              ! corners of the coarse grid element, the numbers of the
              ! edge(-midpoint)s of the element and the number of the midpoint
              ! of the element -- which is defined as the element number itself.

              nod4ids(1) = p_IverticesAtElement(1,iel)
              nod4ids(2) = p_IedgesAtElement(1,iel)
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(4,iel)
              CALL fgmvwrite_cell_type('quad 4',4,nod4ids)

            CASE DEFAULT
              CALL output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            END SELECT
            
          END DO
          
          DO iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            DO i=1,UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(i,iel) .EQ. 0) EXIT
            END DO
            
            ! We have i-1 vertices on that element -- so what is it?
            SELECT CASE (i-1)
              
            CASE (3)
              ! Triangle.
              
              ! Element "NEL+1"
              nod3ids(1) = p_IverticesAtElement(2,iel)
              nod3ids(2) = p_IedgesAtElement(2,iel)
              nod3ids(3) = p_IedgesAtElement(1,iel)
              CALL fgmvwrite_cell_type('tri 3',3,nod3ids)
              
              ! Element "NEL+2"
              nod3ids(1) = p_IverticesAtElement(1,iel)
              nod3ids(2) = p_IedgesAtElement(1,iel)
              nod3ids(3) = p_IedgesAtElement(3,iel)
              CALL fgmvwrite_cell_type('tri 3',3,nod3ids)
              
              ! Element "NEL+3"
              nod3ids(1) = p_IverticesAtElement(3,iel)
              nod3ids(2) = p_IedgesAtElement(3,iel)
              nod3ids(3) = p_IedgesAtElement(2,iel)
              CALL fgmvwrite_cell_type('tri 3',3,nod3ids)
              
            CASE (4)
              ! Quad
              !
              ! Element "NEL+1"
              nod4ids(1) = p_IverticesAtElement(2,iel)
              nod4ids(2) = p_IedgesAtElement(2,iel)
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(1,iel)
              CALL fgmvwrite_cell_type('quad 4',4,nod4ids)
              
              ! Element "NEL+2"
              nod4ids(1) = p_IverticesAtElement(3,iel)
              nod4ids(2) = p_IedgesAtElement(3,iel)
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(2,iel)
              CALL fgmvwrite_cell_type('quad 4',4,nod4ids)
              
              ! Element "NEL+3"
              nod4ids(1) = p_IverticesAtElement(4,iel)
              nod4ids(2) = p_IedgesAtElement(4,iel)
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(3,iel)
              CALL fgmvwrite_cell_type('quad 4',4,nod4ids)

            CASE DEFAULT
              CALL output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            END SELECT
            
          END DO
          
        CASE (NDIM3D)
          
          CALL output_line ('GMV export for 1x refined mesh in 3D'//&
              ' not implemented!', OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
          CALL sys_halt()
          
        END SELECT
        
      END IF
      
    END IF


    !----------------------------------------------------
    ! Write material names -- if specified
    !

    IF (rexport%saltFileMaterials .NE. "") THEN
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!

      CALL fgmvwrite_material_fromfile(TRIM(rexport%saltFileMaterials))
      
    ELSE
      
      ! Cell material
      IF (ASSOCIATED(rexport%ScellMaterials)) THEN
        ! GMV only supports <= 1000 materials!
        CALL fgmvwrite_material_header(MIN(1000,INT(SIZE(rexport%ScellMaterials))),0)
        DO i=1,MIN(1000,SIZE(rexport%ScellMaterials))
          ! GMV supports only <= 8 characters and does not allow spaces
          ! in the material name. We replace all invalid spaces by "_".
          CALL fgmvwrite_material_name(&
              sys_charreplace(TRIM(rexport%ScellMaterials(i)),' ','_'))
        END DO
        
        IF (rexport%hIcellMaterial .NE. ST_NOHANDLE) THEN
          ! Write a list of material id's. For every cell, we specify its
          ! material by the material number.
          CALL storage_getbase_int (rexport%hIcellMaterial,p_Idata)
          DO i=1,SIZE(p_Idata)
            CALL fgmvwrite_material_ids(INT(p_Idata(i)), 0)
          END DO
        END IF
      END IF
      
      ! Vertex materials; coincide with cell materials if not specified.
      IF (ASSOCIATED(rexport%SvertexMaterials)) THEN
        ! GMV only supports <= 1000 materials!
        CALL fgmvwrite_material_header(MIN(1000,INT(SIZE(rexport%SvertexMaterials))),1)
        DO i=1,MIN(1000,SIZE(rexport%SvertexMaterials))
          ! GMV supports only <= 8 characters and does not allow spaces
          ! in the material name. We replace all invalid spaces by "_".
          CALL fgmvwrite_material_name(&
              sys_charreplace(TRIM(rexport%SvertexMaterials(i)),' ','_'))
        END DO
        
        IF (rexport%hIvertexMaterial .NE. ST_NOHANDLE) THEN
          ! Write a list of material id's. For every vertex, we specify its
          ! material by the material number.
          CALL storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
          DO i=1,SIZE(p_Idata)
            CALL fgmvwrite_material_ids(INT(p_Idata(i)), 1)
          END DO
        END IF
      ELSE
        IF (ASSOCIATED(rexport%ScellMaterials)) THEN
          ! GMV only supports <= 1000 materials!
          CALL fgmvwrite_material_header(MIN(1000,INT(SIZE(rexport%ScellMaterials))),1)
          DO i=1,MIN(1000,SIZE(rexport%ScellMaterials))
            ! GMV supports only <= 8 characters and does not allow spaces
            ! in the material name. We replace all invalid spaces by "_".
            CALL fgmvwrite_material_name(&
                sys_charreplace(TRIM(rexport%ScellMaterials(i)),' ','_'), 1)
          END DO
        END IF
        
        IF (rexport%hIvertexMaterial .NE. ST_NOHANDLE) THEN
          ! Write a list of material id's. For every vertex, we specify its
          ! material by the material number.
          CALL storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
          DO i=1,SIZE(p_Idata)
            CALL fgmvwrite_material_ids(INT(p_Idata(i)), 1)
          END DO
        END IF
      END IF
    END IF
    
    
    IF (rexport%nvariables .GT. 0) THEN
      !----------------------------------------------------
      ! Write a velocity field -- if there is one
      !
      ! To write the velocity field, we have to look it up in the
      ! set of variables.
      ! Search for the X-velocity and write it out.
      ! If there is no X-velocity, there is no velocity at all.
      ! There may be a velocity field given in vertices or in cells,
      ! so we have to search twice!
      !
      
      ! Look for cell based velocity.
      DO i=1,rexport%nvariables
        IF ((IAND(rexport%p_IvariableSpec(i),UCD_VAR_XVELOCITY) .NE. 0) .AND. &
            (rexport%p_IvariableBase(i) .EQ. UCD_BASE_ELEMENT)) THEN
          
          CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)

          ! Allocate temporal memory
          ALLOCATE(X(rexport%ncells), Y(rexport%ncells), Z(rexport%ncells))
          
          ! Don't be confused! ivt=number of cell, as we are in the 
          ! 'cell-oriented' case here!!!
          DO ivt=1,rexport%ncells
            X(ivt) = p_Ddata(ivt)
          END DO
          
          ! Find the Y-velocity
          DO j=1,rexport%nvariables
            IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_YVELOCITY) .NE. 0) .AND. &
                (rexport%p_IvariableBase(j) .EQ. UCD_BASE_ELEMENT)) THEN
              
              ! Found it. Write it out.
              CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              DO ivt=1,rexport%ncells
                Y(ivt) = p_Ddata(ivt)
              END DO
              
              EXIT
            END IF
          END DO
          IF (j .GT. rexport%nvariables) THEN
            ! Not found. Write out 0's instead.
            DO ivt=1,rexport%ncells
              Y(ivt) = 0.0E0
            END DO
          END IF
          
          ! Find the Z-velocity
          DO j=1,rexport%nvariables
            IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_ZVELOCITY) .NE. 0) .AND. &
                (rexport%p_IvariableBase(j) .EQ. UCD_BASE_ELEMENT)) THEN
              
              ! Found it. Write it out.
              CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              DO ivt=1,rexport%ncells
                Z(ivt) = p_Ddata(ivt)
              END DO
              
              EXIT
            END IF
          END DO
          IF (j .GT. rexport%nvariables) THEN
            ! Not found. Write out 0's instead.
            DO ivt=1,rexport%ncells
              Z(ivt) = 0.0E0
            END DO
          END IF
          
          ! Write cellbased velocities
          CALL fgmvwrite_velocity_data(0,X,Y,Z)

          ! Deallocate temporal memory
          DEALLOCATE(X,Y,Z)

        END IF
      END DO
      
      ! Look for vertex based velocity.
      DO i=1,rexport%nvariables
        IF ((IAND(rexport%p_IvariableSpec(i),UCD_VAR_XVELOCITY) .NE. 0) .AND. &
            (rexport%p_IvariableBase(i) .EQ. UCD_BASE_VERTEX)) THEN
          
          CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)

          ! Allocate temporal memory
          ALLOCATE(X(rexport%nvertices), Y(rexport%nvertices), Z(rexport%nvertices))

          DO ivt=1,rexport%nvertices
            X(ivt) = p_Ddata(ivt)
          END DO
          
          ! Find the Y-velocity
          DO j=1,rexport%nvariables
            IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_YVELOCITY) .NE. 0) .AND. &
                (rexport%p_IvariableBase(j) .EQ. UCD_BASE_VERTEX)) THEN
              
              ! Found it. Write it out.
              CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              DO ivt=1,rexport%nvertices
                Y(ivt) = p_Ddata(ivt)
              END DO
              
              EXIT
            END IF
          END DO
          IF (j .GT. rexport%nvariables) THEN
            ! Not found. Write out 0's instead.
            DO ivt=1,rexport%nvertices
              Y(ivt) = 0.0E0
            END DO
          END IF
          
          ! Find the Z-velocity
          DO j=1,rexport%nvariables
            IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_ZVELOCITY) .NE. 0) .AND. &
                (rexport%p_IvariableBase(j) .EQ. UCD_BASE_VERTEX)) THEN
              
              ! Found it. Write it out.
              CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              DO ivt=1,rexport%nvertices
                Z(ivt) = p_Ddata(ivt)
              END DO
              
              EXIT
            END IF
          END DO
          IF (j .GT. rexport%nvariables) THEN
            ! Not found. Write out 0's instead.
            DO ivt=1,rexport%nvertices
              Z(ivt) = 0.0E0
            END DO
          END IF
          
          ! Write cellbased velocities
          CALL fgmvwrite_velocity_data(1,X,Y,Z)

          ! Deallocate temporal memory
          DEALLOCATE(X,Y,Z)

        END IF
      END DO
      
      !----------------------------------------------------
      ! Write all variables which are not velocities
      DO i=1,rexport%nvariables
        IF (IAND(rexport%p_IvariableSpec(i), &
            UCD_VAR_XVELOCITY+UCD_VAR_YVELOCITY+UCD_VAR_ZVELOCITY) .EQ. 0) THEN
          
          CALL fgmvwrite_variable_header()
          
          IF (rexport%p_IvariableBase(i) .EQ. UCD_BASE_ELEMENT) THEN
            ! Cell based variable
            CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)

            ! Allocate temporal memory
            ALLOCATE(VAR(rexport%ncells))
            DO ivt=1,rexport%ncells
              VAR(ivt) = p_Ddata(ivt)
            END DO
            
            CALL fgmvwrite_variable_name_data(0, TRIM(rexport%p_SvariableNames(i)), VAR)

            ! Deallocate temporal memory
            DEALLOCATE(VAR)

          ELSE
            ! Vertex based variable
            CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
            
            ! Allocate temporal memory
            ALLOCATE(VAR(rexport%nvertices))
            DO ivt=1,rexport%nvertices
              VAR(ivt) = p_Ddata(ivt)
            END DO

            CALL fgmvwrite_variable_name_data(1, TRIM(rexport%p_SvariableNames(i)), VAR)

            ! Deallocate temporal memory
            DEALLOCATE(VAR)
          END IF
          
          CALL fgmvwrite_variable_endvars()
         
        END IF
        
      END DO ! i
      
    END IF


    !----------------------------------------------------
    ! Write polygon data

    IF (rexport%saltFilePolygons .NE. "") THEN
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!
      
      CALL fgmvwrite_polygons_fromfile(TRIM(rexport%saltFilePolygons))
      
    ELSE
      
      IF (ASSOCIATED(rexport%p_Hpolygons)) THEN
        
        ! At least one polygon
        CALL fgmvwrite_polygons_header()
        
        ! Materials
        CALL storage_getbase_int (rexport%hpolygonMaterial,p_Idata)
        
        ! Write all polygons.
        DO i=1,rexport%npolygons
          
          ! Coordinates        
          CALL storage_getbase_double2D (rexport%p_Hpolygons(i),p_Ddata2D)
          
          ! Allocate temporal memory
          ALLOCATE(X(UBOUND(p_Ddata2D,2)), Y(UBOUND(p_Ddata2D,2)), Z(UBOUND(p_Ddata2D,2)))
          
          ! Either we have 2D or 3D coordinates. 
          SELECT CASE(UBOUND(p_Ddata2D,1))
            
          CASE (NDIM2D)
            DO j=1,UBOUND(p_Ddata2D,2)
              X(j) = REAL(p_Ddata2D(1,j))
              Y(j) = REAL(p_Ddata2D(2,j))
              Z(j) = 0.0E0
            END DO
            
          CASE (NDIM3D)
            DO j=1,UBOUND(p_Ddata2D,2)
              X(j) = REAL(p_Ddata2D(1,j))
              Y(j) = REAL(p_Ddata2D(2,j))
              Z(j) = REAL(p_Ddata2D(3,j))
            END DO
            
          CASE DEFAULT
            CALL output_line ('Invalid spatial dimensions for polygon output!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
          END SELECT
          
          nverts = UBOUND(p_Ddata2D,2)
          matnum = p_Idata(i)
          
          CALL fgmvwrite_polygons_data(nverts, matnum, X, Y, Z)
          
          ! Deallocate temporal memory
          DEALLOCATE(X,Y,Z)
          
        END DO
        
        CALL fgmvwrite_polygons_endpoly()
        
      END IF
    END IF

    !----------------------------------------------------
    ! Write tracer coordinates and data
    IF (rexport%ntracers .NE. 0) THEN
            
      CALL storage_getbase_double2d(rexport%htracers,p_Ddata2D)
      
      ! Allocate temporal memory
      ALLOCATE(X(rexport%ntracers), Y(rexport%ntracers), Z(rexport%ntracers))

      SELECT CASE(UBOUND(p_Ddata,1))

      CASE(NDIM1D)
        DO j=1,rexport%ntracers
          X(j) = REAL(p_Ddata2D(j,1))
          Y(j) = 0.0E0
          Z(j) = 0.0E0
        END DO

      CASE(NDIM2D)
        DO j=1,rexport%ntracers
          X(j) = REAL(p_Ddata2D(j,1))
          Y(j) = REAL(p_Ddata2D(j,2))
          Z(j) = 0.0E0
        END DO

      CASE(NDIM3D)
        DO j=1,rexport%ntracers
          X(j) = REAL(p_Ddata2D(j,1))
          Y(j) = REAL(p_Ddata2D(j,2))
          Z(j) = REAL(p_Ddata2D(j,3))
        END DO
      END SELECT
      
      CALL fgmvwrite_tracers_header(INT(rexport%ntracers), X, Y, Z)
      
      ! Deallocate temporal memory
      DEALLOCATE(X,Y,Z)

!!$      
!!$      
!!$      ! Write tracer variables if specified
!!$      DO i=1,rexport%ntracerVariables
!!$        WRITE (mfile,'(A32)') rexport%StracerVariable(i)
!!$        
!!$        CALL storage_getbase_double (rexport%p_HtracerVariables(i), p_Ddata)
!!$        DO j=1,rexport%ntracers
!!$          WRITE (mfile,rexport%sdataFormat) p_Ddata(j)
!!$        END DO
!!$      END DO


      CALL fgmvwrite_tracers_endtrace()
      
    END IF

    !----------------------------------------------------
    ! Finally close the GMV file, finish
    
    CALL fgmvwrite_closefile()
    
    END SUBROUTINE

    !****************************************************************
    
    SUBROUTINE ucd_writeAVS (rexport)

    ! Specific output routine for AVS output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startAVS
    ! subroutine before.
    
    ! The export structure with all information
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
    ! local variables
    INTEGER :: mfile,i,j,k, num_ndata, num_cdata, ncells, nverts, &
        lenTemp, ivt1, ivt2
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords, p_DvertexRefined, &
        p_Ddata2D
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER, DIMENSION(:), POINTER :: p_IcellMaterial
    REAL(DP), DIMENSION(1:3) :: Dvert
    CHARACTER(LEN=SYS_STRLEN) :: sdl
    TYPE(t_ucdRefine) :: rrefine
      
      ! Get file unit and export format
      mfile = rexport%iunit
      sdl = rexport%sdataFormat
      
      ! Calculate grid refinement
      CALL ucd_refine(rrefine, rexport%p_rtriangulation, rexport%cflags)
      
      ! Get refined vertices
      IF (rrefine%h_DvertexCoords .NE. ST_NOHANDLE) THEN
        CALL storage_getbase_double2d (rrefine%h_DvertexCoords, p_DvertexRefined)
      ELSE
        p_DvertexRefined => NULL()
      END IF

      ! If the mesh is to be refined, then we take the refined elements,
      ! otherwise we use the elements from the triangulation
      IF (rrefine%h_IverticesAtElement .NE. ST_NOHANDLE) THEN
        CALL storage_getbase_int2d (rrefine%h_IverticesAtElement, &
            p_IverticesAtElement)
        ncells = rrefine%ncells
      ELSE
        CALL storage_getbase_int2d (rexport%p_rtriangulation%h_IverticesAtElement,&
            p_IverticesAtElement)
        ncells = rexport%p_rtriangulation%NEL
      END IF

      ! Get corner vertices
      CALL storage_getbase_double2d (rexport%p_Rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
      nverts = rexport%p_rtriangulation%NVT + rrefine%nvertices
      
      ! Get cell materials, if specified
      IF(rexport%hIcellMaterial .NE. ST_NOHANDLE) THEN
        CALL storage_getbase_int (rexport%hIcellMaterial, p_IcellMaterial)
      ELSE
        p_IcellMaterial => NULL()
      END IF
      
      ! First we need to count how many vertex-based and how many element-based
      ! variables we have.
      num_ndata = 0
      num_cdata = 0
      DO i = 1, rexport%nvariables
      
        IF (rexport%p_IvariableBase(i) .EQ. UCD_BASE_ELEMENT) THEN
          num_cdata = num_cdata + 1
        ELSE
          num_ndata = num_ndata + 1
        END IF
        
      END DO
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write the AVS header
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      WRITE(mfile, '(5I10)') nverts, ncells, num_ndata, num_cdata, 0
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write vertice coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      IF (UBOUND(p_DvertexCoords,1) .EQ. 3) THEN
        ! 3D coordinates
        
        ! Write corner vertices
        DO i=1, rexport%p_Rtriangulation%NVT
          WRITE(mfile, '(I10,3E16.7)') i, p_DvertexCoords(1:3, i)
        END DO
        
        ! Write refined vertices
        j = rexport%p_Rtriangulation%NVT
        DO i=1, rrefine%nvertices
          WRITE(mfile, '(I10,3E16.7)') j+i, p_DvertexRefined(1:3, i)
        END DO
        
      ELSE
        ! 2D coordinates
        
        ! Write corner vertices
        DO i=1, rexport%p_Rtriangulation%NVT
          WRITE(mfile, '(I10,3E16.7)') i, p_DvertexCoords(1:2, i), 0.0_DP
        END DO

        ! Write refined vertices
        j = rexport%p_Rtriangulation%NVT
        DO i=1, rrefine%nvertices
          WRITE(mfile, '(I10,3E16.7)') j+i, p_DvertexRefined(1:2, i), 0.0_DP
        END DO
        
      END IF
    
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write elements
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      DO j=1, ncells
      
        ! Count the number of vertices on that element
        DO i=1,UBOUND(p_IverticesAtElement,1)
          IF (p_IverticesAtElement(i,j) .EQ. 0) EXIT
        END DO
        
        ! Do we have cell materials?
        k = 1
        !IF (ASSOCIATED(p_IcellMaterial)) k = p_IcellMaterial(j)
        
        ! We have i-1 vertices on that element -- so what is it?
        SELECT CASE (i-1)
        CASE (3)
          ! Triangle
          WRITE(mfile, '(2I10,A,3I10)') j, k, ' tri', &
              p_IverticesAtElement(1:3,j)
          
        CASE (4)
          ! Quad
          WRITE(mfile, '(2I10,A,4I10)') j, k, ' quad', &
              p_IverticesAtElement(1:4,j)
          
        CASE DEFAULT
          CALL output_line ('Invalid element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeAVS')
          
        END SELECT
        
      END DO
      
      ! We do not need the mesh refinement anymore, so destroy it
      IF (rrefine%h_DvertexCoords .NE. ST_NOHANDLE) THEN
        CALL storage_free(rrefine%h_DvertexCoords)
      END IF
      IF (rrefine%h_IverticesAtElement .NE. ST_NOHANDLE) THEN
        CALL storage_free(rrefine%h_IverticesAtElement)
      END IF
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write node variable info
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      WRITE(mfile, '(I10)', ADVANCE='NO') num_ndata
      DO i=1, num_ndata
        WRITE(mfile, '(A)', ADVANCE='NO') " 1"
      END DO
      WRITE(mfile, '(A)') ""

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write node variable names
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      DO i=1, rexport%nvariables
      
        ! If this is a node variable, write its name
        IF (rexport%p_IvariableBase(i) .EQ. UCD_BASE_VERTEX) THEN
          WRITE(mfile, '(A,A)') &
             sys_charreplace(TRIM(rexport%p_SvariableNames(i)), ' ', '_'), &
             ", nounit"
        END IF
      
      END DO
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write node variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      IF (num_ndata .GT. 0) THEN
        
        ! Write all vertices
        DO i=1, rexport%nvertices
          
          ! Write vertice index
          WRITE(mfile, '(I10)', ADVANCE='NO') i
          
          ! Write all variable values
          DO j=1, rexport%nvariables
            
            IF (rexport%p_IvariableBase(j) .EQ. UCD_BASE_VERTEX) THEN
              CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              WRITE(mfile, sdl, ADVANCE='NO') p_Ddata(i)
            END IF
          
          END DO
          
          ! Write a line break
          WRITE(mfile, '(A)') ""
          
        END DO
        
      END IF
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variable info
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      WRITE(mfile, '(I10)', ADVANCE='NO') num_cdata
      DO i=1, num_cdata
        WRITE(mfile, '(A)', ADVANCE='NO') " 1"
      END DO
      WRITE(mfile, '(A)') ""

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variable names
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      DO i=1, rexport%nvariables
      
        ! If this is a element variable, write its name
        IF (rexport%p_IvariableBase(i) .EQ. UCD_BASE_ELEMENT) THEN
          WRITE(mfile, '(A,A)') &
             sys_charreplace(TRIM(rexport%p_SvariableNames(i)), ' ', '_'), &
             ", nounit"
        END IF
      
      END DO
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      IF (num_cdata .GT. 0) THEN
        
        ! Write all cells
        DO i=1, rexport%ncells
          
          ! Write cell index
          WRITE(mfile, '(I10)', ADVANCE='NO') i
          
          ! Write all variable values
          DO j=1, rexport%nvariables
            
            IF (rexport%p_IvariableBase(j) .EQ. UCD_BASE_ELEMENT) THEN
              CALL storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              WRITE(mfile, sdl, ADVANCE='NO') p_Ddata(i)
            END IF
          
          END DO
          
          ! Write a line break
          WRITE(mfile, '(A)') ""
          
        END DO
        
      END IF
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write comments
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      IF (rexport%ncommentBufSize .GT. 0) THEN
        
        i = 1
        ! Find the end of the first line
        DO j=i,rexport%ncommentBufSize
          IF (rexport%p_Scomments(j) .EQ. NEWLINE) EXIT
        END DO
        
        ! Write out all lines, one after the other
        DO WHILE (j .LE. rexport%ncommentBufSize)
          ! Write a comment identifier
          WRITE(mfile, '(A)', ADVANCE='NO') "# "

          ! Write the line (character i..j-1), continue with the next
          DO k=i,j-1
            WRITE(mfile,'(A)',ADVANCE='NO') rexport%p_Scomments(k)
          END DO
          WRITE(mfile,'(A)')
          
          ! Continue after the NEWLINE character, find the next NEWLINE.
          i = j+1
          DO j=i,rexport%ncommentBufSize
            IF (rexport%p_Scomments(j) .EQ. NEWLINE) EXIT
          END DO
        
        END DO
      
      END IF

    END SUBROUTINE

    !****************************************************************
    
    SUBROUTINE ucd_writeVTK (rexport)

    ! Specific output routine for VTK output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startVTK
    ! subroutine before.
    
    ! The export structure with all information
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
    ! local variables
    INTEGER :: mfile,i,j,k,jy,jz,num_ndata,num_cdata,ncls,nverts,ncells
    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords, p_DvertexRefined, &
        p_Ddata2D
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata, p_Dx, p_Dy, p_Dz
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER, DIMENSION(:), ALLOCATABLE :: p_InumVertsPerCell
    REAL(DP), DIMENSION(1:3) :: Dvert
    CHARACTER(LEN=SYS_STRLEN) :: sdl
    LOGICAL :: bQuadratic, bVec2Sc
    TYPE(t_ucdRefine) :: rrefine
      
      ! Get file unit and export format
      mfile = rexport%iunit
      sdl = rexport%sdataFormat
          
      ! Calculate grid refinement
      CALL ucd_refine(rrefine, rexport%p_rtriangulation, rexport%cflags)
      
      ! Get refined vertices
      IF (rrefine%h_DvertexCoords .NE. ST_NOHANDLE) THEN
        CALL storage_getbase_double2d (rrefine%h_DvertexCoords, p_DvertexRefined)
      ELSE
        p_DvertexRefined => NULL()
      END IF
      
      ! If the mesh is to be refined, then we take the refined elements,
      ! otherwise we use the elements from the triangulation
      IF (rrefine%h_IverticesAtElement .NE. ST_NOHANDLE) THEN
        CALL storage_getbase_int2d (rrefine%h_IverticesAtElement, &
            p_IverticesAtElement)
        ncells = rrefine%ncells
      ELSE
        CALL storage_getbase_int2d (rexport%p_rtriangulation%h_IverticesAtElement,&
            p_IverticesAtElement)
        ncells = rexport%p_rtriangulation%NEL
      END IF

      ! Get corner vertices
      CALL storage_getbase_double2d (rexport%p_rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
      nverts = rexport%p_rtriangulation%NVT + rrefine%nvertices
      
      ! Get edges of elements (needed for quadratic cells)
      IF (rexport%p_rtriangulation%h_IedgesAtElement .NE. ST_NOHANDLE) THEN
        CALL storage_getbase_int2d (rexport%p_rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElement)
      END IF

      ! First we need to count how many vertex-based and how many element-based
      ! variables we have.
      num_ndata = 0
      num_cdata = 0
      DO i = 1, rexport%nvariables
      
        IF (rexport%p_IvariableBase(i) .EQ. UCD_BASE_ELEMENT) THEN
          num_cdata = num_cdata + 1
        ELSE
          num_ndata = num_ndata + 1
        END IF
        
      END DO
      
      ! Should we write vector components as scalars?
      bVec2Sc = (IAND(rexport%cparam, UCD_PARAM_VTK_VECTOR_TO_SCALAR) .NE. 0)

      ! Are we going to write quadratic elements?
      bQuadratic = (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .EQ. 0) .AND. &
                   (IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .AND. &
                   (IAND(rexport%cparam,UCD_PARAM_VTK_USE_QUADRATIC) .NE. 0)
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write VTK header
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      WRITE(mfile, '(A)') "# vtk DataFile Version 2.0"
      WRITE(mfile, '(A)') "Generated by FEATFlow 2.x"
      WRITE(mfile, '(A)') "ASCII"
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write Dataset header
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! By default, our grid is unstructured
      WRITE(mfile, '(A)') "DATASET UNSTRUCTURED_GRID"
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write vertices
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      WRITE(mfile, '(A,I10,A)') "POINTS", nverts, " double"
      
      ! Do we have 2D or 3D points?
      IF (UBOUND(p_DvertexCoords,1) .EQ. 3) THEN
        ! 3D coordinates
        
        ! Write corner vertices
        DO i=1, rexport%p_Rtriangulation%NVT
          WRITE(mfile, '(3E16.7)') p_DvertexCoords(1:3, i)
        END DO
        
        ! Write refined vertices
        DO i=1, rrefine%nvertices
          WRITE(mfile, '(3E16.7)') p_DvertexRefined(1:3, i)
        END DO
        
      ELSE IF (UBOUND(p_DvertexCoords,1) .EQ. 2) THEN
        ! 2D coordinates
        
        ! Write corner vertices
        DO i=1, rexport%p_Rtriangulation%NVT
          WRITE(mfile, '(3E16.7)') p_DvertexCoords(1:2, i), 0.0_DP
        END DO

        ! Write refined vertices
        DO i=1, rrefine%nvertices
          WRITE(mfile, '(3E16.7)') p_DvertexRefined(1:2, i), 0.0_DP
        END DO
        
      ELSE
        ! 1D coordinates
        
        ! Write corner vertices
        DO i=1, rexport%p_Rtriangulation%NVT
          WRITE(mfile, '(3E16.7)') p_DvertexCoords(1, i), 0.0_DP, 0.0_DP
        END DO

        ! Write refined vertices
        DO i=1, rrefine%nvertices
          WRITE(mfile, '(3E16.7)') p_DvertexRefined(1, i), 0.0_DP, 0.0_DP
        END DO

      END IF
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cells
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! First we need to count the number of vertices per cell and the total
      ! number of integers needed for the cell index list
      ncls = 0
      ALLOCATE(p_InumVertsPerCell(ncells))
      DO j = 1, ncells
      
        ! Count the number of vertices on that element
        DO i=1,UBOUND(p_IverticesAtElement,1)
          IF (p_IverticesAtElement(i,j) .EQ. 0) EXIT
        END DO
        
        ! Store number of vertices
        p_InumVertsPerCell(j) = i-1
        
        ! Store the number of vertices for this element + 1
        IF (bQuadratic) THEN
          ncls = ncls + 2*i - 1
        ELSE
          ncls = ncls + i
        END IF
        
      END DO
      
      ! Write cells
      WRITE(mfile, '(A, 2I10)') "CELLS", ncells, ncls

      ! Should we write quadratic cells?
      IF (bQuadratic) THEN
      
        ! Get offset of first edge-midpoint-vertice
        k = rexport%p_rtriangulation%NVT - 1

        DO j=1, ncells
        
          SELECT CASE (p_InumVertsPerCell(j))
          CASE (2)
            ! Quadratic edge
            WRITE(mfile, '(4I10)') 3, (p_IverticesAtElement(1,j)-1), &
                (p_IverticesAtElement(2,j)-1), k + j
                
          CASE (3)
            ! Quadratic triangle
            WRITE(mfile, '(7I10)') 6,(p_IverticesAtElement(1,j)-1), &
                (p_IverticesAtElement(2,j)-1), (p_IverticesAtElement(3,j)-1), &
                k + p_IedgesAtElement(1,j), k + p_IedgesAtElement(2,j), &
                k + p_IedgesAtElement(3,j)
            
          CASE (4)
            ! Quadratic quadrilateral
            WRITE(mfile, '(9I10)') 8, (p_IverticesAtElement(1,j)-1), &
                (p_IverticesAtElement(2,j)-1), (p_IverticesAtElement(3,j)-1), &
                (p_IverticesAtElement(4,j)-1), k + p_IedgesAtElement(1,j), &
                k + p_IedgesAtElement(2,j), k + p_IedgesAtElement(3,j), &
                k + p_IedgesAtElement(4,j)
            
          CASE DEFAULT
            CALL output_line ('Invalid element!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeVTK')
          END SELECT
          
        END DO

        ! Write cell types
        WRITE(mfile, '(A, I10)') "CELL_TYPES", ncells
        DO j=1, ncells
          SELECT CASE (p_InumVertsPerCell(j))
          CASE (2)
            ! quadratic edge
            WRITE(mfile, '(I4)') VTK_QUADRATIC_EDGE 
             
          CASE (3)
            ! quadratic triangle
            WRITE(mfile, '(I4)') VTK_QUADRATIC_TRIANGLE
            
          CASE (4)
            ! quadratic quadrilateral
            WRITE(mfile, '(I4)') VTK_QUADRATIC_QUAD
            
          END SELECT
        END DO

      ELSE

        DO j=1, ncells
        
          SELECT CASE (p_InumVertsPerCell(j))
          CASE (2)
            ! edge
            WRITE(mfile, '(3I10)') 2, (p_IverticesAtElement(1,j)-1), &
                (p_IverticesAtElement(2,j)-1)

          CASE (3)
            ! triangle
            WRITE(mfile, '(4I10)') 3, (p_IverticesAtElement(1,j)-1), &
                (p_IverticesAtElement(2,j)-1), (p_IverticesAtElement(3,j)-1)

          CASE (4)
            ! quadrilateral
            WRITE(mfile, '(5I10)') 4, (p_IverticesAtElement(1,j)-1), &
                (p_IverticesAtElement(2,j)-1), (p_IverticesAtElement(3,j)-1), &
                (p_IverticesAtElement(4,j)-1)
          
          CASE (8)
            ! hexahedron
            WRITE(mfile, '(9I10)') 8, (p_IverticesAtElement(1,j)-1), &
                (p_IverticesAtElement(2,j)-1), (p_IverticesAtElement(3,j)-1), &
                (p_IverticesAtElement(4,j)-1), (p_IverticesAtElement(5,j)-1), &
                (p_IverticesAtElement(6,j)-1), (p_IverticesAtElement(7,j)-1), &
                (p_IverticesAtElement(8,j)-1)
                
          CASE DEFAULT
            CALL output_line ('Invalid element!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeVTK')
          END SELECT
          
        END DO

        ! Write cell types
        WRITE(mfile, '(A, I10)') "CELL_TYPES", ncells
        DO j=1, ncells
          SELECT CASE (p_InumVertsPerCell(j))
          CASE (2)
            ! Edge
            WRITE(mfile, '(I4)') VTK_LINE
            
          CASE (3)
            ! Triangle
            WRITE(mfile, '(I4)') VTK_TRIANGLE
            
          CASE (4)
            ! Quadrilateral
            WRITE(mfile, '(I4)') VTK_QUAD
          
          CASE (8)
            ! Hexahedron
            WRITE(mfile, '(I4)') VTK_HEXAHEDRON
            
          END SELECT
        END DO

      END IF
      
      ! Now we do not need the vertice counts anymore
      DEALLOCATE(p_InumVertsPerCell)
      
      ! We also do not need the mesh refinement anymore, so destroy it
      IF (rrefine%h_DvertexCoords .NE. ST_NOHANDLE) THEN
        CALL storage_free(rrefine%h_DvertexCoords)
      END IF
      IF (rrefine%h_IverticesAtElement .NE. ST_NOHANDLE) THEN
        CALL storage_free(rrefine%h_IverticesAtElement)
      END IF
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write vertex variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Go through all variables
      IF (num_ndata .GT. 0) THEN
      
        WRITE(mfile, '(A,I10)') "POINT_DATA", nverts
      
        ! Loop through all variables
        DO j=1, rexport%nvariables
        
          ! Is it vertex- or element-based?
          IF (rexport%p_IvariableBase(j) .NE. UCD_BASE_VERTEX) CYCLE
            
          ! Is the value scalar?
          IF ((rexport%p_IvariableSpec(j) .NE. UCD_VAR_STANDARD) .AND. &
              (.NOT. bVec2Sc)) CYCLE
            
          ! Print some stuff
          WRITE(mfile, '(A,A,A)') "SCALARS ", &
              sys_charreplace(TRIM(rexport%p_SvariableNames(j)),' ', '_'), &
                " double 1"
          WRITE(mfile, '(A)') "LOOKUP_TABLE default"
          
          ! Go for the data
          CALL storage_getbase_double (rexport%p_Hvariables(j), p_Ddata)
          
          ! Write the data
          DO i=1, UBOUND(p_Ddata,1)
            WRITE(mfile, sdl) p_Ddata(i)
          END DO
            
        END DO
        
        ! Should we write vectors?
        IF (.NOT. bVec2Sc) THEN
        
          ! Loop through all vectors
          DO j=1, rexport%nvectors
          
            ! Is this a vertex-based vector?
            IF(rexport%p_Ivectors(1,j) .NE. 1) CYCLE
        
            p_Dx => NULL()
            p_Dy => NULL()
            p_Dz => NULL()

            ! Go for the X, Y and Z coordinates
            i = rexport%p_Ivectors(2,j)
            CALL storage_getbase_double (rexport%p_Hvariables(i), p_Dx)
            i = rexport%p_Ivectors(3,j)
            IF (i .NE. 0) THEN
              CALL storage_getbase_double (rexport%p_Hvariables(i), p_Dy)
            END IF
            i = rexport%p_Ivectors(4,j)
            IF (i .NE. 0) THEN
              CALL storage_getbase_double (rexport%p_Hvariables(i), p_Dz)
            END IF
            
            ! Make sure we have at least the X-coordinate

            IF (.NOT. ASSOCIATED(p_Dx)) THEN
              CALL output_line ('Error: Variable vector '//&
                  TRIM(sys_siL(j,10))//' does not have X-coordinates!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeVTK')
              
              ! Try next vector
              CYCLE
            END IF
            
            ! Print the stuff
            WRITE(mfile,'(A,A,A)') "VECTORS ", &
                sys_charreplace(TRIM(rexport%p_SvarVecNames(j)), ' ', '_'), &
                  " double"
            
            ! Write the vector
            DO i = 1, UBOUND(p_Dx, 1)
            
              ! Write X-coord
              WRITE(mfile, sdl, ADVANCE='NO') p_Dx(i)
              
              ! Write Y-coord
              IF (ASSOCIATED(p_Dy)) THEN
                WRITE(mfile, sdl, ADVANCE='NO') p_Dy(i)
              ELSE
                WRITE(mfile, sdl, ADVANCE='NO') 0.0_DP
              END IF
              
              ! Write Z-coord
              IF (ASSOCIATED(p_Dz)) THEN
                WRITE(mfile, sdl, ADVANCE='NO') p_Dz(i)
              ELSE
                WRITE(mfile, sdl, ADVANCE='NO') 0.0_DP
              END IF
              
              ! Write a new line
              WRITE(mfile, '(A)') ""
              
            END DO
            
          END DO
           
        END IF
        
      END IF
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Go through all variables
      IF (num_cdata .GT. 0) THEN
      
        WRITE(mfile, '(A,I10)') "CELL_DATA", ncells
      
        ! Loop through all variables
        DO j=1, rexport%nvariables
        
          ! Is it vertex- or element-based?
          IF (rexport%p_IvariableBase(j) .NE. UCD_BASE_ELEMENT) CYCLE
          
          ! Print some stuff
          WRITE(mfile,'(A,A,A)') "SCALARS ", &
              sys_charreplace(TRIM(rexport%p_SvariableNames(j)),' ','_'),&
                " double 1"
          WRITE(mfile, '(A)') "LOOKUP_TABLE default"
          
          ! Go for the data
          CALL storage_getbase_double (rexport%p_Hvariables(j), p_Ddata)
          
          ! Write variable data
          DO i=1, UBOUND(p_Ddata,1)
            WRITE(mfile, sdl) p_Ddata(i)
          END DO
          
        END DO
        
      END IF
      
      ! That's it

    END SUBROUTINE
    
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_moreVariables (rexport)

!<description>
  ! (Re)allocates memory for variables. Increases the number of available
  ! variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    INTEGER(I32), DIMENSION(:), POINTER :: p_IvariableSpec 
    INTEGER, DIMENSION(:), POINTER :: p_Hvariables,p_IvariableBase
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: p_SvariableNames
    INTEGER :: nsize
    
    nsize = 0
    IF (ASSOCIATED(rexport%p_IvariableSpec)) &
      nsize = SIZE(rexport%p_IvariableSpec)

    ALLOCATE(p_IvariableSpec(nsize+16))
    IF (ASSOCIATED(rexport%p_IvariableSpec)) THEN
      p_IvariableSpec(1:nsize) = rexport%p_IvariableSpec
      DEALLOCATE(rexport%p_IvariableSpec)
    END IF
    rexport%p_IvariableSpec => p_IvariableSpec
      
    ALLOCATE(p_IvariableBase(nsize+16))
    IF (ASSOCIATED(rexport%p_IvariableBase)) THEN
      p_IvariableSpec(1:nsize) = rexport%p_IvariableBase
      DEALLOCATE(rexport%p_IvariableBase)
    END IF
    rexport%p_IvariableBase => p_IvariableBase

    ALLOCATE(p_Hvariables(nsize+16))
    IF (ASSOCIATED(rexport%p_Hvariables)) THEN
      p_Hvariables(1:nsize) = rexport%p_Hvariables
      DEALLOCATE(rexport%p_Hvariables)
    END IF
    rexport%p_Hvariables => p_Hvariables

    ALLOCATE(p_SvariableNames(nsize+16))
    IF (ASSOCIATED(rexport%p_SvariableNames)) THEN
      p_SvariableNames(1:nsize) = rexport%p_SvariableNames
      DEALLOCATE(rexport%p_SvariableNames)
    END IF
    rexport%p_SvariableNames => p_SvariableNames

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_moreVectors (rexport)

!<description>
  ! (Re)allocates memory for vector variables. Increases the number of
  ! available variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    INTEGER, DIMENSION(:,:), POINTER :: p_Ivectors
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: p_SvarVecNames
    INTEGER :: nsize
    
    nsize = 0
    IF (ASSOCIATED(rexport%p_Ivectors)) &
      nsize = SIZE(rexport%p_Ivectors)

    ALLOCATE(p_Ivectors(4,nsize+16))
    IF (ASSOCIATED(rexport%p_Ivectors)) THEN
      p_Ivectors(:,1:nsize) = rexport%p_Ivectors(:,:)
      DEALLOCATE(rexport%p_Ivectors)
    END IF
    rexport%p_Ivectors => p_Ivectors
      
    ALLOCATE(p_SvarVecNames(nsize+16))
    IF (ASSOCIATED(rexport%p_SvarVecNames)) THEN
      p_SvarVecNames(1:nsize) = rexport%p_SvarVecNames
      DEALLOCATE(rexport%p_SvarVecNames)
    END IF
    rexport%p_SvarVecNames => p_SvarVecNames

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_morePolygons (rexport)

!<description>
  ! (Re)allocates memory for polygons. Increases the number of available
  ! variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    INTEGER, DIMENSION(:), POINTER :: p_Hpolygons
    INTEGER :: nsize

    nsize = 0
    IF (ASSOCIATED(rexport%p_Hpolygons)) &
      nsize = SIZE(rexport%p_Hpolygons)

    ALLOCATE(p_Hpolygons(nsize+16))
    IF (rexport%hpolygonMaterial .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('ucd_morePolygons','hpolygonMaterial',&
          INT(nsize+16,I32),ST_INT,&
          rexport%hpolygonMaterial,ST_NEWBLOCK_ZERO)
    ELSE
      CALL storage_realloc ('ucd_morePolygons', INT(nsize+16,I32), &
          rexport%hpolygonMaterial, ST_NEWBLOCK_ZERO)
    END IF
    
    IF (ASSOCIATED(rexport%p_Hpolygons)) THEN
      p_Hpolygons(1:SIZE(rexport%p_Hpolygons)) = rexport%p_Hpolygons
      DEALLOCATE(rexport%p_Hpolygons)
    END IF
    rexport%p_Hpolygons => p_Hpolygons

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_moreTracerVariables (rexport)

!<description>
  ! (Re)allocates memory for tracer variables. Increases the number of 
  ! available tracer variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    INTEGER, DIMENSION(:), POINTER :: p_Hvariables 
    INTEGER :: nsize

    nsize = 0
    IF (ASSOCIATED(rexport%p_HtracerVariables)) &
      nsize = SIZE(rexport%p_HtracerVariables)

    ALLOCATE(p_Hvariables(nsize+16))
    IF (ASSOCIATED(rexport%p_HtracerVariables)) THEN
      p_Hvariables(1:nsize) = rexport%p_HtracerVariables
      DEALLOCATE(rexport%p_HtracerVariables)
    END IF
    rexport%p_HtracerVariables => p_Hvariables

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addVariableVertexBased (rexport,sname,cvarSpec, &
      DdataVert, DdataMid, DdataElem)

!<description>
  ! Adds variable data to the ouput file identified by rexport, based
  ! on vertices.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Name of the variable.
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! Specification bitfield for the variable. A combination of the 
  ! UCD_VAR_xxxx flags for special-type variables (like x-/y-velocity).
  ! Standard value=UCD_VAR_STANDARD.
  INTEGER(I32), INTENT(IN) :: cvarSpec
  
  ! DdataVert(I) is the value of the variable in vertex I of the triangulation.
  REAL(DP), DIMENSION(:), INTENT(IN) :: DdataVert

  ! OPTIONAL: DdataMid(I) is the value of the variable in edge midpoint I of 
  ! the triangulation. Must be specified if DdataElem is specified!
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataMid

  ! OPTIONAL: DdataElem(I) is the value of the variable in element midpoint I of 
  ! the triangulation.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataElem
!</input>

!</subroutine>

    ! local varables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    IF (PRESENT(DdataElem) .AND. .NOT. PRESENT(DdataMid)) THEN
      CALL output_line ('Error in the parameters!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableVertexBased')
      CALL sys_halt()
    END IF
    
    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableVertexBased')
      CALL sys_halt()
    END IF
    
    ! Create a new variable. If necessary, increase the size of the buffer.
    IF (.NOT. ASSOCIATED(rexport%p_Hvariables)) THEN
      CALL ucd_moreVariables(rexport)
    END IF
    
    IF (rexport%nvariables .GE. SIZE(rexport%p_Hvariables)) THEN
      CALL ucd_moreVariables(rexport)
    END IF
    
    rexport%nvariables = rexport%nvariables+1
    
    rexport%p_IvariableSpec(rexport%nvariables) = cvarSpec
    
    rexport%p_IvariableBase(rexport%nvariables) = UCD_BASE_VERTEX
    
    rexport%p_SvariableNames(rexport%nvariables) = sname
    
    ! Allocate a new vector for the data
    CALL storage_new ('ucd_addVariableVertexBased','hvariable',&
        INT(rexport%nvertices,I32),ST_DOUBLE,&
        rexport%p_Hvariables(rexport%nvariables),&
        ST_NEWBLOCK_ZERO)
    
    ! Copy the vertex data into that vector
    CALL storage_getbase_double (rexport%p_Hvariables(rexport%nvariables),p_Ddata)
    CALL lalg_copyVectorDble(DdataVert(1:rexport%p_rtriangulation%NVT), &
        p_Ddata(1:rexport%p_rtriangulation%NVT))
    
    ! Copy edge midpoint data if available
    IF ((IAND(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .NE. 0) .OR. &
        (IAND(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .NE. 0) .OR. &
        (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      IF (PRESENT(DdataMid)) THEN
        CALL lalg_copyVectorDble( &
            DdataMid(1:rexport%p_rtriangulation%NMT), &
            p_Ddata(rexport%p_rtriangulation%NVT+1:rexport%p_rtriangulation%NVT+ &
                                                   rexport%p_rtriangulation%NMT))
      ELSE
        CALL output_line ('Warning. No edge midpoint data available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableVertexBased')
      END IF
    END IF
    
    ! Copy element midpoint data if available
    IF ((IAND(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .NE. 0) .OR. &
        (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0)) THEN
      IF (PRESENT(DdataElem)) THEN
        CALL lalg_copyVectorDble( &
            DdataElem(1:rexport%p_rtriangulation%NEL), &
            p_Ddata(rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+1: &
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+ &
                    rexport%p_rtriangulation%NEL))
      ELSE
        CALL output_line ('Warning. No element midpoint data available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableVertexBased')
      END IF
    END IF    

  END SUBROUTINE
  
!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addVarVertBasedVec (rexport, sname, DdataVert_X, &
      DdataVert_Y, DdataVert_Z, DdataMid_X, DdataMid_Y, DdataMid_Z, &
      DdataElem_X, DdataElem_Y, DdataElem_Z)

!<description>
  ! Adds variable vector data to the ouput file identified by rexport, based
  ! on vertices.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Name of the vector.
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! The variable for the X-component of the velocities.
  REAL(DP), DIMENSION(:), INTENT(IN) :: DdataVert_X
  
  ! OPTIONAL: The data array for the Y-component of the velocities.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataVert_Y

  ! OPTIONAL: The data array for the Z-component of the velocities.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataVert_Z

  ! OPTIONAL: DdataMid_X/Y/Z(I) is the X/Y/Z-coordinate of the variable in
  ! edge midpoint I of the triangulation. Must be specified if DdataElem
  ! is specified!
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataMid_X
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataMid_Y
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataMid_Z

  ! OPTIONAL: DdataElem_X/Y/Z(I) is the X/Y/Z-coordinate of the variable in
  ! element midpoint I of the triangulation.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataElem_X
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataElem_Y
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataElem_Z
!</input>

!</subroutine>

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVarVertBasedVec')
      CALL sys_halt()
    END IF
    
    ! Let's create a new vector first
    IF (.NOT. ASSOCIATED(rexport%p_Ivectors)) THEN
      CALL ucd_moreVectors(rexport)
    ENDIF
    
    IF (rexport%nvectors .GE. SIZE(rexport%p_Ivectors, 2)) THEN
      CALL ucd_moreVectors(rexport)
    END IF
    
    ! Increment vector count and store its name
    rexport%nvectors = rexport%nvectors + 1
    rexport%p_SvarVecNames(rexport%nvectors) = sname
    
    ! This vector is vertex-based
    rexport%p_Ivectors(1,rexport%nvectors) = 1
    
    ! There must be an X-component, so add it first
    CALL ucd_addVariableVertexBased(rexport, TRIM(sname)//'_X', &
        UCD_VAR_XVELOCITY, DdataVert_X, DdataMid_X, dDataElem_X)
    
    rexport%p_Ivectors(2,rexport%nvectors) = rexport%nvariables
    
    ! Is there an Y-component?
    IF (PRESENT(DdataVert_Y)) THEN

      CALL ucd_addVariableVertexBased(rexport, TRIM(sname)//'_Y', &
          UCD_VAR_YVELOCITY, DdataVert_Y, DdataMid_Y, dDataElem_Y)
      
      rexport%p_Ivectors(3,rexport%nvectors) = rexport%nvariables
   
    ELSE
      rexport%p_Ivectors(3,rexport%nvectors) = 0
    END IF
    
    ! Is there an Z-component?
    IF (PRESENT(DdataVert_Z)) THEN

      CALL ucd_addVariableVertexBased(rexport, TRIM(sname)//'_Z', &
          UCD_VAR_ZVELOCITY, DdataVert_Z, DdataMid_Z, dDataElem_Z)
      
      rexport%p_Ivectors(4,rexport%nvectors) = rexport%nvariables
   
    ELSE
      rexport%p_Ivectors(4,rexport%nvectors) = 0
    END IF
    
    ! That's it
  
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addVarElemBasedVec (rexport, sname, Ddata_X, Ddata_Y, Ddata_Z)

!<description>
  ! Adds variable vector data to the ouput file identified by rexport, based
  ! on elements.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Name of the vector.
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! The variable for the X-component of the velocities.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Ddata_X
  
  ! OPTIONAL: The data array for the Y-component of the velocities.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Ddata_Y

  ! OPTIONAL: The data array for the Z-component of the velocities.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Ddata_Z

!</input>

!</subroutine>

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVarElemBasedVec')
      CALL sys_halt()
    END IF
    
    ! Let's create a new vector first
    IF (.NOT. ASSOCIATED(rexport%p_Ivectors)) THEN
      CALL ucd_moreVectors(rexport)
    ENDIF
    
    IF (rexport%nvectors .GE. SIZE(rexport%p_Ivectors, 2)) THEN
      CALL ucd_moreVectors(rexport)
    END IF
    
    ! Increment vector count and store its name
    rexport%nvectors = rexport%nvectors + 1
    rexport%p_SvarVecNames(rexport%nvectors) = sname
    
    ! This vector is vertex-based
    rexport%p_Ivectors(1,rexport%nvectors) = 1
    
    ! There must be an X-component, so add it first
    CALL ucd_addVariableElementBased(rexport, TRIM(sname)//'_X', &
        UCD_VAR_XVELOCITY, Ddata_X)
    
    rexport%p_Ivectors(2,rexport%nvectors) = rexport%nvariables
    
    ! Is there an Y-component?
    IF (PRESENT(Ddata_Y)) THEN

      CALL ucd_addVariableElementBased(rexport, TRIM(sname)//'_Y', &
          UCD_VAR_YVELOCITY, Ddata_Y)
      
      rexport%p_Ivectors(3,rexport%nvectors) = rexport%nvariables
   
    ELSE
      rexport%p_Ivectors(3,rexport%nvectors) = 0
    END IF
    
    ! Is there an Z-component?
    IF (PRESENT(Ddata_Z)) THEN

      CALL ucd_addVariableElementBased(rexport, TRIM(sname)//'_Z', &
          UCD_VAR_ZVELOCITY, Ddata_Z)
      
      rexport%p_Ivectors(4,rexport%nvectors) = rexport%nvariables
   
    ELSE
      rexport%p_Ivectors(4,rexport%nvectors) = 0
    END IF
    
    ! That's it
  
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addVariableElementBased (rexport,sname,cvarSpec, &
      Ddata)

!<description>
  ! Adds variable data to the ouput file identified by rexport, based on elements.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Name of the variable.
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! Specification bitfield for the variable. A combination of the 
  ! UCD_VAR_xxxx flags.
  INTEGER(I32), INTENT(IN) :: cvarSpec
  
  ! DdataVert(I) os the value of the variable in element I of the triangulation.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Ddata
!</input>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(DP) :: dv
    INTEGER(PREC_ELEMENTIDX) :: iel

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableElementBased')
      CALL sys_halt()
    END IF
    
    ! Create a new variable. If necessary, increase the size of the buffer.
    IF (.NOT. ASSOCIATED(rexport%p_Hvariables)) THEN
      CALL ucd_moreVariables(rexport)
    END IF

    ! Create a new variable. If necessary, increase the size of the buffer.
    IF (rexport%nvariables .GE. SIZE(rexport%p_Hvariables)) THEN
      CALL ucd_moreVariables(rexport)
    END IF
    rexport%nvariables = rexport%nvariables+1
    rexport%p_IvariableSpec(rexport%nvariables) = cvarSpec
    rexport%p_IvariableBase(rexport%nvariables) = UCD_BASE_ELEMENT
    rexport%p_SvariableNames(rexport%nvariables) = sname
    
    ! Allocate a new vector for the data
    CALL storage_new ('ucd_addVariableVertexBased','hvariable',&
        INT(rexport%ncells,I32),ST_DOUBLE,&
        rexport%p_Hvariables(rexport%nvariables),&
        ST_NEWBLOCK_ZERO)
        
    ! Copy the element data into that vector
    CALL storage_getbase_double (rexport%p_Hvariables(rexport%nvariables),p_Ddata)
    CALL lalg_copyVectorDble(Ddata(1:rexport%p_rtriangulation%NEL), &
        p_Ddata(1:rexport%p_rtriangulation%NEL))
        
    ! On a once refined mesh, generate the data on the sub-cells by
    ! constant interpolation of the data in the large cell.
        
    ! Copy edge midpoint data if available
    IF (IAND(rexport%cflags,UCD_FLAG_ONCEREFINED) .NE. 0) THEN
    
      IF (rexport%p_rtriangulation%NDIM .EQ. NDIM2D) THEN
        ! Implicitly use 2-level ordering to get the numbers of the sub-elements
        DO iel=1,rexport%p_rtriangulation%NEL
          dv = p_Ddata(iel)
          p_Ddata(rexport%p_rtriangulation%NEL+3*(iel-1)+1) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+3*(iel-1)+2) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+3*(iel-1)+3) = dv
        END DO
      ELSE
        DO iel=1,rexport%p_rtriangulation%NEL
          dv = p_Ddata(iel)
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+1) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+2) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+3) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+4) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+5) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+6) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+7) = dv
        END DO
      END IF
      
    END IF
    
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addPolygon (rexport,DpolygonCoords,ipolygonMaterial)

!<description>
  ! Adds a polygon to the output file. A polygon is a set of lines.
  ! DpolygonCoords is a 2D or 3D array with (X,Y)-/(X,Y,Z)-coordinates.
  ! Line i is specified by the start coordinates DpolygonCoords(i)
  ! and the end coordinates DpolygonCoords(i+1).
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! A list of 2D (X,Y) or 3D (X,Y,Z) coordinates specifying the
  ! points which form the polygon.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DpolygonCoords
  
  ! OPTIONAL: A material identifier for the polygon. Whether polygons can have
  ! a material associated depends on the output file format.
  ! If not specified, a default material is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: ipolygonMaterial
!</input>

!</subroutine>

    ! local varables
    INTEGER(I32), DIMENSION(2) :: Ilength
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addPolygon')
      CALL sys_halt()
    END IF
    
    ! Create a new variable. If necessary, increase the size of the buffers.
    IF (.NOT. ASSOCIATED(rexport%p_Hpolygons)) THEN
      CALL ucd_morePolygons(rexport)
    END IF
    
    IF (rexport%npolygons .GE. SIZE(rexport%p_Hpolygons)) THEN
      CALL ucd_morePolygons(rexport)
    END IF
    rexport%npolygons = rexport%npolygons+1
    
    ! Allocate a new vector for the data
    Ilength = UBOUND(DpolygonCoords)
    CALL storage_new2D ('ucd_addPolygon','hpolygon',&
        Ilength,ST_DOUBLE,rexport%p_Hpolygons(rexport%npolygons),&
        ST_NEWBLOCK_NOINIT)
        
    ! Copy the coordinate data into that vector
    CALL storage_getbase_double2d (rexport%p_Hpolygons(rexport%npolygons),p_Ddata)
    p_Ddata = DpolygonCoords
    
    ! Copy the material
    CALL storage_getbase_int (rexport%hpolygonMaterial,p_Idata)
    IF (PRESENT(ipolygonMaterial)) THEN
      p_Idata(rexport%npolygons) = ipolygonMaterial
    ELSE
      ! Default material 0
      p_Idata(rexport%npolygons) = 0
    END IF

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_setTracers (rexport,DtracerCoordinates)

!<description>
  ! This routine adds tracer information to the export structure rexport.
  ! If there is already any tracer data attached to rexport, it will be
  ! deleted.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! A list of 2D (X,Y) or 3D (X,Y,Z) coordinates specifying the
  ! position of all the tracers.
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: DtracerCoordinates
!</input>

!</subroutine>

    ! local varables
    INTEGER(I32), DIMENSION(2) :: Ilength
    REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setTracers')
      CALL sys_halt()
    END IF
    
    ! Create a new variable. If necessary, increase the size of the buffers.
    IF (rexport%ntracers .NE. 0) THEN
      ! There are already tracers and tracer data attached.
      ! Clean up tracer information as we create a new set of tracers.
      CALL ucd_removeTracers (rexport)
    END IF
    
    ! Remember total number of tracers
    rexport%ntracers = UBOUND(DtracerCoordinates,2)
    
    ! Allocate a new vector for the data
    Ilength = UBOUND(DtracerCoordinates)
    CALL storage_new2D ('ucd_setTracers','htracers',&
        Ilength,ST_DOUBLE,rexport%htracers,ST_NEWBLOCK_NOINIT)
        
    ! Copy the coordinate data into that vector
    CALL storage_getbase_double2d (rexport%htracers,p_Ddata)
    p_Ddata = DtracerCoordinates
    
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_removeTracers (rexport)

!<description>
  ! Removes all information about tracers from the structure rexport.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    INTEGER :: i
    
    ! Remove all tracer variables if there are some
    IF (ASSOCIATED(rexport%p_HtracerVariables    )) THEN
      DO i=1,rexport%ntracervariables
        IF (rexport%p_HtracerVariables(i) .NE. ST_NOHANDLE) &
          CALL storage_free (rexport%p_HtracerVariables(i))
      END DO
      DEALLOCATE(rexport%p_HtracerVariables    )
    END IF

    ! Delete the coordinates of the tracers
    IF (rexport%htracers       .NE. ST_NOHANDLE) &
      CALL storage_free(rexport%htracers    )
        
    rexport%ntracers = 0

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addTracerVariable (rexport,sname,Ddata)

!<description>
  ! Adds tracer variable data to the ouput file identified by rexport.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Name of the tracer variable.
  CHARACTER(LEN=*), INTENT(IN) :: sname
  
  ! array [1..#Tracers] of double. For every tracer I, Ddata(I) is the
  ! value of the variable that should be associated to that tracer.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Ddata
!</input>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addTracerVariable')
      CALL sys_halt()
    END IF
    
    IF (rexport%ntracers .LE. 0) THEN
      CALL output_line ('No tracers specified!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addTracerVariable')
      CALL sys_halt()
    END IF
    
    IF (SIZE(Ddata) .LT. rexport%ntracers) THEN
      CALL output_line ('Ddata too small, more tracers than data!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addTracerVariable')
      CALL sys_halt()
    END IF

    ! Create a new variable. If necessary, increase the size of the buffer.
    IF (rexport%ntracerVariables .GE. SIZE(rexport%p_HtracerVariables)) THEN
      CALL ucd_moreTracerVariables(rexport)
    END IF
    rexport%ntracerVariables = rexport%ntracerVariables+1
    rexport%StracerVariable(rexport%ntracerVariables) = sname
    
    ! Allocate a new vector for the data
    CALL storage_new ('ucd_addTracerVariable','hvariable',&
        INT(rexport%ntracers,I32),ST_DOUBLE,&
        rexport%p_HtracerVariables(rexport%ntracerVariables),&
        ST_NEWBLOCK_ZERO)
        
    ! Copy the element data into that vector
    CALL storage_getbase_double (&
        rexport%p_HtracerVariables(rexport%ntracerVariables),p_Ddata)
    CALL lalg_copyVectorDble(Ddata(1:rexport%ntracers), &
        p_Ddata(1:rexport%ntracers))
        
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_setSimulationTime (rexport,dtime,ssimTimeFormat)

!<description>
  ! Specifies the simulation time, the export structure should represent.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Simulation time. Is written to the output file.
  REAL(DP), INTENT(IN) :: dtime
  
  ! OPTIONAL: Fortran format string, e.g. "(F20.5)"
  ! Allows to specify an output format of the simulation time. Whether or
  ! not this is used for writing to the file depends on the type of
  ! output format (GMV, AVS,...).
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ssimTimeFormat
!</input>

!</subroutine>

    CHARACTER(LEN=SYS_STRLEN) :: stext

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setSimulationTime')
      CALL sys_halt()
    END IF
    
    rexport%dsimulationTime = dtime
    
    ! Copy the output format string and overwrite the standard setting.
    IF (PRESENT(ssimTimeFormat)) THEN
      rexport%ssimTimeFormat = ssimTimeFormat
      
      ! Text the output format string -- to be sure it's valid.
      ! If not, a compiler error is thrown here! This is either
      ! a runtime error or simply a message on the screen.
      WRITE(stext,ssimTimeFormat) 0.0_DP
      IF (stext .EQ. "") THEN
        CALL output_line ('Invalid output format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ucd_setSimulationTime')
        CALL sys_halt()
      END IF
    END IF
  
  END SUBROUTINE
  
!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_setOutputNumberFormat (rexport,sformat)

!<description>
  ! Specifies the number format of numbers in the outpuot file.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! Fortran format string, e.g. "(E15.8)"
  ! Specifies an output format of double precision numbers. Whether or
  ! not and where this is used for writing to the file depends on the 
  ! type of output format (GMV, AVS,...).
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sformat
!</input>

!</subroutine>

    CHARACTER(LEN=SYS_STRLEN) :: stext

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setOutputNumberFormat')
      CALL sys_halt()
    END IF
    
    ! Copy the output format string and overwrite the standard setting.
    rexport%sdataFormat = sformat
    
    ! Text the output format string -- to be sure it's valid.
    ! If not, a compiler error is thrown here! This is either
    ! a runtime error or simply a message on the screen.
    WRITE(stext,sformat) 0.0_DP
    IF (stext .EQ. "") THEN
      CALL output_line ('Invalid output format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setOutputNumberFormat')
      CALL sys_halt()
    END IF
  
  END SUBROUTINE
  
!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addCommentLine (rexport,scomment)

!<description>
  ! Adds a comment line to the comment buffer.
  ! All comments are written to a single comment block at the
  ! beginning of the output file.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! The comment to be added to the output.
  CHARACTER(LEN=*), INTENT(IN) :: scomment
!</input>

!</subroutine>

    CHARACTER, DIMENSION(:), POINTER :: p_sbuf
    INTEGER :: i

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addCommentLine')
      CALL sys_halt()
    END IF
    
    ! Is there enough space in the output buffer? If not, reallocate it.
    IF (.NOT. ASSOCIATED(rexport%p_Scomments)) THEN
      ALLOCATE(rexport%p_Scomments(MAX(LEN(scomment)+1,4*SYS_STRLEN)))
    ELSE IF ( rexport%ncommentBufSize+LEN(scomment)+1 .GT. &
              SIZE(rexport%p_Scomments)) THEN
      ALLOCATE(p_sbuf (rexport%ncommentBufSize + MAX(LEN(scomment)+1,4*SYS_STRLEN)))
      p_sbuf(1:rexport%ncommentBufSize) = &
          rexport%p_Scomments(1:rexport%ncommentBufSize)
      DEALLOCATE(rexport%p_Scomments)
      rexport%p_Scomments => p_sbuf
    END IF

    ! Append the string -- character by character
    DO i=1,LEN(scomment)
      rexport%p_Scomments (rexport%ncommentBufSize+i) = scomment (i:i)
    END DO
    
    ! Append a NEWLINE character as line end
    rexport%p_Scomments (rexport%ncommentBufSize+LEN(scomment)+1) = NEWLINE

    rexport%ncommentBufSize = rexport%ncommentBufSize+LEN(scomment)+1
        
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_addParameterList (rexport,rparList)

!<description>
  ! Adds the parameter list rparList to the output file identified by
  ! rexport. All parameters in this parameter list will be written
  ! as comment block to the output file.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!<input>
  ! A parameter list containing all configuration parameters of a simulation.
  TYPE(t_parlist), INTENT(IN) :: rparList
!</input>

!</subroutine>

    CHARACTER, DIMENSION(:), POINTER :: p_sbuf,p_sparams

    IF (rexport%coutputFormat .EQ. UCD_FORMAT_NONE) THEN
      CALL output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addParameterList')
      CALL sys_halt()
    END IF
    
    ! Ask the parameter list to create a large array of characters
    ! containing the configuration.
    CALL parlst_getStringRepresentation (rparlist, p_sparams)
    
    ! May be the buffer is empty...
    IF (.NOT. ASSOCIATED(p_sparams)) RETURN
    
    ! By default, the parameter list uses the same separation character
    ! for all lines as we do. Therefore, we can simply add that large
    ! character block to our buffer.
    
    ! Is there enough space in the output buffer? If not, reallocate it.
    IF (.NOT. ASSOCIATED(rexport%p_Scomments)) THEN
      ALLOCATE(rexport%p_Scomments(MAX(SIZE(p_sparams)+1,4*SYS_STRLEN)))
    ELSE IF ( rexport%ncommentBufSize+SIZE(p_sparams) .GT. &
              SIZE(rexport%p_Scomments)) THEN
      ALLOCATE(p_sbuf (rexport%ncommentBufSize + MAX(SIZE(p_sparams),4*SYS_STRLEN)))
      p_sbuf(1:rexport%ncommentBufSize) = &
          rexport%p_Scomments(1:rexport%ncommentBufSize)
      DEALLOCATE(rexport%p_Scomments)
      rexport%p_Scomments => p_sbuf
    END IF
    
    ! Copy the data
    rexport%p_Scomments ( &
      rexport%ncommentBufSize+1:rexport%ncommentBufSize+SIZE(p_sparams)) = p_sparams

    rexport%ncommentBufSize = rexport%ncommentBufSize+SIZE(p_sparams)
    
    ! Remove data from the heap
    DEALLOCATE(p_sparams)
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE ucd_readGMV (sfilename,rexport,rtriangulation)

!<description>
  ! Reads a GMV file from the hard disc and creates a rexport structure
  ! containing the GMV data. The structure can be used to attach
  ! more data if desired. The filename for the export can be changed 
  ! by ucd_setFilename. Using ucd_write, the data can be written out
  ! to another output file.
  !
  ! If rtriangulation is undefined, the mesh in the GMV is read and
  ! a new triangulation based on the GMV data is created in
  ! rtriangulation.
  ! If rtriangulation is defined, rtriangulation must specify the mesh
  ! that corresponds to the one in the GMV file.
!</description>

!<input>
  ! Name of the GMV file.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
!</input>

!<output>
  ! The export structure where the GMV data is saved to.
  ! The structure can be used to attach more data. The filename for the
  ! export can be changed by ucd_setFilename.
  TYPE(t_ucdExport), INTENT(OUT) :: rexport
!</output>

!<inputoutput>
  ! Triangulation structure.
  ! If this triangulation structure is empty, the mesh is read from the GMV 
  ! and a new triangulation is created in rtriangulation.
  ! If this contains a valid mesh, the mesh is assumed to correspond to
  ! the one in the GMV file and based on that, vector data is read.
  ! The mesh in the GMV file is skipped.
  TYPE(t_triangulation), INTENT(INOUT), TARGET :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    INTEGER :: mfile,ilinelen,ios
    LOGICAL :: biscomment
    CHARACTER(LEN=SYS_STRLEN) :: sline,skey,sdummy
  
    ! Try to open the file
    IF (sfilename .NE. '') THEN
      CALL io_openFileForReading(sfilename, mfile, .TRUE.)
      IF (mfile .LT. 0) THEN
        CALL output_line ('Cannot open file "'//TRIM(sfilename)&
                //'" for writing!',OU_CLASS_ERROR,OU_MODE_STD,'ucd_readGMV')
        CALL sys_halt()
      END IF
    END IF

    !----------------------------------------------------
    ! Read the GMV header
    CALL io_readlinefromfile (mfile, sline, ilinelen, ios)

    ! The triangulation in rexport will be rtriangulation --
    ! regardless whether we build it or have it
    rexport%p_rtriangulation => rtriangulation

    ! Set the other flags in rexport to standard values
    rexport%coutputFormat = UCD_FORMAT_GMV
    rexport%cflags = UCD_FLAG_STANDARD

    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL

    ! Read each line and interpret it
    DO WHILE (ios .EQ. 0)
      CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
      
      ! Get the key word from the line
      CALL sys_tolower (sline)
      READ(sline,*) skey
      
      IF (skey .EQ. 'probtime') THEN
        
        !----------------------------------------------------
        ! Read the simulation time
        READ(sline,*) sdummy,rexport%dsimulationTime
        
      ELSE IF (skey .EQ. 'comments') THEN
      
        !----------------------------------------------------
        ! Read comments
        biscomment = .TRUE.
        DO WHILE ((ios .EQ. 0) .AND. biscomment)
          CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
          CALL sys_tolower (sline)
          READ(sline,'(A)') skey
          biscomment = (skey .NE. 'endcomm')
          IF (biscomment) THEN
            ! Add each comment line to our comment variable in rexport
            CALL ucd_addCommentLine (rexport,sline)
          END IF
        END DO
        
      ELSE IF (skey .EQ. 'nodes') THEN
        
        !----------------------------------------------------
        ! Read triangulation (or ignore it if alread given in rtriangulation)
        CALL read_triangulation (mfile,sline,rtriangulation)
        
        ! NEL/NVT has changed
        rexport%nvertices = rtriangulation%NVT
        rexport%ncells = rtriangulation%NEL
        
      ELSE IF (skey .EQ. 'material') THEN
      
        !----------------------------------------------------
        ! Read materials
        CALL read_materials (mfile,sline,rexport)
        
      ELSE IF (skey .EQ. 'velocity') THEN

        !----------------------------------------------------
        ! Read velocity data
        CALL read_velocity (mfile,sline,rexport)
        
      ELSE IF (skey .EQ. 'variable') THEN
      
        !----------------------------------------------------
        ! Read general variable data
        CALL read_variables (mfile,sline,rexport)
        
      ELSE IF (skey .EQ. 'polygons') THEN
      
        !----------------------------------------------------
        ! Read polygon data
        CALL read_polygondata (mfile,sline,rexport)
        
      ELSE IF (skey .EQ. 'tracers') THEN
      
        !----------------------------------------------------
        ! Read tracer data
        CALL read_tracerdata (mfile,sline,rexport)
        
      ELSE IF (skey .EQ. 'endgmv') THEN
        
        ! GMV finish
        ios = 1
        
      END IF
    END DO
    
    ! Close the file, finish.
    CLOSE (mfile)
    
  CONTAINS
  
    ! ---------------------------------------------------------------
    
    SUBROUTINE read_tracerdata (mfile,scommand,rexport)
    
    ! Reads data about tracers from the GMV file mfile.
    
    ! Handle to the GMV file
    INTEGER, INTENT(IN) :: mfile
    
    ! Command line in the GMV file with information about the tracers
    CHARACTER(LEN=*), INTENT(IN) :: scommand
    
    ! UCD structure where tracer data is saved to.
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
      ! Local variables
      INTEGER :: ilinelen,ios,ntracers
      CHARACTER(LEN=SYS_STRLEN) :: skey,sline
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DtracerCoordinates
      REAL(DP), DIMENSION(:), ALLOCATABLE :: DtracerTemp

      ! Interpret the command line, how many tracers do we have?
      READ (scommand,*) skey,ntracers
      
      ! Allocate memory for X/Y or X/Y/Z coordinates
      ALLOCATE (DtracerCoordinates(rexport%p_rtriangulation%ndim,ntracers))
      ALLOCATE (DtracerTemp(ntracers))
      
      ! Read all X-coordinates
      READ(mfile,*) DtracerCoordinates(1,:)
      
      ! Read all Y-coordinates
      READ(mfile,*) DtracerCoordinates(2,:)
      
      ! Probably read all Z-coordinates -- or ignore them
      IF (rexport%p_rtriangulation%ndim .EQ. NDIM3D) THEN
        READ(mfile,*) DtracerCoordinates(3,:)
      ELSE 
        READ(mfile,*) DtracerTemp(:)
      END IF
      
      ! There must be an 'endtrace' at the end
      CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
      CALL sys_tolower (sline)
      READ(sline,*) skey
      IF (skey .NE. 'endtrace') THEN
        CALL output_line ('Error reading GMV data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_tracerdata')
        CALL sys_halt() 
      END IF
      
      ! Add the tracer data
      CALL ucd_setTracers (rexport,DtracerCoordinates)
      
      DEALLOCATE(DtracerTemp,DtracerCoordinates)

    END SUBROUTINE
    
    ! ---------------------------------------------------------------
    
    SUBROUTINE read_polygondata (mfile,scommand,rexport)
    
    ! Reads data about polygons from the GMV file mfile.
    
    ! Handle to the GMV file
    INTEGER, INTENT(IN) :: mfile
    
    ! Last command line in the GMV file 
    CHARACTER(LEN=*), INTENT(IN) :: scommand
    
    ! UCD structure where data is saved to.
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
      ! Local variables
      INTEGER :: ilinelen,ios,npoints,imaterial
      LOGICAL :: bfinish
      CHARACTER(LEN=SYS_STRLEN) :: skey,sline
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dcoordinates
      REAL(DP), DIMENSION(:), ALLOCATABLE :: Dtemp

      ! Read the next line specifying basic data about the polygon
      CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
      CALL sys_tolower (sline)
      
      READ(sline,*) skey
      bfinish = (skey .NE. 'endtrace') 
      DO WHILE ((ios .EQ. 0) .AND. (.NOT. bfinish))
        
        ! Read material, #points
        WRITE (mfile,*) imaterial,npoints
        
        ! Allocate memory for the polygon
        ALLOCATE (Dcoordinates(rexport%p_rtriangulation%ndim,npoints))
        ALLOCATE (Dtemp(npoints))
        
        ! Read the polygon. All X-coordinates
        READ(mfile,*) Dcoordinates(1,:)
        
        ! Y-coordinates
        READ(mfile,*) Dcoordinates(2,:)
        
        ! Probably Z-coordinates
        IF (rexport%p_rtriangulation%ndim .EQ. 3) THEN
          READ(mfile,*) Dcoordinates(3,:)
        ELSE 
          READ(mfile,*) Dtemp(:)
        END IF
        
        ! Add the polygon
        CALL ucd_addPolygon (rexport,Dcoordinates,imaterial)
        
        ! Deallocate data
        DEALLOCATE (Dtemp,Dcoordinates)
      
        ! Read the next line specifying basic data about the next polygon
        CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
        CALL sys_tolower (sline)
        
        READ(sline,*) skey
        bfinish = (skey .NE. 'endtrace') 
      END DO
      
      IF (ios .NE. 0) THEN
        CALL output_line ('Error reading GMV data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_polygondata')
        CALL sys_halt()
      END IF
      
    END SUBROUTINE
    
    ! ---------------------------------------------------------------
    
    SUBROUTINE read_variables (mfile,scommand,rexport)
    
    ! Reads data about variables from the GMV file mfile.
    
    ! Handle to the GMV file
    INTEGER, INTENT(IN) :: mfile
    
    ! Last command line in the GMV file 
    CHARACTER(LEN=*), INTENT(IN) :: scommand
    
    ! UCD structure where data is saved to.
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
      ! Local variables
      INTEGER :: ilinelen,ios,itype,imaterial
      LOGICAL :: bfinish
      CHARACTER(LEN=SYS_STRLEN) :: skey,sline,sname
      REAL(DP), DIMENSION(:), ALLOCATABLE :: Dcell
      REAL(DP), DIMENSION(:), ALLOCATABLE :: Dnode

      ! Allocate memory for cell- and node-based data
      ALLOCATE(Dcell(rexport%p_rtriangulation%NEL))
      ALLOCATE(Dnode(rexport%p_rtriangulation%NVT))

      ! Read the next line specifying basic data about the variable
      CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
      CALL sys_tolower (sline)
      
      READ(sline,*) skey
      bfinish = (skey .EQ. 'endvars') 
      DO WHILE ((ios .EQ. 0) .AND. (.NOT. bfinish))
        
        ! Read variable name, type
        READ (sline,*) sname,itype
        
        SELECT CASE (itype)
        CASE (0) 
          ! Cell based data. Read and remember
          READ(mfile,*) Dcell
          
          CALL ucd_addVariableElementBased (rexport,sname,UCD_VAR_STANDARD,Dcell)
              
        CASE (1)
          ! Node/Vertex based data.
          READ(mfile,*) Dnode
          
          CALL ucd_addVariableVertexBased (rexport,sname,UCD_VAR_STANDARD,Dnode)
          
        CASE DEFAULT
          CALL output_line ('Error reading GMV data! Unknown variable type!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'read_variables')
          CALL sys_halt()
  
        END SELECT
        
        ! Read the next line specifying basic data about the next variable
        CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
        CALL sys_tolower (sline)
        
        READ(sline,*) skey
        bfinish = (skey .EQ. 'endvars') 
        
      END DO
      
      IF (ios .NE. 0) THEN
        CALL output_line ('Error reading GMV data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_variables')
        CALL sys_halt()
      END IF
      
    END SUBROUTINE

    ! ---------------------------------------------------------------
    
    SUBROUTINE read_materials (mfile,scommand,rexport)
    
    ! Reads data about materials from the GMV file mfile.
    
    ! Handle to the GMV file
    INTEGER, INTENT(IN) :: mfile
    
    ! Last read command line in the GMV file 
    CHARACTER(LEN=*), INTENT(IN) :: scommand
    
    ! UCD structure where tracer data is saved to.
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
      ! Local variables
      INTEGER :: ilinelen,ios,nmats,itype,i
      CHARACTER(LEN=SYS_STRLEN) :: skey,sline
      CHARACTER(LEN=32), DIMENSION(:), ALLOCATABLE :: Snames
      INTEGER, DIMENSION(:), ALLOCATABLE :: Imat

      ! Interpret the command line, how many materials do we have? Type?
      READ (scommand,*) skey,nmats,itype
      
      ! Allocate memory for the variable names
      ALLOCATE(Snames(nmats))
      
      ! Read the material names
      DO i=1,nmats
        READ(mfile,*) Snames(i)
      END DO
      
      ! Allocate memory for vertex/element material classification
      SELECT CASE (itype)
      CASE (0) 
        ! Cell based data. 
        ALLOCATE(Imat(rexport%p_rtriangulation%NEL))
            
      CASE (1)
        ! Node/Vertex based data.
        ALLOCATE(Imat(rexport%p_rtriangulation%NVT))
        
      CASE DEFAULT
        CALL output_line ('Error reading GMV data! Unknown variable type!', &
                OU_CLASS_ERROR,OU_MODE_STD,'read_materials')
        CALL sys_halt()

      END SELECT
      
      ! Read material data of each cell / vertex
      READ(mfile,*) Imat(:)
      
      ! GMV supports only the same material names for both, cells and
      ! vertices. Inform the rexport structure about the material names:
      CALL ucd_setMaterials (rexport,Snames)

      ! Specify the material for vertices / elements
      SELECT CASE (itype)
      CASE (0) 
        ! Cell based data. 
        CALL ucd_setVertexMaterial (rexport,Imat)
            
      CASE (1)
        ! Node/Vertex based data.
        CALL ucd_setVertexMaterial (rexport,Imat)
        
      END SELECT

      ! Deallocate memory      
      DEALLOCATE(Imat,Snames)

    END SUBROUTINE
    
    ! ---------------------------------------------------------------
    
    SUBROUTINE read_velocity(mfile,scommand,rexport)
    
    ! Reads data about velocity from the GMV file mfile.
    
    ! Handle to the GMV file
    INTEGER, INTENT(IN) :: mfile
    
    ! Last command line in the GMV file 
    CHARACTER(LEN=*), INTENT(IN) :: scommand
    
    ! UCD structure where data is saved to.
    TYPE(t_ucdExport), INTENT(INOUT) :: rexport
    
      ! Local variables
      INTEGER :: ilinelen,ios,itype
      LOGICAL :: bfinish
      CHARACTER(LEN=SYS_STRLEN) :: skey,sline,sname
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Ddata
      REAL(DP), DIMENSION(:), ALLOCATABLE :: Dtemp

      ! Type of velocity data? Vertex or element based?
      READ(scommand,*) skey,itype
      
      SELECT CASE (itype)
      CASE (0) 
            
        ! Node based velocity data. Allocate memory.
        ALLOCATE(Ddata(rexport%p_rtriangulation%NEL,rexport%p_rtriangulation%ndim))
        ALLOCATE(Dtemp(rexport%p_rtriangulation%NEL))
        
        ! Read the data. X-velocity
        READ(mfile,*) Ddata(:,1)
        
        ! Y-velocity
        READ(mfile,*) Ddata(:,2)
        
        ! Probably Z-velocity.
        ! Put the data to the rexport structure
        IF (rexport%p_rtriangulation%ndim .EQ. NDIM3D) THEN
          READ(mfile,*) Ddata(:,3)
          CALL ucd_addVarElemBasedVec (rexport, sname, Ddata(:,1), Ddata(:,NDIM3D))
        ELSE 
          READ(mfile,*) Dtemp(:)
          CALL ucd_addVarElemBasedVec (rexport, sname, Ddata(:,1), Ddata(:,NDIM2D))
        END IF
        
        ! Deallocate memory
        DEALLOCATE(Dtemp,Ddata)
            
      CASE (1)
        
        ! Node based velocity data. Allocate memory.
        ALLOCATE(Ddata(rexport%p_rtriangulation%NVT,rexport%p_rtriangulation%ndim))
        ALLOCATE(Dtemp(rexport%p_rtriangulation%NVT))
        
        ! Read the data. X-velocity
        READ(mfile,*) Ddata(:,1)
        
        ! Y-velocity
        READ(mfile,*) Ddata(:,2)
        
        ! Probably Z-velocity.
        ! Put the data to the rexport structure
        IF (rexport%p_rtriangulation%ndim .EQ. NDIM3D) THEN
          READ(mfile,*) Ddata(:,3)
          CALL ucd_addVarVertBasedVec (rexport, sname, Ddata(:,1), Ddata(:,NDIM3D))
        ELSE 
          READ(mfile,*) Dtemp(:)
          CALL ucd_addVarVertBasedVec (rexport, sname, Ddata(:,1), Ddata(:,NDIM2D))
        END IF
        
        ! Deallocate memory
        DEALLOCATE(Dtemp,Ddata)
        
      CASE DEFAULT
        CALL output_line ('Error reading GMV data! Unknown variable type!', &
                OU_CLASS_ERROR,OU_MODE_STD,'read_velocity')
        CALL sys_halt()

      END SELECT
      

      ! Read the next line specifying basic data about the variable
      CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
      CALL sys_tolower (sline)
      
      READ(sline,*) skey
      bfinish = (skey .NE. 'endvars') 
      DO WHILE ((ios .EQ. 0) .AND. (.NOT. bfinish))
        
        ! Read variable name, type
        WRITE (mfile,*) sname,itype
        
        
        ! Read the next line specifying basic data about the next variable
        CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
        CALL sys_tolower (sline)
        
        READ(sline,*) skey
        bfinish = (skey .NE. 'endvars') 
        
      END DO
      
      IF (ios .NE. 0) THEN
        CALL output_line ('Error reading GMV data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_variables')
        CALL sys_halt()
      END IF
      
    END SUBROUTINE

    ! ---------------------------------------------------------------
    
    SUBROUTINE read_triangulation (mfile,scommand,rtriangulation)
    
    ! Reads data about tracers from the GMV file mfile.
    
    ! Handle to the GMV file
    INTEGER, INTENT(IN) :: mfile
    
    ! Last read command line in the GMV file
    CHARACTER(LEN=*), INTENT(IN) :: scommand
    
    ! Triangulation structure. If this contains a valid triangulation,
    ! the mesh data in the GMV file is skipped.
    ! If rtriangulation is empty, a new triangulation with mesh data
    ! from the GMV file is set up.
    TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
    
      ! Local variables
      INTEGER :: ilinelen,ios,ntracers,n,i,nve,ive
      INTEGER(I32) :: nmaxnve
      CHARACTER(LEN=SYS_STRLEN) :: skey,sdummy,sline
      REAL(DP), DIMENSION(:), ALLOCATABLE :: Dx,Dy,Dz
      INTEGER(I32), DIMENSION(:,:), ALLOCATABLE :: IverticesAtElement
      REAL(DP), DIMENSION(:,:), POINTER :: p_Ddata2D
      INTEGER(I32), DIMENSION(:,:), POINTER :: p_Idata2D
      INTEGER(I32), DIMENSION(2) :: Isize

      ! Do we have a "nodes x" or a "nodes fromfile" command?
      READ(scommand,*) skey,sdummy
      
      IF ((rtriangulation%ndim .NE. 0) .AND. (skey .EQ. 'fromfile')) THEN
        CALL output_line ('Cross reference to mesh not supported!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        CALL sys_halt()
      END IF
      
      ! Quit if this is just this line that contains triangulation info
      IF ((rtriangulation%ndim .NE. 0) .AND. (skey .EQ. 'fromfile')) RETURN
      
      ! Get the number of vertices
      READ(scommand,*) skey,n
      
      ! What is this for a command?
      IF (skey .NE. 'nodes') THEN
        CALL output_line ('Unsupported GMV file structure!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        CALL sys_halt()
      END IF
      
      rtriangulation%NVT = n
      
      ! Read point position data
      ALLOCATE(Dx(n),Dy(n),Dz(n))
      READ(mfile,*) Dx(:)
      READ(mfile,*) Dy(:)
      READ(mfile,*) Dz(:)
      
      ! 2D or 3D mesh? It's 3D as soon as one Z-coordinate is <> 0
      rtriangulation%ndim = NDIM2D
      DO i=1,n
        IF (Dz(i) .NE. 0.0_DP) THEN
          rtriangulation%ndim = NDIM3D
          EXIT
        END IF
      END DO
      
      ! Allocate memory for the coordinates with the storage-system
      Isize = (/rtriangulation%ndim,INT(rtriangulation%NVT,I32)/)
      CALL storage_new2D ('tria_read_tri2D', 'DCORVG', Isize, ST_DOUBLE, &
          rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      
      ! Get the pointers to the coordinate array
      CALL storage_getbase_double2D(&
          rtriangulation%h_DvertexCoords,p_Ddata2D)
      
      ! Copy the coordinates
      DO i=1,n
        p_Ddata2D(NDIM1D,i) = Dx(i)
      END DO
      IF (rtriangulation%ndim .GE. NDIM2D) THEN
        DO i=1,n
          p_Ddata2D(NDIM2D,i) = Dy(i)
        END DO
      END IF
      IF (rtriangulation%ndim .GE. NDIM3D) THEN
        DO i=1,n
          p_Ddata2D(NDIM3D,i) = Dz(i)
        END DO
      END IF
      
      ! Deallocate memory
      DEALLOCATE(Dz,Dy,Dx)
      
      ! At next we need KVERT = IverticesAtElement.
      CALL io_readlinefromfile (mfile, sline, ilinelen, ios)
      
      ! Do we have a "nodes x" or a "nodes fromfile" command?
      READ(sline,*) skey,sdummy
      
      IF (skey .EQ. 'fromfile') THEN
        CALL output_line ('Cross reference to mesh not supported!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        CALL sys_halt()
      END IF
      
      ! Get the number of cells
      READ(sline,*) skey,n
      
      ! What is this for a command?
      IF (skey .NE. 'cells') THEN
        CALL output_line ('Unsupported GMV file structure!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        CALL sys_halt()
      END IF
      
      rtriangulation%NEL = n
      
      rtriangulation%InelOfType(:) = 0
      
      ! Allocate memory for connectivity. We support up to TRIA_MAXNVE
      ! vertices per element and don't know a-priori
      ALLOCATE(IverticesAtElement(TRIA_MAXNVE,n))
      
      ! Read the next n cell data tags.
      ! Format: Line 1: "type" NVE. Line 2: Vertex numbers.
      nmaxnve = 0
      DO i=1,n
        IverticesAtElement(:,i) = 0
        READ(mfile,*) sdummy,nve
        
        ! Increase the number of elements of that type
        rtriangulation%InelOfType(nve) = rtriangulation%InelOfType(nve)+1
        
        ! Figure out the size of the first dimension of KVERT
        nmaxnve = MAX(nve,nmaxnve)
        
        ! Read the vertex numbers
        READ(mfile,*) IverticesAtElement(1:nve,i)
      END DO
      
      ! All elements read in. Create the actual IverticesAtElement.

      Isize = (/nmaxnve,INT(rtriangulation%NEL,I32)/)
      CALL storage_new2D ('tria_read_tri2D', 'KVERT', Isize, ST_INT, &
          rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
          
      ! Get the pointer to the IverticesAtElement array and read the array
      CALL storage_getbase_int2D(&
          rtriangulation%h_IverticesAtElement,p_Idata2D)
      
      ! Copy the data.
      DO i=1,n
        p_Idata2D(1:nmaxnve,i) = IverticesAtElement(1:nmaxnve,i)
      END DO
      
      ! We don't need IverticesAtElement anymore...
      DEALLOCATE(IverticesAtElement)
      
      SELECT CASE (rtriangulation%ndim)
      CASE (NDIM2D)
        ! Create some standard mesh information
        CALL tria_genElementsAtVertex2D    (rtriangulation)
        CALL tria_genNeighboursAtElement2D (rtriangulation)
        CALL tria_genEdgesAtElement2D      (rtriangulation)
        CALL tria_genElementsAtEdge2D      (rtriangulation)
        
        ! Reconstruct InodalProperty
        CALL reconstruct_InodalProperty_2D (rtriangulation)
        
        ! Now generate a standard mesh from the raw mesh.
        ! Generate all missing information.
        CALL tria_initStandardMeshFromRaw(rtriangulation)

      CASE DEFAULT
        ! For an (efficient) implementation of reconstructing InodalProperty
        ! in 3D, see [Jens Acker, Diploma Thesis].
        CALL output_line ('Only 2D supported! Cannot reconstruct InodalProperty!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        CALL sys_halt()
      END SELECT

    END SUBROUTINE
    
    SUBROUTINE reconstruct_InodalProperty_2D (rtriangulation)
    
      ! Reconstructs the InodalProperty, IboundaryCpIdx and 
      ! IverticesAtBoundary arrays. Sets NBCT!
      ! DvertexParameterValue is cleared.
      
      ! Triangulation structure. InodalProperty and NBCT are initialised here.
      ! InodalProperty must be allocated and initialised with 0.
      TYPE(t_triangulation), INTENT(INOUT) :: rtriangulation
    
      ! local variables
      INTEGER(I32), DIMENSION(:), POINTER :: p_InodalProperty
      INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
      INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
      INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
      INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IverticesAtBoundary
      INTEGER(PREC_VERTEXIDX), DIMENSION(:), ALLOCATABLE :: IverticesAtBoundary
      INTEGER(I32), DIMENSION(:), POINTER :: p_IboundaryCpIdx
      INTEGER(PREC_VERTEXIDX) :: ivt,ivt2
      INTEGER(PREC_ELEMENTIDX) :: iel,ielidx
      INTEGER :: ive,iveprevious,ivenext,nbct,nve,nvbd,ivbd,ibctidx,icurrentbc
      REAL(DP), DIMENSION(:), POINTER :: p_DvertexParameterValue
      
      ! Allocate memory for the arrays 
      CALL storage_new ('tria_read_tri2D', 'KNPR', &
          INT(rtriangulation%NVT,I32), ST_INT, &
          rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
      
      ! Get the pointer to some arrays
      CALL storage_getbase_int(&
          rtriangulation%h_InodalProperty,p_InodalProperty)
      CALL storage_getbase_int2d(&
          rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
      CALL storage_getbase_int2d(&
          rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
      CALL storage_getbase_int2d(&
          rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
      CALL storage_getbase_int(&
          rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
      CALL storage_getbase_int(&
          rtriangulation%h_IelementsAtVertexIdx,p_IelementsAtVertexIdx)
          
      ! In a first loop find the edges with no neighbour element. These are
      ! boundary edges.
      ! Count the number of boundary vertices.
      nvbd = 0
      DO iel=1,rtriangulation%NEL
        DO ive=1,UBOUND(p_IverticesAtElement,1)
          IF (p_IverticesAtElement(ive,iel) .EQ. 0) EXIT ! Triangle in a quad mesh
          IF (p_IneighboursAtElement (ive,iel) .EQ. 0) THEN
            
            ! Both vertices on the edge are boundary vertices
            ivenext = MOD(ive,UBOUND(p_IverticesAtElement,1))+1
            IF (p_IverticesAtElement(ivenext,iel) .EQ. 0) ivenext = ivenext-1
            
            IF (p_InodalProperty(p_IverticesAtElement(ive,iel)) .NE. -1) THEN
              p_InodalProperty(p_IverticesAtElement(ive,iel))     = -1
              nvbd = nvbd+1
            END IF
            IF (p_InodalProperty(p_IverticesAtElement(ivenext,iel)) .NE. -1) THEN
              p_InodalProperty(p_IverticesAtElement(ivenext,iel)) = -1
              nvbd = nvbd+1
            END IF
          END IF
        END DO
      END DO
      
      rtriangulation%NVBD = nvbd
      
      ! Allocate memory for IverticesAtBoundary.
      CALL storage_new ('tria_generateBasicBoundary', &
          'KVBD', INT(rtriangulation%NVBD,I32), &
          ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)

      CALL storage_getbase_int(&
          rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
      
      ! Now loop through all vertices. Whenever we find a vertex with InodalProperty=-1,
      ! we start going from vertex to vertex until we found all vertices of
      ! that boundary component.
      ! We save all boundary vertices in IverticesAtBoundary.
      nbct = 0
      ivbd = 0
      DO ivt = 1,rtriangulation%NVT
        IF (p_InodalProperty(ivt) .EQ. -1) THEN
          ! New boundary component found
          nbct = nbct+1
          
          ! Go from edge to edge until we hit ivt again. This must happen as the 
          ! boundary component must be closed.
          ivt2 = ivt
          
          ivt2loop: DO
            p_InodalProperty(ivt2) = nbct
            
            ! Save the vertex.
            ivbd = ivbd+1
            p_IverticesAtBoundary(ivbd) = ivt2
            
            ! Check the edges adjacent to the vertex...for this purpose, we must 
            ! find them.
            DO ielidx = p_IelementsAtVertexIdx(ivt2),p_IelementsAtVertexIdx(ivt2+1)-1

              iel = p_IelementsAtVertex(ielidx)
            
              ! Where is the point in the element
              DO ive=1,UBOUND(p_IverticesAtElement,1)
                IF (p_IverticesAtElement(ive,iel) .EQ. ivt2) EXIT 
              END DO
              
              ! Triangle or quad?
              nve = UBOUND(p_IverticesAtElement,1)
              IF (p_IverticesAtElement(nve,iel) .EQ. 0) nve = nve-1

              iveprevious = MOD(ive+nve-2,nve)+1
              ivenext = MOD(ive,nve)+1
              
              ! On that element, check the edge following the vertex if it's a boundary edge.
              IF (p_IneighboursAtElement(ive,iel) .EQ. 0) THEN
              
                ! Yes, that's a boundary edge. does it lead to a vertex that we haven't had?
                IF (p_InodalProperty(p_IverticesAtElement(ivenext,iel)) .EQ. -1) THEN
                  ! Yes! That node is now a boundary node. We continue with it.
                  ivt2 = p_IverticesAtElement(ivenext,iel)
                  CYCLE ivt2loop
                END IF
              
              END IF
                
              ! No, it's not. Check the previous edge
              IF (p_IneighboursAtElement(iveprevious,iel) .EQ. 0) THEN
              
                ! Yes, that's a boundary edge. does it lead to a vertex that we haven't had?
                IF (p_InodalProperty(p_IverticesAtElement(iveprevious,iel)) .EQ. -1) THEN
                  ! Yes! That node is now a boundary node. We continue with it.
                  ivt2 = p_IverticesAtElement(iveprevious,iel)
                  CYCLE ivt2loop
                END IF
              
              END IF
                
            END DO
          
            ! Ok, there is no edge found starting from the vertex ivt2 that leads to 
            ! another vertex on the boundary. That means we found all vertices
            ! on the boundary component!
            ! So we quit the loop here and continue to find the next non-classified
            ! vertex ivt.
            EXIT ivt2loop
            
          END DO ivt2loop
          
        END IF
      END DO
    
      ! nbct is the number of boundary components we found.
      rtriangulation%NBCT = nbct
    
      ! Allocate memory for the boundary component index vector.
      ! Initialise that with zero!
      CALL storage_new ('tria_generateBasicBoundary', &
          'KBCT', INT(rtriangulation%NBCT+1,I32), &
          ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
      CALL storage_getbase_int (&
          rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
      ! Figure out the IboundaryCpIdx array but once checking the BC
      ! number of all vertices on the boundary
      ibctidx = 0
      icurrentbc = 0
      DO ivbd = 1,nvbd
        IF (p_InodalProperty(p_IverticesAtBoundary(ivbd)) .NE. icurrentbc) THEN
          ibctidx = ibctidx+1
          p_IboundaryCpIdx(ibctidx) = ivbd
          icurrentbc = p_InodalProperty(p_IverticesAtBoundary(ivbd))
        END IF
      END DO
      
      p_IboundaryCpIdx(nbct+1) = nvbd+1
      
      ! Allocate memory for  and DvertexParameterValue
      CALL storage_new ('tria_generateBasicBoundary', &
          'DVBDP', INT(rtriangulation%NVBD,I32), &
          ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
      
      ! Get the array where to store boundary parameter values.
      CALL storage_getbase_double (&
          rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)

      ! Clear the array -- we have no parameter values!
      CALL lalg_clearVectorDble (p_DvertexParameterValue)
    
      ! InodalProperty is completely classified -- that's it.
    
    END SUBROUTINE

  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_setFilename(rexport,sfilename)

!<description>
  ! This routine allows to change the filename of an output file.
  ! Whether or not this is possible depends on the file format of the
  ! output file.
  ! Changing the filename is possible as long as no data is written
  ! to the file, e.g. for all file formats that buffer the output
  ! in memory before writing the file.
!</description>

!<input>
  ! New filename for output file.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
!</input>

!<inputoutput>
  ! The ucd export structure that specifies the output.
  TYPE(t_ucdExport), INTENT(INOUT) :: rexport
!</inputoutput>

!</subroutine>

    ! Start the output
    SELECT CASE(rexport%coutputFormat)
    CASE (UCD_FORMAT_GMV,UCD_FORMAT_AVS,UCD_FORMAT_VTK)
      ! Here it's ok to change the filename
      rexport%sfilename = sfilename
    CASE DEFAULT
      CALL output_line ('Cannot change filename', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setFilename')
      CALL sys_halt()
    END SELECT
    
  END SUBROUTINE

!************************************************************************
  
!<subroutine>

  SUBROUTINE ucd_getVariable(rexport,svarName,Ddata,nlength,itype)

!<description>
  ! Reads the data of the variable svarName. If Ddata is specified,
  ! the variable data is written to Ddata. If nlength is specified,
  ! the length of the variable is returned in nlength.
!</description>

!<input>
  ! The ucd export structure containing data.
  TYPE(t_ucdExport), INTENT(IN) :: rexport

  ! Name of the variable whose data should be retrieved.
  CHARACTER(LEN=*), INTENT(IN) :: svarName
!</input>

!<output>
  ! OPTIONAL: Data array. Must be large enough to hold all data.
  ! The data of variable svarName is transferred to Ddata.
  !
  ! If the variable is unknown, Ddata is not changed.
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Ddata
  
  ! OPTIONAL: Length qualifier.
  ! If specified, nlength is set to the number of elements in
  ! the variable svarName, thus to the minimum length Ddata must
  ! have to accept the data of the variable.
  !
  ! If the variable is unknown, -1 is returned.
  INTEGER(I32), INTENT(OUT), OPTIONAL :: nlength
  
  ! OPTIONAL: Type of data. One of the UCD_BASE_xxxx flags.
  ! =UCD_BASE_ELEMENT: element based data.
  ! =UCD_BASE_VERTEX: vertex based data.
  !
  ! Is set to -1 if the variable is unknown.
  INTEGER, INTENT(OUT), OPTIONAL :: itype
!</output>

!</subroutine>

    INTEGER :: i
    CHARACTER(LEN=SYS_NAMELEN) :: sname
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    sname = svarName
    CALL sys_toupper_replace (sname)
    
    ! Find the variable
    DO i=1,SIZE(rexport%p_SvariableNames)
    
      IF (sys_upcase(rexport%p_SvariableNames(i)) .EQ. sname) THEN
        ! Variable found! Return data as desired.
        
        CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
        
        IF (PRESENT(nlength)) nlength = SIZE(p_Ddata)
        IF (PRESENT(itype)) itype = rexport%p_IvariableBase(i)
        IF (PRESENT(Ddata)) THEN
          CALL lalg_copyVectorDble (p_Ddata,Ddata(1:SIZE(p_Ddata)))
        END IF
        
      END IF
    
    END DO
    
    ! Variable not found.
    IF (PRESENT(nlength)) nlength = -1
    IF (PRESENT(itype)) itype = -1
    
  END SUBROUTINE

END MODULE

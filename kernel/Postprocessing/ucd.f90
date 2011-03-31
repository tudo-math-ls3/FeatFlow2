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
!#      -> Allows to specify material names for the material ID`s that
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
!# 23.) ucd_infoVariables
!#      -> Output information about vairables.
!#
!# 24.) ucd_getSimulationTime
!#      -> Get simulation time.
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

module ucd

  use fsystem
  use storage
  use genoutput
  use basicgeometry
  use linearalgebra
  use triangulation
  use paramlist
  use io
  use geometry 
  implicit none
  
  private

!<constants>

!<constantblock description="Output format types">

  ! No output format
  integer, parameter, public :: UCD_FORMAT_NONE = 0
  
  ! GMV file format (ASCII)
  integer, parameter, public :: UCD_FORMAT_GMV  = 1
  
  ! AVS/Express file format (ASCII)
  integer, parameter, public :: UCD_FORMAT_AVS  = 2
  
  ! Visualization Toolkit (VTK) file format (ASCII)
  integer, parameter, public :: UCD_FORMAT_VTK  = 3

  ! GMV file format (binary)
  integer, parameter, public :: UCD_FORMAT_BGMV = 4
  
!</constantblock>

!<constantblock description="Data type flags">

  ! Element-based data
  integer, parameter, public :: UCD_BASE_ELEMENT = 0
  
  ! Vertex-based data
  integer, parameter, public :: UCD_BASE_VERTEX  = 1
  
!</constantblock>

!<constantblock description="Parameter constants for the GMV exporter">

  ! Export vector components as scalars
  integer, parameter, public :: UCD_PARAM_GMV_VECTOR_TO_SCALAR = 2**0

!</constantblock>

!<constantblock description="Parameter constants for the VTK exporter">

  ! Export vector components as scalars
  integer, parameter, public :: UCD_PARAM_VTK_VECTOR_TO_SCALAR = 2**0
  
  ! Use quadratic cell elements when possible
  integer, parameter, public :: UCD_PARAM_VTK_USE_QUADRATIC    = 2**1

!</constantblock>

!<constantblock description="Constants for the cflags variable in ucd_startXXXX. Bitfield.">

  ! Standard flags. Write information in corner vertices only, linear interpolation.
  integer(I32), parameter, public :: UCD_FLAG_STANDARD            = 0
  
  ! Construct edge midpoints and write them as vertices to the output file.
  integer(I32), parameter, public :: UCD_FLAG_USEEDGEMIDPOINTS    = 2**0
  
  ! Construct element midpoints and write them as vertices to the output file.
  integer(I32), parameter, public :: UCD_FLAG_USEELEMENTMIDPOINTS = 2**1
  
  ! Quadratic buld interpolation. Information is given in corner vertices and
  ! edge midpoints. Implies UCD_FLAG_USEEDGEMIDPOINTS.
  integer(I32), parameter, public :: UCD_FLAG_BULBQUADRATIC       = 2**2
  
  ! Output of a linear interpolated solution on a once refined mesh. 
  ! Implies UCD_FLAG_USEEDGEMIDPOINTS and UCD_FLAG_USEELEMENTMIDPOINTS. 
  ! Cannot be used with UCD_FLAG_BULBQUADRATIC.
  integer(I32), parameter, public :: UCD_FLAG_ONCEREFINED         = 2**3

  ! Can be specified additionally to UCD_FLAG_USEEDGEMIDPOINTS and/or
  ! UCD_FLAG_USEELEMENTMIDPOINTS. Prevents a warning if missing nodes are
  ! not specified.
  integer(I32), parameter, public :: UCD_FLAG_IGNOREDEADNODES = 2**4

!</constantblock>

!<constantblock description="Specification flags for variables. Bitfield.">
  ! Standard 'scalar' variable.
  integer(I32), parameter, public :: UCD_VAR_STANDARD             = 0

  ! The variable specifies a velocity component. For some output formats
  ! the velocity field may require some special treatment.
  integer(I32), parameter, public :: UCD_VAR_VELOCITY             = 2**0

  ! The variable specifies the X-component of a vector field.
  ! Cannot be used together with UCD_VAR_YVECTORCOMP or UCD_VAR_ZVECTORCOMP
  integer(I32), parameter, public :: UCD_VAR_XVECTORCOMP          = 2**1
  
  ! The variable specifies the Y-component of a vector field.
  ! Cannot be used together with UCD_VAR_XVECTORCOMP or UCD_VAR_ZVECTORCOMP
  integer(I32), parameter, public :: UCD_VAR_YVECTORCOMP          = 2**2
  
  ! The variable specifies the Z-component of a vector field.
  ! Cannot be used together with UCD_VAR_XVECTORCOMP or UCD_VAR_YVECTORCOMP
  integer(I32), parameter, public :: UCD_VAR_ZVECTORCOMP          = 2**3

  ! The variable specifies the X velocity of a velocity field.
  ! Cannot be used together with UCD_VAR_YVELOCITY or UCD_VAR_ZVELOCITY.
  integer(I32), parameter, public :: UCD_VAR_XVELOCITY            = UCD_VAR_XVECTORCOMP&
                                                                  + UCD_VAR_VELOCITY

  ! The variable specifies the Y velocity of a velocity field.
  ! Cannot be used together with UCD_VAR_XVELOCITY or UCD_VAR_ZVELOCITY.
  integer(I32), parameter, public :: UCD_VAR_YVELOCITY            = UCD_VAR_YVECTORCOMP&
                                                                  + UCD_VAR_VELOCITY

  ! The variable specifies the Z velocity of a velocity field.
  ! Cannot be used together with UCD_VAR_XVELOCITY or UCD_VAR_YVELOCITY.
  integer(I32), parameter, public :: UCD_VAR_ZVELOCITY            = UCD_VAR_ZVECTORCOMP&
                                                                  + UCD_VAR_VELOCITY
  
!</constantblock>

!<constantblock description="Constants for specifying alternative source files. Bitfield.">

  ! Whole triangulation is specified by a different source file.
  integer(I32), parameter, public :: UCD_ASRC_ALL                 = not(0)

  ! Coordinates of mesh points are specified in a different source file.
  integer(I32), parameter, public :: UCD_ASRC_POINTS              = 2**0

  ! Cell information is specified in a different source file.
  integer(I32), parameter, public :: UCD_ASRC_CELLS               = 2**1

  ! Material information is specified in a different source file.
  integer(I32), parameter, public :: UCD_ASRC_MATERIALS           = 2**2

  ! Polygons are specified in a difference source file.
  integer(I32), parameter, public :: UCD_ASRC_POLYGONS            = 2**3

!</constantblock>

!<constantblock description="Cell type constants for VTK exporter">

  ! A single vertex
  integer, parameter, public :: VTK_VERTEX                         =  1
  
  ! A set of vertices
  integer, parameter, public :: VTK_POLY_VERTEX                    =  2
  
  ! A line
  integer, parameter, public :: VTK_LINE                           =  3
  
  ! A line strip
  integer, parameter, public :: VTK_POLY_LINE                      =  4
  
  ! A triangle
  integer, parameter, public :: VTK_TRIANGLE                       =  5
  
  ! A triangle strip
  integer, parameter, public :: VTK_TRIANGLE_STRIP                 =  6
  
  ! A polygon
  integer, parameter, public :: VTK_POLYGON                        =  7
  
  ! A pixel
  integer, parameter, public :: VTK_PIXEL                          =  8
  
  ! A quadrilateral
  integer, parameter, public :: VTK_QUAD                           =  9
  
  ! A tetrahedron
  integer, parameter, public :: VTK_TETRA                          = 10
  
  ! A voxel (cube)
  integer, parameter, public :: VTK_VOXEL                          = 11
  
  ! A hexahedron
  integer, parameter, public :: VTK_HEXAHEDRON                     = 12
  
  ! A wedge
  integer, parameter, public :: VTK_WEDGE                          = 13
  
  ! A pyramid
  integer, parameter, public :: VTK_PYRAMID                        = 14
  
  ! A quadratic edge
  integer, parameter, public :: VTK_QUADRATIC_EDGE                 = 21
  
  ! A quadratic triangle
  integer, parameter, public :: VTK_QUADRATIC_TRIANGLE             = 22
  
  ! A quadratic quadrilateral
  integer, parameter, public :: VTK_QUADRATIC_QUAD                 = 23
  
  ! A quadratic tetrahedron
  integer, parameter, public :: VTK_QUADRATIC_TETRA                = 24
  
  ! A quadratic hexahedron
  integer, parameter, public :: VTK_QUADRATIC_HEXAHEDRON           = 25

!</constantblock>
  
!</constants>

  
!<types>

!<typeblock>
  
  ! UCD export structure. The structure is created by one of the ucd_startXXXX
  ! routines and configured to write output of type XXXX (e.g. GMV or AVS).
  ! After writing out the data, the structure can be released with ucd_release.
  
  type t_ucdExport
    
    ! Output format. One of the UCD_FORMAT_XXXX constants.
    integer :: coutputFormat = UCD_FORMAT_NONE
    
    ! Parameters for output routines. A combination of UCD_PARAM_XXX_YYYY,
    ! where XXX specifies the output format.
    integer :: cparam = 0
    
    ! Name of the output file
    character(LEN=SYS_STRLEN) :: sfilename = ""

    ! For VTK polygons are written into a special file
    character(LEN=SYS_STRLEN) :: sfilepolyvtk = ""
    
    ! IO channel of the file
    integer :: iunit = 0
    
    ! Export flags specifying the output
    integer(I32) :: cflags            = 0
    
    ! Number of currently attached variables
    integer :: nvariables        = 0
    
    ! Number of curently attached polygons
    integer :: npolygons         = 0

    ! Number of curently attached surface triangulations
    integer :: nsurfTri          = 0
    
    ! Number of currently attached tracer data fields
    integer :: ntracerVariables  = 0
    
    ! Number of currently attached variable vectors
    integer :: nvectors = 0

    ! The simulation time. =SYS_INFINITY if no simulation time is specified
    real(DP) :: dsimulationTime  = SYS_INFINITY
    
    ! Format of the simulation time. Fortran format string.
    character(LEN=SYS_STRLEN) :: ssimTimeFormat = "(ES18.8E3)"
    
    ! Format of the output of double-precision numbers. Fortran format string.
    character(LEN=SYS_STRLEN) :: sdataFormat = "(ES18.8E3)"
    
    ! An array containing the names of all the variables
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_SvariableNames => null()
    
    ! An array containing the names of all the variable vectors
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_SvarVecNames => null()
    
    ! Filename of file containing point coordinates.
    ! ""=no alternative source file.
    character(LEN=SYS_STRLEN) :: saltFilePoints = ""

    ! Filename of file containing cell structure. 
    ! ""=no alternative source file.
    character(LEN=SYS_STRLEN) :: saltFileCells = ""

    ! Filename of file containing material structure. 
    ! ""=no alternative source file.
    character(LEN=SYS_STRLEN) :: saltFileMaterials = ""

    ! Filename of file containing polygons.
    ! ""=no alternative source file.
    character(LEN=SYS_STRLEN) :: saltFilePolygons = ""
    
    ! A pointer to the underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation => null()
    
    ! A pointer to an array with specification flags. p_IvariableSpec(I)
    ! is a bitfield for variable I that specifies the type of the variable
    ! and how to handle it.
    integer(I32), dimension(:), pointer :: p_IvariableSpec => null()
    
    ! A pointer to an array that specifies whether a variable is vertex based (1)
    ! or cell based (0).
    integer, dimension(:), pointer :: p_IvariableBase => null()
    
    ! A pointer to an array that specifies the components of a vector variable.
    ! The first dimension of this vector is always 4.
    integer, dimension(:,:), pointer :: p_Ivectors => null()
    
    ! A pointer to a list of handles of double precision pointers. 
    ! p_Hvariables(I) points to the data of variable I.
    integer, dimension(:), pointer :: p_Hvariables => null()
    
    ! A pointer to a list of handles to polygon data.
    ! p_Hpolygon(I) points to a list of (X,Y) or (X,Y,Z) tags
    ! containing the points of a polygon of line segments.
    integer, dimension(:), pointer :: p_Hpolygons => null()

    ! A pointer to a list of handles to surface triangulation vertices.
    ! p_HsurfTris(I) points to a list of (X,Y,Z) tags
    ! containing the points of a the surface triangulation
    integer, dimension(:), pointer :: p_HsurfTris => null()

    ! A pointer to a list of handles to surface triangulation data.
    ! p_HTriangles(I) is a handle an integer array that describes
    ! the connectivity of the vertices in p_HsurfTris(I)
    integer, dimension(:), pointer :: p_HTriangles => null()
    
    ! A pointer to a list of handles to surface triangulation data.
    ! p_HsurfData(I) is a handle an integer array that describes
    ! how many vertices and triangles there are in p_HsurfTris(I)
    integer, dimension(:), pointer :: p_HsurfData => null()
    
    ! A pointer to a list of material identifier for polygonal data.
    ! Element i in this list specifies the material of the polygon.
    integer :: hpolygonMaterial = ST_NOHANDLE
    
    ! A handle to a list of (X,Y) or (X,Y,Z) coordinates of tracers.
    integer :: htracers = ST_NOHANDLE
    
    ! An array containing the names of all tracer variables
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_StracerVariableNames => null()
    
    ! A list of handles to tracer data. Each handle identifies an
    ! "array[1..#tracers] of double", which specifies data for each
    ! tracer. p_StracerVariableNames[i] os the name of the i-th array.
    integer, dimension(:), pointer :: p_HtracerVariables => null()
    
    ! A handle to an array containing for every cell a cell material id.
    ! If not specified, every cell gets a default material id.
    integer :: hIcellMaterial = ST_NOHANDLE

    ! A handle to an array containing for every node/vertex (corners, 
    ! midpoints,...) a material id.
    ! If not specified, every vertex/node gets a default material id.
    integer :: hIvertexMaterial = ST_NOHANDLE
    
    ! Pointer to material names for vertex materials
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_SvertexMaterials => null()

    ! Pointer to material names for cell materials
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_ScellMaterials => null()
    
    ! Current length of the comment buffer
    integer :: ncommentBufSize = 0
    
    ! Pointer to a buffer with comment lines.
    ! This is a character array. the lines are separated by the NEWLINE
    ! character.
    character, dimension(:), pointer :: p_Scomments => null()
    
    ! Status variable: Number of vertices containing data in
    ! vertex based data arrays. All vertex-based variables in p_Hvariables
    ! have this length.
    integer :: nvertices = 0
    
    ! Status variable: Number of cells containing data in
    ! cell based data arrays. All cell-based variables in p_Hvariables
    ! have this length.
    integer :: ncells    = 0

    ! Status variable: Number of tracers.
    integer :: ntracers  = 0
    
  end type
  
  public :: t_ucdExport
  
!</typeblock>

!<typeblock>

  ! This structure is used for mesh refinement which is needed by some
  ! export routines. This structure is not meant to be used outside this
  ! module.

  type t_ucdRefine
  
    ! Number of additional vertices
    integer :: nvertices = 0
    
    ! Number of additional cells
    integer :: ncells = 0
    
    ! Handle to additional vertice array
    integer :: h_DvertexCoords = ST_NOHANDLE
    
    ! Handle to additional cell array
    integer :: h_IverticesAtElement = ST_NOHANDLE

  end type
  
  public :: t_ucdRefine

!</typeblock>

!</types>

  public :: ucd_startGMV
  public :: ucd_startBGMV
  public :: ucd_startAVS
  public :: ucd_startVTK
  public :: ucd_setAlternativeSource
  public :: ucd_addVariableVertexBased, ucd_addVariableElementBased
  public :: ucd_addVarVertBasedVec, ucd_addVarElemBasedVec
  public :: ucd_setVertexMaterial
  public :: ucd_setCellMaterial
  public :: ucd_setMaterials
  public :: ucd_setTracers
  public :: ucd_addTracerVariable
  public :: ucd_removeTracers
  public :: ucd_setSimulationTime
  public :: ucd_addPolygon
  public :: ucd_addCommentLine
  public :: ucd_addParameterList
  public :: ucd_setOutputNumberFormat
  public :: ucd_write
  public :: ucd_release
  public :: ucd_setFilename
  public :: ucd_readGMV
  public :: ucd_getVariable
  public :: ucd_infoVariables
  public :: ucd_getSimulationTime
  public :: ucd_addSurfTri

  interface ucd_addVariableVertexBased
    module procedure ucd_addVariableVertexBased1
    module procedure ucd_addVariableVertexBased2
  end interface

  interface ucd_addVariableElementBased
    module procedure ucd_addVariableElementBased1
    module procedure ucd_addVariableElementBased2
  end interface

  interface ucd_addVarVertBasedVec
    module procedure ucd_addVarVertBasedVec1
    module procedure ucd_addVarVertBasedVec2
  end interface

  interface ucd_addVarElemBasedVec
    module procedure ucd_addVarElemBasedVec1
    module procedure ucd_addVarElemBasedVec2
  end interface

contains

  !************************************************************************

!<subroutine>

  subroutine ucd_refine(rrefine, rtria, cflags)
  
!<description>
  ! This routine calculates the mesh refinement of a given triangulation
  ! and a combination of UCD_FLAG_XXXX flags.
  ! This routine is not meant to be called from outside this module.
!</description>

!<input>
  ! Specification of the underlying triangulation.
  type(t_triangulation), intent(in) :: rtria
  
  ! Bitfield that specifies the output.
  integer(I32), intent(in) :: cflags
!</input>

!<output>
  ! A mesh refinement
  type(t_ucdRefine), intent(out) :: rrefine
!</output>

!</subroutine>

    ! Some local variables
    integer :: i,j,k,ivt,off,numNewElem
    real(DP) :: dx
    real(DP), dimension(:,:), pointer :: p_DvertexCoords, p_DnewVerts
    integer, dimension(:,:), pointer :: p_IvertsAtEdge
    integer, dimension(:,:), pointer :: p_InewVertsAtElement
    integer, dimension(:,:), pointer :: p_IvertsAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(2) :: I_dim
    logical :: bedgeMids, belemMids, brefined
    
    ! Initialize vertice / cell counts
    rrefine%nvertices = 0
    rrefine%ncells = 0
    
    ! Do we need to write edge / element midpoints?
    brefined  = (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)
    belemMids = (iand(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. brefined
    bedgeMids = (iand(cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
                (iand(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. brefined
    
    ! If we use element midpoints or if the mesh is refines, we need to add
    ! the number of elements
    if (belemMids) then
        
      rrefine%nvertices = rtria%NEL
      
    end if
    
    ! Do we need to add edge midpoints?
    if (bedgeMids) then
        
        rrefine%nvertices = rrefine%nvertices + rtria%NMT
    
    end if
    
    ! If we do not allocate any vertices, we can leave this routine here
    if (rrefine%nvertices .eq. 0) return
    
    ! Otherwise allocate the vertices
    if (rtria%ndim .eq. NDIM2D) then
      I_dim(1) = 2
    else
      I_dim(1) = 3
    end if
    I_dim(2) = rrefine%nvertices
    
    call storage_new("ucd_refine", "p_Dvertices", I_dim, ST_DOUBLE, &
        rrefine%h_DvertexCoords, ST_NEWBLOCK_ZERO)

    ! And get a pointer to them
    call storage_getbase_double2D(rrefine%h_DvertexCoords, p_DnewVerts)
    
    ! Get pointers to the triangulation`s arrays
    call storage_getbase_double2D(rtria%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2D(rtria%h_IverticesAtEdge, p_IvertsAtEdge)
    call storage_getbase_int2D(rtria%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2D(rtria%h_IverticesAtElement, p_IvertsAtElement)
    
    ! Calculate the vertices
    off = 0
    if (bedgeMids) then
      do i = 1, rtria%NMT
        do j = 1, ubound(p_DvertexCoords, 1)
          p_DnewVerts(j, i) = 0.5_DP * &
                             (p_DvertexCoords(j, p_IvertsAtEdge(1, i)) + &
                              p_DvertexCoords(j, p_IvertsAtEdge(2, i)))
        end do
      end do
      ! remember we already wrote some vertices
      off = rtria%NMT
    end if
    
    ! Calculate element midpoints?
    if (belemMids) then
    
      ! Go through all elements...
      do i = 1, rtria%NEL
      
        ! ...and all coordinates...
        do j = 1, ubound(p_DvertexCoords, 1)
        
          dx = 0.0_DP
          
          ! ...and through every vertex of the element
          do k = 1, ubound(p_IvertsAtElement, 1)
          
            ivt = p_IvertsAtElement(k, i)
          
            ! Does this element only have k-1 vertices?
            if (ivt .eq. 0) exit
          
            dx = dx + p_DvertexCoords(j, ivt)
            
          end do
          
          ! Store element midpoint
          p_DnewVerts(j, off+i) = dx / real(k-1, DP)
          
        end do
        
      end do
      
    end if
    
    ! Do we also need to refine the elements?
    ! If not, then we can return here
    if (.not. brefined) return
    
    ! We first need to count the number of elements we will create for the
    ! refinement...
    ! Since we work in 2D here, every refinement takes 4 new elements...
    ! TODO: This code may need to be replaced for 3D grids...
    numNewElem = rtria%NEL * 4
    
    ! Allocate elements
    I_dim(1) = 4
    I_dim(2) = numNewElem
    call storage_new("ucd_refine", "p_DnewVertsAtElement", I_dim, &
        ST_INT, rrefine%h_IverticesAtElement, ST_NEWBLOCK_ZERO)

    ! And get a pointer to them
    call storage_getbase_int2D(rrefine%h_IverticesAtElement, p_InewVertsAtElement)

    ! Now go through all elements
    do i = 1, rtria%NEL
    
      ! Count the number of vertices for this element
      do k = 1, ubound(p_IvertsAtElement, 1)
        if (p_IvertsAtElement(k, i) .eq. 0) exit
      end do
      
      ! Now this element has k-1 vertices
      k = k - 1
      
      select case(k)
      case (3)
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

      case (4)
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
        
        ! Get the element`s mid-point offset
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
        
      end select
    
    end do
    
    ! That is it
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ucd_startGMV (rexport,cflags,rtriangulation,sfilename)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the GMV file
  character(LEN=*), intent(in) :: sfilename
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  integer(I32), intent(in) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  type(t_triangulation), intent(in), target :: rtriangulation
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  type(t_ucdExport), intent(out) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(out) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_GMV
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    if ((iand(cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    end if

    if ((iand(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      if (rtriangulation%NDIM .eq. NDIM2D) then
        rexport%ncells = rtriangulation%NEL*4
      else
        rexport%ncells = rtriangulation%NEL*8
      end if
    end if
  
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ucd_startBGMV (rexport,cflags,rtriangulation,sfilename)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the GMV file
  character(LEN=*), intent(in) :: sfilename
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  integer(I32), intent(in) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  type(t_triangulation), intent(in), target :: rtriangulation
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  type(t_ucdExport), intent(out) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(out) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_BGMV
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    if ((iand(cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    end if

    if ((iand(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      if (rtriangulation%NDIM .eq. NDIM2D) then
        rexport%ncells = rtriangulation%NEL*4
      else
        rexport%ncells = rtriangulation%NEL*8
      end if
    end if
  
  end subroutine
  
  !************************************************************************

!<subroutine>

  subroutine ucd_startAVS (rexport,cflags,rtriangulation,sfilename)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the AVS file
  character(LEN=*), intent(in) :: sfilename
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  integer(I32), intent(in) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  type(t_triangulation), intent(in), target :: rtriangulation
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  type(t_ucdExport), intent(out) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(out) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_AVS
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    rexport%p_Scomments => null()
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    if ((iand(cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    end if

    if ((iand(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      if (rtriangulation%NDIM .eq. NDIM2D) then
        rexport%ncells = rtriangulation%NEL*4
      else
        rexport%ncells = rtriangulation%NEL*8
      end if
    end if
  
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ucd_startVTK (rexport,cflags,rtriangulation,sfilename,sfilepoly,cparam)
 
!<description>
  ! Initialises the UCD output to a file sfilename. A UCD export structure
  ! rexport is created that specifies that file. This structure must be
  ! passed to all UCD output routines.
!</description>
 
!<input>
  ! Filename of the VTK file
  character(LEN=*), intent(in) :: sfilename

  ! Filename of the VTKpoly file
  character(LEN=*), optional, intent(in) :: sfilepoly
  
  ! Bitfield that specifies the output. Standard value is UCD_FLAG_STANDARD.
  integer(I32), intent(in) :: cflags
  
  ! Specification of the underlying triangulation. A pointer to this
  ! object is saved until the otput is finished.
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! OPTIONAL: Parameters for the VTK exporter
  integer, optional, intent(in) :: cparam
!</input>
  
!<output>
  ! An UCD export structure which collects information about the output.
  ! Must be passed to all export subroutines.
  type(t_ucdExport), intent(out) :: rexport
!</output>
 
!</subroutine>

    ! Most of the things in rexport is initialised by INTENT(out) with standard
    ! values automatically. We only have to initialise minor things.
    
    rexport%coutputFormat = UCD_FORMAT_VTK
    rexport%cflags = cflags
    rexport%sfilename = sfilename
    rexport%p_rtriangulation => rtriangulation
    
    ! Do we have a poly outputfile?
    if (present(sfilepoly)) then
      rexport%sfilepolyvtk = sfilepoly
    end if
    
    ! Do we have any parameters?
    if (present(cparam)) then
      rexport%cparam = cparam
    end if
    
    ! How many vertices do we have in the trangulation that have to be 
    ! filled with values?
    rexport%nvertices = rtriangulation%NVT
    rexport%ncells = rtriangulation%NEL
    
    if ((iand(cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NMT
    end if

    if ((iand(cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
        (iand(cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      rexport%nvertices = rexport%nvertices + rtriangulation%NEL
      if (rtriangulation%NDIM .eq. NDIM2D) then
        rexport%ncells = rtriangulation%NEL*4
      else
        rexport%ncells = rtriangulation%NEL*8
      end if
    end if
  
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_release (rexport)

!<description>
  ! Releases all memory in the UCD export structure rexport.
!</description>

!<inputoutput>
  ! The ucd export structure that is to be released.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>
    integer :: i
    
    ! Release all memory
    
    if (associated(rexport%p_SvariableNames)) deallocate(rexport%p_SvariableNames)
    if (associated(rexport%p_SvarVecNames)) deallocate(rexport%p_SvarVecNames)
    if (associated(rexport%p_IvariableSpec)) deallocate(rexport%p_IvariableSpec)
    if (associated(rexport%p_IvariableBase)) deallocate(rexport%p_IvariableBase)
    if (associated(rexport%p_Ivectors)) deallocate(rexport%p_Ivectors)
    if (associated(rexport%p_StracerVariableNames)) deallocate(rexport%p_StracerVariableNames)
    if (associated(rexport%p_SvertexMaterials)) deallocate(rexport%p_SvertexMaterials)
    if (associated(rexport%p_ScellMaterials)) deallocate(rexport%p_ScellMaterials)
    if (associated(rexport%p_Scomments))      deallocate(rexport%p_Scomments)
    
    if(rexport%hpolygonMaterial .ne. ST_NOHANDLE) call storage_free(rexport%hpolygonMaterial)
    
    
    
    if (associated(rexport%p_Hvariables)) then
      do i=1,rexport%nvariables
        if (rexport%p_Hvariables(i) .ne. ST_NOHANDLE) &
          call storage_free (rexport%p_Hvariables(i))
      end do
      deallocate(rexport%p_Hvariables   )
    end if
    
    if (associated(rexport%p_Hpolygons    )) then
      do i=1,rexport%npolygons
        if (rexport%p_Hpolygons(i) .ne. ST_NOHANDLE) &
          call storage_free (rexport%p_Hpolygons(i))
      end do
      deallocate(rexport%p_Hpolygons    )
    end if
    

    if (associated(rexport%p_HsurfTris    )) then
      do i=1,rexport%nsurfTri
        if (rexport%p_HsurfTris(i) .ne. ST_NOHANDLE) &
          call storage_free (rexport%p_HsurfTris(i))
          call storage_free (rexport%p_HTriangles(i))
          call storage_free (rexport%p_HsurfData(i))
      end do
      deallocate(rexport%p_HsurfTris)
      deallocate(rexport%p_HTriangles)
      deallocate(rexport%p_HsurfData)
    end if

    if (associated(rexport%p_SvertexMaterials)) deallocate(rexport%p_SvertexMaterials)
    if (associated(rexport%p_ScellMaterials)) deallocate(rexport%p_ScellMaterials)
    
    rexport%nvectors          = 0
    rexport%nvariables        = 0
    rexport%npolygons         = 0
    rexport%nsurfTri          = 0
    
    ! Release all tracer information
    call ucd_removeTracers (rexport)

    ! Do not deallocate the tringulation -- we are not the owner!!! :-)

  end subroutine

  !************************************************************************

!<subroutine>
  
  subroutine ucd_setAlternativeSource (rexport,sfilename,caltFlags)
 
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
  character(LEN=*), intent(in) :: sfilename
  
  ! Bitfield. Combination if UCD_ASRC_xxxx constants that specify which
  ! parts of the mesh are to be found in sfilename. UCD_ASRC_ALL specifies
  ! that the whole mesh is to be found in sfilename.
  integer(I32), intent(in) :: caltFlags
!</input>
  
!<inputoutput>
  ! The ucd export structure that specifies the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setAlternativeSource')
      call sys_halt()
    end if
    
    ! Depending on the flags in caltFlags, put sfilename to the
    ! corresponding filename-strings in rexport that specify an alternative
    ! source for that part of the triangulation.
    !
    ! Point coordinates:
    if (iand(caltFlags,UCD_ASRC_POINTS) .ne. 0) then
      rexport%saltFilePoints = sfilename
    end if

    ! Cell structure
    if (iand(caltFlags,UCD_ASRC_CELLS) .ne. 0) then
      rexport%saltFileCells = sfilename
    end if

    ! Material structure
    if (iand(caltFlags,UCD_ASRC_MATERIALS) .ne. 0) then
      rexport%saltFileMaterials = sfilename
    end if

    ! Polygons
    if (iand(caltFlags,UCD_ASRC_POLYGONS) .ne. 0) then
      rexport%saltFilePolygons = sfilename
    end if
    
    ! The output routine must take care, that rexport%saltFileXXXX is
    ! correctly interpreted!
  
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ucd_setMaterials (rexport,SmaterialsCells,SmaterialsVert)
  
!<description>
  ! This routine allows to specify names for material ID`s.
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
  ! The i-th string specifies a material id of material i.
  character(LEN=SYS_NAMELEN), dimension(:), intent(in) :: SmaterialsCells
  
  ! OPTIONAL: Array with strings for the vertex/node materials.
  ! The i-th string specifies a material id of material i.
  ! If not specified, the same material names will be used for both,
  ! cells and vertices.
  character(LEN=SYS_NAMELEN), dimension(:), intent(in), optional :: SmaterialsVert
!</input>
  
!<inputoutput>
  ! The ucd export structure that specifies the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setMaterials')
      call sys_halt()
    end if
    
    ! Make a copy of the strings. Use ALLOCATE/DEALLOCATE directly.
    if (associated(rexport%p_ScellMaterials)) deallocate(rexport%p_ScellMaterials)
    allocate(rexport%p_ScellMaterials(size(SmaterialsCells)))
    rexport%p_ScellMaterials = SmaterialsCells
    
    if (associated(rexport%p_SvertexMaterials)) deallocate(rexport%p_SvertexMaterials)
    
    if (present(SmaterialsVert)) then
      allocate(rexport%p_SvertexMaterials(size(SmaterialsVert)))
      rexport%p_SvertexMaterials = SmaterialsVert
    else
      allocate(rexport%p_SvertexMaterials(size(SmaterialsCells)))
      rexport%p_SvertexMaterials = SmaterialsCells
    end if
  
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ucd_setCellMaterial (rexport,Imaterials)
  
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
  integer, intent(in), dimension(:) :: Imaterials
!</input>
  
!<inputoutput>
  ! The ucd export structure that specifies the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_Idata
  integer :: NEL

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setCellMaterial')
      call sys_halt()
    end if
    
    NEL = rexport%p_rtriangulation%NEL
    
    if (size(Imaterials) .lt. NEL) then
      call output_line ('Imaterials invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setCellMaterial')
      call sys_halt()
    end if
    
    ! Copy that data and save it to the rexport structure.
    ! Create a new hImaterials handle if it does not exist.
    if (rexport%hIcellMaterial .eq. ST_NOHANDLE) then
      call storage_new ('ucd_setCellMaterial','hIcellMaterial',&
          NEL,ST_INT,rexport%hIcellMaterial,ST_NEWBLOCK_NOINIT)
    end if
    
    call storage_getbase_int (rexport%hIcellMaterial,p_Idata)
    call lalg_copyVectorInt (Imaterials(1:NEL),p_Idata(1:NEL))
  
  end subroutine
  
  !************************************************************************
  
!<subroutine>

  subroutine ucd_setVertexMaterial (rexport,&
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
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Array with as many elements as NVT in the triangulation. For every
  ! vertex i, ImaterialsVert(i) specifies the vertex material id that
  ! should be assigned to that vertex.
  integer, intent(in), dimension(:) :: ImaterialsVert

  ! OPTIONAL: Array with as many elements as NMT in the triangulation. 
  ! For every edge i, ImaterialsMid(i) specifies the material id
  ! that should be assigned to the corresponding edge midpoint.
  integer, intent(in), dimension(:), optional :: ImaterialsMid

  ! OPTIONAL: Array with as many elements as NEL in the triangulation. 
  ! For every element i, ImaterialsElem(i) specifies the material id
  ! that should be assigned to the corresponding element midpoint.
  ! Note: The material of the element midpoint need not to coincide
  !  with the material of the element -- which is specified
  !  in ucd_setCellMaterial!
  integer, intent(in), dimension(:), optional :: ImaterialsElem

!</input>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_Idata
  integer :: NEL
  integer :: NMT
  integer :: NVT

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
      call sys_halt()
    end if
    
    NVT = rexport%p_rtriangulation%NVT
    NMT = rexport%p_rtriangulation%NMT
    NEL = rexport%p_rtriangulation%NEL
    
    if (size(ImaterialsVert) .lt. NVT) then
      call output_line ('ImaterialsVert invalid!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
      call sys_halt()
    end if

    if (present(ImaterialsMid)) then
      if (size(ImaterialsMid) .lt. NMT) then
        call output_line ('ImaterialsMid invalid!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
        call sys_halt()
      end if
    end if

    if (present(ImaterialsElem)) then
      if (size(ImaterialsElem) .lt. NEL) then
        call output_line ('ImaterialsElem invalid!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setVertexMaterial')
        call sys_halt()
      end if
    end if
    
    ! Create a new data array for the vertex materials if necessary.
    ! Fill it with 0, which is the default material id.
    
    ! Create a new hImaterials handle if it does not exist.
    if (rexport%hIvertexMaterial .eq. ST_NOHANDLE) then
      call storage_new ('ucd_setVertexMaterial','hIvertexMaterial',&
          rexport%nvertices,ST_INT,rexport%hIvertexMaterial,ST_NEWBLOCK_ZERO)
    end if

    call storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
    
    ! Copy that data and save it to the rexport structure.
    call lalg_copyVectorInt (ImaterialsVert(1:NVT),p_Idata(1:NVT))

    ! Copy edge midpoint data if available
    if ((iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
        (iand(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
        (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      if (present(ImaterialsMid)) then
        call lalg_copyVectorInt( &
            ImaterialsMid(1:NMT),p_Idata(NVT+1:NVT+NMT))
      end if
    end if
    
    ! Copy element midpoint data if available
    if ((iand(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
        (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      if (present(ImaterialsElem)) then
        call lalg_copyVectorInt( &
            ImaterialsElem(1:NEL),p_Idata(NVT+NMT+1:NVT+NMT+NEL))
      end if
    end if    

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_write (rexport)

!<description>
  ! Writes the output of UCD data into a postprocessig file.
  ! All pending data in rexport is written into the file. 
  ! If the file is identified by a filename, a new file is opened,
  ! data is written to it and the file is closed at the end.
!</description>

!<inputoutput>
  ! The ucd export structure that specifies the output.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    ! Check for special binary file output
    select case(rexport%coutputFormat)
    case (UCD_FORMAT_BGMV)
      call ucd_writeBGMV (rexport)
      return
    end select

    ! If there is a new filename, open the output file.
    if (rexport%sfilename .ne. '') then
      call io_openFileForWriting(rexport%sfilename, rexport%iunit, SYS_REPLACE)
      if (rexport%iunit .lt. 0) then
        call output_line ('Cannot open file "'//trim(rexport%sfilename)&
                //'" for writing!', OU_CLASS_ERROR,OU_MODE_STD,'ucd_write')
        call sys_halt()
      end if
    end if
    
    if (rexport%iunit .eq. 0) then
      call output_line ('Cannot write UCD output: No output channel/filename!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_write')
      call sys_halt()
    end if

    ! Start the output
    select case(rexport%coutputFormat)
    case (UCD_FORMAT_GMV)
      call ucd_writeGMV (rexport)
    case (UCD_FORMAT_AVS)
      call ucd_writeAVS (rexport)
    case (UCD_FORMAT_VTK)
      call ucd_writeVTK (rexport)
    end select
    
    if (rexport%sfilename .ne. '') then
      ! Close the file if it was opened previously.
      close (rexport%iunit)
      rexport%iunit = 0
    end if
    
    ! if vtk output is desired
    ! polygons are written into a seperate file
    if ((rexport%coutputFormat .eq. UCD_FORMAT_VTK).and. &
        ((rexport%npolygons .gt. 0) .or. (rexport%nsurfTri .gt. 0)))then
      call ucd_writeVTKPolygon(rexport)
    end if

    
  contains
    
    !****************************************************************
    
    subroutine ucd_writeGMV (rexport)

    ! Specific output routine for GMV output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startGMV
    ! subroutine before.
    
    ! The export structure with all information
    type(t_ucdExport), intent(inout) :: rexport
    
    ! local variables
    integer :: mfile,i,j,k,icoor,ncomp
    integer :: ivt,ivt1,ivt2,nnodes
    integer :: imt
    integer :: iel
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Idata
    real(DP), dimension(:,:), pointer :: p_DvertexCoords,p_Ddata2D
    integer, dimension(:,:), pointer :: p_IverticesAtEdge 
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    real(DP) :: dx
    logical :: bVec2Sc,bhasVelocity
    
      mfile = rexport%iunit

      call storage_getbase_double2d (rexport%p_Rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
      call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtElement,&
          p_IverticesAtElement)
      
      ! Should we write vector components as scalars?
      bVec2Sc = (iand(rexport%cparam, UCD_PARAM_GMV_VECTOR_TO_SCALAR) .ne. 0)

      !----------------------------------------------------
      ! Write the GMV header
      write(mfile,'(A)') 'gmvinput ascii'
      
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
      
      if (rexport%saltFilePoints .ne. "") then
      
        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disk space!
        
        write(mfile,'(A)')'nodes fromfile "'//trim(rexport%saltFilePoints)//'"'
      
      else

        write(mfile,'(A,1X,I10)') 'nodes',rexport%nvertices
        
        ! Loop through the X/Y/Z coordinates
        
        do icoor = 1,min(ubound(p_DvertexCoords,1),3)
        
          do ivt=1,rexport%p_Rtriangulation%NVT
            write(mfile,rexport%sdataFormat) p_DvertexCoords(icoor,ivt)
          end do

          ! Write coordinates of edge midpoints?
          if ((iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
              (iand(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
              (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
              
            call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                p_IverticesAtEdge)

            ! We construct them by hand.
            ! In a later implementation, one could take the coordinates
            ! of h_DfreecornerCoordinates...
            do imt=1,rexport%p_Rtriangulation%NMT
              ivt1 = p_IverticesAtEdge(1,imt)
              ivt2 = p_IverticesAtEdge(2,imt)
              dx = 0.5_DP*(p_DvertexCoords(icoor,ivt1) + &
                           p_DvertexCoords(icoor,ivt2))
              write(mfile,rexport%sdataFormat) dx
            end do
              
          end if

          ! Write coordinates of element midpoints?
          if ((iand(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
              (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
              
            call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                p_IverticesAtEdge)

            ! We construct them by hand.
            ! In a later implementation, one could take the coordinates
            ! of h_DfreecornerCoordinates...
            do iel=1,rexport%p_Rtriangulation%NEL
              
              dx = 0.0_DP
              
              do i=1,ubound(p_IverticesAtElement,1)
                ivt = p_IverticesAtElement(i,iel)
                if (ivt .ne. 0) then
                  dx = dx + p_DvertexCoords(icoor,ivt)
                else
                  ! We have only (i-1) vertices in that element; 
                  ! happens e.g. in triangles that are mixed into a quad mesh.
                  ! We stop here.
                  exit
                end if
              end do
              
              ! If all vertices of the element are touched, there is i=NVE+1.
              ! Divide by the number of vertices to get the coordinate of the 
              ! midpoint of the element.
              dx = dx / real(i-1,DP)
              
              write(mfile,rexport%sdataFormat) dx
            end do
              
          end if
        
        end do ! icoor
        
        ! If there are not enough coordinates, we must add 0`s -- as
        ! GMV always expects 3D data.

        do icoor = ubound(p_DvertexCoords,1)+1 , 3
        
          do ivt=1,rexport%p_Rtriangulation%NVT
            write(mfile,rexport%sdataFormat) 0.0_DP
          end do

          ! Write coordinates of edge midpoints?
          if ((iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
              (iand(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
              (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
              
            do imt=1,rexport%p_Rtriangulation%NMT
              write(mfile,rexport%sdataFormat) 0.0_DP
            end do
              
          end if

          ! Write coordinates of element midpoints?
          if ((iand(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
              (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
              
            do iel=1,rexport%p_Rtriangulation%NEL
              write(mfile,rexport%sdataFormat) 0.0_DP
            end do
              
          end if
        
        end do ! icoor
        
      end if
      
      ! Mesh connectivity / Cells:
      
      if (rexport%saltFileCells .ne. "") then
      
        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disc space!
        
        write(mfile,'(A)')'cells fromfile "'//trim(rexport%saltFileCells)//'"'
      
      else
      
        ! Write the connectivity to the mesh - i.e. the cells.

        write(mfile,'(A,1X,I10)') 'cells',rexport%ncells
        
        if (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .eq. 0) then
        
          select case (rexport%p_rtriangulation%ndim)
          
          case (NDIM1D)
        
            ! Standard mesh.
            do iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              do i=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(i,iel) .eq. 0) exit
              end do
              
              ! We have i-1 vertices on that element -- so what is it?
              select case (i-1)
              case (2)
                ! Line in 1D
                write(mfile,'(A)') 'line 2'
                write(mfile,'(2I8)') p_IverticesAtElement(1:2,iel)
              
              case DEFAULT
                call output_line ('Invalid element!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              end select
                
            end do
            
          case (NDIM2D)

            ! Standard mesh.
            do iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              do i=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(i,iel) .eq. 0) exit
              end do
              
              ! We have i-1 vertices on that element -- so what is it?
              select case (i-1)
              
              case (3)
                ! Triangle
                write(mfile,'(A)') 'tri 3'
                write(mfile,'(3I8)') p_IverticesAtElement(1:3,iel)
                
              case (4)
                ! Quad
                write(mfile,'(A)')'quad 4'
                write(mfile,'(4I8)') p_IverticesAtElement(1:4,iel)
                
              case DEFAULT
                call output_line ('Invalid element!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              end select
                
            end do

          case (NDIM3D)
          
            ! Standard mesh.
            do iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              do i=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(i,iel) .eq. 0) exit
              end do
              
              ! We have i-1 vertices on that element -- so what is it?
              select case (i-1)
              
              case (4)
                ! Tetrahedron
                write(mfile,'(A)') 'ptet4 4'
                write(mfile,'(4I8)') p_IverticesAtElement(1:4,iel)
                
              case (5)
                ! Pyramid
                write(mfile,'(A)') 'ppyrmd5 5'
                write(mfile,'(5I8)') p_IverticesAtElement(1:5,iel)

              case (6)
                ! Prism
                write(mfile,'(A)') 'pprism6 6'
                write(mfile,'(6I8)') p_IverticesAtElement(1:6,iel)

              case (8)
                ! Hexahedron
                write(mfile,'(A)')'phex8 8'
                write(mfile,'(8I8)') p_IverticesAtElement(1:8,iel)
                
              case DEFAULT
                call output_line ('Invalid element!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              end select
                
            end do

          end select          
            
        else

          ! 1x refined mesh
          
          ! The edge numbers give the numbers of the edge midpoints and thus
          ! the numbers of the new vertices on the once refined mesh.
          call storage_getbase_int2d (rexport%p_Rtriangulation%h_IedgesAtElement,&
              p_IedgesAtElement)
              
          select case (rexport%p_rtriangulation%ndim)
          
          case (NDIM1D)
              
            do iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              do i=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(i,iel) .eq. 0) exit
              end do
              
              ! We have i-1 vertices on that element -- so what is it?
              select case (i-1)
              case (2)
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
                write(mfile,'(A)') 'line 2'
                write(mfile,'(2I8)') p_IverticesAtElement(1,iel),iel+rexport%p_rtriangulation%NVT

              end select
              
            end do

            do iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              do i=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(i,iel) .eq. 0) exit
              end do
              
              ! We have i-1 vertices on that element -- so what is it?
              select case (i-1)
              case (2)
                ! Line in 1D.
                !
                ! Element "NEL+1"
                write(mfile,'(A)') 'line 2'
                write(mfile,'(2I8)') iel+rexport%p_rtriangulation%NVT,p_IverticesAtElement(2,iel)
                
              end select
              
            end do

          case (NDIM2D)
              
            do iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              do i=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(i,iel) .eq. 0) exit
              end do
              
              ! We have i-1 vertices on that element -- so what is it?
              select case (i-1)

              case (3)
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
                    
                write(mfile,'(A)') 'tri 3'
                write(mfile,'(3I8)') p_IedgesAtElement(1:3,iel)+rexport%p_rtriangulation%NVT
                
              case (4)
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
                
                write(mfile,'(A)')'quad 4'
                write(mfile,'(4I8)') &
                    p_IverticesAtElement(1,iel),&
                    p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT,&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(4,iel)+rexport%p_rtriangulation%NVT
                
              end select
              
            end do

            do iel = 1,rexport%p_rtriangulation%NEL
            
              ! Count the number of vertices on that element
              do i=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(i,iel) .eq. 0) exit
              end do
              
              ! We have i-1 vertices on that element -- so what is it?
              select case (i-1)
                
              case (3)
                ! Triangle.
                    
                ! Element "NEL+1"
                write(mfile,'(A)') 'tri 3'
                write(mfile,'(3I8)') p_IverticesAtElement(2,iel), &
                    p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT,&
                    p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT

                ! Element "NEL+2"
                write(mfile,'(A)') 'tri 3'
                write(mfile,'(3I8)') p_IverticesAtElement(1,iel), &
                    p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT,&
                    p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT

                ! Element "NEL+3"
                write(mfile,'(A)') 'tri 3'
                write(mfile,'(3I8)') p_IverticesAtElement(3,iel), &
                    p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT,&
                    p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT
                
              case (4)
                ! Quad
                !
                ! Element "NEL+1"
                write(mfile,'(A)')'quad 4'
                write(mfile,'(4I8)') &
                    p_IverticesAtElement(2,iel),&
                    p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT,&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT

                ! Element "NEL+2"
                write(mfile,'(A)')'quad 4'
                write(mfile,'(4I8)') &
                    p_IverticesAtElement(3,iel),&
                    p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT,&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT

                ! Element "NEL+3"
                write(mfile,'(A)')'quad 4'
                write(mfile,'(4I8)') &
                    p_IverticesAtElement(4,iel),&
                    p_IedgesAtElement(4,iel)+rexport%p_rtriangulation%NVT,&
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel, &
                    p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT
                
              end select
              
            end do
            
          case (NDIM3D)

            call output_line ('GMV export for 1x refined mesh in 3D'//&
                ' not implemented!', OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
            call sys_halt()
          
          end select
          
        end if
        
      end if

      !----------------------------------------------------
      ! Write material names -- if specified
      !

      if (rexport%saltFileMaterials .ne. "") then

        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disc space!
        
        write(mfile,'(A)')'material fromfile "'//trim(rexport%saltFileMaterials)//'"'
        
      else

        ! Cell material
        if (associated(rexport%p_ScellMaterials) .and. &
            (rexport%hIcellMaterial .ne. ST_NOHANDLE)) then
          ! GMV only supports <= 1000 materials!
          write(mfile,'(A,1X,I10,I10)') 'material',min(1000,size(rexport%p_ScellMaterials)),0
          do i=1,min(1000,size(rexport%p_ScellMaterials))
            ! GMV supports only <= 8 characters and does not allow spaces
            ! in the material name. We replace all invalid spaces by "_".
            write(mfile,'(A8)') &
                sys_charreplace(trim(rexport%p_ScellMaterials(i)),' ','_')
          end do
          
          if (rexport%hIcellMaterial .ne. ST_NOHANDLE) then
            ! Write a list of material ID`s. For every cell, we specify its
            ! material by the material number.
            call storage_getbase_int (rexport%hIcellMaterial,p_Idata)
            do i=1,size(p_Idata)
              write(mfile,'(I4)') p_Idata(i)
            end do
          end if
        end if
        
        ! Vertex materials; coincide with cell materials if not specified.
        if (associated(rexport%p_SvertexMaterials) .and. &
            (rexport%hIvertexMaterial .ne. ST_NOHANDLE)) then
          ! GMV only supports <= 1000 materials!
          write(mfile,'(A,1X,I10,I10)') 'material',min(1000,size(rexport%p_SvertexMaterials)),1
          do i=1,min(1000,size(rexport%p_SvertexMaterials))
            ! GMV supports only <= 8 characters and does not allow spaces
            ! in the material name. We replace all invalid spaces by "_".
            write(mfile,'(A8)') &
                sys_charreplace(trim(rexport%p_SvertexMaterials(i)),' ','_')
          end do
          
          if (rexport%hIvertexMaterial .ne. ST_NOHANDLE) then
            ! Write a list of material ID`s. For every vertex, we specify its
            ! material by the material number.
            call storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
            do i=1,size(p_Idata)
              write(mfile,'(I4)') p_Idata(i)
            end do
          end if
        else
          if (associated(rexport%p_ScellMaterials) .and. &
              (rexport%hIvertexMaterial .ne. ST_NOHANDLE)) then
            ! GMV only supports <= 1000 materials!
            write(mfile,'(A,1X,I10,I10)') 'material',min(1000,size(rexport%p_ScellMaterials)),1
            do i=1,min(1000,size(rexport%p_ScellMaterials))
              ! GMV supports only <= 8 characters and does not allow spaces
              ! in the material name. We replace all invalid spaces by "_".
              write(mfile,'(A8)') &
                  sys_charreplace(trim(rexport%p_ScellMaterials(i)),' ','_')
            end do
          end if
          
          if (rexport%hIvertexMaterial .ne. ST_NOHANDLE) then
            ! Write a list of material ID`s. For every vertex, we specify its
            ! material by the material number.
            call storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
            do i=1,size(p_Idata)
              write(mfile,'(I4)') p_Idata(i)
            end do
          end if
        end if
      end if
      
      if (rexport%nvariables .gt. 0) then
        
        !----------------------------------------------------
        ! Write all variables which are not components of a vector field
        do j=1,rexport%nvariables

          ! Is the value scalar?
          if ((rexport%p_IvariableSpec(j) .ne. UCD_VAR_STANDARD) .and. &
              (.not. bVec2Sc)) cycle
          
          write (mfile,'(A)') 'variable'
          
          if (rexport%p_IvariableBase(j) .eq. UCD_BASE_ELEMENT) then
            ! Cell based variable
            write (mfile,'(A,1X,I5)') trim(rexport%p_SvariableNames(j)),0
            call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
            ivt1 = rexport%ncells
          else
            ! Vertex based variable
            write (mfile,'(A,1X,I5)') trim(rexport%p_SvariableNames(j)),1
            call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
            ivt1 = rexport%nvertices
          end if
          
          do ivt=1,ivt1
            write (mfile,rexport%sdataFormat) p_Ddata(ivt)
          end do
          
          write (mfile,'(A)') 'endvars'      
          
        end do ! j

        ! Should we write vectors?
        if (.not. bVec2Sc) then
        
          ! Check if we have an explicit velocity vector. Otherwise,
          ! the first vector is used as velocity vector.
          bhasVelocity = .false.
          do j=1, rexport%nvectors
            if (iand(rexport%p_IvariableSpec(rexport%p_Ivectors(2,j)),&
                UCD_VAR_VELOCITY) .ne. 0) bhasVelocity=.true.
          end do
          
          ! Loop through all vectors
          do j=1, rexport%nvectors

            ! Is this a vertex-based vector?
            if(rexport%p_Ivectors(1,j) .ne. 1) cycle
            
            ! Go for the X coordinate
            i = rexport%p_Ivectors(2,j); ncomp = 1
            call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
            
            ! Make sure we have at least the X-coordinate
            if (.not. associated(p_Ddata)) then
              call output_line ('Error: Variable vector '//&
                  trim(sys_siL(j,10))//' does not have X-coordinates!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              
              ! Try next vector
              cycle
            end if

            ! Go for the Y and Z coordinate
            if (rexport%p_Ivectors(3,j) .ne. 0) ncomp = ncomp+1
            if (rexport%p_Ivectors(4,j) .ne. 0) ncomp = ncomp+1

            ! Is this vector the velocity field? Or do we have no
            ! explicit velocity vector, then use the first vectors.
            if ((.not. bhasVelocity .and. j.eq.1) .or.&
                iand(rexport%p_IvariableSpec(i),UCD_VAR_VELOCITY) .ne. 0) then
              write (mfile,'(A)') 'velocity 1'
            else
              write (mfile,'(A)') 'vectors'
              write (mfile,'(A,1X,I5,1X,I5,1X,I5)')&
                  sys_charreplace(trim(rexport%p_SvarVecNames(j)), ' ', '_'), &
                  1, ncomp, 1
              do k = 2,4
                i = rexport%p_Ivectors(k,j)
                if (i .ne. 0) write(mfile,'(A,1X)', ADVANCE='NO')&
                    trim(rexport%p_SvariableNames(i))
              end do
              write(mfile, '(A)') ""
            end if

            ! Write X coordinate
            do ivt=1,rexport%nvertices
              write (mfile,rexport%sdataFormat) p_Ddata(ivt)
            end do

            ! Write Y coordinate (if available)
            i = rexport%p_Ivectors(3,j)
            if (i .ne. 0) then
              call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
              do ivt=1,rexport%nvertices
                write (mfile,rexport%sdataFormat) p_Ddata(ivt)
              end do
            end if
            
            ! Write Z coordinate (if available)
            i = rexport%p_Ivectors(4,j)
            if (i .ne. 0) then
              call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
              do ivt=1,rexport%nvertices
                write (mfile,rexport%sdataFormat) p_Ddata(ivt)
              end do
            end if
            
            ! Write 'endvect' tag for vectors
            i = rexport%p_Ivectors(2,j)
            if ((bhasVelocity .or. j.ne.1) .and.&
                iand(rexport%p_IvariableSpec(i), UCD_VAR_VELOCITY) .eq. 0)&
                write (mfile,'(A)') 'endvect'
            
          end do ! j

          ! Loop through all vectors
          do j=1, rexport%nvectors

            ! Is this an element-based vector?
            if(rexport%p_Ivectors(1,j) .ne. 0) cycle
            
            ! Go for the X coordinate
            i = rexport%p_Ivectors(2,j); ncomp = 1
            call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
            
            ! Make sure we have at least the X-coordinate
            if (.not. associated(p_Ddata)) then
              call output_line ('Error: Variable vector '//&
                  trim(sys_siL(j,10))//' does not have X-coordinates!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
              
              ! Try next vector
              cycle
            end if

            ! Go for the Y and Z coordinate
            if (rexport%p_Ivectors(3,j) .ne. 0) ncomp = ncomp+1
            if (rexport%p_Ivectors(4,j) .ne. 0) ncomp = ncomp+1

            ! Is this vector the velocity field? Or do we have no
            ! explicit velocity vector, then use the first vectors.
            if ((.not. bhasVelocity .and. j.eq.1) .or.&
                iand(rexport%p_IvariableSpec(i),UCD_VAR_VELOCITY) .ne. 0) then
              write (mfile,'(A)') 'velocity 0'
            else
              write (mfile,'(A)') 'vectors'
              write (mfile,'(A,1X,I5,1X,I5,1X,I5)')&
                  sys_charreplace(trim(rexport%p_SvarVecNames(j)), ' ', '_'), &
                  0, ncomp, 1
              do k = 2,4
                i = rexport%p_Ivectors(k,j)
                if (i .ne. 0) write(mfile,'(A,1X)', ADVANCE='NO')&
                    trim(rexport%p_SvariableNames(i))
              end do
              write(mfile, '(A)') ""
            end if

            ! Write X coordinate
            ! Do not be confused! ivt=number of cell, as we are in the 
            ! 'cell-oriented' case here!!!
            do ivt=1,rexport%ncells
              write (mfile,rexport%sdataFormat) p_Ddata(ivt)
            end do

            ! Write Y coordinate (if available)
            i = rexport%p_Ivectors(3,j)
            if (i .ne. 0) then
              call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
              ! Do not be confused! ivt=number of cell, as we are in the 
              ! 'cell-oriented' case here!!!
              do ivt=1,rexport%ncells
                write (mfile,rexport%sdataFormat) p_Ddata(ivt)
              end do
            end if
            
            ! Write Z coordinate (if available)
            i = rexport%p_Ivectors(4,j)
            if (i .ne. 0) then
              call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
              ! Do not be confused! ivt=number of cell, as we are in the 
              ! 'cell-oriented' case here!!!
              do ivt=1,rexport%ncells
                write (mfile,rexport%sdataFormat) p_Ddata(ivt)
              end do
            end if
            
            ! Write 'endvect' tag for vectors
            i = rexport%p_Ivectors(2,j)
            if ((bhasVelocity .or. j.ne.1) .and.&
                iand(rexport%p_IvariableSpec(i), UCD_VAR_VELOCITY) .eq. 0)&
                write (mfile,'(A)') 'endvect'

          end do ! j

        end if

      end if

      !----------------------------------------------------
      ! Write polygon data

      if (rexport%saltFilePolygons .ne. "") then
      
        ! Write only a reference to the alternative source file
        ! to the GMV. Saves disc space!
        
        write(mfile,'(A)')'polygons fromfile "'//trim(rexport%saltFilePolygons)//'"'
      
      else
        
        if (associated(rexport%p_Hpolygons)) then
          
          ! At least one polygon
          write (mfile,'(A)') 'polygons'
          
          ! Materials
          call storage_getbase_int (rexport%hpolygonMaterial,p_Idata)
          
          ! Write all polygons.
          do i=1,rexport%npolygons
            
            ! Coordinates        
            call storage_getbase_double2D (rexport%p_Hpolygons(i),p_Ddata2D)
            
            ! Write material, #points
            write (mfile,'(2I10)') p_Idata(i),ubound(p_Ddata2D,2)
            
            ! Either we have 2D or 3D coordinates. 
            ! Write coordinates of the points forming the line segments 
            ! of the polygon
            ! First all X-, then all Y- and at the end all Z-coordinates -- or 0.0.
            do k=1,NDIM3D
              if (ubound(p_Ddata2D,1) .ge. k) then
                do j=1,ubound(p_Ddata2D,2)
                  write (mfile,'(E15.7)') p_Ddata2D(k,j)
                end do
              else
                do j=1,ubound(p_Ddata2D,2)
                  write (mfile,'(E15.7)') 0.0_DP
                end do
              end if
            end do
            
          end do
          
          write (mfile,'(A)') 'endpoly'
          
        end if

      end if

      !----------------------------------------------------
      ! Write tracer coordinates and data
      if (rexport%ntracers .ne. 0) then
        write(mfile,'(A,1X,I10)') 'tracers',rexport%ntracers
        
        call storage_getbase_double2d(rexport%htracers,p_Ddata2D)
        
        ! First write all X-coordinates, then Y-coordinates, 
        ! then Z-coordinates -- if specified.
        do i=1,3
          if (i .le. ubound(p_Ddata2d,1)) then
            do j=1,rexport%ntracers
              write (mfile,rexport%sdataFormat) p_Ddata2D(i,j)
            end do
          else
            do j=1,rexport%ntracers
              write (mfile,rexport%sdataFormat) 0.0_DP
            end do
          end if
        end do
        
        ! Write tracer variables if specified
        do i=1,rexport%ntracerVariables
          write (mfile,'(A32)') rexport%p_StracerVariableNames(i)
          
          call storage_getbase_double (rexport%p_HtracerVariables(i), p_Ddata)
          do j=1,rexport%ntracers
            write (mfile,rexport%sdataFormat) p_Ddata(j)
          end do
        end do
        
        write (mfile,'(A)') 'endtrace'
        
      end if

      !----------------------------------------------------
      ! Simulation time
      if (rexport%dsimulationTime .ne. SYS_INFINITY) then
        write(mfile,'(A)',ADVANCE='NO')     'probtime '
        write(mfile,rexport%ssimTimeFormat) rexport%dsimulationTime
      end if

      !----------------------------------------------------
      ! Write all comments
      if (rexport%ncommentBufSize .gt. 0) then
        write(mfile,'(A)') 'comments'
        
        i = 1
        ! Find the end of the first line
        do j=i,rexport%ncommentBufSize
          if (rexport%p_Scomments(j) .eq. NEWLINE) exit
        end do
        
        ! Write out all lines, one after the other
        do while (j .le. rexport%ncommentBufSize)
          ! Write the line (character i..j-1), continue with the next
          do k=i,j-1
            write(mfile,'(A)',ADVANCE='NO') rexport%p_Scomments(k)
          end do
          write(mfile,'(A)')
          
          ! Continue after the NEWLINE character, find the next NEWLINE.
          i = j+1
          do j=i,rexport%ncommentBufSize
            if (rexport%p_Scomments(j) .eq. NEWLINE) exit
          end do
        
        end do
        
        write(mfile,'(A)') 'endcomm'
      end if

      !----------------------------------------------------
      ! Write information about generating code
      write(mfile,'(A)') 'codename Featflow'
      write(mfile,'(A)') 'codever  2.0'
      
      !----------------------------------------------------
      ! Finally write the GMV footer, finish
      
      write(mfile,'(A)') 'endgmv'
      
    end subroutine

    !****************************************************************
    
    subroutine ucd_writeBGMV (rexport)

    ! Specific output routine for GMV output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startBGMV
    ! subroutine before.
    
    ! The export structure with all information
    type(t_ucdExport), intent(inout) :: rexport
    
    ! local variables
    integer :: i,j,k,icoor
    integer :: ivt,ivt1,ivt2,nnodes,nvt,nverts,matnum
    integer :: imt
    integer :: iel
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Idata
    real(DP), dimension(:,:), pointer :: p_DvertexCoords,p_Ddata2D
    integer, dimension(:,:), pointer :: p_IverticesAtEdge 
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement

    integer :: nod2ids(2), nod3ids(3), nod4ids(4), nod5ids(5), nod6ids(6), nod8ids(8)
    real, dimension(:), allocatable :: X,Y,Z,VAR
    real dx, dy,dz

    ! Open file for binary output
    call fgmvwrite_openfile(rexport%sfilename)
    
    call storage_getbase_double2d (rexport%p_Rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    
    !----------------------------------------------------
    ! Simulation time
    if (rexport%dsimulationTime .ne. SYS_INFINITY) then
      call fgmvwrite_probtime(rexport%dsimulationTime)
    end if

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
    
    if (rexport%saltFilePoints .ne. "") then
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!
      
      call fgmvwrite_nodes_fromfile(trim(rexport%saltFilePoints), int(rexport%nvertices))
      
    else

      ! Allocate temporal memory
      allocate(X(rexport%nvertices), Y(rexport%nvertices), Z(rexport%nvertices))
      
      select case(rexport%p_Rtriangulation%ndim)
      
      case(NDIM1D)
        do ivt=1,rexport%p_Rtriangulation%NVT
          X(ivt) =  real(p_DvertexCoords(1,ivt))
          Y(ivt) = 0.0E0
          Z(ivt) = 0.0E0
        end do

        ! Store number of vertives already processed
        nvt = rexport%p_Rtriangulation%NVT

        ! Write coordinates of edge midpoints?
        if ((iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
          
          call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          do imt=1,rexport%p_Rtriangulation%NMT
            ivt1 = p_IverticesAtEdge(1,imt)
            ivt2 = p_IverticesAtEdge(2,imt)
            X(nvt+imt) = 0.5*(real(p_DvertexCoords(1,ivt1)) + &
                              real(p_DvertexCoords(1,ivt2)))
            Y(nvt+imt) = 0.0E0
            Z(nvt+imt) = 0.0E0
          end do
          
          ! Store number of vertives already processed
          nvt = rexport%p_Rtriangulation%NVT + rexport%p_Rtriangulation%NMT
        end if

        ! Write coordinates of element midpoints?
        if ((iand(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
          
          call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
              p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          do iel=1,rexport%p_Rtriangulation%NEL
            
            dx = 0.0E0
            
            do i=1,ubound(p_IverticesAtElement,1)
              ivt = p_IverticesAtElement(i,iel)
              if (ivt .ne. 0) then
                dx = dx + real(p_DvertexCoords(1,ivt))
              else
                ! We have only (i-1) vertices in that element; 
                ! happens e.g. in triangles that are mixed into a quad mesh.
                ! We stop here.
                exit
              end if
            end do
            
            ! If all vertices of the element are touched, there is i=NVE+1.
            ! Divide by the number of vertices to get the coordinate of the 
            ! midpoint of the element.
            dx = dx / real(i-1)
            
            X(nvt+iel) = dx
            Y(nvt+iel) = 0.0E0
            Z(nvt+iel) = 0.0E0
          end do
          
        end if


      case(NDIM2D)
        do ivt=1,rexport%p_Rtriangulation%NVT
          X(ivt) =  real(p_DvertexCoords(1,ivt))
          Y(ivt) =  real(p_DvertexCoords(2,ivt))
          Z(ivt) = 0.0E0
        end do

        ! Store number of vertives already processed
        nvt = rexport%p_Rtriangulation%NVT

        ! Write coordinates of edge midpoints?
        if ((iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
          
          call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          do imt=1,rexport%p_Rtriangulation%NMT
            ivt1 = p_IverticesAtEdge(1,imt)
            ivt2 = p_IverticesAtEdge(2,imt)
            X(nvt+imt) = 0.5*(real(p_DvertexCoords(1,ivt1)) + &
                              real(p_DvertexCoords(1,ivt2)))
            Y(nvt+imt) = 0.5*(real(p_DvertexCoords(2,ivt1)) + &
                              real(p_DvertexCoords(2,ivt2)))
            Z(nvt+imt) = 0.0E0
          end do
          
          ! Store number of vertives already processed
          nvt = rexport%p_Rtriangulation%NVT + rexport%p_Rtriangulation%NMT
        end if

        ! Write coordinates of element midpoints?
        if ((iand(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
          
          call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          do iel=1,rexport%p_Rtriangulation%NEL
            
            dx = 0.0E0
            dy = 0.0E0
            
            do i=1,ubound(p_IverticesAtElement,1)
              ivt = p_IverticesAtElement(i,iel)
              if (ivt .ne. 0) then
                dx = dx + real(p_DvertexCoords(1,ivt))
                dy = dy + real(p_DvertexCoords(2,ivt))
              else
                ! We have only (i-1) vertices in that element; 
                ! happens e.g. in triangles that are mixed into a quad mesh.
                ! We stop here.
                exit
              end if
            end do
            
            ! If all vertices of the element are touched, there is i=NVE+1.
            ! Divide by the number of vertices to get the coordinate of the 
            ! midpoint of the element.
            dx = dx / real(i-1)
            dy = dy / real(i-1)
            
            X(nvt+iel) = dx
            Y(nvt+iel) = dy
            Z(nvt+iel) = 0.0E0
          end do
          
        end if


      case(NDIM3D)
        do ivt=1,rexport%p_Rtriangulation%NVT
          X(ivt) =  real(p_DvertexCoords(1,ivt))
          Y(ivt) =  real(p_DvertexCoords(2,ivt))
          Z(ivt) =  real(p_DvertexCoords(3,ivt))
        end do

        ! Store number of vertives already processed
        nvt = rexport%p_Rtriangulation%NVT

        ! Write coordinates of edge midpoints?
        if ((iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
          
          call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          do imt=1,rexport%p_Rtriangulation%NMT
            ivt1 = p_IverticesAtEdge(1,imt)
            ivt2 = p_IverticesAtEdge(2,imt)
            X(nvt+imt) = 0.5*(real(p_DvertexCoords(1,ivt1)) + &
                              real(p_DvertexCoords(1,ivt2)))
            Y(nvt+imt) = 0.5*(real(p_DvertexCoords(2,ivt1)) + &
                              real(p_DvertexCoords(2,ivt2)))
            Z(nvt+imt) = 0.5*(real(p_DvertexCoords(3,ivt1)) + &
                              real(p_DvertexCoords(3,ivt2)))
          end do

          ! Store number of vertives already processed
          nvt = rexport%p_Rtriangulation%NVT + rexport%p_Rtriangulation%NMT
        end if

        ! Write coordinates of element midpoints?
        if ((iand(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
            (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
          
          call storage_getbase_int2d (rexport%p_Rtriangulation%h_IverticesAtEdge,&
                                      p_IverticesAtEdge)
          
          ! We construct them by hand.
          ! In a later implementation, one could take the coordinates
          ! of h_DfreecornerCoordinates...
          do iel=1,rexport%p_Rtriangulation%NEL
            
            dx = 0.0E0
            dy = 0.0E0
            dz = 0.0E0
            
            do i=1,ubound(p_IverticesAtElement,1)
              ivt = p_IverticesAtElement(i,iel)
              if (ivt .ne. 0) then
                dx = dx + real(p_DvertexCoords(1,ivt))
                dy = dy + real(p_DvertexCoords(2,ivt))
                dz = dz + real(p_DvertexCoords(3,ivt))
              else
                ! We have only (i-1) vertices in that element; 
                ! happens e.g. in triangles that are mixed into a quad mesh.
                ! We stop here.
                exit
              end if
            end do
            
            ! If all vertices of the element are touched, there is i=NVE+1.
            ! Divide by the number of vertices to get the coordinate of the 
            ! midpoint of the element.
            dx = dx / real(i-1)
            dy = dy / real(i-1)
            dz = dz / real(i-1)
            
            X(nvt+iel) = dx
            Y(nvt+iel) = dy
            Z(nvt+iel) = dz
          end do
          
        end if
        
      end select
      
      ! Write nodes to file
      call fgmvwrite_node_data(rexport%nvertices, X, Y, Z)
      
      ! Deallocate temporal memory
      deallocate(X, Y, Z)

    end if

    
    ! Mesh connectivity / Cells:
    
    if (rexport%saltFileCells .ne. "") then
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!
      
      call fgmvwrite_cells_fromfile(trim(rexport%saltFileCells), int(rexport%ncells))
      
    else

      ! Write the connectivity to the mesh - i.e. the cells.
      call fgmvwrite_cell_header(int(rexport%ncells))
      
      if (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .eq. 0) then
        
        select case (rexport%p_rtriangulation%ndim)
          
        case (NDIM1D)
          
          ! Standard mesh.
          do iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            do i=1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(i,iel) .eq. 0) exit
            end do
            
            ! We have i-1 vertices on that element -- so what is it?
            select case (i-1)
            case (2)
              ! Line in 1D
              nod2ids(1) = int(p_IverticesAtElement(1,iel))
              nod2ids(2) = int(p_IverticesAtElement(2,iel))
              call fgmvwrite_cell_type('line 2',2,nod2ids)
              
            case DEFAULT
              call output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            end select
            
          end do
          
        case (NDIM2D)
          
          ! Standard mesh.
          do iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            do i=1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(i,iel) .eq. 0) exit
            end do
            
            ! We have i-1 vertices on that element -- so what is it?
            select case (i-1)
              
            case (3)
              ! Triangle
              nod3ids(1) = p_IverticesAtElement(1,iel)
              nod3ids(2) = p_IverticesAtElement(2,iel)
              nod3ids(3) = p_IverticesAtElement(3,iel)
              call fgmvwrite_cell_type('tri 3',3,nod3ids)
              
            case (4)
              ! Quad
              nod4ids(1) = p_IverticesAtElement(1,iel)
              nod4ids(2) = p_IverticesAtElement(2,iel)
              nod4ids(3) = p_IverticesAtElement(3,iel)
              nod4ids(4) = p_IverticesAtElement(4,iel)
              call fgmvwrite_cell_type('quad 4',4,nod4ids)
              
            case DEFAULT
              call output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            end select
            
          end do
          
        case (NDIM3D)
          
          ! Standard mesh.
          do iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            do i=1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(i,iel) .eq. 0) exit
            end do
            
            ! We have i-1 vertices on that element -- so what is it?
            select case (i-1)
              
            case (4)
              ! Tetrahedron
              nod4ids(1) = p_IverticesAtElement(1,iel)
              nod4ids(2) = p_IverticesAtElement(2,iel)
              nod4ids(3) = p_IverticesAtElement(3,iel)
              nod4ids(4) = p_IverticesAtElement(4,iel)
              call fgmvwrite_cell_type('ptet4 4',4,nod4ids)

            case (5)
              ! Pyramid
              nod5ids(1) = p_IverticesAtElement(1,iel)
              nod5ids(2) = p_IverticesAtElement(2,iel)
              nod5ids(3) = p_IverticesAtElement(3,iel)
              nod5ids(4) = p_IverticesAtElement(4,iel)
              nod5ids(5) = p_IverticesAtElement(5,iel)
              call fgmvwrite_cell_type('ppyrmd5 5',5,nod5ids)
              
            case(6)
              ! Prism
              nod6ids(1) = p_IverticesAtElement(1,iel)
              nod6ids(2) = p_IverticesAtElement(2,iel)
              nod6ids(3) = p_IverticesAtElement(3,iel)
              nod6ids(4) = p_IverticesAtElement(4,iel)
              nod6ids(5) = p_IverticesAtElement(5,iel)
              nod6ids(6) = p_IverticesAtElement(6,iel)
              call fgmvwrite_cell_type('pprism6 6',6,nod6ids)

            case (8)
              ! Hexahedron
              nod8ids(1) = p_IverticesAtElement(1,iel)
              nod8ids(2) = p_IverticesAtElement(2,iel)
              nod8ids(3) = p_IverticesAtElement(3,iel)
              nod8ids(4) = p_IverticesAtElement(4,iel)
              nod8ids(5) = p_IverticesAtElement(5,iel)
              nod8ids(6) = p_IverticesAtElement(6,iel)
              nod8ids(7) = p_IverticesAtElement(7,iel)
              nod8ids(8) = p_IverticesAtElement(8,iel)
              call fgmvwrite_cell_type('phex8 8',8,nod8ids)
              
            case DEFAULT
              call output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            end select
            
          end do
          
        end select
        
      else

        ! 1x refined mesh
        
        ! The edge numbers give the numbers of the edge midpoints and thus
        ! the numbers of the new vertices on the once refined mesh.
        call storage_getbase_int2d (rexport%p_Rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElement)
        
        select case (rexport%p_rtriangulation%ndim)
          
        case (NDIM1D)
          
          do iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            do i=1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(i,iel) .eq. 0) exit
            end do
            
            ! We have i-1 vertices on that element -- so what is it?
            select case (i-1)
            case (2)
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
              nod2ids(2) = iel+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('line 2',2,nod2ids)

            case DEFAULT
              call output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')             
            end select
            
          end do
          
          do iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            do i=1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(i,iel) .eq. 0) exit
            end do
            
            ! We have i-1 vertices on that element -- so what is it?
            select case (i-1)
            case (2)
              ! Line in 1D.
              !
              ! Element "NEL+1"
              nod2ids(1) = iel+rexport%p_rtriangulation%NVT
              nod2ids(2) = p_IverticesAtElement(2,iel)
              call fgmvwrite_cell_type('line 2',2,nod2ids)

            case DEFAULT
              call output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')   
            end select
            
          end do
          
        case (NDIM2D)
          
          do iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            do i=1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(i,iel) .eq. 0) exit
            end do
            
            ! We have i-1 vertices on that element -- so what is it?
            select case (i-1)
              
            case (3)
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
              
              nod3ids(1) = p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT
              nod3ids(2) = p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT
              nod3ids(3) = p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('tri 3',3,nod3ids)
              
            case (4)
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
              nod4ids(2) = p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(4,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('quad 4',4,nod4ids)

            case DEFAULT
              call output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            end select
            
          end do
          
          do iel = 1,rexport%p_rtriangulation%NEL
            
            ! Count the number of vertices on that element
            do i=1,ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(i,iel) .eq. 0) exit
            end do
            
            ! We have i-1 vertices on that element -- so what is it?
            select case (i-1)
              
            case (3)
              ! Triangle.
              
              ! Element "NEL+1"
              nod3ids(1) = p_IverticesAtElement(2,iel)
              nod3ids(2) = p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT
              nod3ids(3) = p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('tri 3',3,nod3ids)
              
              ! Element "NEL+2"
              nod3ids(1) = p_IverticesAtElement(1,iel)
              nod3ids(2) = p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT
              nod3ids(3) = p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('tri 3',3,nod3ids)
              
              ! Element "NEL+3"
              nod3ids(1) = p_IverticesAtElement(3,iel)
              nod3ids(2) = p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT
              nod3ids(3) = p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('tri 3',3,nod3ids)
              
            case (4)
              ! Quad
              !
              ! Element "NEL+1"
              nod4ids(1) = p_IverticesAtElement(2,iel)
              nod4ids(2) = p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(1,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('quad 4',4,nod4ids)
              
              ! Element "NEL+2"
              nod4ids(1) = p_IverticesAtElement(3,iel)
              nod4ids(2) = p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(2,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('quad 4',4,nod4ids)
              
              ! Element "NEL+3"
              nod4ids(1) = p_IverticesAtElement(4,iel)
              nod4ids(2) = p_IedgesAtElement(4,iel)+rexport%p_rtriangulation%NVT
              nod4ids(3) = rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+iel
              nod4ids(4) = p_IedgesAtElement(3,iel)+rexport%p_rtriangulation%NVT
              call fgmvwrite_cell_type('quad 4',4,nod4ids)

            case DEFAULT
              call output_line ('Invalid element!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
            end select
            
          end do
          
        case (NDIM3D)
          
          call output_line ('GMV export for 1x refined mesh in 3D'//&
              ' not implemented!', OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeGMV')
          call sys_halt()
          
        end select
        
      end if
      
    end if


    !----------------------------------------------------
    ! Write material names -- if specified
    !

    if (rexport%saltFileMaterials .ne. "") then
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!

      call fgmvwrite_material_fromfile(trim(rexport%saltFileMaterials))
      
    else
      
      ! Cell material
      if (associated(rexport%p_ScellMaterials)) then
        ! GMV only supports <= 1000 materials!
        call fgmvwrite_material_header(min(1000,int(size(rexport%p_ScellMaterials))),0)
        do i=1,min(1000,size(rexport%p_ScellMaterials))
          ! GMV supports only <= 8 characters and does not allow spaces
          ! in the material name. We replace all invalid spaces by "_".
          call fgmvwrite_material_name(&
              sys_charreplace(trim(rexport%p_ScellMaterials(i)),' ','_'))
        end do
        
        if (rexport%hIcellMaterial .ne. ST_NOHANDLE) then
          ! Write a list of material ID`s. For every cell, we specify its
          ! material by the material number.
          call storage_getbase_int (rexport%hIcellMaterial,p_Idata)
          do i=1,size(p_Idata)
            call fgmvwrite_material_ids(int(p_Idata(i)), 0)
          end do
        end if
      end if
      
      ! Vertex materials; coincide with cell materials if not specified.
      if (associated(rexport%p_SvertexMaterials)) then
        ! GMV only supports <= 1000 materials!
        call fgmvwrite_material_header(min(1000,int(size(rexport%p_SvertexMaterials))),1)
        do i=1,min(1000,size(rexport%p_SvertexMaterials))
          ! GMV supports only <= 8 characters and does not allow spaces
          ! in the material name. We replace all invalid spaces by "_".
          call fgmvwrite_material_name(&
              sys_charreplace(trim(rexport%p_SvertexMaterials(i)),' ','_'))
        end do
        
        if (rexport%hIvertexMaterial .ne. ST_NOHANDLE) then
          ! Write a list of material ID`s. For every vertex, we specify its
          ! material by the material number.
          call storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
          do i=1,size(p_Idata)
            call fgmvwrite_material_ids(int(p_Idata(i)), 1)
          end do
        end if
      else
        if (associated(rexport%p_ScellMaterials)) then
          ! GMV only supports <= 1000 materials!
          call fgmvwrite_material_header(min(1000,int(size(rexport%p_ScellMaterials))),1)
          do i=1,min(1000,size(rexport%p_ScellMaterials))
            ! GMV supports only <= 8 characters and does not allow spaces
            ! in the material name. We replace all invalid spaces by "_".
            call fgmvwrite_material_name(&
                sys_charreplace(trim(rexport%p_ScellMaterials(i)),' ','_'))
          end do
        end if
        
        if (rexport%hIvertexMaterial .ne. ST_NOHANDLE) then
          ! Write a list of material ID`s. For every vertex, we specify its
          ! material by the material number.
          call storage_getbase_int (rexport%hIvertexMaterial,p_Idata)
          do i=1,size(p_Idata)
            call fgmvwrite_material_ids(int(p_Idata(i)), 1)
          end do
        end if
      end if
    end if
    
    
    if (rexport%nvariables .gt. 0) then
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
      do i=1,rexport%nvariables
        if ((iand(rexport%p_IvariableSpec(i),UCD_VAR_XVELOCITY) .ne. 0) .and. &
            (rexport%p_IvariableBase(i) .eq. UCD_BASE_ELEMENT)) then
          
          call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)

          ! Allocate temporal memory
          allocate(X(rexport%ncells), Y(rexport%ncells), Z(rexport%ncells))
          
          ! Do not be confused! ivt=number of cell, as we are in the 
          ! 'cell-oriented' case here!!!
          do ivt=1,rexport%ncells
            X(ivt) = p_Ddata(ivt)
          end do
          
          ! Find the Y-velocity
          do j=1,rexport%nvariables
            if ((iand(rexport%p_IvariableSpec(j),UCD_VAR_YVELOCITY) .ne. 0) .and. &
                (rexport%p_IvariableBase(j) .eq. UCD_BASE_ELEMENT)) then
              
              ! Found it. Write it out.
              call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              do ivt=1,rexport%ncells
                Y(ivt) = p_Ddata(ivt)
              end do
              
              exit
            end if
          end do
          if (j .gt. rexport%nvariables) then
            ! Not found. Write out 0`s instead.
            do ivt=1,rexport%ncells
              Y(ivt) = 0.0E0
            end do
          end if
          
          ! Find the Z-velocity
          do j=1,rexport%nvariables
            if ((iand(rexport%p_IvariableSpec(j),UCD_VAR_ZVELOCITY) .ne. 0) .and. &
                (rexport%p_IvariableBase(j) .eq. UCD_BASE_ELEMENT)) then
              
              ! Found it. Write it out.
              call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              do ivt=1,rexport%ncells
                Z(ivt) = p_Ddata(ivt)
              end do
              
              exit
            end if
          end do
          if (j .gt. rexport%nvariables) then
            ! Not found. Write out 0`s instead.
            do ivt=1,rexport%ncells
              Z(ivt) = 0.0E0
            end do
          end if
          
          ! Write cellbased velocities
          call fgmvwrite_velocity_data(0,X,Y,Z)

          ! Deallocate temporal memory
          deallocate(X,Y,Z)

        end if
      end do
      
      ! Look for vertex based velocity.
      do i=1,rexport%nvariables
        if ((iand(rexport%p_IvariableSpec(i),UCD_VAR_XVELOCITY) .ne. 0) .and. &
            (rexport%p_IvariableBase(i) .eq. UCD_BASE_VERTEX)) then
          
          call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)

          ! Allocate temporal memory
          allocate(X(rexport%nvertices), Y(rexport%nvertices), Z(rexport%nvertices))

          do ivt=1,rexport%nvertices
            X(ivt) = p_Ddata(ivt)
          end do
          
          ! Find the Y-velocity
          do j=1,rexport%nvariables
            if ((iand(rexport%p_IvariableSpec(j),UCD_VAR_YVELOCITY) .ne. 0) .and. &
                (rexport%p_IvariableBase(j) .eq. UCD_BASE_VERTEX)) then
              
              ! Found it. Write it out.
              call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              do ivt=1,rexport%nvertices
                Y(ivt) = p_Ddata(ivt)
              end do
              
              exit
            end if
          end do
          if (j .gt. rexport%nvariables) then
            ! Not found. Write out 0`s instead.
            do ivt=1,rexport%nvertices
              Y(ivt) = 0.0E0
            end do
          end if
          
          ! Find the Z-velocity
          do j=1,rexport%nvariables
            if ((iand(rexport%p_IvariableSpec(j),UCD_VAR_ZVELOCITY) .ne. 0) .and. &
                (rexport%p_IvariableBase(j) .eq. UCD_BASE_VERTEX)) then
              
              ! Found it. Write it out.
              call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              do ivt=1,rexport%nvertices
                Z(ivt) = p_Ddata(ivt)
              end do
              
              exit
            end if
          end do
          if (j .gt. rexport%nvariables) then
            ! Not found. Write out 0`s instead.
            do ivt=1,rexport%nvertices
              Z(ivt) = 0.0E0
            end do
          end if
          
          ! Write cellbased velocities
          call fgmvwrite_velocity_data(1,X,Y,Z)

          ! Deallocate temporal memory
          deallocate(X,Y,Z)

        end if
      end do
      
      !----------------------------------------------------
      ! Write all variables which are not velocities
      do i=1,rexport%nvariables
        if (iand(rexport%p_IvariableSpec(i), UCD_VAR_VELOCITY) .eq. 0) then
          
          call fgmvwrite_variable_header()
          
          if (rexport%p_IvariableBase(i) .eq. UCD_BASE_ELEMENT) then
            ! Cell based variable
            call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)

            ! Allocate temporal memory
            allocate(VAR(rexport%ncells))
            do ivt=1,rexport%ncells
              VAR(ivt) = p_Ddata(ivt)
            end do
            
            call fgmvwrite_variable_name_data(0, trim(rexport%p_SvariableNames(i)), VAR)

            ! Deallocate temporal memory
            deallocate(VAR)

          else
            ! Vertex based variable
            call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
            
            ! Allocate temporal memory
            allocate(VAR(rexport%nvertices))
            do ivt=1,rexport%nvertices
              VAR(ivt) = p_Ddata(ivt)
            end do

            call fgmvwrite_variable_name_data(1, trim(rexport%p_SvariableNames(i)), VAR)

            ! Deallocate temporal memory
            deallocate(VAR)
          end if
          
          call fgmvwrite_variable_endvars()
         
        end if
        
      end do ! i
      
    end if


    !----------------------------------------------------
    ! Write polygon data

    if (rexport%saltFilePolygons .ne. "") then
      
      ! Write only a reference to the alternative source file
      ! to the GMV. Saves disc space!
      
      call fgmvwrite_polygons_fromfile(trim(rexport%saltFilePolygons))
      
    else
      
      if (associated(rexport%p_Hpolygons)) then
        
        ! At least one polygon
        call fgmvwrite_polygons_header()
        
        ! Materials
        call storage_getbase_int (rexport%hpolygonMaterial,p_Idata)
        
        ! Write all polygons.
        do i=1,rexport%npolygons
          
          ! Coordinates        
          call storage_getbase_double2D (rexport%p_Hpolygons(i),p_Ddata2D)
          
          ! Allocate temporal memory
          allocate(X(ubound(p_Ddata2D,2)), Y(ubound(p_Ddata2D,2)), Z(ubound(p_Ddata2D,2)))
          
          ! Either we have 2D or 3D coordinates. 
          select case(ubound(p_Ddata2D,1))
            
          case (NDIM2D)
            do j=1,ubound(p_Ddata2D,2)
              X(j) = real(p_Ddata2D(1,j))
              Y(j) = real(p_Ddata2D(2,j))
              Z(j) = 0.0E0
            end do
            
          case (NDIM3D)
            do j=1,ubound(p_Ddata2D,2)
              X(j) = real(p_Ddata2D(1,j))
              Y(j) = real(p_Ddata2D(2,j))
              Z(j) = real(p_Ddata2D(3,j))
            end do
            
          case DEFAULT
            call output_line ('Invalid spatial dimensions for polygon output!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeBGMV')
          end select
          
          nverts = ubound(p_Ddata2D,2)
          matnum = p_Idata(i)
          
          call fgmvwrite_polygons_data(nverts, matnum, X, Y, Z)
          
          ! Deallocate temporal memory
          deallocate(X,Y,Z)
          
        end do
        
        call fgmvwrite_polygons_endpoly()
        
      end if
    end if

    !----------------------------------------------------
    ! Write tracer coordinates and data
    if (rexport%ntracers .ne. 0) then
            
      call storage_getbase_double2d(rexport%htracers,p_Ddata2D)
      
      ! Allocate temporal memory
      allocate(X(rexport%ntracers), Y(rexport%ntracers), Z(rexport%ntracers))

      select case(ubound(p_Ddata2D,1))

      case(NDIM1D)
        do j=1,rexport%ntracers
          X(j) = real(p_Ddata2D(1,j))
          Y(j) = 0.0E0
          Z(j) = 0.0E0
        end do

      case(NDIM2D)
        do j=1,rexport%ntracers
          X(j) = real(p_Ddata2D(1,j))
          Y(j) = real(p_Ddata2D(2,j))
          Z(j) = 0.0E0
        end do

      case(NDIM3D)
        do j=1,rexport%ntracers
          X(j) = real(p_Ddata2D(1,j))
          Y(j) = real(p_Ddata2D(2,j))
          Z(j) = real(p_Ddata2D(3,j))
        end do
      end select
      
      call fgmvwrite_tracers_header(int(rexport%ntracers), X, Y, Z)
      
      ! Deallocate temporal memory
      deallocate(X,Y,Z)

!!$      
!!$      
!!$      ! Write tracer variables if specified
!!$      DO i=1,rexport%ntracerVariables
!!$        WRITE (mfile,'(A32)') rexport%p_StracerVariableNames(i)
!!$        
!!$        CALL storage_getbase_double (rexport%p_HtracerVariables(i), p_Ddata)
!!$        DO j=1,rexport%ntracers
!!$          WRITE (mfile,rexport%sdataFormat) p_Ddata(j)
!!$        END DO
!!$      END DO


      call fgmvwrite_tracers_endtrace()
      
    end if

    !----------------------------------------------------
    ! Finally close the GMV file, finish
    
    call fgmvwrite_closefile()
    
    end subroutine

    !****************************************************************
    
    subroutine ucd_writeAVS (rexport)

    ! Specific output routine for AVS output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startAVS
    ! subroutine before.
    
    ! The export structure with all information
    type(t_ucdExport), intent(inout) :: rexport
    
    ! local variables
    integer :: mfile,i,j,k, num_ndata, num_cdata, ncells, nverts, &
        lenTemp, ivt1, ivt2
    real(DP), dimension(:,:), pointer :: p_DvertexCoords, p_DvertexRefined, &
        p_Ddata2D
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_IcellMaterial
    real(DP), dimension(1:3) :: Dvert
    character(LEN=SYS_STRLEN) :: sdl
    type(t_ucdRefine) :: rrefine
      
      ! Get file unit and export format
      mfile = rexport%iunit
      sdl = rexport%sdataFormat
      
      ! Calculate grid refinement
      call ucd_refine(rrefine, rexport%p_rtriangulation, rexport%cflags)
      
      ! Get refined vertices
      if (rrefine%h_DvertexCoords .ne. ST_NOHANDLE) then
        call storage_getbase_double2d (rrefine%h_DvertexCoords, p_DvertexRefined)
      else
        p_DvertexRefined => null()
      end if

      ! If the mesh is to be refined, then we take the refined elements,
      ! otherwise we use the elements from the triangulation
      if (rrefine%h_IverticesAtElement .ne. ST_NOHANDLE) then
        call storage_getbase_int2d (rrefine%h_IverticesAtElement, &
            p_IverticesAtElement)
        ncells = rrefine%ncells
      else
        call storage_getbase_int2d (rexport%p_rtriangulation%h_IverticesAtElement,&
            p_IverticesAtElement)
        ncells = rexport%p_rtriangulation%NEL
      end if

      ! Get corner vertices
      call storage_getbase_double2d (rexport%p_Rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
      nverts = rexport%p_rtriangulation%NVT + rrefine%nvertices
      
      ! Get cell materials, if specified
      if(rexport%hIcellMaterial .ne. ST_NOHANDLE) then
        call storage_getbase_int (rexport%hIcellMaterial, p_IcellMaterial)
      else
        p_IcellMaterial => null()
      end if
      
      ! First we need to count how many vertex-based and how many element-based
      ! variables we have.
      num_ndata = 0
      num_cdata = 0
      do i = 1, rexport%nvariables
      
        if (rexport%p_IvariableBase(i) .eq. UCD_BASE_ELEMENT) then
          num_cdata = num_cdata + 1
        else
          num_ndata = num_ndata + 1
        end if
        
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write the AVS header
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      write(mfile, '(5I10)') nverts, ncells, num_ndata, num_cdata, 0
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write vertice coordinates
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if (ubound(p_DvertexCoords,1) .eq. 3) then
        ! 3D coordinates
        
        ! Write corner vertices
        do i=1, rexport%p_Rtriangulation%NVT
          write(mfile, '(I10,3E16.7)') i, p_DvertexCoords(1:3, i)
        end do
        
        ! Write refined vertices
        j = rexport%p_Rtriangulation%NVT
        do i=1, rrefine%nvertices
          write(mfile, '(I10,3E16.7)') j+i, p_DvertexRefined(1:3, i)
        end do
        
      else
        ! 2D coordinates
        
        ! Write corner vertices
        do i=1, rexport%p_Rtriangulation%NVT
          write(mfile, '(I10,3E16.7)') i, p_DvertexCoords(1:2, i), 0.0_DP
        end do

        ! Write refined vertices
        j = rexport%p_Rtriangulation%NVT
        do i=1, rrefine%nvertices
          write(mfile, '(I10,3E16.7)') j+i, p_DvertexRefined(1:2, i), 0.0_DP
        end do
        
      end if
    
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write elements
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      do j=1, ncells
      
        ! Count the number of vertices on that element
        do i=1,ubound(p_IverticesAtElement,1)
          if (p_IverticesAtElement(i,j) .eq. 0) exit
        end do
        
        ! Do we have cell materials?
        k = 1
        !IF (ASSOCIATED(p_IcellMaterial)) k = p_IcellMaterial(j)
        
        ! We have i-1 vertices on that element -- so what is it?
        select case (i-1)
        case (3)
          ! Triangle
          write(mfile, '(2I10,A,3I10)') j, k, ' tri', &
              p_IverticesAtElement(1:3,j)
          
        case (4)
          ! Quad
          write(mfile, '(2I10,A,4I10)') j, k, ' quad', &
              p_IverticesAtElement(1:4,j)
          
        case DEFAULT
          call output_line ('Invalid element!',&
              OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeAVS')
          
        end select
        
      end do
      
      ! We do not need the mesh refinement anymore, so destroy it
      if (rrefine%h_DvertexCoords .ne. ST_NOHANDLE) then
        call storage_free(rrefine%h_DvertexCoords)
      end if
      if (rrefine%h_IverticesAtElement .ne. ST_NOHANDLE) then
        call storage_free(rrefine%h_IverticesAtElement)
      end if
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write node variable info
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      write(mfile, '(I10)', ADVANCE='NO') num_ndata
      do i=1, num_ndata
        write(mfile, '(A)', ADVANCE='NO') " 1"
      end do
      write(mfile, '(A)') ""

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write node variable names
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      do i=1, rexport%nvariables
      
        ! If this is a node variable, write its name
        if (rexport%p_IvariableBase(i) .eq. UCD_BASE_VERTEX) then
          write(mfile, '(A,A)') &
             sys_charreplace(trim(rexport%p_SvariableNames(i)), ' ', '_'), &
             ", nounit"
        end if
      
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write node variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if (num_ndata .gt. 0) then
        
        ! Write all vertices
        do i=1, rexport%nvertices
          
          ! Write vertice index
          write(mfile, '(I10)', ADVANCE='NO') i
          
          ! Write all variable values
          do j=1, rexport%nvariables
            
            if (rexport%p_IvariableBase(j) .eq. UCD_BASE_VERTEX) then
              call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              write(mfile, sdl, ADVANCE='NO') p_Ddata(i)
            end if
          
          end do
          
          ! Write a line break
          write(mfile, '(A)') ""
          
        end do
        
      end if
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variable info
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      write(mfile, '(I10)', ADVANCE='NO') num_cdata
      do i=1, num_cdata
        write(mfile, '(A)', ADVANCE='NO') " 1"
      end do
      write(mfile, '(A)') ""

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variable names
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      do i=1, rexport%nvariables
      
        ! If this is a element variable, write its name
        if (rexport%p_IvariableBase(i) .eq. UCD_BASE_ELEMENT) then
          write(mfile, '(A,A)') &
             sys_charreplace(trim(rexport%p_SvariableNames(i)), ' ', '_'), &
             ", nounit"
        end if
      
      end do
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if (num_cdata .gt. 0) then
        
        ! Write all cells
        do i=1, rexport%ncells
          
          ! Write cell index
          write(mfile, '(I10)', ADVANCE='NO') i
          
          ! Write all variable values
          do j=1, rexport%nvariables
            
            if (rexport%p_IvariableBase(j) .eq. UCD_BASE_ELEMENT) then
              call storage_getbase_double (rexport%p_Hvariables(j),p_Ddata)
              write(mfile, sdl, ADVANCE='NO') p_Ddata(i)
            end if
          
          end do
          
          ! Write a line break
          write(mfile, '(A)') ""
          
        end do
        
      end if
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write comments
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      if (rexport%ncommentBufSize .gt. 0) then
        
        i = 1
        ! Find the end of the first line
        do j=i,rexport%ncommentBufSize
          if (rexport%p_Scomments(j) .eq. NEWLINE) exit
        end do
        
        ! Write out all lines, one after the other
        do while (j .le. rexport%ncommentBufSize)
          ! Write a comment identifier
          write(mfile, '(A)', ADVANCE='NO') "# "

          ! Write the line (character i..j-1), continue with the next
          do k=i,j-1
            write(mfile,'(A)',ADVANCE='NO') rexport%p_Scomments(k)
          end do
          write(mfile,'(A)')
          
          ! Continue after the NEWLINE character, find the next NEWLINE.
          i = j+1
          do j=i,rexport%ncommentBufSize
            if (rexport%p_Scomments(j) .eq. NEWLINE) exit
          end do
        
        end do
      
      end if

    end subroutine

    !****************************************************************
    
    subroutine ucd_writeVTKPolygon (rexport)
    ! Specific output routine for VTK polygon output. 
    
    ! The export structure with all information
    type(t_ucdExport), intent(inout) :: rexport
    integer :: mfile,i,j,k,jy,jz,ipoints,num_ndata,num_cdata,ncls,nverts,ncells,ipart
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    integer, dimension(:,:), pointer :: p_Itriangles
    integer, dimension(:), pointer :: p_IsurfData
    character(LEN=SYS_STRLEN) :: sdl
    character(LEN=SYS_STRLEN) :: sfilepoly
    integer, dimension(GEOM_STANDARD_VCOUNT) :: iconnect    
    integer :: ioffset1,ioffset2,ilengthPolylist,ioffset,ipolygonsTotal
    ! Get file name
    sfilepoly = rexport%sfilepolyvtk
    
    ! If there is a new filename, open the output file.
    if (sfilepoly .ne. '') then !		SFILENAME	'./gmv/u.vtk
      call io_openFileForWriting(sfilepoly, rexport%iunit, SYS_REPLACE)
      if (rexport%iunit .lt. 0) then
        call output_line ('Cannot open file "'//trim(sfilepoly)&
                //'" for writing!', OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeVTKPolygon')
        call sys_halt()
      end if
    end if

    if (rexport%iunit .eq. 0) then
      call output_line ('Cannot write UCD output: No output channel/filename!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_write')
      call sys_halt()
    end if
    
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Write VTK header
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    write(rexport%iunit, '(A)') "# vtk DataFile Version 2.0"
    write(rexport%iunit, '(A)') "Generated by FEATFlow 2.x"
    write(rexport%iunit, '(A)') "ASCII"

    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Write Dataset header
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! we write polydata
    write(rexport%iunit, '(A)') "DATASET POLYDATA"
    
    ! calculate the number of total points
    ipoints = GEOM_STANDARD_VCOUNT * rexport%npolygons

    do ipart=1,rexport%nsurfTri
      call storage_getbase_int(rexport%p_HsurfData(ipart),p_IsurfData)    
      ipoints = ipoints + p_IsurfData(1)
    end do
    write(rexport%iunit, '(A,I10,A)') "POINTS", ipoints, " double"
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Write vertices of simple polygons
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    do ipart=1,rexport%npolygons
      call storage_getbase_double2D (rexport%p_Hpolygons(ipart),p_Ddata2D)
      do j=1,GEOM_STANDARD_VCOUNT
         write(rexport%iunit, '(3E15.7)') p_Ddata2D(1, j),p_Ddata2D(2, j), 0.0_DP
      end do    
    end do
    
    do j=1,GEOM_STANDARD_VCOUNT
      iconnect(j)=j-1
    end do
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Write vertices of surface triangulations
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ioffset = GEOM_STANDARD_VCOUNT * rexport%npolygons
    do ipart=1,rexport%nsurfTri
      call storage_getbase_double2D(rexport%p_HsurfTris(ipart),p_Ddata2D)
      call storage_getbase_int(rexport%p_HsurfData(ipart),p_IsurfData)   
      do j=1,p_IsurfData(1)
         write(rexport%iunit, '(3E15.7)') p_Ddata2D(1, j),p_Ddata2D(2, j), p_Ddata2D(3, j)
      end do    
    end do
        
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Write list of simple polygons
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ipolygonsTotal=rexport%npolygons
    ilengthPolylist=rexport%npolygons*(GEOM_STANDARD_VCOUNT+1)
    do i=1,rexport%nsurfTri
      call storage_getbase_int(rexport%p_HsurfData(i),p_IsurfData)   
      ipolygonsTotal  = ipolygonsTotal+p_IsurfData(2) 
      ilengthPolylist=ilengthPolylist + p_IsurfData(2)*4
    end do

    write(rexport%iunit, '(A,2I10)') "POLYGONS", ipolygonsTotal, ilengthPolylist
    do ipart=1,rexport%npolygons
        write(rexport%iunit, '(65I10)') GEOM_STANDARD_VCOUNT, iconnect(:)
        do j=1,GEOM_STANDARD_VCOUNT
          iconnect(j)=iconnect(j)+GEOM_STANDARD_VCOUNT
        end do
    end do
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Write list of triangles
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    do ipart=1,rexport%nsurfTri
        call storage_getbase_int2d(rexport%p_Htriangles(ipart),p_Itriangles)
        call storage_getbase_int(rexport%p_HsurfData(ipart),p_IsurfData)       
        do j=1,p_IsurfData(2)   
        write(rexport%iunit, '(4I10)') 3,p_Itriangles(1,j)+ioffset,&
            p_Itriangles(2,j)+ioffset ,p_Itriangles(3,j)+ioffset  
        end do
        ioffset=ioffset+p_IsurfData(1)   
    end do
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Write list of lines
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    write(rexport%iunit,'(A,2I10)')"LINES",rexport%npolygons,rexport%npolygons*2+rexport%npolygons
    do ipart=1,rexport%npolygons
        ioffset1=(ipart-1)*GEOM_STANDARD_VCOUNT
        ioffset2=ioffset1+GEOM_STANDARD_VCOUNT/2 
        write(rexport%iunit, '(3I10)') 2, ioffset1, ioffset2
    end do
    
    ! first we write the polygons and lines then the triangles
    if (sfilepoly .ne. '') then
      ! Close the file if it was opened previously.
      close (rexport%iunit)
      rexport%iunit = 0
    end if


    end subroutine
    
    !****************************************************************    
    
    subroutine ucd_writeVTK (rexport)

    ! Specific output routine for VTK output. Writes the whole
    ! structure to the output channel that was opened in the ucd_startVTK
    ! subroutine before.
    
    ! The export structure with all information
    type(t_ucdExport), intent(inout) :: rexport
    
    ! local variables
    integer :: mfile,i,j,k,jy,jz,num_ndata,num_cdata,ncls,nverts,ncells
    real(DP), dimension(:,:), pointer :: p_DvertexCoords, p_DvertexRefined, &
        p_Ddata2D
    real(DP), dimension(:), pointer :: p_Ddata, p_Dx, p_Dy, p_Dz
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:), allocatable :: p_InumVertsPerCell
    real(DP), dimension(1:3) :: Dvert
    character(LEN=SYS_STRLEN) :: sdl
    logical :: bQuadratic, bVec2Sc
    type(t_ucdRefine) :: rrefine
      
      ! Get file unit and export format
      mfile = rexport%iunit
      sdl = rexport%sdataFormat
          
      ! Calculate grid refinement
      call ucd_refine(rrefine, rexport%p_rtriangulation, rexport%cflags)
      
      ! Get refined vertices
      if (rrefine%h_DvertexCoords .ne. ST_NOHANDLE) then
        call storage_getbase_double2d (rrefine%h_DvertexCoords, p_DvertexRefined)
      else
        p_DvertexRefined => null()
      end if
      
      ! If the mesh is to be refined, then we take the refined elements,
      ! otherwise we use the elements from the triangulation
      if (rrefine%h_IverticesAtElement .ne. ST_NOHANDLE) then
        call storage_getbase_int2d (rrefine%h_IverticesAtElement, &
            p_IverticesAtElement)
        ncells = rrefine%ncells
      else
        call storage_getbase_int2d (rexport%p_rtriangulation%h_IverticesAtElement,&
            p_IverticesAtElement)
        ncells = rexport%p_rtriangulation%NEL
      end if

      ! Get corner vertices
      call storage_getbase_double2d (rexport%p_rtriangulation%h_DvertexCoords,&
          p_DvertexCoords)
      nverts = rexport%p_rtriangulation%NVT + rrefine%nvertices
      
      ! Get edges of elements (needed for quadratic cells)
      if (rexport%p_rtriangulation%h_IedgesAtElement .ne. ST_NOHANDLE) then
        call storage_getbase_int2d (rexport%p_rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElement)
      end if

      ! First we need to count how many vertex-based and how many element-based
      ! variables we have.
      num_ndata = 0
      num_cdata = 0
      do i = 1, rexport%nvariables
      
        if (rexport%p_IvariableBase(i) .eq. UCD_BASE_ELEMENT) then
          num_cdata = num_cdata + 1
        else
          num_ndata = num_ndata + 1
        end if
        
      end do
      
      ! Should we write vector components as scalars?
      bVec2Sc = (iand(rexport%cparam, UCD_PARAM_VTK_VECTOR_TO_SCALAR) .ne. 0)

      ! Are we going to write quadratic elements?
      bQuadratic = (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .eq. 0) .and. &
                   (iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .and. &
                   (iand(rexport%cparam,UCD_PARAM_VTK_USE_QUADRATIC) .ne. 0)
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write VTK header
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      write(mfile, '(A)') "# vtk DataFile Version 2.0"
      write(mfile, '(A)') "Generated by FEATFlow 2.x"
      write(mfile, '(A)') "ASCII"
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write Dataset header
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! By default, our grid is unstructured
      write(mfile, '(A)') "DATASET UNSTRUCTURED_GRID"
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write vertices
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      write(mfile, '(A,I10,A)') "POINTS", nverts, " double"
      
      ! Do we have 2D or 3D points?
      if (ubound(p_DvertexCoords,1) .eq. 3) then
        ! 3D coordinates
        
        ! Write corner vertices
        do i=1, rexport%p_Rtriangulation%NVT
          write(mfile, '(3E16.7)') p_DvertexCoords(1:3, i)
        end do
        
        ! Write refined vertices
        do i=1, rrefine%nvertices
          write(mfile, '(3E16.7)') p_DvertexRefined(1:3, i)
        end do
        
      else if (ubound(p_DvertexCoords,1) .eq. 2) then
        ! 2D coordinates
        
        ! Write corner vertices
        do i=1, rexport%p_Rtriangulation%NVT
          write(mfile, '(3E16.7)') p_DvertexCoords(1:2, i), 0.0_DP
        end do

        ! Write refined vertices
        do i=1, rrefine%nvertices
          write(mfile, '(3E16.7)') p_DvertexRefined(1:2, i), 0.0_DP
        end do
        
      else
        ! 1D coordinates
        
        ! Write corner vertices
        do i=1, rexport%p_Rtriangulation%NVT
          write(mfile, '(3E16.7)') p_DvertexCoords(1, i), 0.0_DP, 0.0_DP
        end do

        ! Write refined vertices
        do i=1, rrefine%nvertices
          write(mfile, '(3E16.7)') p_DvertexRefined(1, i), 0.0_DP, 0.0_DP
        end do

      end if
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cells
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! First we need to count the number of vertices per cell and the total
      ! number of integers needed for the cell index list
      ncls = 0
      allocate(p_InumVertsPerCell(ncells))
      do j = 1, ncells
      
        ! Count the number of vertices on that element
        do i=1,ubound(p_IverticesAtElement,1)
          if (p_IverticesAtElement(i,j) .eq. 0) exit
        end do
        
        ! Store number of vertices
        p_InumVertsPerCell(j) = i-1
        
        ! Store the number of vertices for this element + 1
        if (bQuadratic) then
          ncls = ncls + 2*i - 1
        else
          ncls = ncls + i
        end if
        
      end do
      
      ! Write cells
      write(mfile, '(A, 2I10)') "CELLS", ncells, ncls

      ! Should we write quadratic cells?
      if (bQuadratic) then
      
        ! Get offset of first edge-midpoint-vertice
        k = rexport%p_rtriangulation%NVT - 1

        do j=1, ncells
        
          select case (p_InumVertsPerCell(j))
          case (2)
            ! Quadratic edge
            write(mfile, '(4I10)') 3, p_IverticesAtElement(1,j)-1, &
                p_IverticesAtElement(2,j)-1, k + j
                
          case (3)
            ! Quadratic triangle
            write(mfile, '(7I10)') 6, p_IverticesAtElement(1,j)-1, &
                p_IverticesAtElement(2,j)-1, p_IverticesAtElement(3,j)-1, &
                k + p_IedgesAtElement(1,j), k + p_IedgesAtElement(2,j), &
                k + p_IedgesAtElement(3,j)
            
          case (4)
            ! Quadratic quadrilateral
            write(mfile, '(9I10)') 8, p_IverticesAtElement(1,j)-1, &
                p_IverticesAtElement(2,j)-1, p_IverticesAtElement(3,j)-1, &
                p_IverticesAtElement(4,j)-1, k + p_IedgesAtElement(1,j), &
                k + p_IedgesAtElement(2,j), k + p_IedgesAtElement(3,j), &
                k + p_IedgesAtElement(4,j)
            
          case DEFAULT
            call output_line ('Invalid element!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeVTK')
          end select
          
        end do

        ! Write cell types
        write(mfile, '(A, I10)') "CELL_TYPES", ncells
        do j=1, ncells
          select case (p_InumVertsPerCell(j))
          case (2)
            ! quadratic edge
            write(mfile, '(I4)') VTK_QUADRATIC_EDGE 
             
          case (3)
            ! quadratic triangle
            write(mfile, '(I4)') VTK_QUADRATIC_TRIANGLE
            
          case (4)
            ! quadratic quadrilateral
            write(mfile, '(I4)') VTK_QUADRATIC_QUAD
            
          end select
        end do

      else

        do j=1, ncells
        
          select case (p_InumVertsPerCell(j))
          case (2)
            ! edge
            write(mfile, '(3I10)') 2, p_IverticesAtElement(1,j)-1, &
                p_IverticesAtElement(2,j)-1

          case (3)
            ! triangle
            write(mfile, '(4I10)') 3, p_IverticesAtElement(1,j)-1, &
                p_IverticesAtElement(2,j)-1, p_IverticesAtElement(3,j)-1

          case (4)
            ! quadrilateral
            write(mfile, '(5I10)') 4, p_IverticesAtElement(1,j)-1, &
                p_IverticesAtElement(2,j)-1, p_IverticesAtElement(3,j)-1, &
                p_IverticesAtElement(4,j)-1
          
          case (8)
            ! hexahedron
            write(mfile, '(9I10)') 8, p_IverticesAtElement(1,j)-1, &
                p_IverticesAtElement(2,j)-1, p_IverticesAtElement(3,j)-1, &
                p_IverticesAtElement(4,j)-1, p_IverticesAtElement(5,j)-1, &
                p_IverticesAtElement(6,j)-1, p_IverticesAtElement(7,j)-1, &
                p_IverticesAtElement(8,j)-1
                
          case DEFAULT
            call output_line ('Invalid element!',&
                OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeVTK')
          end select
          
        end do

        ! Write cell types
        write(mfile, '(A, I10)') "CELL_TYPES", ncells
        do j=1, ncells
          select case (p_InumVertsPerCell(j))
          case (2)
            ! Edge
            write(mfile, '(I4)') VTK_LINE
            
          case (3)
            ! Triangle
            write(mfile, '(I4)') VTK_TRIANGLE
            
          case (4)
            ! Quadrilateral
            write(mfile, '(I4)') VTK_QUAD
          
          case (8)
            ! Hexahedron
            write(mfile, '(I4)') VTK_HEXAHEDRON
            
          end select
        end do

      end if
      
      ! Now we do not need the vertice counts anymore
      deallocate(p_InumVertsPerCell)
      
      ! We also do not need the mesh refinement anymore, so destroy it
      if (rrefine%h_DvertexCoords .ne. ST_NOHANDLE) then
        call storage_free(rrefine%h_DvertexCoords)
      end if
      if (rrefine%h_IverticesAtElement .ne. ST_NOHANDLE) then
        call storage_free(rrefine%h_IverticesAtElement)
      end if
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write vertex variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Go through all variables
      if (num_ndata .gt. 0) then
      
        write(mfile, '(A,I10)') "POINT_DATA", nverts
      
        ! Loop through all variables
        do j=1, rexport%nvariables
        
          ! Is it vertex- or element-based?
          if (rexport%p_IvariableBase(j) .ne. UCD_BASE_VERTEX) cycle

          ! Is the value scalar?
          if ((rexport%p_IvariableSpec(j) .ne. UCD_VAR_STANDARD) .and. &
              (.not. bVec2Sc)) cycle
            
          ! Print some stuff
          write(mfile, '(A,A,A)') "SCALARS ", &
              sys_charreplace(trim(rexport%p_SvariableNames(j)),' ', '_'), &
                " double 1"
          write(mfile, '(A)') "LOOKUP_TABLE default"
          
          ! Go for the data
          call storage_getbase_double (rexport%p_Hvariables(j), p_Ddata)
          
          ! Write the data
          do i=1, ubound(p_Ddata,1)
            write(mfile, sdl) p_Ddata(i)
          end do
            
        end do
        
        ! Should we write vectors?
        if (.not. bVec2Sc) then
        
          ! Loop through all vectors
          do j=1, rexport%nvectors
          
            ! Is this a vertex-based vector?
            if(rexport%p_Ivectors(1,j) .ne. 1) cycle
        
            p_Dx => null()
            p_Dy => null()
            p_Dz => null()

            ! Go for the X, Y and Z coordinates
            i = rexport%p_Ivectors(2,j)
            call storage_getbase_double (rexport%p_Hvariables(i), p_Dx)
            i = rexport%p_Ivectors(3,j)
            if (i .ne. 0) then
              call storage_getbase_double (rexport%p_Hvariables(i), p_Dy)
            end if
            i = rexport%p_Ivectors(4,j)
            if (i .ne. 0) then
              call storage_getbase_double (rexport%p_Hvariables(i), p_Dz)
            end if
            
            ! Make sure we have at least the X-coordinate
            if (.not. associated(p_Dx)) then
              call output_line ('Error: Variable vector '//&
                  trim(sys_siL(j,10))//' does not have X-coordinates!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'ucd_writeVTK')
              
              ! Try next vector
              cycle
            end if
            
            ! Print the stuff
            write(mfile,'(A,A,A)') "VECTORS ", &
                sys_charreplace(trim(rexport%p_SvarVecNames(j)), ' ', '_'), &
                  " double"
            
            ! Write the vector
            do i = 1, ubound(p_Dx, 1)
            
              ! Write X-coord
              write(mfile, sdl, ADVANCE='NO') p_Dx(i)
              
              ! Write Y-coord
              if (associated(p_Dy)) then
                write(mfile, sdl, ADVANCE='NO') p_Dy(i)
              else
                write(mfile, sdl, ADVANCE='NO') 0.0_DP
              end if
              
              ! Write Z-coord
              if (associated(p_Dz)) then
                write(mfile, sdl, ADVANCE='NO') p_Dz(i)
              else
                write(mfile, sdl, ADVANCE='NO') 0.0_DP
              end if
              
              ! Write a new line
              write(mfile, '(A)') ""
              
            end do
            
          end do
           
        end if
        
      end if
      
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Write cell variables
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Go through all variables
      if (num_cdata .gt. 0) then
      
        write(mfile, '(A,I10)') "CELL_DATA", ncells
      
        ! Loop through all variables
        do j=1, rexport%nvariables
        
          ! Is it vertex- or element-based?
          if (rexport%p_IvariableBase(j) .ne. UCD_BASE_ELEMENT) cycle
          
          ! Print some stuff
          write(mfile,'(A,A,A)') "SCALARS ", &
              sys_charreplace(trim(rexport%p_SvariableNames(j)),' ','_'),&
                " double 1"
          write(mfile, '(A)') "LOOKUP_TABLE default"
          
          ! Go for the data
          call storage_getbase_double (rexport%p_Hvariables(j), p_Ddata)
          
          ! Write variable data
          do i=1, ubound(p_Ddata,1)
            write(mfile, sdl) p_Ddata(i)
          end do
          
        end do
        
      end if
      
      ! That is it
stop
    end subroutine
    
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_moreVariables (rexport)

!<description>
  ! (Re)allocates memory for variables. Increases the number of available
  ! variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    integer(I32), dimension(:), pointer :: p_IvariableSpec
    integer, dimension(:), pointer :: p_Hvariables,p_IvariableBase
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_SvariableNames
    integer :: nsize
    
    nsize = 0
    if (associated(rexport%p_IvariableSpec)) &
      nsize = size(rexport%p_IvariableSpec)

    allocate(p_IvariableSpec(nsize+16))
    if (associated(rexport%p_IvariableSpec)) then
      p_IvariableSpec(1:nsize) = rexport%p_IvariableSpec
      deallocate(rexport%p_IvariableSpec)
    end if
    rexport%p_IvariableSpec => p_IvariableSpec
      
    allocate(p_IvariableBase(nsize+16))
    if (associated(rexport%p_IvariableBase)) then
      p_IvariableBase(1:nsize) = rexport%p_IvariableBase
      deallocate(rexport%p_IvariableBase)
    end if
    rexport%p_IvariableBase => p_IvariableBase

    allocate(p_Hvariables(nsize+16))
    if (associated(rexport%p_Hvariables)) then
      p_Hvariables(1:nsize) = rexport%p_Hvariables
      deallocate(rexport%p_Hvariables)
    end if
    rexport%p_Hvariables => p_Hvariables

    allocate(p_SvariableNames(nsize+16))
    p_SvariableNames = ''
    if (associated(rexport%p_SvariableNames)) then
      p_SvariableNames(1:nsize) = rexport%p_SvariableNames
      deallocate(rexport%p_SvariableNames)
    end if
    rexport%p_SvariableNames => p_SvariableNames

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_moreVectors (rexport)

!<description>
  ! (Re)allocates memory for vector variables. Increases the number of
  ! available variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    integer, dimension(:,:), pointer :: p_Ivectors
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_SvarVecNames
    integer :: nsize
    
    nsize = 0
    if (associated(rexport%p_Ivectors)) &
      nsize = size(rexport%p_Ivectors)

    allocate(p_Ivectors(4,nsize+16))
    if (associated(rexport%p_Ivectors)) then
      p_Ivectors(:,1:nsize) = rexport%p_Ivectors(:,:)
      deallocate(rexport%p_Ivectors)
    end if
    rexport%p_Ivectors => p_Ivectors
      
    allocate(p_SvarVecNames(nsize+16))
    p_SvarVecNames = ''
    if (associated(rexport%p_SvarVecNames)) then
      p_SvarVecNames(1:nsize) = rexport%p_SvarVecNames
      deallocate(rexport%p_SvarVecNames)
    end if
    rexport%p_SvarVecNames => p_SvarVecNames

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_morePolygons (rexport)

!<description>
  ! (Re)allocates memory for polygons. Increases the number of available
  ! variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    integer, dimension(:), pointer :: p_Hpolygons
    integer :: nsize

    nsize = 0
    if (associated(rexport%p_Hpolygons)) &
      nsize = size(rexport%p_Hpolygons)

    allocate(p_Hpolygons(nsize+16))
    if (rexport%hpolygonMaterial .eq. ST_NOHANDLE) then
      call storage_new ('ucd_morePolygons','hpolygonMaterial',&
          nsize+16,ST_INT,&
          rexport%hpolygonMaterial,ST_NEWBLOCK_ZERO)
    else
      call storage_realloc ('ucd_morePolygons', nsize+16, &
          rexport%hpolygonMaterial, ST_NEWBLOCK_ZERO)
    end if
    
    if (associated(rexport%p_Hpolygons)) then
      p_Hpolygons(1:size(rexport%p_Hpolygons)) = rexport%p_Hpolygons
      deallocate(rexport%p_Hpolygons)
    end if
    rexport%p_Hpolygons => p_Hpolygons

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_moreTracerVariables (rexport)

!<description>
  ! (Re)allocates memory for tracer variables. Increases the number of 
  ! available tracer variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    integer, dimension(:), pointer :: p_Hvariables
    character(LEN=SYS_NAMELEN), dimension(:), pointer :: p_StracerVariableNames
    integer :: nsize

    nsize = 0
    if (associated(rexport%p_HtracerVariables)) &
      nsize = size(rexport%p_HtracerVariables)

    allocate(p_Hvariables(nsize+16))
    if (associated(rexport%p_HtracerVariables)) then
      p_Hvariables(1:nsize) = rexport%p_HtracerVariables
      deallocate(rexport%p_HtracerVariables)
    end if
    rexport%p_HtracerVariables => p_Hvariables
    
    allocate(p_StracerVariableNames(nsize+16))
    if (associated(rexport%p_StracerVariableNames)) then
      p_StracerVariableNames(1:nsize) = rexport%p_StracerVariableNames
      deallocate(rexport%p_StracerVariableNames)
    end if
    rexport%p_StracerVariableNames => p_StracerVariableNames

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVariableVertexBased1 (rexport, sname, &
      DdataVert, DdataMid, DdataElem, cvarSpec)

!<description>
  ! Adds variable data to the ouput file identified by rexport, based
  ! on vertices.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the variable.
  character(LEN=*), intent(in) :: sname
  
  ! DdataVert(I) is the value of the variable in vertex I of the triangulation.
  real(DP), dimension(:), intent(in) :: DdataVert

  ! OPTIONAL: DdataMid(I) is the value of the variable in edge midpoint I of 
  ! the triangulation. Must be specified if DdataElem is specified!
  real(DP), dimension(:), intent(in), optional :: DdataMid

  ! OPTIONAL: DdataElem(I) is the value of the variable in element midpoint I of 
  ! the triangulation.
  real(DP), dimension(:), intent(in), optional :: DdataElem

  ! OPTIONAL: Specification bitfield for the variable. A combination of the 
  ! UCD_VAR_xxxx flags for special-type variables (like x-/y-velocity).
  ! If not specified, UCD_VAR_STANDARD is used.
  integer(I32), intent(in), optional :: cvarSpec
!</input>

!</subroutine>

    if (present(cvarSpec)) then
      call ucd_addVariableVertexBased (rexport, sname, cvarSpec, &
          DdataVert, DdataMid, DdataElem)
    else
      call ucd_addVariableVertexBased (rexport, sname, UCD_VAR_STANDARD, &
          DdataVert, DdataMid, DdataElem)
    end if

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVariableVertexBased2 (rexport, sname, cvarSpec, &
      DdataVert, DdataMid, DdataElem)

!<description>
  ! Adds variable data to the ouput file identified by rexport, based
  ! on vertices.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the variable.
  character(LEN=*), intent(in) :: sname
  
  ! Specification bitfield for the variable. A combination of the 
  ! UCD_VAR_xxxx flags for special-type variables (like x-/y-velocity).
  ! Standard value=UCD_VAR_STANDARD.
  integer(I32), intent(in) :: cvarSpec
  
  ! DdataVert(I) is the value of the variable in vertex I of the triangulation.
  real(DP), dimension(:), intent(in) :: DdataVert

  ! OPTIONAL: DdataMid(I) is the value of the variable in edge midpoint I of 
  ! the triangulation. Must be specified if DdataElem is specified!
  real(DP), dimension(:), intent(in), optional :: DdataMid

  ! OPTIONAL: DdataElem(I) is the value of the variable in element midpoint I of 
  ! the triangulation.
  real(DP), dimension(:), intent(in), optional :: DdataElem
!</input>

!</subroutine>

    ! local varables
    real(DP), dimension(:), pointer :: p_Ddata

    if (present(DdataElem) .and. .not. present(DdataMid)) then
      call output_line ('Error in the parameters!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableVertexBased2')
      call sys_halt()
    end if
    
    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableVertexBased2')
      call sys_halt()
    end if
    
    ! Create a new variable. If necessary, increase the size of the buffer.
    if (.not. associated(rexport%p_Hvariables)) then
      call ucd_moreVariables(rexport)
    end if
    
    if (rexport%nvariables .ge. size(rexport%p_Hvariables)) then
      call ucd_moreVariables(rexport)
    end if
    
    rexport%nvariables = rexport%nvariables+1
    
    rexport%p_IvariableSpec(rexport%nvariables) = cvarSpec
    
    rexport%p_IvariableBase(rexport%nvariables) = UCD_BASE_VERTEX
    
    rexport%p_SvariableNames(rexport%nvariables) = sname
    
    ! Allocate a new vector for the data
    call storage_new ('ucd_addVariableVertexBased2','hvariable',&
        rexport%nvertices,ST_DOUBLE,&
        rexport%p_Hvariables(rexport%nvariables),&
        ST_NEWBLOCK_ZERO)
    
    ! Copy the vertex data into that vector
    call storage_getbase_double (rexport%p_Hvariables(rexport%nvariables),p_Ddata)
    call lalg_copyVectorDble(DdataVert(1:rexport%p_rtriangulation%NVT), &
        p_Ddata(1:rexport%p_rtriangulation%NVT))
    
    ! Copy edge midpoint data if available
    if ((iand(rexport%cflags,UCD_FLAG_BULBQUADRATIC) .ne. 0) .or. &
        (iand(rexport%cflags,UCD_FLAG_USEEDGEMIDPOINTS) .ne. 0) .or. &
        (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      if (present(DdataMid)) then
        call lalg_copyVectorDble( &
            DdataMid(1:rexport%p_rtriangulation%NMT), &
            p_Ddata(rexport%p_rtriangulation%NVT+1:rexport%p_rtriangulation%NVT+ &
                                                   rexport%p_rtriangulation%NMT))
      else if (iand(rexport%cflags,UCD_FLAG_IGNOREDEADNODES) .eq. 0) then
        call output_line ('Warning. No edge midpoint data available!',&
            OU_CLASS_WARNING,OU_MODE_STD,'ucd_addVariableVertexBased2')
      end if
    end if
    
    ! Copy element midpoint data if available
    if ((iand(rexport%cflags,UCD_FLAG_USEELEMENTMIDPOINTS) .ne. 0) .or. &
        (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0)) then
      if (present(DdataElem)) then
        call lalg_copyVectorDble( &
            DdataElem(1:rexport%p_rtriangulation%NEL), &
            p_Ddata(rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+1: &
                    rexport%p_rtriangulation%NVT+rexport%p_rtriangulation%NMT+ &
                    rexport%p_rtriangulation%NEL))
      else if (iand(rexport%cflags,UCD_FLAG_IGNOREDEADNODES) .eq. 0) then
        call output_line ('Warning. No element midpoint data available!',&
            OU_CLASS_WARNING,OU_MODE_STD,'ucd_addVariableVertexBased2')
      end if
    end if    

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVarVertBasedVec1 (rexport, sname, DdataVert_X, &
      DdataVert_Y, DdataVert_Z, DdataMid_X, DdataMid_Y, DdataMid_Z, &
      DdataElem_X, DdataElem_Y, DdataElem_Z, cvarSpec)
  
!<description>
  ! Adds variable vector data to the ouput file identified by rexport, based
  ! on vertices.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the vector.
  character(LEN=*), intent(in) :: sname
  
  ! The variable for the X-component of the velocities.
  real(DP), dimension(:), intent(in) :: DdataVert_X
  
  ! OPTIONAL: The data array for the Y-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: DdataVert_Y

  ! OPTIONAL: The data array for the Z-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: DdataVert_Z

  ! OPTIONAL: DdataMid_X/Y/Z(I) is the X/Y/Z-coordinate of the variable in
  ! edge midpoint I of the triangulation. Must be specified if DdataElem
  ! is specified!
  real(DP), dimension(:), intent(in), optional :: DdataMid_X
  real(DP), dimension(:), intent(in), optional :: DdataMid_Y
  real(DP), dimension(:), intent(in), optional :: DdataMid_Z

  ! OPTIONAL: DdataElem_X/Y/Z(I) is the X/Y/Z-coordinate of the variable in
  ! element midpoint I of the triangulation.
  real(DP), dimension(:), intent(in), optional :: DdataElem_X
  real(DP), dimension(:), intent(in), optional :: DdataElem_Y
  real(DP), dimension(:), intent(in), optional :: DdataElem_Z

  ! OPTIONAL: Specification bitfield for the variable. A combination
  ! of the UCD_VAR_xxxx flags.
  ! If not specified, UCD_VAR_STANDARD is used.
  integer(I32), intent(in), optional :: cvarSpec
!</input>

!</subroutine>

    if (present(cvarSpec)) then
      call ucd_addVarVertBasedVec2(rexport, sname, cvarSpec,&
          DdataVert_X, DdataVert_Y, DdataVert_Z,&
          DdataMid_X, DdataMid_Y, DdataMid_Z, &
          DdataElem_X, DdataElem_Y, DdataElem_Z)
    else
      call ucd_addVarVertBasedVec2(rexport, sname, UCD_VAR_STANDARD,&
          DdataVert_X, DdataVert_Y, DdataVert_Z,&
          DdataMid_X, DdataMid_Y, DdataMid_Z, &
          DdataElem_X, DdataElem_Y, DdataElem_Z)
    end if

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVarVertBasedVec2 (rexport, sname, cvarSpec,&
      DdataVert_X, DdataVert_Y, DdataVert_Z, DdataMid_X,&
      DdataMid_Y, DdataMid_Z, DdataElem_X, DdataElem_Y, DdataElem_Z)

!<description>
  ! Adds variable vector data to the ouput file identified by rexport, based
  ! on vertices.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the vector.
  character(LEN=*), intent(in) :: sname
  
  ! Specification bitfield for the variable.
  ! Standard value=UCD_VAR_STANDARD.
  integer(I32), intent(in) :: cvarSpec

  ! The variable for the X-component of the velocities.
  real(DP), dimension(:), intent(in) :: DdataVert_X
  
  ! OPTIONAL: The data array for the Y-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: DdataVert_Y

  ! OPTIONAL: The data array for the Z-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: DdataVert_Z

  ! OPTIONAL: DdataMid_X/Y/Z(I) is the X/Y/Z-coordinate of the variable in
  ! edge midpoint I of the triangulation. Must be specified if DdataElem
  ! is specified!
  real(DP), dimension(:), intent(in), optional :: DdataMid_X
  real(DP), dimension(:), intent(in), optional :: DdataMid_Y
  real(DP), dimension(:), intent(in), optional :: DdataMid_Z

  ! OPTIONAL: DdataElem_X/Y/Z(I) is the X/Y/Z-coordinate of the variable in
  ! element midpoint I of the triangulation.
  real(DP), dimension(:), intent(in), optional :: DdataElem_X
  real(DP), dimension(:), intent(in), optional :: DdataElem_Y
  real(DP), dimension(:), intent(in), optional :: DdataElem_Z
!</input>

!</subroutine>

    ! local variables
    integer(i32) :: ivarSpec

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVarVertBasedVec')
      call sys_halt()
    end if
    
    ! Let us create a new vector first
    if (.not. associated(rexport%p_Ivectors)) then
      call ucd_moreVectors(rexport)
    endif
    
    if (rexport%nvectors .ge. size(rexport%p_Ivectors, 2)) then
      call ucd_moreVectors(rexport)
    end if
    
    ! Increment vector count and store its name
    rexport%nvectors = rexport%nvectors + 1
    rexport%p_SvarVecNames(rexport%nvectors) = sname
    
    ! This vector is vertex-based
    rexport%p_Ivectors(1,rexport%nvectors) = 1
    
    ! There must be an X-component, so add it first
    ivarSpec = ior(cvarSpec, UCD_VAR_XVECTORCOMP)
    call ucd_addVariableVertexBased(rexport, trim(sname)//'_X', &
        ivarSpec, DdataVert_X, DdataMid_X, dDataElem_X)
    
    rexport%p_Ivectors(2,rexport%nvectors) = rexport%nvariables
    
    ! Is there an Y-component?
    if (present(DdataVert_Y)) then

      ivarSpec = ior(cvarSpec, UCD_VAR_YVECTORCOMP)
      call ucd_addVariableVertexBased(rexport, trim(sname)//'_Y', &
          ivarSpec, DdataVert_Y, DdataMid_Y, dDataElem_Y)
      
      rexport%p_Ivectors(3,rexport%nvectors) = rexport%nvariables
   
    else
      rexport%p_Ivectors(3,rexport%nvectors) = 0
    end if
    
    ! Is there an Z-component?
    if (present(DdataVert_Z)) then

      ivarSpec = ior(cvarSpec, UCD_VAR_ZVECTORCOMP)
      call ucd_addVariableVertexBased(rexport, trim(sname)//'_Z', &
          ivarSpec, DdataVert_Z, DdataMid_Z, dDataElem_Z)
      
      rexport%p_Ivectors(4,rexport%nvectors) = rexport%nvariables
   
    else
      rexport%p_Ivectors(4,rexport%nvectors) = 0
    end if
    
    ! That is it
  
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVarElemBasedVec1 (rexport, sname,&
      Ddata_X, Ddata_Y, Ddata_Z, cvarSpec)

!<description>
  ! Adds variable vector data to the ouput file identified by rexport, based
  ! on elements.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the vector.
  character(LEN=*), intent(in) :: sname
  
  ! The variable for the X-component of the velocities.
  real(DP), dimension(:), intent(in) :: Ddata_X
  
  ! OPTIONAL: The data array for the Y-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: Ddata_Y

  ! OPTIONAL: The data array for the Z-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: Ddata_Z

  ! OPTIONAL: Specification bitfield for the variable. A combination
  ! of the UCD_VAR_xxxx flags.
  ! If not specified, UCD_VAR_STANDARD is used.
  integer(I32), intent(in), optional :: cvarSpec
!</input>

!</subroutine>

    if (present(cvarSpec)) then
      call ucd_addVarElemBasedVec2 (rexport, sname, cvarSpec,&
          Ddata_X, Ddata_Y, Ddata_Z)
    else
      call ucd_addVarElemBasedVec2 (rexport, sname, UCD_VAR_STANDARD,&
          Ddata_X, Ddata_Y, Ddata_Z)
    end if

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVarElemBasedVec2 (rexport, sname, cvarSpec,&
      Ddata_X, Ddata_Y, Ddata_Z)

!<description>
  ! Adds variable vector data to the ouput file identified by rexport, based
  ! on elements.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the vector.
  character(LEN=*), intent(in) :: sname
  
  ! The variable for the X-component of the velocities.
  real(DP), dimension(:), intent(in) :: Ddata_X
  
  ! OPTIONAL: The data array for the Y-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: Ddata_Y

  ! OPTIONAL: The data array for the Z-component of the velocities.
  real(DP), dimension(:), intent(in), optional :: Ddata_Z

  ! Specification bitfield for the variable.
  ! Standard value=UCD_VAR_STANDARD.
  integer(I32), intent(in) :: cvarSpec
!</input>

!</subroutine>

    ! local variable
    integer(i32) :: ivarSpec

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVarElemBasedVec2')
      call sys_halt()
    end if
    
    ! Let us create a new vector first
    if (.not. associated(rexport%p_Ivectors)) then
      call ucd_moreVectors(rexport)
    endif
    
    if (rexport%nvectors .ge. size(rexport%p_Ivectors, 2)) then
      call ucd_moreVectors(rexport)
    end if
    
    ! Increment vector count and store its name
    rexport%nvectors = rexport%nvectors + 1
    rexport%p_SvarVecNames(rexport%nvectors) = sname
    
    ! This vector is vertex-based
    rexport%p_Ivectors(1,rexport%nvectors) = 1
    
    ! There must be an X-component, so add it first
    ivarSpec = ior(cvarSpec, UCD_VAR_XVECTORCOMP)
    call ucd_addVariableElementBased(rexport, trim(sname)//'_X', &
        ivarSpec, Ddata_X)
    
    rexport%p_Ivectors(2,rexport%nvectors) = rexport%nvariables
    
    ! Is there an Y-component?
    if (present(Ddata_Y)) then

      ivarSpec = ior(cvarSpec, UCD_VAR_YVECTORCOMP)
      call ucd_addVariableElementBased(rexport, trim(sname)//'_Y', &
          ivarSpec, Ddata_Y)
      
      rexport%p_Ivectors(3,rexport%nvectors) = rexport%nvariables
   
    else
      rexport%p_Ivectors(3,rexport%nvectors) = 0
    end if
    
    ! Is there an Z-component?
    if (present(Ddata_Z)) then

      ivarSpec = ior(cvarSpec, UCD_VAR_ZVECTORCOMP)
      call ucd_addVariableElementBased(rexport, trim(sname)//'_Z', &
          ivarSpec, Ddata_Z)
      
      rexport%p_Ivectors(4,rexport%nvectors) = rexport%nvariables
   
    else
      rexport%p_Ivectors(4,rexport%nvectors) = 0
    end if
    
    ! That is it
  
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVariableElementBased1 (rexport, sname, Ddata, cvarSpec)

!<description>
  ! Adds variable data to the ouput file identified by rexport, based on elements.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the variable.
  character(LEN=*), intent(in) :: sname
  
  ! DdataVert(I) os the value of the variable in element I of the triangulation.
  real(DP), dimension(:), intent(in) :: Ddata

  ! OPTIONAL: Specification bitfield for the variable. A combination
  ! of the UCD_VAR_xxxx flags.
  ! If not specified, UCD_VAR_STANDARD is used.
  integer(I32), intent(in), optional :: cvarSpec
!</input>

!</subroutine>

    if (present(cvarSpec)) then
      call ucd_addVariableElementBased2 (rexport, sname, cvarSpec, Ddata)
    else
      call ucd_addVariableElementBased2 (rexport, sname, UCD_VAR_STANDARD, Ddata)
    end if

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addVariableElementBased2 (rexport, sname, cvarSpec, Ddata)

!<description>
  ! Adds variable data to the ouput file identified by rexport, based on elements.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the variable.
  character(LEN=*), intent(in) :: sname
  
  ! Specification bitfield for the variable. A combination of the 
  ! UCD_VAR_xxxx flags.
  integer(I32), intent(in) :: cvarSpec
  
  ! DdataVert(I) os the value of the variable in element I of the triangulation.
  real(DP), dimension(:), intent(in) :: Ddata
!</input>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP) :: dv
    integer :: iel

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addVariableElementBased2')
      call sys_halt()
    end if
    
    ! Create a new variable. If necessary, increase the size of the buffer.
    if (.not. associated(rexport%p_Hvariables)) then
      call ucd_moreVariables(rexport)
    end if

    ! Create a new variable. If necessary, increase the size of the buffer.
    if (rexport%nvariables .ge. size(rexport%p_Hvariables)) then
      call ucd_moreVariables(rexport)
    end if
    rexport%nvariables = rexport%nvariables+1
    rexport%p_IvariableSpec(rexport%nvariables) = cvarSpec
    rexport%p_IvariableBase(rexport%nvariables) = UCD_BASE_ELEMENT
    rexport%p_SvariableNames(rexport%nvariables) = sname
    
    ! Allocate a new vector for the data
    call storage_new ('ucd_addVariableVertexBased2','hvariable',&
        rexport%ncells,ST_DOUBLE,&
        rexport%p_Hvariables(rexport%nvariables),&
        ST_NEWBLOCK_ZERO)
        
    ! Copy the element data into that vector
    call storage_getbase_double (rexport%p_Hvariables(rexport%nvariables),p_Ddata)
    call lalg_copyVectorDble(Ddata(1:rexport%p_rtriangulation%NEL), &
        p_Ddata(1:rexport%p_rtriangulation%NEL))
        
    ! On a once refined mesh, generate the data on the sub-cells by
    ! constant interpolation of the data in the large cell.
        
    ! Copy edge midpoint data if available
    if (iand(rexport%cflags,UCD_FLAG_ONCEREFINED) .ne. 0) then
    
      if (rexport%p_rtriangulation%NDIM .eq. NDIM2D) then
        ! Implicitly use 2-level ordering to get the numbers of the sub-elements
        do iel=1,rexport%p_rtriangulation%NEL
          dv = p_Ddata(iel)
          p_Ddata(rexport%p_rtriangulation%NEL+3*(iel-1)+1) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+3*(iel-1)+2) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+3*(iel-1)+3) = dv
        end do
      else
        do iel=1,rexport%p_rtriangulation%NEL
          dv = p_Ddata(iel)
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+1) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+2) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+3) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+4) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+5) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+6) = dv
          p_Ddata(rexport%p_rtriangulation%NEL+7*(iel-1)+7) = dv
        end do
      end if
      
    end if
    
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_moreSurfTri(rexport)

!<description>
  ! (Re)allocates memory for surface triangulations. Increases the number of available
  ! variables in the export structure.
!</description>

!<inputoutput>
  ! The ucd export structure which is to be tackled.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    integer :: nsize
    integer, dimension(:), pointer :: p_HsurfTris
    integer, dimension(:), pointer :: p_HTriangles
    integer, dimension(:), pointer :: p_HsurfData
    

!    ! A pointer to a list of handles to surface triangulation vertices.
!    ! p_HsurfTris(I) points to a list of (X,Y,Z) tags
!    ! containing the points of a the surface triangulation
!    integer, dimension(:), pointer :: p_HsurfTris => null()
!
!    ! A pointer to a list of handles to surface triangulation data.
!    ! p_HTriangles(I) is a handle an integer array that describes
!    ! the connectivity of the vertices in p_HsurfTris(I)
!    integer, dimension(:), pointer :: p_HTriangles => null()
!    
!    ! A pointer to a list of handles to surface triangulation data.
!    ! p_HsurfData(I) is a handle an integer array that describes
!    ! how many vertices and triangles there are in p_HsurfTris(I)
!    integer, dimension(:), pointer :: p_HsurfData => null()
    



    nsize = 0
    if (associated(rexport%p_HsurfTris)) &
      nsize = size(rexport%p_HsurfTris)

    allocate(p_HsurfTris(nsize+16))
    allocate(p_HTriangles(nsize+16))
    allocate(p_HsurfData(nsize+16))
    
    if (associated(rexport%p_HsurfTris)) then
      p_HsurfTris(1:size(rexport%p_HsurfTris)) = rexport%p_HsurfTris
      p_HTriangles(1:size(rexport%p_HTriangles)) = rexport%p_HTriangles
      p_HsurfData(1:size(rexport%p_HsurfData)) = rexport%p_HsurfData
      deallocate(rexport%p_HsurfTris)
      deallocate(rexport%p_HTriangles)
      deallocate(rexport%p_HsurfData)
    end if
    rexport%p_HsurfTris => p_HsurfTris
    rexport%p_HTriangles => p_HTriangles
    rexport%p_HsurfData => p_HsurfData

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addSurfTri(rexport,DpolygonCoords,Itriangles,Idata)

!<description>
  ! add a surface triangulation to the output file
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! A list of 3D (X,Y,Z) coordinates specifying the
  ! points which form the surface triangulation
  real(DP), dimension(:,:), intent(in) :: DpolygonCoords
  ! array describing triangle connectivity
  integer, dimension(:,:), intent(in)  :: Itriangles
  ! Idata(1)=number vertices in triangulation
  ! Idata(2)=number triangles
  integer, dimension(2), intent(in)    :: Idata
!</input>

!</subroutine>
    ! local varables
    integer, dimension(2) :: Ilength
    integer :: isize
    real(DP), dimension(:,:), pointer :: p_Ddata
    integer, dimension(:,:), pointer :: p_ItriangleData
    integer, dimension(:), pointer :: p_Idata
    
    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addSurfTri')
      call sys_halt()
    end if
    
    ! Create a new variable. If necessary, increase the size of the buffers.
    if (.not. associated(rexport%p_HsurfTris)) then
      call ucd_moreSurfTri(rexport)
    end if
    
    if (rexport%nsurfTri .ge. size(rexport%p_HsurfTris)) then
      call ucd_moreSurfTri(rexport)
    end if
    rexport%nsurfTri = rexport%nsurfTri+1
    
    ! Allocate a new vector for the data
    Ilength = ubound(DpolygonCoords)
    call storage_new ('ucd_addSurfTri','p_HsurfTris',&
        Ilength,ST_DOUBLE,rexport%p_HsurfTris(rexport%nsurfTri),&
        ST_NEWBLOCK_NOINIT)

    ! Allocate a new vector for the data
    Ilength = ubound(Itriangles)
    call storage_new ('ucd_addSurfTri','p_HTriangles',&
        Ilength,ST_INT,rexport%p_HTriangles(rexport%nsurfTri),&
        ST_NEWBLOCK_NOINIT)

    isize=2
    call storage_new ('ucd_addSurfTri','p_HsurfData',&
        isize,ST_INT,rexport%p_HsurfData(rexport%nsurfTri),&
        ST_NEWBLOCK_NOINIT)

        
    ! Copy the coordinate data into that vector
    call storage_getbase_double2d (rexport%p_HsurfTris(rexport%nsurfTri),p_Ddata)
    p_Ddata = DpolygonCoords

    ! Copy the coordinate data into that vector
    call storage_getbase_int2d (rexport%p_HTriangles(rexport%nsurfTri),p_ItriangleData)
    p_ItriangleData = Itriangles

    ! Copy the coordinate data into that vector
    call storage_getbase_int(rexport%p_HsurfData(rexport%nsurfTri),p_Idata)
    p_Idata = Idata

  end subroutine ucd_addSurfTri
  
  !************************************************************************  
  
!<subroutine>

  subroutine ucd_addPolygon (rexport,DpolygonCoords,ipolygonMaterial)

!<description>
  ! Adds a polygon to the output file. A polygon is a set of lines.
  ! DpolygonCoords is a 2D or 3D array with (X,Y)-/(X,Y,Z)-coordinates.
  ! Line i is specified by the start coordinates DpolygonCoords(i)
  ! and the end coordinates DpolygonCoords(i+1).
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! A list of 2D (X,Y) or 3D (X,Y,Z) coordinates specifying the
  ! points which form the polygon.
  real(DP), dimension(:,:), intent(in) :: DpolygonCoords
  
  ! OPTIONAL: A material identifier for the polygon. Whether polygons can have
  ! a material associated depends on the output file format.
  ! If not specified, a default material is assumed.
  integer, intent(in), optional :: ipolygonMaterial
!</input>

!</subroutine>

    ! local varables
    integer, dimension(2) :: Ilength
    real(DP), dimension(:,:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_Idata

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addPolygon')
      call sys_halt()
    end if
    
    ! Create a new variable. If necessary, increase the size of the buffers.
    if (.not. associated(rexport%p_Hpolygons)) then
      call ucd_morePolygons(rexport)
    end if
    
    if (rexport%npolygons .ge. size(rexport%p_Hpolygons)) then
      call ucd_morePolygons(rexport)
    end if
    rexport%npolygons = rexport%npolygons+1
    
    ! Allocate a new vector for the data
    Ilength = ubound(DpolygonCoords)
    call storage_new ('ucd_addPolygon','hpolygon',&
        Ilength,ST_DOUBLE,rexport%p_Hpolygons(rexport%npolygons),&
        ST_NEWBLOCK_NOINIT)
        
    ! Copy the coordinate data into that vector
    call storage_getbase_double2d (rexport%p_Hpolygons(rexport%npolygons),p_Ddata)
    p_Ddata = DpolygonCoords
    
    ! Copy the material
    call storage_getbase_int (rexport%hpolygonMaterial,p_Idata)
    if (present(ipolygonMaterial)) then
      p_Idata(rexport%npolygons) = ipolygonMaterial
    else
      ! Default material 0
      p_Idata(rexport%npolygons) = 0
    end if

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_setTracers (rexport,DtracerCoordinates)

!<description>
  ! This routine adds tracer information to the export structure rexport.
  ! If there is already any tracer data attached to rexport, it will be
  ! deleted.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! A list of 2D (X,Y) or 3D (X,Y,Z) coordinates specifying the
  ! position of all the tracers.
  real(DP), dimension(:,:), intent(in) :: DtracerCoordinates
!</input>

!</subroutine>

    ! local varables
    integer, dimension(2) :: Ilength
    real(DP), dimension(:,:), pointer :: p_Ddata

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setTracers')
      call sys_halt()
    end if
    
    ! Create a new variable. If necessary, increase the size of the buffers.
    if (rexport%ntracers .ne. 0) then
      ! There are already tracers and tracer data attached.
      ! Clean up tracer information as we create a new set of tracers.
      call ucd_removeTracers (rexport)
    end if
    
    ! Remember total number of tracers
    rexport%ntracers = ubound(DtracerCoordinates,2)

    ! Allocate a new vector for the data
    Ilength = ubound(DtracerCoordinates)
    call storage_new ('ucd_setTracers','htracers',&
        Ilength,ST_DOUBLE,rexport%htracers,ST_NEWBLOCK_NOINIT)
        
    ! Copy the coordinate data into that vector
    call storage_getbase_double2d (rexport%htracers,p_Ddata)
    p_Ddata = DtracerCoordinates
    
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_removeTracers (rexport)

!<description>
  ! Removes all information about tracers from the structure rexport.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    integer :: i
    
    ! Remove all tracer variables if there are some
    if (associated(rexport%p_HtracerVariables    )) then
      do i=1,rexport%ntracervariables
        if (rexport%p_HtracerVariables(i) .ne. ST_NOHANDLE) &
          call storage_free (rexport%p_HtracerVariables(i))
      end do
      deallocate(rexport%p_HtracerVariables    )
    end if

    ! Delete the coordinates of the tracers
    if (rexport%htracers .ne. ST_NOHANDLE) &
      call storage_free(rexport%htracers    )
        
    rexport%ntracers = 0

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_addTracerVariable (rexport,sname,Ddata)

!<description>
  ! Adds tracer variable data to the ouput file identified by rexport.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Name of the tracer variable.
  character(LEN=*), intent(in) :: sname
  
  ! array [1..#Tracers] of double. For every tracer I, Ddata(I) is the
  ! value of the variable that should be associated to that tracer.
  real(DP), dimension(:), intent(in) :: Ddata
!</input>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addTracerVariable')
      call sys_halt()
    end if
    
    if (rexport%ntracers .le. 0) then
      call output_line ('No tracers specified!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addTracerVariable')
      call sys_halt()
    end if
    
    if (size(Ddata) .lt. rexport%ntracers) then
      call output_line ('Ddata too small, more tracers than data!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addTracerVariable')
      call sys_halt()
    end if

    ! Create a new variable. If necessary, increase the size of the buffer.
    if (.not. associated(rexport%p_HtracerVariables)) then
      call ucd_moreTracerVariables(rexport)
    end if

    if (rexport%ntracerVariables .ge. size(rexport%p_HtracerVariables)) then
      call ucd_moreTracerVariables(rexport)
    end if
    
    rexport%ntracerVariables = rexport%ntracerVariables+1
    
    rexport%p_StracerVariableNames(rexport%ntracerVariables) = sname
    
    ! Allocate a new vector for the data
    call storage_new ('ucd_addTracerVariable','hvariable',&
        rexport%ntracers,ST_DOUBLE,&
        rexport%p_HtracerVariables(rexport%ntracerVariables),&
        ST_NEWBLOCK_ZERO)
        
    ! Copy the element data into that vector
    call storage_getbase_double (&
        rexport%p_HtracerVariables(rexport%ntracerVariables),p_Ddata)
    call lalg_copyVectorDble(Ddata(1:rexport%ntracers), &
        p_Ddata(1:rexport%ntracers))
        
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_setSimulationTime (rexport,dtime,ssimTimeFormat)

!<description>
  ! Specifies the simulation time, the export structure should represent.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Simulation time. Is written to the output file.
  real(DP), intent(in) :: dtime
  
  ! OPTIONAL: Fortran format string, e.g. "(F20.5)"
  ! Allows to specify an output format of the simulation time. Whether or
  ! not this is used for writing to the file depends on the type of
  ! output format (GMV, AVS,...).
  character(LEN=*), intent(in), optional :: ssimTimeFormat
!</input>

!</subroutine>

    character(LEN=SYS_STRLEN) :: stext

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setSimulationTime')
      call sys_halt()
    end if
    
    rexport%dsimulationTime = dtime
    
    ! Copy the output format string and overwrite the standard setting.
    if (present(ssimTimeFormat)) then
      rexport%ssimTimeFormat = ssimTimeFormat
      
      ! Text the output format string -- to be sure it is valid.
      ! If not, a compiler error is thrown here! This is either
      ! a runtime error or simply a message on the screen.
      write(stext,ssimTimeFormat) 0.0_DP
      if (stext .eq. "") then
        call output_line ('Invalid output format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'ucd_setSimulationTime')
        call sys_halt()
      end if
    end if
  
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_getSimulationTime (rexport,dtime)

!<description>
  ! Gets the simulation time from the export structure.
!</description>

!<input>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(in) :: rexport
!</input>

!<output>
  ! Simulation time.
  real(DP), intent(out) :: dtime
!</output>

!</subroutine>

    dtime=rexport%dsimulationTime
  
  end subroutine
  
  !************************************************************************
  
!<subroutine>

  subroutine ucd_setOutputNumberFormat (rexport,sformat)

!<description>
  ! Specifies the number format of numbers in the output file.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! Fortran format string, e.g. "(E15.8)"
  ! Specifies an output format of double precision numbers. Whether or
  ! not and where this is used for writing to the file depends on the 
  ! type of output format (GMV, AVS,...).
  character(LEN=*), intent(in), optional :: sformat
!</input>

!</subroutine>

    character(LEN=SYS_STRLEN) :: stext

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setOutputNumberFormat')
      call sys_halt()
    end if
    
    ! Copy the output format string and overwrite the standard setting.
    rexport%sdataFormat = sformat
    
    ! Text the output format string -- to be sure it is valid.
    ! If not, a compiler error is thrown here! This is either
    ! a runtime error or simply a message on the screen.
    write(stext,sformat) 0.0_DP
    if (stext .eq. "") then
      call output_line ('Invalid output format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_setOutputNumberFormat')
      call sys_halt()
    end if
  
  end subroutine
  
  !************************************************************************
  
!<subroutine>

  subroutine ucd_addCommentLine (rexport,scomment)

!<description>
  ! Adds a comment line to the comment buffer.
  ! All comments are written to a single comment block at the
  ! beginning of the output file.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! The comment to be added to the output.
  character(LEN=*), intent(in) :: scomment
!</input>

!</subroutine>

    character, dimension(:), pointer :: p_sbuf
    integer :: i

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addCommentLine')
      call sys_halt()
    end if
    
    ! Is there enough space in the output buffer? If not, reallocate it.
    if (.not. associated(rexport%p_Scomments)) then
      allocate(rexport%p_Scomments(max(len(scomment)+1,4*SYS_STRLEN)))
    else if ( rexport%ncommentBufSize+len(scomment)+1 .gt. &
              size(rexport%p_Scomments)) then
      allocate(p_sbuf (rexport%ncommentBufSize + max(len(scomment)+1,4*SYS_STRLEN)))
      p_sbuf(1:rexport%ncommentBufSize) = &
          rexport%p_Scomments(1:rexport%ncommentBufSize)
      deallocate(rexport%p_Scomments)
      rexport%p_Scomments => p_sbuf
    end if

    ! Append the string -- character by character
    do i=1,len(scomment)
      rexport%p_Scomments (rexport%ncommentBufSize+i) = scomment (i:i)
    end do
    
    ! Append a NEWLINE character as line end
    rexport%p_Scomments (rexport%ncommentBufSize+len(scomment)+1) = NEWLINE

    rexport%ncommentBufSize = rexport%ncommentBufSize+len(scomment)+1
        
  end subroutine

!************************************************************************
  
!<subroutine>

  subroutine ucd_addParameterList (rexport,rparList)

!<description>
  ! Adds the parameter list rparList to the output file identified by
  ! rexport. All parameters in this parameter list will be written
  ! as comment block to the output file.
!</description>

!<inputoutput>
  ! The ucd export structure identifying the output file.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!<input>
  ! A parameter list containing all configuration parameters of a simulation.
  type(t_parlist), intent(in) :: rparList
!</input>

!</subroutine>

    character, dimension(:), pointer :: p_sbuf,p_sparams

    if (rexport%coutputFormat .eq. UCD_FORMAT_NONE) then
      call output_line ('Export structure not initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'ucd_addParameterList')
      call sys_halt()
    end if
    
    ! Ask the parameter list to create a large array of characters
    ! containing the configuration.
    call parlst_getStringRepresentation (rparlist, p_sparams)
    
    ! May be the buffer is empty...
    if (.not. associated(p_sparams)) return
    
    ! By default, the parameter list uses the same separation character
    ! for all lines as we do. Therefore, we can simply add that large
    ! character block to our buffer.
    
    ! Is there enough space in the output buffer? If not, reallocate it.
    if (.not. associated(rexport%p_Scomments)) then
      allocate(rexport%p_Scomments(max(size(p_sparams)+1,4*SYS_STRLEN)))
    else if ( rexport%ncommentBufSize+size(p_sparams) .gt. &
              size(rexport%p_Scomments)) then
      allocate(p_sbuf (rexport%ncommentBufSize + max(size(p_sparams),4*SYS_STRLEN)))
      p_sbuf(1:rexport%ncommentBufSize) = &
          rexport%p_Scomments(1:rexport%ncommentBufSize)
      deallocate(rexport%p_Scomments)
      rexport%p_Scomments => p_sbuf
    end if
    
    ! Copy the data
    rexport%p_Scomments ( &
      rexport%ncommentBufSize+1:rexport%ncommentBufSize+size(p_sparams)) = p_sparams

    rexport%ncommentBufSize = rexport%ncommentBufSize+size(p_sparams)
    
    ! Remove data from the heap
    deallocate(p_sparams)
    
  end subroutine

  !************************************************************************

!<subroutine>

  subroutine ucd_readGMV (sfilename,rexport,rtriangulation)

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
  character(LEN=*), intent(in) :: sfilename
!</input>

!<output>
  ! The export structure where the GMV data is saved to.
  ! The structure can be used to attach more data. The filename for the
  ! export can be changed by ucd_setFilename.
  type(t_ucdExport), intent(out) :: rexport
!</output>

!<inputoutput>
  ! Triangulation structure.
  ! If this triangulation structure is empty, the mesh is read from the GMV 
  ! and a new triangulation is created in rtriangulation.
  ! If this contains a valid mesh, the mesh is assumed to correspond to
  ! the one in the GMV file and based on that, vector data is read.
  ! The mesh in the GMV file is skipped.
  type(t_triangulation), intent(inout), target :: rtriangulation
!</inputoutput>

!</subroutine>
  
    ! local variables
    integer :: mfile,ilinelen,ios
    logical :: biscomment
    character(LEN=SYS_STRLEN) :: sline,skey,sdummy
  
    ! Try to open the file
    if (sfilename .ne. '') then
      call io_openFileForReading(sfilename, mfile, .true.)
      if (mfile .lt. 0) then
        call output_line ('Cannot open file "'//trim(sfilename)&
                //'" for writing!',OU_CLASS_ERROR,OU_MODE_STD,'ucd_readGMV')
        call sys_halt()
      end if
    end if

    !----------------------------------------------------
    ! Read the GMV header
    call io_readlinefromfile (mfile, sline, ilinelen, ios)

    ! The triangulation in rexport will be rtriangulation --
    ! regardless whether we build it or have it
    rexport%p_rtriangulation => rtriangulation

    ! Set the other flags in rexport to standard values
    rexport%coutputFormat = UCD_FORMAT_GMV
    rexport%cflags = UCD_FLAG_STANDARD

    rexport%nvertices = rexport%p_rtriangulation%NVT
    rexport%ncells = rexport%p_rtriangulation%NEL

    ! Read each line and interpret it
    do while (ios .eq. 0)
      call io_readlinefromfile (mfile, sline, ilinelen, ios)
      
      ! Get the key word from the line
      call sys_tolower (sline)
      read(sline,*) skey
      
      if (trim(adjustl(skey)) .eq. "probtime") then
        
        !----------------------------------------------------
        ! Read the simulation time
        read(sline,*) sdummy,rexport%dsimulationTime
        
      elseif (trim(adjustl(skey)) .eq. "comments") then
      
        !----------------------------------------------------
        ! Read comments
        biscomment = .true.
        do while ((ios .eq. 0) .and. biscomment)
          call io_readlinefromfile (mfile, sline, ilinelen, ios)
          call sys_tolower (sline)
          read(sline,'(A)') skey
          biscomment = (trim(adjustl(skey)) .ne. 'endcomm')
          if (biscomment) then
            ! Add each comment line to our comment variable in rexport
            call ucd_addCommentLine (rexport,sline)
          end if
        end do
        
      elseif (trim(adjustl(skey)) .eq. "nodes") then
        
        !----------------------------------------------------
        ! Read triangulation (or ignore it if already given in rtriangulation)
        if (rexport%p_rtriangulation%ndim .eq. 0)&
            call read_triangulation (mfile,sline,rexport%p_rtriangulation)
        
        ! NEL/NVT has changed
        rexport%nvertices = rexport%p_rtriangulation%NVT
        rexport%ncells = rexport%p_rtriangulation%NEL
        
      elseif (trim(adjustl(skey)) .eq. "material") then
      
        !----------------------------------------------------
        ! Read materials
        call read_materials (mfile,sline,rexport)
        
      elseif (trim(adjustl(skey)) .eq. "velocity") then

        !----------------------------------------------------
        ! Read velocity data
        call read_velocity (mfile,sline,rexport)
        
      elseif (trim(adjustl(skey)) .eq. "variable") then
      
        !----------------------------------------------------
        ! Read general variable data
        call read_variables (mfile,sline,rexport)
        
      elseif (trim(adjustl(skey)) .eq. "polygons") then
      
        !----------------------------------------------------
        ! Read polygon data
        call read_polygondata (mfile,sline,rexport)
        
      elseif (trim(adjustl(skey)) .eq. "tracers") then
      
        !----------------------------------------------------
        ! Read tracer data
        call read_tracerdata (mfile,sline,rexport)
        
      elseif (trim(adjustl(skey)) .eq. "endgmv") then
        
        ! GMV finish
        ios = 1
        
      end if
    end do
    
    ! Close the file, finish.
    close (mfile)
    
  contains
  
    ! ---------------------------------------------------------------
    
    subroutine read_tracerdata (mfile,scommand,rexport)
    
    ! Reads data about tracers from the GMV file mfile.
    
    ! Handle to the GMV file
    integer, intent(in) :: mfile
    
    ! Command line in the GMV file with information about the tracers
    character(LEN=*), intent(in) :: scommand
    
    ! UCD structure where tracer data is saved to.
    type(t_ucdExport), intent(inout) :: rexport
    
      ! Local variables
      integer :: ilinelen,ios,ntracers
      character(LEN=SYS_STRLEN) :: skey,sline
      real(DP), dimension(:,:), allocatable :: DtracerCoordinates
      real(DP), dimension(:), allocatable :: DtracerTemp

      ! Interpret the command line, how many tracers do we have?
      read (scommand,*) skey,ntracers
      
      ! Allocate memory for X/Y or X/Y/Z coordinates
      allocate (DtracerCoordinates(rexport%p_rtriangulation%ndim,ntracers))
      allocate (DtracerTemp(ntracers))
      
      ! Read all X-coordinates
      read(mfile,*) DtracerCoordinates(1,:)
      
      ! Read all Y-coordinates
      read(mfile,*) DtracerCoordinates(2,:)
      
      ! Probably read all Z-coordinates -- or ignore them
      if (rexport%p_rtriangulation%ndim .eq. NDIM3D) then
        read(mfile,*) DtracerCoordinates(3,:)
      else 
        read(mfile,*) DtracerTemp(:)
      end if
      
      ! There must be an 'endtrace' at the end
      call io_readlinefromfile (mfile, sline, ilinelen, ios)
      call sys_tolower (sline)
      read(sline,*) skey
      if (trim(adjustl(skey)) .ne. 'endtrace') then
        call output_line ('Error reading GMV data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_tracerdata')
        call sys_halt() 
      end if
      
      ! Add the tracer data
      call ucd_setTracers (rexport,DtracerCoordinates)
      
      deallocate(DtracerTemp,DtracerCoordinates)

    end subroutine
    
    ! ---------------------------------------------------------------
    
    subroutine read_polygondata (mfile,scommand,rexport)
    
    ! Reads data about polygons from the GMV file mfile.
    
    ! Handle to the GMV file
    integer, intent(in) :: mfile
    
    ! Last command line in the GMV file 
    character(LEN=*), intent(in) :: scommand
    
    ! UCD structure where data is saved to.
    type(t_ucdExport), intent(inout) :: rexport
    
      ! Local variables
      integer :: ilinelen,ios,npoints,imaterial
      logical :: bfinish
      character(LEN=SYS_STRLEN) :: skey,sline
      real(DP), dimension(:,:), allocatable :: Dcoordinates
      real(DP), dimension(:), allocatable :: Dtemp

      ! Read the next line specifying basic data about the polygon
      call io_readlinefromfile (mfile, sline, ilinelen, ios)
      call sys_tolower (sline)
      
      read(sline,*) skey
      bfinish = (trim(adjustl(skey)) .ne. 'endtrace') 
      do while ((ios .eq. 0) .and. (.not. bfinish))
        
        ! Read material, #points
        write (mfile,*) imaterial,npoints
        
        ! Allocate memory for the polygon
        allocate (Dcoordinates(rexport%p_rtriangulation%ndim,npoints))
        allocate (Dtemp(npoints))
        
        ! Read the polygon. All X-coordinates
        read(mfile,*) Dcoordinates(1,:)
        
        ! Y-coordinates
        read(mfile,*) Dcoordinates(2,:)
        
        ! Probably Z-coordinates
        if (rexport%p_rtriangulation%ndim .eq. 3) then
          read(mfile,*) Dcoordinates(3,:)
        else 
          read(mfile,*) Dtemp(:)
        end if
        
        ! Add the polygon
        call ucd_addPolygon (rexport,Dcoordinates,imaterial)
        
        ! Deallocate data
        deallocate (Dtemp,Dcoordinates)
      
        ! Read the next line specifying basic data about the next polygon
        call io_readlinefromfile (mfile, sline, ilinelen, ios)
        call sys_tolower (sline)
        
        read(sline,*) skey
        bfinish = (trim(adjustl(skey)) .ne. 'endtrace') 
      end do
      
      if (ios .ne. 0) then
        call output_line ('Error reading GMV data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_polygondata')
        call sys_halt()
      end if
      
    end subroutine
    
    ! ---------------------------------------------------------------
    
    subroutine read_variables (mfile,scommand,rexport)
    
    ! Reads data about variables from the GMV file mfile.
    
    ! Handle to the GMV file
    integer, intent(in) :: mfile
    
    ! Last command line in the GMV file 
    character(LEN=*), intent(in) :: scommand
    
    ! UCD structure where data is saved to.
    type(t_ucdExport), intent(inout) :: rexport
    
      ! Local variables
      integer :: iel,ivt
      integer :: ilinelen,ios,itype,imaterial
      logical :: bfinish
      character(LEN=SYS_STRLEN) :: skey,sline,sname
      real(DP), dimension(:), allocatable :: Dcell
      real(DP), dimension(:), allocatable :: Dnode

      ! Allocate memory for cell- and node-based data
      allocate(Dcell(rexport%p_rtriangulation%NEL))
      allocate(Dnode(rexport%p_rtriangulation%NVT))

      ! Read the next line specifying basic data about the variable
      call io_readlinefromfile (mfile, sline, ilinelen, ios)
      call sys_tolower (sline)
      
      read(sline,*) skey
      bfinish = (trim(adjustl(skey)) .eq. 'endvars') 
      do while ((ios .eq. 0) .and. (.not. bfinish))
        
        ! Read variable name, type
        read (sline,*) sname,itype

        select case (itype)
        case (0) 
          ! Cell based data. Read and remember
          do iel = 1, rexport%p_rtriangulation%NEL
            read(mfile,FMT=rexport%sdataFormat) Dcell(iel)
          end do
          
          call ucd_addVariableElementBased (rexport,sname,UCD_VAR_STANDARD,Dcell)
              
        case (1)
          ! Node/Vertex based data.
          do ivt = 1, rexport%p_rtriangulation%NVT
            read(mfile,FMT=rexport%sdataFormat) Dnode(ivt)
          end do
          
          call ucd_addVariableVertexBased (rexport,sname,UCD_VAR_STANDARD,Dnode)
          
        case DEFAULT
          call output_line ('Error reading GMV data! Unknown variable type!', &
                  OU_CLASS_ERROR,OU_MODE_STD,'read_variables')
          call sys_halt()
  
        end select
        
        ! Read the next line specifying basic data about the next variable
        call io_readlinefromfile (mfile, sline, ilinelen, ios)
        call sys_tolower (sline)
        
        read(sline,*) skey
        bfinish = (trim(adjustl(skey)) .eq. 'endvars') 
        
      end do
      
      if (ios .ne. 0) then
        call output_line ('Error reading GMV data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_variables')
        call sys_halt()
      end if
      
    end subroutine

    ! ---------------------------------------------------------------
    
    subroutine read_materials (mfile,scommand,rexport)
    
    ! Reads data about materials from the GMV file mfile.
    
    ! Handle to the GMV file
    integer, intent(in) :: mfile
    
    ! Last read command line in the GMV file 
    character(LEN=*), intent(in) :: scommand
    
    ! UCD structure where tracer data is saved to.
    type(t_ucdExport), intent(inout) :: rexport
    
      ! Local variables
      integer :: ilinelen,ios,nmats,itype,i
      character(LEN=SYS_STRLEN) :: skey,sline
      character(LEN=32), dimension(:), allocatable :: Snames
      integer, dimension(:), allocatable :: Imat

      ! Interpret the command line, how many materials do we have? Type?
      read (scommand,*) skey,nmats,itype
      
      ! Allocate memory for the variable names
      allocate(Snames(nmats))
      
      ! Read the material names
      do i=1,nmats
        read(mfile,*) Snames(i)
      end do
      
      ! Allocate memory for vertex/element material classification
      select case (itype)
      case (0) 
        ! Cell based data. 
        allocate(Imat(rexport%p_rtriangulation%NEL))
            
      case (1)
        ! Node/Vertex based data.
        allocate(Imat(rexport%p_rtriangulation%NVT))
        
      case DEFAULT
        call output_line ('Error reading GMV data! Unknown variable type!', &
                OU_CLASS_ERROR,OU_MODE_STD,'read_materials')
        call sys_halt()

      end select
      
      ! Read material data of each cell / vertex
      read(mfile,*) Imat(:)
      
      ! GMV supports only the same material names for both, cells and
      ! vertices. Inform the rexport structure about the material names:
      call ucd_setMaterials (rexport,Snames)

      ! Specify the material for vertices / elements
      select case (itype)
      case (0) 
        ! Cell based data. 
        call ucd_setVertexMaterial (rexport,Imat)
            
      case (1)
        ! Node/Vertex based data.
        call ucd_setVertexMaterial (rexport,Imat)
        
      end select

      ! Deallocate memory      
      deallocate(Imat,Snames)

    end subroutine
    
    ! ---------------------------------------------------------------
    
    subroutine read_velocity(mfile,scommand,rexport)
    
    ! Reads data about velocity from the GMV file mfile.
    
    ! Handle to the GMV file
    integer, intent(in) :: mfile
    
    ! Last command line in the GMV file 
    character(LEN=*), intent(in) :: scommand
    
    ! UCD structure where data is saved to.
    type(t_ucdExport), intent(inout) :: rexport
    
      ! Local variables
      integer :: ilinelen,ios,itype
      logical :: bfinish
      character(LEN=SYS_STRLEN) :: sline,sname
      real(DP), dimension(:,:), allocatable :: Ddata

      ! Type of velocity data? Vertex or element based?
      read(scommand,*) sname,itype

      select case (itype)
      case (0) 
            
        ! Element based velocity data. Allocate memory.
        allocate(Ddata(rexport%p_rtriangulation%NEL,rexport%p_rtriangulation%ndim))
        
        ! Read the X-, Y- and Z-velocity
        read(mfile,*) Ddata(:,1)
        read(mfile,*) Ddata(:,2)
        read(mfile,*) Ddata(:,3)
        
        select case(rexport%p_rtriangulation%ndim)
        case (NDIM1D)
          ! Put the data to the rexport structure
          call ucd_addVarElemBasedVec (rexport, sname, Ddata(:,1))
        
        case (NDIM2D)
          ! Put the data to the rexport structure
          call ucd_addVarElemBasedVec (rexport, sname, Ddata(:,1), Ddata(:,2))

        case (NDIM3D)
          ! Put the data to the rexport structure
          call ucd_addVarElemBasedVec (rexport, sname, Ddata(:,1), Ddata(:,2), Ddata(:,3))
        end select
        
        ! Deallocate memory
        deallocate(Ddata)
            
      case (1)
        
        ! Node based velocity data. Allocate memory.
        allocate(Ddata(rexport%p_rtriangulation%NVT,rexport%p_rtriangulation%ndim))
        
        ! Read the X-, Y- and Z-velocity
        read(mfile,*) Ddata(:,1)
        read(mfile,*) Ddata(:,2)
        read(mfile,*) Ddata(:,3)
        
        select case(rexport%p_rtriangulation%ndim)
        case (NDIM1D)
          ! Put the data to the rexport structure
          call ucd_addVarVertBasedVec (rexport, sname, Ddata(:,1))
        
        case (NDIM2D)
          ! Put the data to the rexport structure
          call ucd_addVarVertBasedVec (rexport, sname, Ddata(:,1), Ddata(:,2))

        case (NDIM3D)
          ! Put the data to the rexport structure
          call ucd_addVarVertBasedVec (rexport, sname, Ddata(:,1), Ddata(:,2), Ddata(:,3))
        end select
        
        ! Deallocate memory
        deallocate(Ddata)
        
      case DEFAULT
        call output_line ('Error reading GMV data! Unknown variable type!', &
                OU_CLASS_ERROR,OU_MODE_STD,'read_velocity')
        call sys_halt()

      end select

    end subroutine

    ! ---------------------------------------------------------------
    
    subroutine read_triangulation (mfile,scommand,rtriangulation)
    
    ! Reads triangulation data from the GMV file mfile.
    
    ! Handle to the GMV file
    integer, intent(in) :: mfile
    
    ! Last read command line in the GMV file
    character(LEN=*), intent(in) :: scommand
    
    ! Triangulation structure. If this contains a valid triangulation,
    ! the mesh data in the GMV file is skipped.
    ! If rtriangulation is empty, a new triangulation with mesh data
    ! from the GMV file is set up.
    type(t_triangulation), intent(inout) :: rtriangulation
    
      ! Local variables
      integer :: ilinelen,ios,ntracers,n,i,nve,ive
      integer :: nmaxnve
      character(LEN=SYS_STRLEN) :: skey,sdummy,sline
      real(DP), dimension(:), allocatable :: Dx,Dy,Dz
      integer, dimension(:,:), allocatable :: IverticesAtElement
      real(DP), dimension(:,:), pointer :: p_Ddata2D
      integer, dimension(:,:), pointer :: p_Idata2D
      integer, dimension(2) :: Isize

      ! Do we have a "nodes x" or a "nodes fromfile" command?
      read(scommand,*) skey,sdummy
      
      if ((rtriangulation%ndim .ne. 0) .and. (skey .eq. 'fromfile')) then
        call output_line ('Cross reference to mesh not supported!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        call sys_halt()
      end if
      
      ! Quit if this is just this line that contains triangulation info
      if ((rtriangulation%ndim .ne. 0) .and. (skey .eq. 'fromfile')) return

      ! Get the number of vertices
      read(scommand,*) skey,n
      
      ! What is this for a command?
      if (skey .ne. 'nodes') then
        call output_line ('Unsupported GMV file structure!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        call sys_halt()
      end if
      
      rtriangulation%NVT = n
      
      ! Read point position data
      allocate(Dx(n),Dy(n),Dz(n))
      read(mfile,*) Dx(:)
      read(mfile,*) Dy(:)
      read(mfile,*) Dz(:)
      
      ! 1D or 2D mesh? It is 2D as soon as one Y-coordinate is <> 0
      rtriangulation%ndim = NDIM1D
      do i=1,n
        if (Dy(i) .ne. 0.0_DP) then
          rtriangulation%ndim = NDIM2D
          exit
        end if
      end do

      ! 1D/2D or 3D mesh? It is 3D as soon as one Z-coordinate is <> 0
      do i=1,n
        if (Dz(i) .ne. 0.0_DP) then
          rtriangulation%ndim = NDIM3D
          exit
        end if
      end do
      
      ! Allocate memory for the coordinates with the storage-system
      Isize = (/rtriangulation%ndim,rtriangulation%NVT/)
      call storage_new ('read_triangulation', 'DCORVG', Isize, ST_DOUBLE, &
          rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
      
      ! Get the pointers to the coordinate array
      call storage_getbase_double2D(&
          rtriangulation%h_DvertexCoords,p_Ddata2D)
      
      ! Copy the coordinates
      do i=1,n
        p_Ddata2D(NDIM1D,i) = Dx(i)
      end do
      if (rtriangulation%ndim .ge. NDIM2D) then
        do i=1,n
          p_Ddata2D(NDIM2D,i) = Dy(i)
        end do
      end if
      if (rtriangulation%ndim .ge. NDIM3D) then
        do i=1,n
          p_Ddata2D(NDIM3D,i) = Dz(i)
        end do
      end if
      
      ! Deallocate memory
      deallocate(Dz,Dy,Dx)
      
      ! At next we need KVERT = IverticesAtElement.
      call io_readlinefromfile (mfile, sline, ilinelen, ios)
      
      ! Do we have a "nodes x" or a "nodes fromfile" command?
      read(sline,*) skey,sdummy
      
      if (skey .eq. 'fromfile') then
        call output_line ('Cross reference to mesh not supported!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        call sys_halt()
      end if
      
      ! Get the number of cells
      read(sline,*) skey,n
      
      ! What is this for a command?
      if (skey .ne. 'cells') then
        call output_line ('Unsupported GMV file structure!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        call sys_halt()
      end if
      
      rtriangulation%NEL = n
      
      rtriangulation%InelOfType(:) = 0
      
      ! Allocate memory for connectivity. We support up to TRIA_MAXNVE
      ! vertices per element and do not know a-priori
      allocate(IverticesAtElement(TRIA_MAXNVE,n))
      
      ! Read the next n cell data tags.
      ! Format: Line 1: "type" NVE. Line 2: Vertex numbers.
      nmaxnve = 0
      do i=1,n
        IverticesAtElement(:,i) = 0
        read(mfile,*) sdummy,nve
        
        ! Increase the number of elements of that type
        rtriangulation%InelOfType(nve) = rtriangulation%InelOfType(nve)+1
        
        ! Figure out the size of the first dimension of KVERT
        nmaxnve = max(nve,nmaxnve)
        
        ! Read the vertex numbers
        read(mfile,*) IverticesAtElement(1:nve,i)
      end do
      
      rtriangulation%NNVE = nmaxnve
      
      ! All elements read in. Create the actual IverticesAtElement.

      Isize = (/nmaxnve,rtriangulation%NEL/)
      call storage_new ('read_triangulation', 'KVERT', Isize, ST_INT, &
          rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
          
      ! Get the pointer to the IverticesAtElement array and read the array
      call storage_getbase_int2D(&
          rtriangulation%h_IverticesAtElement,p_Idata2D)
      
      ! Copy the data.
      do i=1,n
        p_Idata2D(1:nmaxnve,i) = IverticesAtElement(1:nmaxnve,i)
      end do
      
      ! We do not need IverticesAtElement anymore...
      deallocate(IverticesAtElement)

      select case (rtriangulation%ndim)
      case (NDIM1D)
        ! Create some standard mesh information
        call tria_genElementsAtVertex1D2D  (rtriangulation)
        call tria_genNeighboursAtElement1D (rtriangulation)

        ! Reconstruct InodalProperty
        call reconstruct_InodalProperty_1D (rtriangulation)

        ! Now generate a standard mesh from the raw mesh.
        ! Generate all missing information.
        call tria_initStandardMeshFromRaw(rtriangulation)

      case (NDIM2D)
        ! Create some standard mesh information
        call tria_genElementsAtVertex1D2D  (rtriangulation)
        call tria_genNeighboursAtElement2D (rtriangulation)
        call tria_genEdgesAtElement2D      (rtriangulation)
        call tria_genElementsAtEdge2D      (rtriangulation)
        
        ! Reconstruct InodalProperty
        call reconstruct_InodalProperty_2D (rtriangulation)
        
        ! Now generate a standard mesh from the raw mesh.
        ! Generate all missing information.
        call tria_initStandardMeshFromRaw(rtriangulation)

      case DEFAULT
        ! For an (efficient) implementation of reconstructing InodalProperty
        ! in 3D, see [Jens Acker, Diploma Thesis].
        call output_line ('Only 2D supported! Cannot reconstruct InodalProperty!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'read_triangulation')
        call sys_halt()
      end select

    end subroutine
    
    subroutine reconstruct_InodalProperty_1D (rtriangulation)

      ! Reconstructs the InodalProperty, IboundaryCpIdx and 
      ! IverticesAtBoundary arrays. Sets NBCT!

      ! Triangulation structure. InodalProperty and NBCT are initialised here.
      ! InodalProperty must be allocated and initialised with 0.
      type(t_triangulation), intent(inout) :: rtriangulation

      ! local variables
      integer, dimension(:), pointer :: p_InodalProperty
      integer, dimension(:,:), pointer :: p_IverticesAtElement
      integer, dimension(:,:), pointer :: p_IneighboursAtElement
      integer, dimension(:), pointer :: p_IverticesAtBoundary
      integer, dimension(:), allocatable :: IverticesAtBoundary
      integer, dimension(:), pointer :: p_IboundaryCpIdx
      integer :: ivt,iel,ive,nbct,nvbd,ivbd,ibctidx,icurrentbc

      ! Allocate memory for the arrays 
      call storage_new ('reconstruct_InodalProperty_1D', 'KNPR', &
          rtriangulation%NVT, ST_INT, &
          rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)

      ! Get the pointer to some arrays
      call storage_getbase_int(&
          rtriangulation%h_InodalProperty,p_InodalProperty)
      call storage_getbase_int2d(&
          rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
      call storage_getbase_int2d(&
          rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)

      ! In a first loop find the nodes with no neighbour element. These are
      ! boundary edges.
      ! Count the number of boundary vertices.
      nvbd = 0
      do iel=1,rtriangulation%NEL
        do ive=1,ubound(p_IverticesAtElement,1)
          if (p_IneighboursAtElement (ive,iel) .eq. 0) then
            if (p_InodalProperty(p_IverticesAtElement(ive,iel)) .ne. -1) then
              p_InodalProperty(p_IverticesAtElement(ive,iel))     = -1
              nvbd = nvbd+1
            end if
          end if
        end do
      end do

      ! There must be two boundary components 
      ! - the interval start and end point
      if (nvbd .ne. 2) then
        call output_line ('Triangulation structure is invalid: NVBD does not match 2!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'reconstruct_InodalProperty_1D')
        call sys_halt()
      else
        rtriangulation%NVBD = nvbd
      end if

      ! Allocate memory for IverticesAtBoundary.
      call storage_new ('reconstruct_InodalProperty_1D', &
          'KVBD', rtriangulation%NVBD, &
          ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)

      call storage_getbase_int(&
          rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)

      ! Now loop through all vertices. Whenever we find a vertex with InodalProperty=-1,
      ! we start going from vertex to vertex until we found all vertices of
      ! that boundary component.
      ! We save all boundary vertices in IverticesAtBoundary.
      nbct = 0
      ivbd = 0
      do ivt = 1,rtriangulation%NVT
        if (p_InodalProperty(ivt) .eq. -1) then
          ! New boundary component found
          nbct = nbct+1

          ! Save the boundary component.
          p_InodalProperty(ivt) = nbct

          ! Save the vertex.
          ivbd = ivbd+1
          p_IverticesAtBoundary(ivbd) = ivt
        end if
      end do

      ! nbct is the number of boundary components we found.
      rtriangulation%NBCT = nbct

      ! Allocate memory for the boundary component index vector.
      ! Initialise that with zero!
      call storage_new ('reconstruct_InodalProperty_1D', &
          'KBCT', rtriangulation%NBCT+1, &
          ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
      call storage_getbase_int (&
          rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)

      ! Figure out the IboundaryCpIdx array but once checking the BC
      ! number of all vertices on the boundary
      ibctidx = 0
      icurrentbc = 0
      do ivbd = 1,nvbd
        if (p_InodalProperty(p_IverticesAtBoundary(ivbd)) .ne. icurrentbc) then
          ibctidx = ibctidx+1
          p_IboundaryCpIdx(ibctidx) = ivbd
          icurrentbc = p_InodalProperty(p_IverticesAtBoundary(ivbd))
        end if
      end do
      
      ! InodalProperty is completely classified -- that is it.

    end subroutine

    subroutine reconstruct_InodalProperty_2D (rtriangulation)
    
      ! Reconstructs the InodalProperty, IboundaryCpIdx and 
      ! IverticesAtBoundary arrays. Sets NBCT!
      ! DvertexParameterValue is cleared.
      
      ! Triangulation structure. InodalProperty and NBCT are initialised here.
      ! InodalProperty must be allocated and initialised with 0.
      type(t_triangulation), intent(inout) :: rtriangulation
    
      ! local variables
      integer, dimension(:), pointer :: p_InodalProperty
      integer, dimension(:,:), pointer :: p_IverticesAtElement
      integer, dimension(:,:), pointer :: p_IedgesAtElement
      integer, dimension(:,:), pointer :: p_IneighboursAtElement
      integer, dimension(:), pointer :: p_IelementsAtVertex
      integer, dimension(:), pointer :: p_IelementsAtVertexIdx
      integer, dimension(:), pointer :: p_IverticesAtBoundary
      integer, dimension(:), allocatable :: IverticesAtBoundary
      integer, dimension(:), pointer :: p_IboundaryCpIdx
      integer :: ivt,ivt2
      integer :: iel,ielidx
      integer :: ive,iveprevious,ivenext,nbct,nve,nvbd,ivbd,ibctidx,icurrentbc
      real(DP), dimension(:), pointer :: p_DvertexParameterValue
      
      ! Allocate memory for the arrays 
      call storage_new ('reconstruct_InodalProperty_2D', 'KNPR', &
          rtriangulation%NVT, ST_INT, &
          rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
      
      ! Get the pointer to some arrays
      call storage_getbase_int(&
          rtriangulation%h_InodalProperty,p_InodalProperty)
      call storage_getbase_int2d(&
          rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
      call storage_getbase_int2d(&
          rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
      call storage_getbase_int2d(&
          rtriangulation%h_IneighboursAtElement,p_IneighboursAtElement)
      call storage_getbase_int(&
          rtriangulation%h_IelementsAtVertex,p_IelementsAtVertex)
      call storage_getbase_int(&
          rtriangulation%h_IelementsAtVertexIdx,p_IelementsAtVertexIdx)
          
      ! In a first loop find the edges with no neighbour element. These are
      ! boundary edges.
      ! Count the number of boundary vertices.
      nvbd = 0
      do iel=1,rtriangulation%NEL
        do ive=1,ubound(p_IverticesAtElement,1)
          if (p_IverticesAtElement(ive,iel) .eq. 0) exit ! Triangle in a quad mesh
          if (p_IneighboursAtElement (ive,iel) .eq. 0) then
            
            ! Both vertices on the edge are boundary vertices
            ivenext = mod(ive,ubound(p_IverticesAtElement,1))+1
            if (p_IverticesAtElement(ivenext,iel) .eq. 0) ivenext = ivenext-1
            
            if (p_InodalProperty(p_IverticesAtElement(ive,iel)) .ne. -1) then
              p_InodalProperty(p_IverticesAtElement(ive,iel))     = -1
              nvbd = nvbd+1
            end if
            if (p_InodalProperty(p_IverticesAtElement(ivenext,iel)) .ne. -1) then
              p_InodalProperty(p_IverticesAtElement(ivenext,iel)) = -1
              nvbd = nvbd+1
            end if
          end if
        end do
      end do
      
      rtriangulation%NVBD = nvbd
      
      ! Allocate memory for IverticesAtBoundary.
      call storage_new ('reconstruct_InodalProperty_2D', &
          'KVBD', rtriangulation%NVBD, &
          ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)

      call storage_getbase_int(&
          rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
      
      ! Now loop through all vertices. Whenever we find a vertex with InodalProperty=-1,
      ! we start going from vertex to vertex until we found all vertices of
      ! that boundary component.
      ! We save all boundary vertices in IverticesAtBoundary.
      nbct = 0
      ivbd = 0
      do ivt = 1,rtriangulation%NVT
        if (p_InodalProperty(ivt) .eq. -1) then
          ! New boundary component found
          nbct = nbct+1
          
          ! Go from edge to edge until we hit ivt again. This must happen as the 
          ! boundary component must be closed.
          ivt2 = ivt
          
          ivt2loop: do
            p_InodalProperty(ivt2) = nbct
            
            ! Save the vertex.
            ivbd = ivbd+1
            p_IverticesAtBoundary(ivbd) = ivt2
            
            ! Check the edges adjacent to the vertex...for this purpose, we must 
            ! find them.
            do ielidx = p_IelementsAtVertexIdx(ivt2),p_IelementsAtVertexIdx(ivt2+1)-1

              iel = p_IelementsAtVertex(ielidx)
            
              ! Where is the point in the element
              do ive=1,ubound(p_IverticesAtElement,1)
                if (p_IverticesAtElement(ive,iel) .eq. ivt2) exit 
              end do
              
              ! Triangle or quad?
              nve = ubound(p_IverticesAtElement,1)
              if (p_IverticesAtElement(nve,iel) .eq. 0) nve = nve-1

              iveprevious = mod(ive+nve-2,nve)+1
              ivenext = mod(ive,nve)+1
              
              ! On that element, check the edge following the vertex if it is a boundary edge.
              if (p_IneighboursAtElement(ive,iel) .eq. 0) then
              
                ! Yes, that is a boundary edge. does it lead to a vertex that we have not had?
                if (p_InodalProperty(p_IverticesAtElement(ivenext,iel)) .eq. -1) then
                  ! Yes! That node is now a boundary node. We continue with it.
                  ivt2 = p_IverticesAtElement(ivenext,iel)
                  cycle ivt2loop
                end if
              
              end if
                
              ! No, it is not. Check the previous edge
              if (p_IneighboursAtElement(iveprevious,iel) .eq. 0) then
              
                ! Yes, that is a boundary edge. does it lead to a vertex that we have not had?
                if (p_InodalProperty(p_IverticesAtElement(iveprevious,iel)) .eq. -1) then
                  ! Yes! That node is now a boundary node. We continue with it.
                  ivt2 = p_IverticesAtElement(iveprevious,iel)
                  cycle ivt2loop
                end if
              
              end if
                
            end do
          
            ! Ok, there is no edge found starting from the vertex ivt2 that leads to 
            ! another vertex on the boundary. That means we found all vertices
            ! on the boundary component!
            ! So we quit the loop here and continue to find the next non-classified
            ! vertex ivt.
            exit ivt2loop
            
          end do ivt2loop
          
        end if
      end do
    
      ! nbct is the number of boundary components we found.
      rtriangulation%NBCT = nbct
    
      ! Allocate memory for the boundary component index vector.
      ! Initialise that with zero!
      call storage_new ('reconstruct_InodalProperty_2D', &
          'KBCT', rtriangulation%NBCT+1, &
          ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
      call storage_getbase_int (&
          rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
      ! Figure out the IboundaryCpIdx array but once checking the BC
      ! number of all vertices on the boundary
      ibctidx = 0
      icurrentbc = 0
      do ivbd = 1,nvbd
        if (p_InodalProperty(p_IverticesAtBoundary(ivbd)) .ne. icurrentbc) then
          ibctidx = ibctidx+1
          p_IboundaryCpIdx(ibctidx) = ivbd
          icurrentbc = p_InodalProperty(p_IverticesAtBoundary(ivbd))
        end if
      end do
      
      p_IboundaryCpIdx(nbct+1) = nvbd+1
      
      ! Allocate memory for  and DvertexParameterValue
      call storage_new ('reconstruct_InodalProperty_2D', &
          'DVBDP', rtriangulation%NVBD, &
          ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
      
      ! Get the array where to store boundary parameter values.
      call storage_getbase_double (&
          rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)

      ! Clear the array -- we have no parameter values!
      call lalg_clearVectorDble (p_DvertexParameterValue)
    
      ! InodalProperty is completely classified -- that is it.
    
    end subroutine

  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_setFilename(rexport,sfilename)

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
  character(LEN=*), intent(in) :: sfilename
!</input>

!<inputoutput>
  ! The ucd export structure that specifies the output.
  type(t_ucdExport), intent(inout) :: rexport
!</inputoutput>

!</subroutine>

    ! Start the output
    select case(rexport%coutputFormat)
    case (UCD_FORMAT_GMV,UCD_FORMAT_AVS,UCD_FORMAT_VTK)
      ! Here it is ok to change the filename
      rexport%sfilename = sfilename
    case DEFAULT
      call output_line ('Cannot change filename', &
                        OU_CLASS_ERROR,OU_MODE_STD,'ucd_setFilename')
      call sys_halt()
    end select
    
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_getVariable(rexport,svarName,Ddata,nlength,itype)

!<description>
  ! Reads the data of the variable svarName. If Ddata is specified,
  ! the variable data is written to Ddata. If nlength is specified,
  ! the length of the variable is returned in nlength.
!</description>

!<input>
  ! The ucd export structure containing data.
  type(t_ucdExport), intent(in) :: rexport

  ! Name of the variable whose data should be retrieved.
  character(LEN=*), intent(in) :: svarName
!</input>

!<output>
  ! OPTIONAL: Data array. Must be large enough to hold all data.
  ! The data of variable svarName is transferred to Ddata.
  !
  ! If the variable is unknown, Ddata is not changed.
  real(DP), dimension(:), intent(inout), optional :: Ddata
  
  ! OPTIONAL: Length qualifier.
  ! If specified, nlength is set to the number of elements in
  ! the variable svarName, thus to the minimum length Ddata must
  ! have to accept the data of the variable.
  !
  ! If the variable is unknown, -1 is returned.
  integer, intent(out), optional :: nlength
  
  ! OPTIONAL: Type of data. One of the UCD_BASE_xxxx flags.
  ! =UCD_BASE_ELEMENT: element based data.
  ! =UCD_BASE_VERTEX: vertex based data.
  !
  ! Is set to -1 if the variable is unknown.
  integer, intent(out), optional :: itype
!</output>

!</subroutine>

    integer :: i
    character(LEN=SYS_NAMELEN) :: sname
    real(DP), dimension(:), pointer :: p_Ddata
    
    sname = svarName
    call sys_toupper_replace (sname)
    
    ! Find the variable
    do i=1,size(rexport%p_SvariableNames)
      if (sys_upcase(rexport%p_SvariableNames(i)) .eq. sname) then
        ! Variable found! Return data as desired.
        
        call storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
        
        if (present(nlength)) nlength = size(p_Ddata)
        if (present(itype)) itype = rexport%p_IvariableBase(i)
        if (present(Ddata)) then
          call lalg_copyVectorDble (p_Ddata,Ddata(1:size(p_Ddata)))
        end if
        
        ! That is it, return.
        return
      end if
    
    end do
    
    ! Variable not found.
    if (present(nlength)) nlength = -1
    if (present(itype)) itype = -1
    
  end subroutine

  !************************************************************************
  
!<subroutine>

  subroutine ucd_infoVariables(rexport)

!<description>
  ! Output list of variables present in rexport.
!</description>

!<input>
  ! The ucd export structure containing data.
  type(t_ucdExport), intent(in) :: rexport
!</input>
!</subroutine>

    integer :: i

    ! Loop over all variables
    write(*,fmt='(A)') 'Names of variables'
    do i = 1, size(rexport%p_SvariableNames)
      write(*,fmt='(I3,1X,A)') i, rexport%p_SvariableNames(i)
    end do
    write(*,*)
    
    ! Loop over all variable vectors
    write(*,fmt='(A)') 'Names of variable vectors'
    do i = 1, size(rexport%p_SvarVecNames)
      write(*,fmt='(I3,1X,A)') i, rexport%p_SvarVecNames(i)
    end do
    write(*,*)
    
  end subroutine ucd_infoVariables

end module ucd

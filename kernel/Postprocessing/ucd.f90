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
!#  2.) ucd_setAlternativeSource
!#      -> Allows to specify an alternative source file for parts of a
!#         triangulation or the whole triangulation. (Saves disc space
!#         if a sequence of output files are written.)
!#
!#  3.) ucd_addVariableVertexBased / ucd_addVariableElementBased
!#      -> Add a vertex- or element-based data vector to the output.
!#
!#  4.) ucd_setVertexMaterial
!#      -> Allows to specify for every vertex/node (corner vertex, edge 
!#         midpoint,...) a material id.
!#
!#  5.) ucd_setCellMaterial
!#      -> Allows to specify for every cell a material id.
!#
!#  6.) ucd_setMaterials
!#      -> Allows to specify material names for the material id's that
!#         are assigned to cells/vertices in ucd_setVertexMaterial and
!#         ucd_setCellMaterial, respectively.
!#
!#  7.) ucd_setTracers
!#      -> Specify tracers by (X,Y) or (X,Y,Z) coordinates
!#
!#  8.) ucd_addTracerVariable
!#      -> Adds a tracer variable, i.e. a named value for every tracer.
!#
!#  9.) ucd_removeTracers
!#      -> Removes all tracer information.
!#
!# 10.) ucd_setSimulationTime
!#      -> Specify simulation time.
!#
!# 11.) ucd_addPolygon
!#      -> Adds a point-array as polygon of line-segments to the output.
!#
!# 12.) ucd_addCommentLine
!#      -> Add a comment line.
!#
!# 13.) ucd_addParameterList
!#      -> Add a configuration parameter list as comments to the output.
!#
!# 14.) ucd_setOutputNumberFormat
!#      -> Allows to specify a custom format of double precision numbers 
!#         in the output file.
!#
!# 15.) ucd_write
!#      -> Writes the output file to the hard disc, closes the file.
!#
!# 16.) ucd_release
!#      -> Releases a UCD export structure that was created by ucd_startXXXX.
!#
!# To write a file (GMV, AVS, MATLAB or what else), one has to use the
!# following sequence of commands:
!#
!# a) ucd_startGMV      - Create an output structure for a file
!#
!# b) ucd_addVariableVertexBased  
!#                      - Add a variable consisting of values in vertices
!#    ucd_addVariableElementBased 
!#                      - Add a variable consisting of values in elements
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

!</constantblock>

!</constants>

  
!<types>

!<typeblock>
  
  ! UCD export structure. The structure is created by one of the ucd_startXXXX
  ! routines and configured to write output of type XXXX (e.g. GMV or AVS).
  ! After writing out the data, the structure can be released with ucd_release.
  
  TYPE t_ucdExport
  
    PRIVATE
    
    ! Output format. 0=undefined, 1=GMV, 2=AVS, 3=MATLAB
    INTEGER :: coutputFormat = 0
    
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
    
    ! The simulation time. =SYS_INFINITY if no simulation time is specified
    REAL(DP) :: dsimulationTime  = SYS_INFINITY
    
    ! Format of the simulation time. Fortran format string.
    CHARACTER(LEN=SYS_STRLEN) :: ssimTimeFormat = "(F20.5)"
    
    ! Format of the output of double-precision numbers. Fortran format string.
    CHARACTER(LEN=SYS_STRLEN) :: sdataFormat = "(E15.8)"
    
    
    ! An array containing the names of all the variables
    CHARACTER(LEN=SYS_NAMELEN), DIMENSION(:), POINTER :: p_SvariableNames => NULL()
    
    ! Filename of file containing point coordinates.
    ! ""=no alternative source file.
    CHARACTER(LEN=SYS_STRLEN) :: saltFilePoints = ""

    ! Filename of file containing cell structure. 
    ! ""=no alternative source file.
    CHARACTER(LEN=SYS_STRLEN) :: saltFileCells = ""
    
    ! A pointer to the underlying triangulation
    TYPE(t_triangulation), POINTER :: p_rtriangulation => NULL()
    
    ! A pointer to an array with specification flags. p_IvariableSpec(I)
    ! is a bitfield for variable I that specifies the type of the variable
    ! and how to handle it.
    INTEGER(I32), DIMENSION(:), POINTER :: p_IvariableSpec => NULL()
    
    ! A pointer to an array that specifies whether a variable is vertex based (1)
    ! or cell based (0).
    INTEGER, DIMENSION(:), POINTER :: p_IvariableBase => NULL()
    
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

!</types>

  PRIVATE :: ucd_moreVariables, ucd_morePolygons, ucd_moreTracerVariables

CONTAINS

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
    
    rexport%coutputFormat = 1
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

    IF (ASSOCIATED(rexport%p_SvariableNames )) DEALLOCATE(rexport%p_SvariableNames )
    IF (ASSOCIATED(rexport%p_IvariableSpec)) DEALLOCATE(rexport%p_IvariableSpec)
    IF (ASSOCIATED(rexport%p_IvariableBase)) DEALLOCATE(rexport%p_IvariableBase)
    
    IF (ASSOCIATED(rexport%p_Hvariables   )) THEN
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_setAlternativeSource: Export structure not initialised!'
      STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_setMaterials: Export structure not initialised!'
      STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_setCellMaterial: Export structure not initialised!'
      STOP
    END IF
    
    NEL = rexport%p_rtriangulation%NEL
    
    IF (SIZE(Imaterials) .LT. NEL) THEN
      PRINT *,'ucd_setCellMaterial error: Imaterials invalid!'
      STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_setVertexMaterial: Export structure not initialised!'
      STOP
    END IF
    
    NVT = rexport%p_rtriangulation%NVT
    NMT = rexport%p_rtriangulation%NMT
    NEL = rexport%p_rtriangulation%NEL
    
    IF (SIZE(ImaterialsVert) .LT. NVT) THEN
      PRINT *,'ucd_setVertexMaterial error: ImaterialsVert invalid!'
      STOP
    END IF

    IF (PRESENT(ImaterialsMid)) THEN
      IF (SIZE(ImaterialsMid) .LT. NMT) THEN
        PRINT *,'ucd_setVertexMaterial error: ImaterialsMid invalid!'
        STOP
      END IF
    END IF

    IF (PRESENT(ImaterialsElem)) THEN
      IF (SIZE(ImaterialsElem) .LT. NEL) THEN
        PRINT *,'ucd_setVertexMaterial error: ImaterialsElem invalid!'
        STOP
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

    ! If there is a new filename, open the output file.
    IF (rexport%sfilename .NE. '') THEN
      CALL io_openFileForWriting(rexport%sfilename, rexport%iunit, SYS_REPLACE)
      IF (rexport%iunit .LT. 0) THEN
        PRINT *,'ucd_write: Cannot open file "'//TRIM(rexport%sfilename)&
                //'" for writing!'
        STOP
      END IF
    END IF

    IF (rexport%coutputFormat .EQ. 1) THEN
      ! Start the GMV output
      CALL ucd_writeGMV (rexport)
    END IF
    
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
        DO j=i,rexport%ncommentBufSize
          IF (rexport%p_Scomments(j) .EQ. NEWLINE) EXIT
        END DO
        
        DO WHILE (j .LE. rexport%ncommentBufSize)
          ! Write the line, continue with the next
          DO k=i,j-1
            WRITE(mfile,'(A)',ADVANCE='NO') rexport%p_Scomments(k)
          END DO
          WRITE(mfile,'(A)')
          
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
              PRINT *,'Invalid element!'
            END SELECT
              
          END DO
            
        ELSE

          ! 1x refined mesh
          
          ! The edge numbers give the numbers of the edge midpoints and thus
          ! the numbers of the new vertices on the once refined mesh.
          CALL storage_getbase_int2d (rexport%p_Rtriangulation%h_IedgesAtElement,&
              p_IedgesAtElement)
              
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
          
        END IF
        
      END IF

      !----------------------------------------------------
      ! Write material names -- if specified
      !
      ! Cell material
      IF (ASSOCIATED(rexport%ScellMaterials)) THEN
        ! GMV only supports <= 1000 materials!
        WRITE(mfile,'(A,I10)') 'material ',MIN(1000,SIZE(rexport%ScellMaterials)),0
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
      IF (ASSOCIATED(rexport%SvertexMaterials)) THEN
        ! GMV only supports <= 1000 materials!
        WRITE(mfile,'(A,I10)') 'material ',MIN(1000,SIZE(rexport%SvertexMaterials)),1
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
        IF (ASSOCIATED(rexport%ScellMaterials)) THEN
          ! GMV only supports <= 1000 materials!
          WRITE(mfile,'(A,I10)') 'material ',MIN(1000,SIZE(rexport%ScellMaterials)),1
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
              (rexport%p_IvariableBase(i) .EQ. 0)) THEN
              
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
                  (rexport%p_IvariableBase(j) .EQ. 0)) THEN
                  
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
                  (rexport%p_IvariableBase(j) .EQ. 0)) THEN
                  
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
              (rexport%p_IvariableBase(i) .EQ. 1)) THEN
              
            ! Found it. Write it out.
            WRITE (mfile,'(A)') 'velocity 1'
            
            CALL storage_getbase_double (rexport%p_Hvariables(i),p_Ddata)
            DO ivt=1,rexport%nvertices
              WRITE (MFILE,rexport%sdataFormat) p_Ddata(ivt)
            END DO
            
            ! Find the Y-velocity
            DO j=1,rexport%nvariables
              IF ((IAND(rexport%p_IvariableSpec(j),UCD_VAR_YVELOCITY) .NE. 0) .AND. &
                  (rexport%p_IvariableBase(j) .EQ. 1)) THEN
                  
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
                  (rexport%p_IvariableBase(j) .EQ. 0)) THEN
                  
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
            
            IF (rexport%p_IvariableBase(i) .EQ. 0) THEN
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
      IF (ASSOCIATED(rexport%p_Hpolygons)) THEN
      
        ! At least one polygon
        WRITE (MFILE,'(A)') 'polygons'

        ! Coordinates        
        CALL storage_getbase_double2D (rexport%p_Hpolygons(i),p_Ddata2D)
        
        ! Materials
        CALL storage_getbase_int (rexport%hpolygonMaterial,p_Idata)

        ! Write all polygons.
        ! Either we have 2D or 3D coordinates...
        IF (UBOUND(p_Ddata2D,1) .EQ. 2) THEN
          DO i=1,SIZE(rexport%p_Hpolygons)
            
            ! Write material, #points
            WRITE (MFILE,'(2I10)',ADVANCE='NO') p_Idata(i),UBOUND(p_Ddata2D,2)
            
            ! Write coordinates of the points forming the line segments 
            ! of the polygon
            DO j=1,UBOUND(p_Ddata2D,2)
              WRITE (MFILE,'(3E15.8)',ADVANCE='NO') &
                p_Ddata2D(1,j),p_Ddata2D(2,j),0.0_DP
            END DO
                
          END DO
        ELSE
          DO i=1,SIZE(rexport%p_Hpolygons)
            
            ! Write material, #points
            WRITE (MFILE,'(2I10)',ADVANCE='NO') p_Idata(i),UBOUND(p_Ddata2D,2)
            
            ! Write coordinates of the points forming the line segments 
            ! of the polygon
            DO j=1,UBOUND(p_Ddata2D,2)
              WRITE (MFILE,'(3E15.8)',ADVANCE='NO') &
                p_Ddata2D(1,j),p_Ddata2D(2,j),p_Ddata2D(3,j)
            END DO
                
          END DO
        END IF
        
        WRITE (mfile,'(A)') 'endpoly'
      
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
    CALL storage_realloc ('ucd_morePolygons', INT(nsize+16,I32), &
        rexport%hpolygonMaterial, ST_NEWBLOCK_ZERO)
    
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
  
  ! DdataVert(I) os the value of the variable in vertex I of the triangulation.
  REAL(DP), DIMENSION(:), INTENT(IN) :: DdataVert

  ! OPTIONAL: DdataMid(I) os the value of the variable in edge midpoint I of 
  ! the triangulation. Must be specified if DdataElem is specified!
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataMid

  ! OPTIONAL: DdataElem(I) os the value of the variable in element midpoint I of 
  ! the triangulation.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DdataElem
!</input>

!</subroutine>

    ! local varables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    IF (PRESENT(DdataMid) .AND. .NOT. PRESENT(DdataElem)) THEN
      PRINT *,'ucd_addVariableVertexBased: Error in the parameters!'
      STOP
    END IF
    
    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_addVariableVertexBased: Export structure not initialised!'
      STOP
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
    
    rexport%p_IvariableBase(rexport%nvariables) = 1
    
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
        PRINT *,'ucd_addVariableVertexBased: Warning. No edge midpoint data available!'
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
        PRINT *,'ucd_addVariableVertexBased: Warning. No element midpoint '//&
                'data available!'
      END IF
    END IF    

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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_addVariableElementBased: Export structure not initialised!'
      STOP
    END IF
    
    ! Create a new variable. If necessary, increase the size of the buffer.
    IF (rexport%nvariables .GE. SIZE(rexport%p_Hvariables)) THEN
      CALL ucd_moreVariables(rexport)
    END IF
    rexport%nvariables = rexport%nvariables+1
    rexport%p_IvariableSpec(rexport%nvariables) = cvarSpec
    rexport%p_IvariableBase(rexport%nvariables) = 0
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_addPolygon: Export structure not initialised!'
      STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_setTracers: Export structure not initialised!'
      STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_addTracerVariable: Export structure not initialised!'
      STOP
    END IF
    
    IF (rexport%ntracers .LE. 0) THEN
      PRINT *,'ucd_addTracerVariable: No tracers specified!'
      STOP
    END IF
    
    IF (SIZE(Ddata) .LT. rexport%ntracers) THEN
      PRINT *,'ucd_addTracerVariable: Ddata too small, more tracers than data!'
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_addTracerVariable: Export structure not initialised!'
      STOP
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
        PRINT *,'ucd_setSimulationTime: Invalid output format!'
        STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_setOutputNumberFormat: Export structure not initialised!'
      STOP
    END IF
    
    ! Copy the output format string and overwrite the standard setting.
    rexport%sdataFormat = sformat
    
    ! Text the output format string -- to be sure it's valid.
    ! If not, a compiler error is thrown here! This is either
    ! a runtime error or simply a message on the screen.
    WRITE(stext,sformat) 0.0_DP
    IF (stext .EQ. "") THEN
      PRINT *,'ucd_setOutputNumberFormat: Invalid output format!'
      STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_addCommentLine: Export structure not initialised!'
      STOP
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

    IF (rexport%coutputFormat .EQ. 0) THEN
      PRINT *,'ucd_addParameterList: Export structure not initialised!'
      STOP
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

END MODULE

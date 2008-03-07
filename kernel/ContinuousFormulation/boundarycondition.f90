!##############################################################################
!# ****************************************************************************
!# <name> boundarycondition </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module defines the analytical representation of boundary conditions.
!# In detainl, the following types of boundary conditions are supported:
!# - Neumann boundary conditions
!# - Dirichlet boundary conditions
!# - Pressure drop boundary conditions
!# The structure here only describe the boundary conditions. There are no
!# routines in this module to impose them into a solution/defect vector.
!# For imposing boundary conditions into vectors, the boundary conditions
!# first have to be 'discretised' using the module 'DiscreteBC'!
!#
!# IMPORTANT NOTE:
!# THIS MODULE IS OUTDATED AND NOT USED ANYMORE! IT'S STILL CONTENT OF THE
!# FEAT2 LIBRARY TO ALLOW A USER TO DEFINE ANALYTIC BONUDARY CONDITIONS ON
!# APPLICATION LEVEL. THE KERNAL ROUTINES HOWEVER DON'T USE ANALYTICAL
!# BOUNDARY CONDITIONS ANYMORE: BOUNDARY CONDITIONS ARE ALWAYS DIRECTLY
!# ASSEMBLED INTO THEIR DISCRETE COUNTERPART SUCH THAT THEY CAN BE IMPLEMENTED
!# INTO MATRICES/VECTORS!
!#
!# Neumann boundary conditions are usually imposed as 'do nothing' boundary
!# conditions and therefore not realised in a special structure. This
!# means, the module here realises all 'other' types of boundary conditions
!# not being Neumann. Neumann is standard, and all other types of BC's
!# must be imposed.
!#
!# Method description \\
!# ------------------ \\
!# The t_boundaryConditions structure collects all boundary conditions to a
!# domain. It contains three lists of 'boundary condition regions'
!# (structure t_bcRegion), for
!# - boundary conditions bound to segments on the real boundary
!# - boundary conditions anywhere on the real boundary
!# - boundary conditions, imposed by fictitious boundary objects.
!#
!# A 'boundary condition region' object t_bcRegion realises a boundary
!# condition. It contains information about the type of the boundary conditions
!# and where they are, on the real boundary as well as on fictitious boundary
!# objects. Furthermore, it contains information about the solution component
!# the boundary condition refers to (e.g. X-velocity, Y-velocity,...)
!#
!# A boundary condition can be deactivated by setting the cbcRegionType flag
!# in the boundary region to BC_DONOTHING.
!#
!# BC's on real boundary components in 2D \\
!# -------------------------------------- \\
!# Boundary conditions on real boundary components are described via a boundary
!# region struture as defined in boundary.f90. This contains basically a minimum
!# and maximum parameter value on the boundary as well as s number of a
!# boundary component. Optionally, the boundary region may be specified to
!# be attached to a real boundary segment, e.g. a full circle or a line on
!# the boundary.
!#
!# BC's on fictitious boundary components \\
!# -------------------------------------- \\
!# The concept of (analytic) fictitious boundary components is somehow
!# different to the concept of standard bondary conditions.
!# An analytic fictitious boundary component is a 'virtual object' in the
!# domain that is associated to a special boundary condition. Each FB-object
!# therefore covers a 'group of cells' inside the domain.
!# Two basic principles drive this idea:\\
!#
!# a) The fictitious boundary object' is not necessary one single object like
!#    a circle or a square. It's can even be a group of objects, all sharing
!#    the same boundary condition (e.g. 10000 flying balls, all described
!#    by a Dirichlet-type fictitious boundary object).\\
!#
!# b) All objects are described in 'analytical fashion' without any grid in the
!#    background. During the discretisation, the 'boundary conditions' that are
!#    associated to a fictitious boundary object (group) are 'discretised'
!#    to get a discrete, mesh- and solution dependent representation.
!#    By using filter techniques, this data can be imposed to the PDE, e.g.
!#    into the right hand side, the solution or into the matrix.\\
!#
!# Fictitious boundary objects are identified to the application by 
!# t_fictBoundaryRegion structures. The content of the structure is more or less
!# application specific; it simply contains a couple of information that allow
!# the evaluation routines to identify a fictitious boundary object.
!#
!#
!# Linear and nonlinear boundary conditions
!# ----------------------------------------
!# There are two types of boundary conditions: Linear (like Dirichlet) and
!# nonlinear ones (like Slip boundary conditions). Linear boundary conditions
!# are standard boundary conditions that are realised by discretising and
!# implementing them into solution/rhs/defect vectors during the linear
!# iteration. Nonlinear boundary conditions are implemented usually in a
!# nonlinear loop and need a special implementation there, which might
!# be different from when a linear system is solved.
!# 
!#
!# The following routines can be found here:\\
!#
!# 1.) bcond_initBC
!#     -> Initialises a boundary condition structure
!#
!# 2.) bcond_doneBC
!#     -> Cleans up a boundary condition structure, releases memory
!#
!# 3.) bcond_newBC
!#     -> Adds a general boundary condition for the real boundary to a 
!#        boundary-condition object; normally used only internally.
!#
!# 4.) bcond_getBCRegion
!#     -> Finds the index of the first boundary condition region containing a point.
!#
!# 5.) bcond_newDirichletBConRealBD
!#     -> Adds a Dirichlet boundary condition the real boundary to a 
!#        boundary-condition object
!#
!# 6.) bcond_newPressureDropBConRealBD
!#     -> Adds a pressure-drop boundary condition the real boundary to a 
!#        boundary-condition object
!#
!# 7.) bcond_newSlipBConRealBD
!#     -> Adds (nonlinear) slip boundary conditions on the real boundary
!#        to a boundary-condition object
!#
!# 8.) bcond_newFeastMirrorBConRealBD
!#     -> Adds a Dirichlet boundary condition the real boundary to a 
!#        boundary-condition object
!#
!# </purpose>
!##############################################################################

MODULE boundarycondition

  USE fsystem
  USE boundary
  USE fictitiousboundary
  
  IMPLICIT NONE

!<constants>

!<constantblock description="The type identifier for (linear) boundary conditions">

  ! Do-nothing boundary conditions (Neumann)
  INTEGER, PARAMETER :: BC_DONOTHING      = 0

  ! Dirichlet boundary conditions.
  ! Dirichlet boundary conditions are always specified for exactly one
  ! equation. t_bcRegion\%nequations is set =1 and t_bcRegion\%Iequations(1)
  ! identifies the number of the equation, the boundary conditions refer to.
  INTEGER, PARAMETER :: BC_DIRICHLET      = 1
  
  ! Robin boundary conditions
  INTEGER, PARAMETER :: BC_ROBIN          = 2
  
  ! Pressure-drop boundary conditions
  INTEGER, PARAMETER :: BC_PRESSUREDROP   = 3
  
  ! Flux boundary conditions
  INTEGER, PARAMETER :: BC_FLUX           = 4
  
  ! FEAST mirror boundary for domain decomposition
  INTEGER, PARAMETER :: BC_FEASTMIRROR    = 5

!</constantblock>

!<constantblock description="The type identifier for nonlinar boundary conditions">

  ! Slip boundary condition
  INTEGER, PARAMETER :: BC_SLIP           = 100

!</constantblock>

!<constantblock description="The type identifier for boundary regions">

  ! The boundary segment is unspecified.
  INTEGER, PARAMETER :: BC_RTYPE_UNDEFINED = 0
  
  ! The boundary dition region corresponds to a specific region on
  ! the real boundary.
  INTEGER, PARAMETER :: BC_RTYPE_REAL      = 1

  ! The boundary condition region 'lives' on the real boundary
  ! but does not correspond to a specific segment ('free' real
  ! boundar condition).
  INTEGER, PARAMETER :: BC_RTYPE_FREE      = 2
  
  ! The boundary condition region corresponds to a fictitious boundary
  ! object.
  INTEGER, PARAMETER :: BC_RTYPE_FBCOBJECT = 3
  
!</constantblock>

!<constantblock>
  
  ! Default blocksize for allocating new structures in the
  ! p_RregionsSpecific / p_RregionsFree list - if the list is full.
  INTEGER, PARAMETER :: BC_LISTBLOCKSIZE = 10
  
  ! Maximum number of equations that are supported simultaneously by
  ! boundary conditions.
  INTEGER, PARAMETER :: BC_MAXEQUATIONS = 16
  
!</constantblock>

!</constants>

!<types>
!<typeblock>
  
  ! The following structure realises a 'region' on the boundary where
  ! boundary conditions are defined. This is realised using the
  ! 'boundary region' structure from the 'boundary' module.
  ! The structure contains  user defined integer, double precision and
  ! string tag, which can be set by the application to arbitrary,
  ! problem dependent values (e.g. identifier tags for callback routines).
  
  TYPE t_bcRegion
  
    ! Type of BC region, this structure defines. One of the BC_RTYPE_xxxx
    ! constants.
    INTEGER :: cbcRegionType = BC_RTYPE_UNDEFINED
    
    ! Type of boundary conditions, this structure describes.
    ! One of the BC_xxxx constants.
    INTEGER :: ctype = BC_DONOTHING
    
    ! Number of equations, this BC refers to. This is usually =1 but can be
    ! higher if the BC couples multiple equations (like $u_x + u_y = c$ or so).
    ! Normally, this is the length of the Iequations list below, but the 
    ! meaning might be different for special type BC's.
    INTEGER :: nequations = 0
    
    ! A list of up to BC_MAXEQUATIONS numbers identifying the equations,
    ! a boundary condition refers to. The use of this is defined
    ! as follows:
    !
    ! a) A simple Dirichlet-boundary condition for the X-velocity on real
    !   boundary is identified by "nequations=1" + "Iequations=[1]".
    !
    ! b) A simple Dirichlet-boundary condition for the Y-velocity on real
    !   boundary is identified by "nequations=1" + "Iequations=[2]".
    !
    ! c) "Pressure drop" boundary conditions on real boundary that must 
    !   modify two "velocity" components 1(=x), 2(=y) are identified by
    !   "nequations=2" + "Iequations=[1 2]".
    !
    ! d) A simple Dirichlet boundary conditions for X- and Y-velocity
    !   for a fictitious boundary region is identified by "nequations=2"
    !   + "Iequations=[1 2]" (1=x, 2=y component)
    !
    ! e) A boundary condition like "u_x + u_y = const" on real boundary
    !   may be identified by "nequations=2" + "Iequations=[1 2]"
    !   (let's see if that's really the case if such a thing is implemented...)
    !
    ! Basically, this is a list of the blocks in a block solution vector
    ! that are affected by a boundary condition. The list here is actually
    ! bounadry-condition specific, i.e. another BC than Dirichlet can use
    ! the list in a different way to remember which equations take part
    ! on a special BC.
    INTEGER, DIMENSION(BC_MAXEQUATIONS) :: Iequations = 0
    
    ! Whether or not this boundary region is a 'static' region.
    ! If TRUE, the discretisation routine will discretise this boundary
    !  region only once and then never touch it anymore when called again.
    ! If FALSE, the discretisation routine assumes that the boundary
    !  region might have changed from one call to the other (e.g. in a
    !  time dependent simulation, the position or values might have changed) 
    !  and will therefore always rebuild the information for this boundary
    !  region when being called again.
    LOGICAL :: bisStatic = .FALSE.
    
    ! User defined tag to identify what to evaluate on the boundary.
    ! =0: undefined.
    ! This tag can be e.g. an identifier that tells the application whether
    ! to evaluate a constant, a parabolic profile or an expression
    ! on the boundary. Not used by the framework internally.
    INTEGER      :: ibdrexprtype = 0
    
    ! user defined integer tag
    INTEGER(I32) :: itag = 0
    
    ! User defined double precision tag
    REAL(DP) :: dtag = 0.0_DP
    
    ! User defined string tag; usually set to the name of an expression to
    ! evaluate in this boundary condition region.
    CHARACTER(LEN=SYS_STRLEN) :: stag = ''
     
    ! Definition of a region on the boundary where the boundary conditions
    ! 'live'. This can be either a 'specific' region (i.e. a region
    ! connected with a boundary segment) in case the boundary region
    ! corresponds to a boundary segment, or it can be a 'free' boundary
    ! region not corresponding to any special segment on the boundary.
    ! For boundary conditions of fictitious boundary objects,
    ! this structure is undefined!
    TYPE (t_boundaryRegion) :: rboundaryRegion
    
    ! Definition of a fictitious boundary region. If boundary conditions
    ! in a fictitious boundary region are discretised, this structure
    ! allows the application to identify the boundary region.
    ! For boundary conditions on the real boundary, this structure
    ! is undefined.
    TYPE (t_fictBoundaryRegion) :: rfictBoundaryRegion
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! Definition of a boundary condition. The structure basically contains
  ! a type-identifier specifying the type of boundary conditions
  ! the structure describes as well as a list of 'boundary region'
  ! segments that describe the position of the boundary condition
  ! on tghe boundary.
  
  TYPE t_boundaryConditions
  
    ! Pointer to the domain that is connected with this boundary 
    ! condition
    TYPE(t_boundary), POINTER :: p_rboundary => NULL()
    
    ! Number of regions in the list of boundary condition regions,
    ! corresponding to boundary regions.
    ! On each region, the boundary condition of type ctype
    ! is specified.
    INTEGER :: iregionCount = 0
    
    ! A list of t_BCregion structures that define the position
    ! of boundary condition ctype on the boundary.
    ! The first iregionCount elements in the array are valid,
    ! the rest may be undefined.
    TYPE(t_bcRegion), DIMENSION(:), POINTER :: p_Rregions => NULL()
    
    ! Number of regions in the list of boundary condition regions,
    ! not corresponding to specific boundary regions.
    ! On each region, the boundary condition of type ctype
    ! is specified.
    !INTEGER :: iregionCountFree = 0
    
    ! A list of t_BCregion structures that define the position
    ! of boundary condition ctype on the boundary.
    ! The first iregionCountFree elements in the array are valid,
    ! the rest may be undefined.
    !TYPE(t_bcRegion), DIMENSION(:), POINTER :: p_RregionsFree => NULL()
    
    ! Number of regions in the list of boundary condition regions,
    ! corresponding to fictitious boundary objects.
    ! On each region, the boundary condition of type ctype
    ! is specified.
    INTEGER :: iregionCountFBC = 0
    
    ! A list of t_BCregion structures that define the position
    ! of boundary condition ctype on the boundary.
    ! The first iregionCountFBC elements in the array are valid,
    ! the rest may be undefined.
    TYPE(t_bcRegion), DIMENSION(:), POINTER :: p_RregionsFBC => NULL()
    
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcond_initBC (p_rboundaryConditions,rboundary,ibcRegionsCount)
  
!<description>
  ! This routine initialises a boundary condition structure.
  ! If p_rboundaryConditions is NULL(), a new structure will be created. 
  ! Otherwise, the existing structure is recreated/updated.
!</description>

!<input>
  ! The domain which is to be connected to the boundary conditions.
  TYPE(t_boundary), INTENT(IN), TARGET :: rboundary
  
  ! OPTIONAL: The initial size of the lists saving boundary conditions.
  ! When adding boundary conditions to the rboundaryConditions structure,
  ! if there's not enough space, the lists saving the boundary conditions
  ! are dynamically increased (in terms of BC_LISTBLOCKSIZE). 
  INTEGER, INTENT(IN), OPTIONAL :: ibcRegionsCount
!</input>

!<output>
  ! The structure to be initialised.
  TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
!</output>

!</subroutine>

  ! local variables
  INTEGER ibcCount

  ! Do we have a structure?
  IF (.NOT. ASSOCIATED(p_rboundaryConditions)) THEN
    ALLOCATE(p_rboundaryConditions)
  ELSE
    ! Release the old structure without removing it from the heap.
    CALL bcond_doneBC(p_rboundaryConditions,.TRUE.)
  END IF

  ! The 'default constructor' does most of the necessary work, as
  ! rboundaryConditions is assumed as 'intent=out'. We only have to make
  ! the connection to the domain.
  
  p_rboundaryConditions%p_rboundary => rboundary
  
  ! Allocate memory for boundary condition lists
  ibcCount = BC_LISTBLOCKSIZE
  IF (PRESENT(ibcRegionsCount)) ibcCount = MAX(1,ibcRegionsCount)
  
  ALLOCATE(p_rboundaryConditions%p_Rregions(ibcCount))
  !ALLOCATE(p_rboundaryConditions%p_RregionsFree(ibcCount))
  ALLOCATE(p_rboundaryConditions%p_RregionsFBC(ibcCount))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE bcond_doneBC (p_rboundaryConditions,bkeepStructure)
  
!<description>
  ! This routine cleans up a boundary condition structure. All reserved
  ! memory is released.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the structure p_rboundaryConditions itself is not 
  ! released from memory. If set to FALSE or not existent (the usual setting), 
  ! the structure p_rboundaryConditions will also be removed from the heap after 
  ! cleaning up.
  LOGICAL, INTENT(IN), OPTIONAL :: bkeepStructure
!</input>

!<inputoutput>
  ! The structure to be cleaned up..
  TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
!</inputoutput>

!</subroutine>

  IF (.NOT. ASSOCIATED(p_rboundaryConditions)) RETURN

  ! Clean up
  DEALLOCATE(p_rboundaryConditions%p_RregionsFBC)
  !DEALLOCATE(p_rboundaryConditions%p_RregionsFree)
  DEALLOCATE(p_rboundaryConditions%p_Rregions)
  p_rboundaryConditions%iregionCountFBC = 0
  !p_rboundaryConditions%iregionCountFree = 0
  p_rboundaryConditions%iregionCount = 0
  p_rboundaryConditions%p_rboundary => NULL()

  ! Deallocate the structure (if we are allowed to), finish.
  IF (.NOT. PRESENT(bkeepStructure)) THEN
    DEALLOCATE(p_rboundaryConditions)
  ELSE
    IF (.NOT. bkeepStructure) DEALLOCATE(p_rboundaryConditions)
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcond_newBC (ctype,cbdtype,rboundaryConditions,&
                          p_rbcRegion, rboundaryRegion, rfictBoundaryRegion)
  
!<description>
  ! Adds a general boundary condition region to the boundary condition structure.. 
  ! A pointer to the region is returned
  ! in p_rbcRegion, the caller can add some user-defined information there
  ! if necessary (e.g. the name, a tag or something else).
  ! The boundary conditions are assumed to be on the 'real' boundary -
  ! fictitious-boundary boundary conditions are supported by another routine.
!</description>

!<input>
  ! The type of boundary conditions.
  ! This is a BC_xxxx flag.
  INTEGER, INTENT(IN) :: ctype

  ! Specifies the type of the boundary where to add boundary conditions.
  ! This is a BC_RTYPE_xxxx flag.
  INTEGER, INTENT(IN) :: cbdtype

  ! OPTIONAL: A boundary-region object, describing the position 
  ! on the boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  ! This parameter must be present for boundary conditions on the real
  ! boundary but can be omitted when adding a fictitious boundary condition.
  TYPE(t_boundaryRegion), INTENT(IN), OPTIONAL :: rboundaryRegion

  ! OPTIONAL: A fictitious-boundary-region object, describing the
  ! fictitious boundary region.
  ! A copy of this is added to the rboundaryConditions structure.
  ! This parameter must be present for boundary conditions on the fictitious
  ! boundary but can be omitted when adding boundary condition on the real 
  ! boundary.
  TYPE(t_fictBoundaryRegion), INTENT(IN), OPTIONAL :: rfictBoundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  TYPE(t_boundaryConditions), INTENT(INOUT), TARGET :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make more specific
  ! modifications to this.
  TYPE(t_bcRegion), OPTIONAL, POINTER :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  TYPE(t_bcRegion), DIMENSION(:), POINTER :: p_Rregion,p_temp
  TYPE(t_bcRegion), POINTER :: p_rbcRegionLocal
  INTEGER, POINTER :: ifull

  ! Select the type of boundary where to add:
  SELECT CASE (cbdtype)
  CASE (BC_RTYPE_REAL,BC_RTYPE_FREE)
    p_Rregion => rboundaryConditions%p_Rregions
    ifull => rboundaryConditions%iregionCount
  !CASE (BC_RTYPE_FREE)
  !  p_Rregion => rboundaryConditions%p_RregionsFree
  !  ifull => rboundaryConditions%iregionCountFree
  CASE (BC_RTYPE_FBCOBJECT)
    p_Rregion => rboundaryConditions%p_RregionsFBC
    ifull => rboundaryConditions%iregionCountFBC
  CASE DEFAULT
    PRINT *,'Not implemented boundary condition.'
    CALL sys_halt()
  END SELECT

  ! Space left, or do we have to reallocate?
  IF (ifull .GE. SIZE(p_Rregion)) THEN
    ALLOCATE(p_temp(ifull+BC_LISTBLOCKSIZE))
    p_temp(1:SIZE(p_Rregion)) = p_Rregion(:)
    DEALLOCATE(p_Rregion)
    p_Rregion => p_temp
    
    SELECT CASE (cbdtype)
    CASE (BC_RTYPE_REAL,BC_RTYPE_FREE)
      rboundaryConditions%p_Rregions => p_Rregion
    !CASE (BC_RTYPE_FREE)
    !  rboundaryConditions%p_RregionsFree => p_Rregion
    CASE (BC_RTYPE_FBCOBJECT)
      rboundaryConditions%p_RregionsFBC => p_Rregion
    END SELECT
  END IF
  
  ! Add the region
  ifull = ifull + 1
  p_rbcRegionLocal => p_Rregion(ifull)

  ! Initialise the structure
  p_rbcRegionLocal%cbcRegionType = cbdtype
  p_rbcRegionLocal%ctype = ctype
  
  IF (PRESENT(rboundaryRegion)) THEN
    p_rbcRegionLocal%rboundaryRegion = rboundaryRegion
  ELSE
    IF (cbdtype .NE. BC_RTYPE_FBCOBJECT) THEN
      PRINT *,'bcond_newBC: Boundary not specified'
    END IF
  END IF

  IF (PRESENT(rfictBoundaryRegion)) THEN
    p_rbcRegionLocal%rfictBoundaryRegion = rfictBoundaryRegion
  ELSE
    IF (cbdtype .EQ. BC_RTYPE_FBCOBJECT) THEN
      PRINT *,'bcond_newBC: Fictitious Boundary not specified'
    END IF
  END IF
  
  ! If p_rbcRegion is given, return the pointer
  IF (PRESENT(p_rbcRegion)) THEN
    p_rbcRegion => p_rbcRegionLocal
  END IF
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<function>

  INTEGER FUNCTION bcond_getBCRegion (rboundaryConditions,&
                                      iboundCompIdx,dparam,cparType,istartIndex)
  
!<description>
  ! The tupel (iboundCompIdx,dparam) specifies a point on the real 2D boundary.
  ! The routine tries to find the index of the first bonudary condition region
  ! that contains this point. This index is returned in iindex.
  ! Afterwards, the caller can access that boundary condition region
  ! using rboundaryConditions%p_Rregions(iindexBC).
!</description>

!<input>
  ! A structure containing all boundary condition regions on the real boundary.
  TYPE(t_boundaryConditions), INTENT(INOUT), TARGET :: rboundaryConditions

  ! The number of the boundary component of the point.
  INTEGER, INTENT(IN) :: iboundCompIdx

  ! The parameter value of the point to be checked.
  REAL(DP), INTENT(IN) :: dparam
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: cparType

  ! OPTIONAL: Start index in rboundaryConditions%p_Rregions(:) where to start
  ! the search for the point identified by (iboundCompIndex,dparam).
  ! If not specified, 1 is assumed, thus the routine will return the first
  ! boundary condition region that contains the point.
  ! If specified, the routine will search for the point in the boundary
  ! condition regions rboundaryConditions%p_Rregions(istartIndex:).
  ! This allows to skip some regions: E.g. if the caller assumes a point
  ! to be in multiple regions and one region is found, the caller can
  ! specify the number of the next region here where to continue the search.
  INTEGER, INTENT(IN), OPTIONAL :: istartIndex
!</input>

!<result>
  ! Number of the boundary region in rboundaryConditions%p_Rregions(:) that
  ! contains the point.
  ! =0 if no boundary condition region containing the point was found.
!</result>

!</function>

    REAL(DP) :: dparValue
    INTEGER :: cactParType, iindexBC
    
    cactParType = BDR_PAR_01
    IF (PRESENT(cparType)) THEN
      cactParType = cparType
    END IF

    ! Initialise iindexBC by istartIndex and use it as a counter.
    iindexBC = 1
    IF (PRESENT(istartIndex)) THEN
      IF (istartIndex .GT. 0) THEN
        iindexBC = istartIndex
      END IF
    END IF

    ! Loop through the boundary condition regions    
    DO WHILE (iindexBC .LE. rboundaryConditions%iregionCount) 
    
      ! Convert the parameter value to the correct parametrisation if necessary
      dparValue = boundary_convertParameter(rboundaryConditions%p_rboundary, &
          iboundCompIdx, dparam, &
          rboundaryConditions%p_Rregions(iindexBC)%rboundaryRegion%cparType, &
          cactParType)
          
      ! Check if it's inside of the region
      SELECT CASE (rboundaryConditions%p_Rregions(iindexBC)%cbcRegionType)
      CASE (BC_RTYPE_REAL,BC_RTYPE_FREE)

        IF (boundary_isInRegion ( &
            rboundaryConditions%p_Rregions(iindexBC)%rboundaryRegion,&
            iboundCompIdx,dparValue)) THEN
          ! Leave the subroutine. iindexBC has the index the caller is searching for.
          bcond_getBCRegion = iindexBC
          RETURN
        END IF

      END SELECT
      
      iindexBC = iindexBC+1
      
    END DO
              
    ! No region was found. Return 0.
    bcond_getBCRegion = 0

  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcond_newDirichletBConRealBD (rboundaryConditions,iequation,&
                                           rboundaryRegion,p_rbcRegion)
  
!<description>
  ! Adds a Dirichlet boundary condition region to the boundary condition 
  ! structure. A pointer to the region is returned in p_rbcRegion, the caller 
  ! can add some user-defined information there if necessary 
  ! (e.g. the name, a tag or something else).
  ! The boundary conditions are assumed to be on the 'real' boundary -
  ! fictitious-boundary boundary conditions are supported by another routine.
!</description>

!<input>
  ! An identifier for the equation, this boundary condition refers to.
  ! >= 1. 1=first equation (e.g. X-velocity), 2=2nd equation (e.g. 
  ! Y-velocity), etc.
  INTEGER, INTENT(IN) :: iequation

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  TYPE(t_boundaryRegion), INTENT(IN) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  TYPE(t_boundaryConditions), INTENT(INOUT), TARGET :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make 
  ! more specific modifications to this if necessary.
  TYPE(t_bcRegion), OPTIONAL, POINTER :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  TYPE(t_bcRegion), POINTER :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  IF (rboundaryRegion%iboundSegIdx .NE. 0) THEN
    CALL bcond_newBC (BC_DIRICHLET,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  ELSE
    CALL bcond_newBC (BC_DIRICHLET,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  END IF

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = 1
  p_rbcReg%Iequations(1) = iequation

  ! Eventually, return the pointer
  IF (PRESENT(p_rbcRegion)) p_rbcRegion => p_rbcReg

  END SUBROUTINE
      
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcond_newFeastMirrorBConRealBD (rboundaryConditions,iequation,&
                                             rboundaryRegion,p_rbcRegion)
  
!<description>
  ! Adds a FEAST mirror boundary condition region to the boundary condition 
  ! structure. A pointer to the region is returned in p_rbcRegion, the caller 
  ! can add some user-defined information there if necessary 
  ! (e.g. the name, a tag or something else).
  ! The boundary conditions are assumed to be on the 'real' boundary -
  ! fictitious-boundary boundary conditions are supported by another routine.
!</description>

!<input>
  ! An identifier for the equation, this boundary condition refers to.
  ! >= 1. 1=first equation (e.g. X-velocity), 2=2nd equation (e.g. 
  ! Y-velocity), etc.
  INTEGER, INTENT(IN) :: iequation

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  TYPE(t_boundaryRegion), INTENT(IN) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  TYPE(t_boundaryConditions), INTENT(INOUT), TARGET :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make 
  ! more specific modifications to this if necessary.
  TYPE(t_bcRegion), OPTIONAL, POINTER :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  TYPE(t_bcRegion), POINTER :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  IF (rboundaryRegion%iboundSegIdx .NE. 0) THEN
    CALL bcond_newBC (BC_FEASTMIRROR,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  ELSE
    CALL bcond_newBC (BC_FEASTMIRROR,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  END IF

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = 1
  p_rbcReg%Iequations(1) = iequation

  ! Eventually, return the pointer
  IF (PRESENT(p_rbcRegion)) p_rbcRegion => p_rbcReg

  END SUBROUTINE
      
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcond_newPressureDropBConRealBD (rboundaryConditions,IvelEqns,&
                                              rboundaryRegion,p_rbcRegion)
  
!<description>
  ! Adds a Pressure-Drop boundary condition region to the boundary condition 
  ! structure. A pointer to the region is returned in p_rbcRegion, the caller 
  ! can add some user-defined information there if necessary 
  ! (e.g. the name, a tag or something else).
!</description>

!<input>
  ! A list of identifiers for the velocity equations, this boundary condition
  ! modifies. Usually (1,2) for X- and Y-velocity.
  INTEGER, DIMENSION(:), INTENT(IN) :: IvelEqns

  ! A boundary-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  TYPE(t_boundaryRegion), INTENT(IN) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  TYPE(t_boundaryConditions), INTENT(INOUT), TARGET :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make more specific
  ! modifications to this.
  TYPE(t_bcRegion), OPTIONAL, POINTER :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  TYPE(t_bcRegion), POINTER :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  IF (rboundaryRegion%iboundSegIdx .NE. 0) THEN
    CALL bcond_newBC (BC_PRESSUREDROP,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  ELSE
    CALL bcond_newBC (BC_PRESSUREDROP,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  END IF

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = SIZE(IvelEqns)
  p_rbcReg%Iequations(1:SIZE(IvelEqns)) = IvelEqns

  ! Eventually, return the pointer
  IF (PRESENT(p_rbcRegion)) p_rbcRegion => p_rbcReg

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcond_newSlipBConRealBD (rboundaryConditions,IvelEqns,&
                                      rboundaryRegion,p_rbcRegion)
  
!<description>
  ! Adds a nonlinear Slip boundary condition region to the boundary condition 
  ! structure. A pointer to the region is returned in p_rbcRegion, the caller 
  ! can add some user-defined information there if necessary 
  ! (e.g. the name, a tag or something else).
!</description>

!<input>
  ! A list of identifiers for the velocity equations, this boundary condition
  ! modifies. Usually (1,2) for X- and Y-velocity.
  INTEGER, DIMENSION(:), INTENT(IN) :: IvelEqns

  ! A boundary-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  TYPE(t_boundaryRegion), INTENT(IN) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  TYPE(t_boundaryConditions), INTENT(INOUT), TARGET :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make more specific
  ! modifications to this.
  TYPE(t_bcRegion), OPTIONAL, POINTER :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  TYPE(t_bcRegion), POINTER :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  IF (rboundaryRegion%iboundSegIdx .NE. 0) THEN
    CALL bcond_newBC (BC_SLIP,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  ELSE
    CALL bcond_newBC (BC_SLIP,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  END IF

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = SIZE(IvelEqns)
  p_rbcReg%Iequations(1:SIZE(IvelEqns)) = IvelEqns

  ! Eventually, return the pointer
  IF (PRESENT(p_rbcRegion)) p_rbcRegion => p_rbcReg

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE bcond_newDirichletBConFictBD (rboundaryConditions,Iequations,&
                                           rboundaryRegion,p_rbcRegion)
  
!<description>
  ! Adds a Dirichlet boundary condition region for a fictitious boundary 
  ! component to the boundary condition structure. 
  ! A pointer to the region is returned in p_rfbcRegion, the caller 
  ! can add some user-defined information there if necessary 
  ! (e.g. the name, a tag or something else).
!</description>

!<input>
  ! An array of identifiers for the equations, this boundary condition 
  ! refers to. Example: Iequations = [1 2] for X-velocity-component (1) and
  ! Y-velocity component (2).
  INTEGER, DIMENSION(:), INTENT(IN) :: Iequations

  ! A fictitious-boundary-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  TYPE(t_fictBoundaryRegion), INTENT(IN) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  TYPE(t_boundaryConditions), INTENT(INOUT), TARGET :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary condition region. The caller can 
  ! make more specific modifications to this if necessary.
  TYPE(t_bcRegion), OPTIONAL, POINTER :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  TYPE(t_bcRegion), POINTER :: p_rbcReg

  ! Add a general boundary condition region for a fictitiouos boundary object.
  CALL bcond_newBC (BC_DIRICHLET,BC_RTYPE_FBCOBJECT,rboundaryConditions,&
                    p_rbcReg,rfictBoundaryRegion=rboundaryRegion)

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = SIZE(Iequations)
  p_rbcReg%Iequations(1:SIZE(Iequations)) = Iequations

  ! Eventually, return the pointer
  IF (PRESENT(p_rbcRegion)) p_rbcRegion => p_rbcReg

  END SUBROUTINE
      
END MODULE

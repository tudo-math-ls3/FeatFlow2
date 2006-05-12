!##############################################################################
!# ****************************************************************************
!# <name> scalarbc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module defines the analytical representation of boundary conditions
!# for 1D scalar equations. In detain, the following types of boundary
!# conditions are supported:
!# - Neumann boundary conditions
!# - Dirichlet boundary conditions
!# The structure here only describe the boundary conditions. There are no
!# routines in this module to impose them into a solution/defect vector.
!# For imposing boundary conditions into vectors, the boundary conditions
!# first have to be 'discretised' using the module 'DiscreteBC'!
!#
!# Neumann boundary conditions are usually imposed as 'do nothing' boundary
!# conditions and therefore not realised in a special structure. This
!# means, the module here realises all 'other' types of boundary conditions
!# not being Neumann. Neumann is standard, and all other types of BC's
!# must be imposed.
!#
!# Method description
!# ------------------
!# The t_boundaryConditions structure collects all boundary conditions to a
!# domain. It contains three lists of 'boundary condition regions'
!# (structure t_bcRegion), for
!# - boundary conditions bound to segments on the real boundary
!# - boundary conditions anywhere on the real boundary
!# - boundary conditions, imposed by fictitious boundary objects.
!# A 'boundary condition region' object t_bcRegion realises a boundary
!# condition. It contains information about the type of the boundary conditions
!# and where they are, on the real boundary as well as on fictitious boundary
!# objects.
!# A boundary condition can be deactivated by setting the cbcRegionType flag
!# in the boundary region to BC_DONOTHING.
!#
!# The following routines can be found here:
!#
!# 1.) scbc_initScalarBC
!#     -> Initialises a boundary condition structure
!#
!# 2.) scbc_doneScalarBC
!#     -> Cleans up a boundary condition structure, releases memory
!#
!# 3.) scbc_newBConRealBD
!#     -> Adds a boundary condition for the real boundary to a 
!#        boundary-condition object
!#
!# </purpose>
!##############################################################################

MODULE scalarbc

  USE fsystem
  USE boundary
  
  IMPLICIT NONE

!<constants>

!<constantblock description="The type identifier for boundary conditions">

  ! Do-nothing boundary conditions (Neumann)
  INTEGER, PARAMETER :: BC_DONOTHING      = 0

  ! Dirichlet boundary conditions
  INTEGER, PARAMETER :: BC_DIRICHLET      = 1
  
  ! Robin boundary conditions
  INTEGER, PARAMETER :: BC_ROBIN          = 2
  
  ! Flux boundary conditions
  INTEGER, PARAMETER :: BC_FLUX           = 3

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
  
  ! The boundary condition region corresponds to a fictitious boundar
  ! object.
  INTEGER, PARAMETER :: BC_RTYPE_FBCOBJECT = 3
  
!</constantblock>

!<constantblock>
  
  ! Default blocksize for allocating new structures in the
  ! p_RregionsSpecific / p_RregionsFree list - if the list is full.
  INTEGER, PARAMETER :: BC_LISTBLOCKSIZE = 10
  
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
    
    ! Whether or not this boundary region is a 'static' region.
    ! If TRUE, the discretisation routine will discretise this boundary
    !  region only once and then never touch it anymore when called again.
    ! If FALSE, the discretisation routine assumes that the boundary
    !  region might have changed from one call to the other (e.g. in a
    !  time dependent simulation, the position or values might have changed) 
    !  and will therefore always rebuild the information for this boundary
    !  region when being called again.
    LOGICAL :: bisStatic = .FALSE.
    
    ! user defined integer tag
    INTEGER(I32) :: itag = 0
    
    ! uesr defined double precision tag
    REAL(DP) :: dtag = 0.0_DP
    
    ! user defined string tag
    CHARACTER(LEN=SYS_STRLEN) :: stag = ''
     
    ! Definition of a region on the boundary where the boundary conditions
    ! 'live'. This can be either a 'specific' region (i.e. a region
    ! connected with a boundary segment) in case the boundary region
    ! corresponds to a boundary segment, or it can be a 'free' boundary
    ! region not corresponding to any special segment on the boundary.
    ! For boundary conditions of fictitious boundary objects,
    ! the structure is undefined!
    TYPE (t_boundaryRegion) :: rboundaryRegion
    
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
    TYPE(t_boundary), POINTER :: rdomain => NULL()
    
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

  SUBROUTINE scbc_initScalarBC (rboundaryConditions,rdomain,ibcRegionsCount)
  
!<description>
  ! This routine initialises a boundary condition structure.
!</description>

!<input>
  ! The domain which is to be connected to the boundary conditions.
  TYPE(t_boundary), INTENT(IN), TARGET :: rdomain
  
  ! OPTIONAL: The initial size of the lists saving boundary conditions.
  ! When adding boundary conditions to the rboundaryConditions structure,
  ! if there's not enough space, the lists saving the boundary conditions
  ! are dynamically increased (in terms of BC_LISTBLOCKSIZE). 
  INTEGER, INTENT(IN), OPTIONAL :: ibcRegionsCount
!</input>

!<output>
  ! The structure to be initialised.
  TYPE(t_boundaryConditions), INTENT(OUT) :: rboundaryConditions
!</output>

!</subroutine>

  ! local variables
  INTEGER ibcCount

  ! The 'default constructor' does most of the necessary work, as
  ! rboundaryConditions is assumed as 'intent=out'. We only have to make
  ! the connection to the domain.
  
  rboundaryConditions%rdomain => rdomain
  
  ! Allocate memory for boundary condition lists
  ibcCount = BC_LISTBLOCKSIZE
  IF (PRESENT(ibcRegionsCount)) ibcCount = MAX(1,ibcRegionsCount)
  
  ALLOCATE(rboundaryConditions%p_Rregions(ibcCount))
  !ALLOCATE(rboundaryConditions%p_RregionsFree(ibcCount))
  ALLOCATE(rboundaryConditions%p_RregionsFBC(ibcCount))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE scbc_doneScalarBC (rboundaryConditions)
  
!<description>
  ! This routine cleans up a boundary condition structure. All reserved
  ! memory is released.
!</description>

!<inputoutput>
  ! The structure to be initialised.
  TYPE(t_boundaryConditions), INTENT(INOUT) :: rboundaryConditions
!</inputoutput>

!</subroutine>

  DEALLOCATE(rboundaryConditions%p_RregionsFBC)
  !DEALLOCATE(rboundaryConditions%p_RregionsFree)
  DEALLOCATE(rboundaryConditions%p_Rregions)
  rboundaryConditions%iregionCountFBC = 0
  !rboundaryConditions%iregionCountFree = 0
  rboundaryConditions%iregionCount = 0
  rboundaryConditions%rdomain => NULL()

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE scbc_newBConRealBD (ctype,cbdtype,rboundaryConditions,rboundaryRegion,p_rbcRegion)
  
!<description>
  ! Adds a boundary condition region to the boundary condition structure,
  ! configured for Dirichlet boundary. A pointer to the region is returned
  ! in p_rbcRegion, the caller cann add some user-defined information there
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
  ! A pointer to the added boundary region. The caller can make more specific
  ! modifications to this.
  TYPE(t_bcRegion), POINTER :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  TYPE(t_bcRegion), DIMENSION(:), POINTER :: p_region,p_temp
  INTEGER, POINTER :: ifull

  ! Select the type of boundary where to add:
  SELECT CASE (cbdtype)
  CASE (BC_RTYPE_REAL,BC_RTYPE_FREE)
    p_region => rboundaryConditions%p_Rregions
    ifull => rboundaryConditions%iregionCount
  !CASE (BC_RTYPE_FREE)
  !  p_region => rboundaryConditions%p_RregionsFree
  !  ifull => rboundaryConditions%iregionCountFree
  CASE (BC_RTYPE_FBCOBJECT)
    PRINT *,'Fictitious boundary conditions not supported by scbc_newBCrealBD!'
    STOP
  CASE DEFAULT
    PRINT *,'Not implemented boundary condition.'
    STOP
  END SELECT

  ! Space left, or do we have to reallocate?
  IF (ifull .GE. SIZE(p_region)) THEN
    ALLOCATE(p_temp(ifull+BC_LISTBLOCKSIZE))
    p_temp = p_region
    DEALLOCATE(p_region)
    p_region => p_temp
    
    SELECT CASE (cbdtype)
    CASE (BC_RTYPE_REAL,BC_RTYPE_FREE)
      rboundaryConditions%p_Rregions => p_region
    !CASE (BC_RTYPE_FREE)
    !  rboundaryConditions%p_RregionsFree => p_region
    CASE (BC_RTYPE_FBCOBJECT)
      rboundaryConditions%p_RregionsFBC => p_region
    END SELECT
  END IF
  
  ! Add the region
  ifull = ifull + 1
  p_rbcRegion => p_region(ifull)

  ! Initialise the structure
  p_rbcRegion%cbcRegionType = cbdtype
  p_rbcRegion%ctype = ctype
  p_rbcRegion%rboundaryRegion = rboundaryRegion
  
  END SUBROUTINE
      
END MODULE

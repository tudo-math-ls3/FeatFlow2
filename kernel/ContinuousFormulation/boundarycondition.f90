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
!# THIS MODULE IS OUTDATED AND NOT USED ANYMORE! IT IS STILL CONTENT OF THE
!# FEAT2 LIBRARY TO ALLOW A USER TO DEFINE ANALYTIC BONUDARY CONDITIONS ON
!# APPLICATION LEVEL. THE KERNEL ROUTINES HOWEVER DO NOT USE ANALYTICAL
!# BOUNDARY CONDITIONS ANYMORE: BOUNDARY CONDITIONS ARE ALWAYS DIRECTLY
!# ASSEMBLED INTO THEIR DISCRETE COUNTERPART SUCH THAT THEY CAN BE IMPLEMENTED
!# INTO MATRICES/VECTORS!
!#
!# Neumann boundary conditions are usually imposed as 'do nothing' boundary
!# conditions and therefore not realised in a special structure. This
!# means, the module here realises all 'other' types of boundary conditions
!# not being Neumann. Neumann is standard, and all other types of BC`s
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
!# BC`s on real boundary components in 2D \\
!# -------------------------------------- \\
!# Boundary conditions on real boundary components are described via a boundary
!# region struture as defined in boundary.f90. This contains basically a minimum
!# and maximum parameter value on the boundary as well as s number of a
!# boundary component. Optionally, the boundary region may be specified to
!# be attached to a real boundary segment, e.g. a full circle or a line on
!# the boundary.
!#
!# BC`s on fictitious boundary components \\
!# -------------------------------------- \\
!# The concept of (analytic) fictitious boundary components is somehow
!# different to the concept of standard bondary conditions.
!# An analytic fictitious boundary component is a 'virtual object' in the
!# domain that is associated to a special boundary condition. Each FB-object
!# therefore covers a 'group of cells' inside the domain.
!# Two basic principles drive this idea:\\
!#
!# a) The fictitious boundary object is not necessarily one single object like
!#    a circle or a square. It can even be a group of objects, all sharing
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

module boundarycondition

!$use omp_lib
  use fsystem
  use boundary
  use fictitiousboundary
  
  implicit none

  private
  public :: t_bcRegion
  public :: t_boundaryConditions
  public :: bcond_initBC
  public :: bcond_doneBC
  public :: bcond_newBC
  public :: bcond_getBCRegion
  public :: bcond_newDirichletBConRealBD
  public :: bcond_newFeastMirrorBConRealBD
  public :: bcond_newPressureDropBConRealBD
  public :: bcond_newSlipBConRealBD
  public :: bcond_newDirichletBConFictBD

!<constants>

!<constantblock description="The type identifier for (linear) boundary conditions">

  ! Do-nothing boundary conditions (Neumann)
  integer, parameter, public :: BC_DONOTHING      = 0

  ! Dirichlet boundary conditions.
  ! Dirichlet boundary conditions are always specified for exactly one
  ! equation. t_bcRegion\%nequations is set =1 and t_bcRegion\%Iequations(1)
  ! identifies the number of the equation, the boundary conditions refer to.
  integer, parameter, public :: BC_DIRICHLET      = 1
  
  ! Robin boundary conditions
  integer, parameter, public :: BC_ROBIN          = 2
  
  ! Pressure-drop boundary conditions
  integer, parameter, public :: BC_PRESSUREDROP   = 3
  
  ! Flux boundary conditions
  integer, parameter, public :: BC_FLUX           = 4
  
  ! FEAST mirror boundary for domain decomposition
  integer, parameter, public :: BC_FEASTMIRROR    = 5

!</constantblock>

!<constantblock description="The type identifier for nonlinar boundary conditions">

  ! Slip boundary condition
  integer, parameter, public :: BC_SLIP           = 100

!</constantblock>

!<constantblock description="The type identifier for boundary regions">

  ! The boundary segment is unspecified.
  integer, parameter, public :: BC_RTYPE_UNDEFINED = 0
  
  ! The boundary dition region corresponds to a specific region on
  ! the real boundary.
  integer, parameter, public :: BC_RTYPE_REAL      = 1

  ! The boundary condition region 'lives' on the real boundary
  ! but does not correspond to a specific segment ('free' real
  ! boundar condition).
  integer, parameter, public :: BC_RTYPE_FREE      = 2
  
  ! The boundary condition region corresponds to a fictitious boundary
  ! object.
  integer, parameter, public :: BC_RTYPE_FBCOBJECT = 3
  
!</constantblock>

!<constantblock>
  
  ! Default blocksize for allocating new structures in the
  ! p_RregionsSpecific / p_RregionsFree list - if the list is full.
  integer, parameter, public :: BC_LISTBLOCKSIZE = 10
  
  ! Maximum number of equations that are supported simultaneously by
  ! boundary conditions.
  integer, parameter, public :: BC_MAXEQUATIONS = 16
  
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
  
  type t_bcRegion
  
    ! Type of BC region, this structure defines. One of the BC_RTYPE_xxxx
    ! constants.
    integer :: cbcRegionType = BC_RTYPE_UNDEFINED
    
    ! Type of boundary conditions, this structure describes.
    ! One of the BC_xxxx constants.
    integer :: ctype = BC_DONOTHING
    
    ! Number of equations, this BC refers to. This is usually =1 but can be
    ! higher if the BC couples multiple equations (like <tex>$u_x + u_y = c$</tex> or so).
    ! Normally, this is the length of the Iequations list below, but the
    ! meaning might be different for special type BC`s.
    integer :: nequations = 0
    
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
    !   (let us see if that is really the case if such a thing is implemented...)
    !
    ! Basically, this is a list of the blocks in a block solution vector
    ! that are affected by a boundary condition. The list here is actually
    ! bounadry-condition specific, i.e. another BC than Dirichlet can use
    ! the list in a different way to remember which equations take part
    ! on a special BC.
    integer, dimension(BC_MAXEQUATIONS) :: Iequations = 0
    
    ! Whether or not this boundary region is a 'static' region.
    ! If TRUE, the discretisation routine will discretise this boundary
    !  region only once and then never touch it anymore when called again.
    ! If FALSE, the discretisation routine assumes that the boundary
    !  region might have changed from one call to the other (e.g. in a
    !  time dependent simulation, the position or values might have changed)
    !  and will therefore always rebuild the information for this boundary
    !  region when being called again.
    logical :: bisStatic = .false.
    
    ! User defined tag to identify what to evaluate on the boundary.
    ! =0: undefined.
    ! This tag can be e.g. an identifier that tells the application whether
    ! to evaluate a constant, a parabolic profile or an expression
    ! on the boundary. Not used by the framework internally.
    integer      :: ibdrexprtype = 0
    
    ! user defined integer tag
    integer(I32) :: itag = 0
    
    ! User defined double precision tag
    real(DP) :: dtag = 0.0_DP
    
    ! User defined string tag; usually set to the name of an expression to
    ! evaluate in this boundary condition region.
    character(LEN=SYS_STRLEN) :: stag = ''
     
    ! Definition of a region on the boundary where the boundary conditions
    ! 'live'. This can be either a 'specific' region (i.e. a region
    ! connected with a boundary segment) in case the boundary region
    ! corresponds to a boundary segment, or it can be a 'free' boundary
    ! region not corresponding to any special segment on the boundary.
    ! For boundary conditions of fictitious boundary objects,
    ! this structure is undefined!
    type (t_boundaryRegion) :: rboundaryRegion
    
    ! Definition of a fictitious boundary region. If boundary conditions
    ! in a fictitious boundary region are discretised, this structure
    ! allows the application to identify the boundary region.
    ! For boundary conditions on the real boundary, this structure
    ! is undefined.
    type (t_fictBoundaryRegion) :: rfictBoundaryRegion
    
  end type
  
!</typeblock>

!<typeblock>
  
  ! Definition of a boundary condition. The structure basically contains
  ! a type-identifier specifying the type of boundary conditions
  ! the structure describes as well as a list of 'boundary region'
  ! segments that describe the position of the boundary condition
  ! on tghe boundary.
  
  type t_boundaryConditions
  
    ! Pointer to the domain that is connected with this boundary
    ! condition
    type(t_boundary), pointer :: p_rboundary => null()
    
    ! Number of regions in the list of boundary condition regions,
    ! corresponding to boundary regions.
    ! On each region, the boundary condition of type ctype
    ! is specified.
    integer :: iregionCount = 0
    
    ! A list of t_BCregion structures that define the position
    ! of boundary condition ctype on the boundary.
    ! The first iregionCount elements in the array are valid,
    ! the rest may be undefined.
    type(t_bcRegion), dimension(:), pointer :: p_Rregions => null()
    
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
    integer :: iregionCountFBC = 0
    
    ! A list of t_BCregion structures that define the position
    ! of boundary condition ctype on the boundary.
    ! The first iregionCountFBC elements in the array are valid,
    ! the rest may be undefined.
    type(t_bcRegion), dimension(:), pointer :: p_RregionsFBC => null()
    
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine bcond_initBC (p_rboundaryConditions,rboundary,ibcRegionsCount)
  
!<description>
  ! This routine initialises a boundary condition structure.
  ! If p_rboundaryConditions is NULL(), a new structure will be created.
  ! Otherwise, the existing structure is recreated/updated.
!</description>

!<input>
  ! The domain which is to be connected to the boundary conditions.
  type(t_boundary), intent(in), target :: rboundary
  
  ! OPTIONAL: The initial size of the lists saving boundary conditions.
  ! When adding boundary conditions to the rboundaryConditions structure,
  ! if there is not enough space, the lists saving the boundary conditions
  ! are dynamically increased (in terms of BC_LISTBLOCKSIZE).
  integer, intent(in), optional :: ibcRegionsCount
!</input>

!<output>
  ! The structure to be initialised.
  type(t_boundaryConditions), pointer :: p_rboundaryConditions
!</output>

!</subroutine>

  ! local variables
  integer ibcCount

  ! Do we have a structure?
  if (.not. associated(p_rboundaryConditions)) then
    allocate(p_rboundaryConditions)
  else
    ! Release the old structure without removing it from the heap.
    call bcond_doneBC(p_rboundaryConditions,.true.)
  end if

  ! The 'default constructor' does most of the necessary work, as
  ! rboundaryConditions is assumed as 'intent=out'. We only have to make
  ! the connection to the domain.
  
  p_rboundaryConditions%p_rboundary => rboundary
  
  ! Allocate memory for boundary condition lists
  ibcCount = BC_LISTBLOCKSIZE
  if (present(ibcRegionsCount)) ibcCount = max(1,ibcRegionsCount)
  
  allocate(p_rboundaryConditions%p_Rregions(ibcCount))
  !ALLOCATE(p_rboundaryConditions%p_RregionsFree(ibcCount))
  allocate(p_rboundaryConditions%p_RregionsFBC(ibcCount))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcond_doneBC (p_rboundaryConditions,bkeepStructure)
  
!<description>
  ! This routine cleans up a boundary condition structure. All reserved
  ! memory is released.
!</description>

!<input>
  ! OPTIONAL: If set to TRUE, the structure p_rboundaryConditions itself is not
  ! released from memory. If set to FALSE or not existent (the usual setting),
  ! the structure p_rboundaryConditions will also be removed from the heap after
  ! cleaning up.
  logical, intent(in), optional :: bkeepStructure
!</input>

!<inputoutput>
  ! The structure to be cleaned up..
  type(t_boundaryConditions), pointer :: p_rboundaryConditions
!</inputoutput>

!</subroutine>

  if (.not. associated(p_rboundaryConditions)) return

  ! Clean up
  deallocate(p_rboundaryConditions%p_RregionsFBC)
  !DEALLOCATE(p_rboundaryConditions%p_RregionsFree)
  deallocate(p_rboundaryConditions%p_Rregions)
  p_rboundaryConditions%iregionCountFBC = 0
  !p_rboundaryConditions%iregionCountFree = 0
  p_rboundaryConditions%iregionCount = 0
  p_rboundaryConditions%p_rboundary => null()

  ! Deallocate the structure (if we are allowed to), finish.
  if (.not. present(bkeepStructure)) then
    deallocate(p_rboundaryConditions)
  else
    if (.not. bkeepStructure) deallocate(p_rboundaryConditions)
  end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcond_newBC (ctype,cbdtype,rboundaryConditions,&
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
  integer, intent(in) :: ctype

  ! Specifies the type of the boundary where to add boundary conditions.
  ! This is a BC_RTYPE_xxxx flag.
  integer, intent(in) :: cbdtype

  ! OPTIONAL: A boundary-region object, describing the position
  ! on the boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  ! This parameter must be present for boundary conditions on the real
  ! boundary but can be omitted when adding a fictitious boundary condition.
  type(t_boundaryRegion), intent(in), optional :: rboundaryRegion

  ! OPTIONAL: A fictitious-boundary-region object, describing the
  ! fictitious boundary region.
  ! A copy of this is added to the rboundaryConditions structure.
  ! This parameter must be present for boundary conditions on the fictitious
  ! boundary but can be omitted when adding boundary condition on the real
  ! boundary.
  type(t_fictBoundaryRegion), intent(in), optional :: rfictBoundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  type(t_boundaryConditions), intent(inout), target :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make more specific
  ! modifications to this.
  type(t_bcRegion), optional, pointer :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  type(t_bcRegion), dimension(:), pointer :: p_Rregion,p_temp
  type(t_bcRegion), pointer :: p_rbcRegionLocal
  integer, pointer :: ifull

  ! Select the type of boundary where to add:
  select case (cbdtype)
  case (BC_RTYPE_REAL,BC_RTYPE_FREE)
    p_Rregion => rboundaryConditions%p_Rregions
    ifull => rboundaryConditions%iregionCount
  !CASE (BC_RTYPE_FREE)
  !  p_Rregion => rboundaryConditions%p_RregionsFree
  !  ifull => rboundaryConditions%iregionCountFree
  case (BC_RTYPE_FBCOBJECT)
    p_Rregion => rboundaryConditions%p_RregionsFBC
    ifull => rboundaryConditions%iregionCountFBC
  case DEFAULT
    print *,'Not implemented boundary condition.'
    call sys_halt()
  end select

  ! Space left, or do we have to reallocate?
  if (ifull .ge. size(p_Rregion)) then
    allocate(p_temp(ifull+BC_LISTBLOCKSIZE))
    p_temp(1:size(p_Rregion)) = p_Rregion(:)
    deallocate(p_Rregion)
    p_Rregion => p_temp
    
    select case (cbdtype)
    case (BC_RTYPE_REAL,BC_RTYPE_FREE)
      rboundaryConditions%p_Rregions => p_Rregion
    !CASE (BC_RTYPE_FREE)
    !  rboundaryConditions%p_RregionsFree => p_Rregion
    case (BC_RTYPE_FBCOBJECT)
      rboundaryConditions%p_RregionsFBC => p_Rregion
    end select
  end if
  
  ! Add the region
  ifull = ifull + 1
  p_rbcRegionLocal => p_Rregion(ifull)

  ! Initialise the structure
  p_rbcRegionLocal%cbcRegionType = cbdtype
  p_rbcRegionLocal%ctype = ctype
  
  if (present(rboundaryRegion)) then
    p_rbcRegionLocal%rboundaryRegion = rboundaryRegion
  else
    if (cbdtype .ne. BC_RTYPE_FBCOBJECT) then
      print *,'bcond_newBC: Boundary not specified'
    end if
  end if

  if (present(rfictBoundaryRegion)) then
    p_rbcRegionLocal%rfictBoundaryRegion = rfictBoundaryRegion
  else
    if (cbdtype .eq. BC_RTYPE_FBCOBJECT) then
      print *,'bcond_newBC: Fictitious Boundary not specified'
    end if
  end if
  
  ! If p_rbcRegion is given, return the pointer
  if (present(p_rbcRegion)) then
    p_rbcRegion => p_rbcRegionLocal
  end if
  
  end subroutine

  ! ***************************************************************************
  
!<function>

  integer function bcond_getBCRegion (rboundaryConditions,&
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
  type(t_boundaryConditions), intent(inout), target :: rboundaryConditions

  ! The number of the boundary component of the point.
  integer, intent(in) :: iboundCompIdx

  ! The parameter value of the point to be checked.
  real(DP), intent(in) :: dparam
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(in), optional :: cparType

  ! OPTIONAL: Start index in rboundaryConditions%p_Rregions(:) where to start
  ! the search for the point identified by (iboundCompIndex,dparam).
  ! If not specified, 1 is assumed, thus the routine will return the first
  ! boundary condition region that contains the point.
  ! If specified, the routine will search for the point in the boundary
  ! condition regions rboundaryConditions%p_Rregions(istartIndex:).
  ! This allows to skip some regions: E.g. if the caller assumes a point
  ! to be in multiple regions and one region is found, the caller can
  ! specify the number of the next region here where to continue the search.
  integer, intent(in), optional :: istartIndex
!</input>

!<result>
  ! Number of the boundary region in rboundaryConditions%p_Rregions(:) that
  ! contains the point.
  ! =0 if no boundary condition region containing the point was found.
!</result>

!</function>

    real(DP) :: dparValue
    integer :: cactParType, iindexBC
    
    cactParType = BDR_PAR_01
    if (present(cparType)) then
      cactParType = cparType
    end if

    ! Initialise iindexBC by istartIndex and use it as a counter.
    iindexBC = 1
    if (present(istartIndex)) then
      if (istartIndex .gt. 0) then
        iindexBC = istartIndex
      end if
    end if

    ! Loop through the boundary condition regions
    do while (iindexBC .le. rboundaryConditions%iregionCount)
    
      ! Convert the parameter value to the correct parametrisation if necessary
      dparValue = boundary_convertParameter(rboundaryConditions%p_rboundary, &
          iboundCompIdx, dparam, &
          rboundaryConditions%p_Rregions(iindexBC)%rboundaryRegion%cparType, &
          cactParType)
          
      ! Check if it is inside of the region
      select case (rboundaryConditions%p_Rregions(iindexBC)%cbcRegionType)
      case (BC_RTYPE_REAL,BC_RTYPE_FREE)

        if (boundary_isInRegion ( &
            rboundaryConditions%p_Rregions(iindexBC)%rboundaryRegion,&
            iboundCompIdx,dparValue)) then
          ! Leave the subroutine. iindexBC has the index the caller is searching for.
          bcond_getBCRegion = iindexBC
          return
        end if

      end select
      
      iindexBC = iindexBC+1
      
    end do
              
    ! No region was found. Return 0.
    bcond_getBCRegion = 0

  end function

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcond_newDirichletBConRealBD (rboundaryConditions,iequation,&
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
  integer, intent(in) :: iequation

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  type(t_boundaryConditions), intent(inout), target :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make
  ! more specific modifications to this if necessary.
  type(t_bcRegion), optional, pointer :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  type(t_bcRegion), pointer :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  if (rboundaryRegion%iboundSegIdx .ne. 0) then
    call bcond_newBC (BC_DIRICHLET,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  else
    call bcond_newBC (BC_DIRICHLET,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  end if

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = 1
  p_rbcReg%Iequations(1) = iequation

  ! Eventually, return the pointer
  if (present(p_rbcRegion)) p_rbcRegion => p_rbcReg

  end subroutine
      
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcond_newFeastMirrorBConRealBD (rboundaryConditions,iequation,&
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
  integer, intent(in) :: iequation

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  type(t_boundaryConditions), intent(inout), target :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make
  ! more specific modifications to this if necessary.
  type(t_bcRegion), optional, pointer :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  type(t_bcRegion), pointer :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  if (rboundaryRegion%iboundSegIdx .ne. 0) then
    call bcond_newBC (BC_FEASTMIRROR,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  else
    call bcond_newBC (BC_FEASTMIRROR,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  end if

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = 1
  p_rbcReg%Iequations(1) = iequation

  ! Eventually, return the pointer
  if (present(p_rbcRegion)) p_rbcRegion => p_rbcReg

  end subroutine
      
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcond_newPressureDropBConRealBD (rboundaryConditions,IvelEqns,&
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
  integer, dimension(:), intent(in) :: IvelEqns

  ! A boundary-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  type(t_boundaryConditions), intent(inout), target :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make more specific
  ! modifications to this.
  type(t_bcRegion), optional, pointer :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  type(t_bcRegion), pointer :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  if (rboundaryRegion%iboundSegIdx .ne. 0) then
    call bcond_newBC (BC_PRESSUREDROP,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  else
    call bcond_newBC (BC_PRESSUREDROP,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  end if

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = size(IvelEqns)
  p_rbcReg%Iequations(1:size(IvelEqns)) = IvelEqns

  ! Eventually, return the pointer
  if (present(p_rbcRegion)) p_rbcRegion => p_rbcReg

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcond_newSlipBConRealBD (rboundaryConditions,IvelEqns,&
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
  integer, dimension(:), intent(in) :: IvelEqns

  ! A boundary-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  type(t_boundaryConditions), intent(inout), target :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary region. The caller can make more specific
  ! modifications to this.
  type(t_bcRegion), optional, pointer :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  type(t_bcRegion), pointer :: p_rbcReg

  ! Add a general boundary condition region.
  ! Add a 'real' boundary condition region if rboundaryRegion%iboundSegIdx
  ! identifies that the boundary region belongs to a real boundary segment.
  ! Add a 'free' boundary condition region if the rboundaryRegion does not belong
  ! to any real segment (rboundaryRegion%iboundSegIdx=0).
  if (rboundaryRegion%iboundSegIdx .ne. 0) then
    call bcond_newBC (BC_SLIP,BC_RTYPE_REAL,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  else
    call bcond_newBC (BC_SLIP,BC_RTYPE_FREE,rboundaryConditions,&
                      p_rbcReg,rboundaryRegion)
  end if

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = size(IvelEqns)
  p_rbcReg%Iequations(1:size(IvelEqns)) = IvelEqns

  ! Eventually, return the pointer
  if (present(p_rbcRegion)) p_rbcRegion => p_rbcReg

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcond_newDirichletBConFictBD (rboundaryConditions,Iequations,&
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
  integer, dimension(:), intent(in) :: Iequations

  ! A fictitious-boundary-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  ! A copy of this is added to the rboundaryConditions structure.
  type(t_fictBoundaryRegion), intent(in) :: rboundaryRegion
!</input>

!<inputoutput>
  ! The structure where the boundary condition region is to be added.
  type(t_boundaryConditions), intent(inout), target :: rboundaryConditions
!</inputoutput>

!<output>
  ! OPTIONAL: A pointer to the added boundary condition region. The caller can
  ! make more specific modifications to this if necessary.
  type(t_bcRegion), optional, pointer :: p_rbcRegion
!</output>

!</subroutine>

  ! local variables
  type(t_bcRegion), pointer :: p_rbcReg

  ! Add a general boundary condition region for a fictitiouos boundary object.
  call bcond_newBC (BC_DIRICHLET,BC_RTYPE_FBCOBJECT,rboundaryConditions,&
                    p_rbcReg,rfictBoundaryRegion=rboundaryRegion)

  ! Modify the structure, impose additionally needed information
  ! for Dirichlet boundary - in detail, set the equation where the BC
  ! applies to.
  p_rbcReg%nequations = size(Iequations)
  p_rbcReg%Iequations(1:size(Iequations)) = Iequations

  ! Eventually, return the pointer
  if (present(p_rbcRegion)) p_rbcRegion => p_rbcReg

  end subroutine
      
end module

!##############################################################################
!# ****************************************************************************
!# <name> discretebc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements discrete boundary conditions. Discrete boundary
!# conditions are an element-dependent way to represent analytical boundary
!# conditions. After "discretising", the boundary condition can quickly
!# be implemented into a vector. Therefore, one can see discrete boundary
!# conditions also as a kind of 'precalculated' boundary conditions. This is
!# exploited e.g. in the filter approach, where a filter routine 'applies' a
!# discrete BC to a vector.
!#
!# The following types of discrete boundary conditions are available:
!#
!# - Dirichlet boundary conditions are typically represented as a
!#   list of DOF`s where a value must be prescribed, plus the value that
!#   must be described in that special DOF.
!#
!# - Pressure drop conditions consist of a list of DOF`s and a modifier
!#   how to modify the actual DOF`s.
!#
!# - Slip boundary conditions consists like Dirichlet boundary conditions
!#   of a list of DOF`s, wher zero must be implemented into the defect
!#   vector. Slip boundary conditions are nonlinear boundary conditions:
!#   In case of a linear system, the corresponding DOF`s are treated
!#   as Dirichlet. The actual BC is implemented during a nonlinear
!#   loop by modifying the appearing defect vector.
!#
!# </purpose>
!##############################################################################

module discretebc

!$ use omp_lib
  use fsystem
  use storage
  use dofmapping

  implicit none

  private

!<constants>

!<constantblock description="The type identifier for discrete (linear) boundary conditions">

  ! undefined discrete BC`s
  integer, parameter, public :: DISCBC_TPUNDEFINED    = 0

  ! Discrete Dirichlet boundary conditions
  integer, parameter, public :: DISCBC_TPDIRICHLET    = 1

  ! Discrete FEAST mirror boundary conditions
  integer, parameter, public :: DISCBC_TPFEASTMIRROR  = 2

!</constantblock>

!<constantblock description="The type identifier for discrete nonlinear boundary conditions">

  ! Discrete pressure drop boundary conditions
  integer, parameter, public :: DISCBC_TPPRESSUREDROP = 100

  ! Discrete slip boundary conditions
  integer, parameter, public :: DISCBC_TPSLIP         = 101

!</constantblock>

!<constantblock description="Type identifiers for the callback routine during discretisation of BC`s">

  ! Calculate the function value in a point on the boundary
  integer, parameter, public :: DISCBC_NEEDFUNC         = 0

  ! Calculate the function value in a point on an edge on the boundary
  integer, parameter, public :: DISCBC_NEEDFUNCMID      = 1

  ! Calculate the x- and y-derivative in a point on the boundary
  integer, parameter, public :: DISCBC_NEEDDERIV        = 2

  ! Calculate the integral mean value over an edge on the boundary
  integer, parameter, public :: DISCBC_NEEDINTMEAN      = 3

  ! Calculate the normal stress in a point. For flow-like problems,
  ! this corresponds to a prescribed pressure in pressure-drop problems.
  integer, parameter, public :: DISCBC_NEEDNORMALSTRESS = 4

!</constantblock>

!<constantblock>

  ! Default blocksize for allocating new structures in the p_RdiscBCList
  ! list - if the list is full.
  integer, parameter, public :: DISCBC_LISTBLOCKSIZE = 10

!</constantblock>
!</constants>

!<types>

!<typeblock>

  ! This structure describes the typical way, Dirichlet boundary conditions
  ! can be discretised. This is done by two arrays: one array is a list of all
  ! DOF`s that refer do Dirichlet nodes. The second array refers to the value
  ! that must be imposed in this DOF.
  ! The variable icomponent describes the number of the component/equation
  ! in the PDE that must be treated that way.

  type t_discreteBCDirichlet

    ! The component of the equation, this discrete BC is specified for
    ! (e.g. 1=X-velocity, 2=Y-velocity or similar)
    integer :: icomponent = 0
    
    ! Number of equations in the vector. Used to assure that the BC structure
    ! is not applied to a wrong vector.
    integer :: NEQ = 0

    ! Number of Dirichlet nodes; may be different from the length of the array!
    integer :: nDOF = 0

    ! Handle to array with all DOF`s that refer to Dirichlet nodes
    !   array [1..*] of integer
    integer :: h_IdirichletDOFs   = ST_NOHANDLE

    ! Handle to array with the Dirichlet value that should be imposed in these nodes
    !   array [1..*] of double
    integer :: h_DdirichletValues = ST_NOHANDLE

  end type

  public :: t_discreteBCDirichlet

!</typeblock>

!<typeblock>

  ! This structure contains the information to implement slip boundary
  ! conditions in case of a Navier-Stokes solver.
  ! Slip boundary conditions are handled the same way as Dirichlet-
  ! zero boundary conditions, but have a different treatment due the nonlinear
  ! loop. We have one array is a list of all DOF`s that refer do slip boundary
  ! nodes.
  ! The variable ncomponent describes the number ov velocity components/
  ! equations in the PDE. Icomponents(1..ncomponents) on the other hand
  ! is a list of numbers of the velocity equations in the PDE.
  ! Normally there is Icomponents(1)=1=X-velocity, Icomponents(1)=2=Y-velocity,
  ! probably Icomponents(3)=3=Z-velocity,

  type t_discreteBCSlip

    ! Number of velocity components that take part on the slip
    ! boundary conditions.
    integer :: ncomponents = 0

    ! List of all velocity components in the PDE that take part on
    ! the slip boundary conditions.
    ! (e.g. 1=X-velocity, 2=Y-velocity or similar)
    integer, dimension(:), pointer :: Icomponents => null()

    ! Number of equations in the vector. Used to assure that the BC structure
    ! is not applied to a wrong vector.
    integer :: NEQ = 0

    ! Number of Dirichlet nodes; may be different from the length of the array!
    integer :: nDOF = 0

    ! Handle to array with all DOF`s that refer to Dirichlet nodes
    !   array [1..*] of integer
    integer :: h_IslipDOFs = ST_NOHANDLE

    ! Handle to an array that contains the normal vectors of the
    ! boundary edges.
    !   DnormalVectors = array [1..ncomponents,1..*] of double with
    !   DnormalVectors(1)=X-component,
    !   DnormalVectors(2)=Y-component,...
    ! of the normal vector.
    integer :: h_DnormalVectors = ST_NOHANDLE

  end type

  public :: t_discreteBCSlip

!</typeblock>

!<typeblock>

  ! This structure describes the way, pressure drop boundary conditions
  ! can be discretised. This is done by two arrays: one array is a list of all
  ! (velocity) DOF`s. The second array specifies for each of these DOF`s
  ! the (non-zero) value, the pressure is multiplied with.
  ! The variable icomponent describes the number of the component/equation
  ! in the PDE that must be treated that way.

  type t_discreteBCpressureDrop

    ! Number of "velocity" components that must be modified when implementing
    ! pressure drop boundary conditions.
    integer :: ncomponents = 0

    ! The components of the solutions/RHS vectors that should be modified
    ! by pressure drop boundary conditions.
    ! Each of the 1..ncomponents entries in the vector specifies a component
    ! in the solution vector that is modified (e.g. 1=X-velocity, 2=Y-velocity
    ! or similar)
    integer, dimension(:), pointer :: Icomponents => null()

    ! Number of equations in the vector. Used to assure that the BC structure
    ! is not applied to a wrong vector.
    integer :: NEQ = 0

    ! Number of DOF`s in the arrays below; may be different from the length of
    ! the array!
    integer :: nDOF = 0

    ! Handle to array with all velocity DOF`s on the boundary that must be
    ! modified.
    !   array [1..*] of integer
    integer :: h_IpressureDropDOFs = ST_NOHANDLE

    ! Handle to array with additive content that must be added to the DOF`s
    ! in the h_IpressureDropDOFs array.
    !   array [1..NDIM2D,1..*] of double
    integer :: h_Dmodifier = ST_NOHANDLE

  end type

  public :: t_discreteBCpressureDrop

!</typeblock>

!<typeblock>

  ! This structure describes the way, FEAST mirror boundary conditions
  ! can be discretised.
  ! The variable icomponent describes the number of the component/equation
  ! in the PDE that must be treated that way.
  ! h_ImirrorBCs specifies a bitfield that defines for every DOF if it is
  ! a FEAST mirror BC DOF or not.

  type t_discreteBCFeastMirror

    ! The component of the equation, this discrete BC is specified for
    ! (e.g. 1=X-velocity, 2=Y-velocity or similar)
    integer  :: icomponent        = 0

    ! Coarsening level. Used when adding additional contributions to the
    ! matrix/vector. Usually = 0. Must be increased for every level coarser
    ! than the maximum one in a mesh hierarchy.
    real(DP) :: icoarseningLevel  = 0

    ! =0: Modify matrix and defect vectors.
    ! =1: Modify matrix, treat defect vectors as Dirichlet.
    integer :: isubtype           = 0

    ! Handle to a list of all DOF`s in the FEAST mirror boundary region.
    ! The list is sorted for increasing DOF numbers.
    integer :: h_ImirrorDOFs   = ST_NOHANDLE

    ! Handle to a list of all DOF`s in the FEAST mirror boundary region
    ! plus the start- and endpoint
    ! The list is sorted for increasing DOF numbers.
    integer :: h_ImirrorDOFsClosed   = ST_NOHANDLE

  end type

  public :: t_discreteBCFeastMirror

!</typeblock>


!<typeblock>

  ! This describes the basic structure for discrete boundary conditions.
  ! A type identifier decides on which boundary conditions this structure
  ! describes. Depending on the type, one of the information blocks
  ! is filled with data about the discrete BC`s.

  type t_discreteBCEntry

    ! The type identifier. Identifies the type of discrete BC`s, this
    ! structure describes.
    integer                             :: itype = DISCBC_TPUNDEFINED

    ! Structure for discrete Dirichlet BC`s.
    ! Only valid if itype=DISCBC_TPDIRICHLET.
    type(t_discreteBCDirichlet)         :: rdirichletBCs

    ! Structure for discrete pressure drop BC`s.
    ! Only valid if itype=DISCBC_TPPRESSUREDROP.
    type(t_discreteBCpressureDrop)      :: rpressuredropBCs

    ! Structure for discrete Slip BC`s.
    ! Only valid if itype=DISCBC_TPSLIP.
    type(t_discreteBCSlip)              :: rslipBCs

    ! Structure for discrete FEAST mirror BC`s.
    ! Only valid if itype=DISCBC_TPFEASTMIRROR.
    type(t_discreteBCFeastMirror)         :: rfeastMirrorBCs

  end type

  public :: t_discreteBCEntry

!</typeblock>

!<typeblock>

  ! The main structure for discrete boundary conditions.
  ! This is just an array of t_discreteBCEntry structures, each describing
  ! a single discrete boundary condition (so to speak, a segment on the
  ! boundary discretised in a special way).
  type t_discreteBC

    ! Total number of allocated t_discreteBCEntry structures in p_RdiscBCList.
    integer :: inumEntriesAlloc = 0

    ! Total number of used t_discreteBCEntry structures in p_RdiscBCList.
    integer :: inumEntriesUsed = 0

    ! An array of t_discreteBCEntry structures. Each structure describes
    ! one discrete boundary condition - so one part of the boundary discretised
    ! in a special, discretisation-dependent way.
    type(t_discreteBCEntry), dimension(:), pointer :: p_RdiscBCList => null()

  end type

  public :: t_discreteBC

!</typeblock>

!</types>

end module

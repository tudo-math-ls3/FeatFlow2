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
!# conditionsalso as a kind of 'precalculated' boundary conditions. This is 
!# exploited e.g. in the filter approach, where a filter routine 'applies' a 
!# discrete BC to a vector. 
!#
!# The following types of discrete boundary conditions are available:
!#
!# - Dirichlet boundary conditions are typically represented as a
!#   list of DOF's where a value must be prescribed, plus the value that
!#   must be described in that special DOF.
!#
!# - Pressure drop conditions consist of a list of DOF's and a modifier
!#   how to modify the actual DOF's.
!#   
!# </purpose>
!##############################################################################

MODULE discretebc

  USE fsystem
  USE boundarycondition
  USE dofmapping
  
  IMPLICIT NONE

!<constants>

!<constantblock description="General constants concerning filters">

  ! A standard length for arrays holding a set of discretised BC's
  INTEGER, PARAMETER :: DISCBC_MAXDISCBC         =  32

!</constantblock>

!<constantblock description="The type identifier for discrete boundary conditions">

  ! undefined discrete BC's
  INTEGER, PARAMETER :: DISCBC_TPUNDEFINED    = 0

  ! Discrete Dirichlet boundary conditions
  INTEGER, PARAMETER :: DISCBC_TPDIRICHLET    = 1
  
  ! Discrete pressure drop boundary conditions
  INTEGER, PARAMETER :: DISCBC_TPPRESSUREDROP = 2

!</constantblock>

!<constantblock description="Type identifiers for the callback routine during discretisation of BC's">
  
  ! Calculate the function value in a point on the boundary
  INTEGER, PARAMETER :: DISCBC_NEEDFUNC         = 0

  ! Calculate the function value in a point on an edge on the boundary
  INTEGER, PARAMETER :: DISCBC_NEEDFUNCMID      = 1
  
  ! Calculate the x- and y-derivative in a point on the boundary
  INTEGER, PARAMETER :: DISCBC_NEEDDERIV        = 2

  ! Calculate the integral mean value over an edge on the boundary
  INTEGER, PARAMETER :: DISCBC_NEEDINTMEAN      = 3

  ! Calculate the normal stress in a point. For flow-like problems,
  ! this corresponds to a prescribed pressure in pressure-drop problems.
  INTEGER, PARAMETER :: DISCBC_NEEDNORMALSTRESS = 4

!</constantblock>
  
!</constants>

!<types>
  
!<typeblock>
  
  ! This structure describes the typical way, Dirichlet boundary conditions
  ! can be discretised. This is done by two arrays: one array is a list of all
  ! DOF's that refer do Dirichlet nodes. The second array refers to the value
  ! that must be imposed in this DOF.
  ! The variable icomponent describes the number of the component/equation
  ! in the PDE that must be treated that way.
  
  TYPE t_discreteBCDirichlet
    
    ! The component of the equation, this discrete BC is specified for
    ! (e.g. 1=X-velocity, 2=Y-velocity or similar)
    INTEGER                            :: icomponent        = 0
    
    ! Number of Dirichlet nodes; may be different from the length of the array!
    INTEGER(PREC_DOFIDX)               :: nDOF              = 0
    
    ! Handle to array with all DOF's that refer to Dirichlet nodes
    !   array [1..*] of integer
    INTEGER :: h_IdirichletDOFs   = ST_NOHANDLE
    
    ! Handle to array with the Dirichlet value that should be imposed in these nodes
    !   array [1..*] of double
    INTEGER :: h_DdirichletValues = ST_NOHANDLE
    
  END TYPE
  
!</typeblock>
  
!<typeblock>
  
  ! This structure describes the way, pressure drop boundary conditions
  ! can be discretised. This is done by two arrays: one array is a list of all
  ! (velocity) DOF's. The second array specifies for each of these DOF's
  ! the (non-zero) value, the pressure is multiplied with.
  ! The variable icomponent describes the number of the component/equation
  ! in the PDE that must be treated that way.
  
  TYPE t_discreteBCpressureDrop
    
    ! Number of "velocity" components that must be modified when implementing
    ! pressure drop boundary conditions.
    INTEGER :: ncomponents = 0
    
    ! The components of the solutions/RHS vectors that should be modified
    ! by pressure drop boundary conditions.
    ! Each of the 1..ncomponents entries in the vector specifies a component
    ! in the solution vector that is modified (e.g. 1=X-velocity, 2=Y-velocity 
    ! or similar)
    INTEGER, DIMENSION(SPDISC_MAXEQUATIONS) :: Icomponents        = 0
    
    ! Number of DOF's in the arrays below; may be different from the length of 
    ! the array!
    INTEGER(PREC_DOFIDX)               :: nDOF              = 0
    
    ! Handle to array with all velocity DOF's on the boundary that must be 
    ! modified.
    !   array [1..*] of integer
    INTEGER :: h_IpressureDropDOFs   = ST_NOHANDLE
    
    ! Handle to array with additive content that must be added to the DOF's
    ! in the h_IpressureDropDOFs array.
    !   array [1..NDIM2D,1..*] of double
    INTEGER :: h_Dmodifier = ST_NOHANDLE
    
  END TYPE
  
!</typeblock>
  
!<typeblock>
  
  ! This describes the basic structure for discrete boundary conditions.
  ! A type identifier decides on which boundary conditions this structure
  ! describes. Depending on the type, one of the information blocks
  ! is filled with data about the discrete BC's.
  
  TYPE t_discreteBCEntry
    
    ! The type identifier. Identifies the type of discrete BC's, this
    ! structure describes.
    INTEGER                             :: itype = DISCBC_TPUNDEFINED
    
    ! Pointer to the boundary-condition structure defining the analytic
    ! boundary conditions - or NULL(), if these BC does not arise
    ! from analytic boundary conditions.
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions => NULL()
    
    ! Structure for discrete Dirichlet BC's.
    ! Only valid if itype=DISCBC_TPDIRICHLET.
    TYPE(t_discreteBCDirichlet)         :: rdirichletBCs
    
    ! Structure for discrete pressure drop BC's.
    ! Only valid if itype=DISCBC_TPPRESSUREDROP.
    TYPE(t_discreteBCpressureDrop)      :: rpressuredropBCs
    
  END TYPE
  
!</typeblock>

!<typeblock>

  ! The main structure for discrete boundary conditions.
  ! This is just an array of t_discreteBCEntry structures, each describing
  ! a single discrete boundary condition (so to speak, a segment on the
  ! boundary discretised in a special way).
  TYPE t_discreteBC
  
    ! An array of t_discreteBCEntry structures. Each structure describes
    ! one discrete boundary condition - so one part of the boundary discretised
    ! in a special, discretisation-dependent way.
    TYPE(t_discreteBCEntry), DIMENSION(:), POINTER :: p_RdiscBCList => NULL()
  
  END TYPE

!</typeblock>

!</types>
  
END MODULE

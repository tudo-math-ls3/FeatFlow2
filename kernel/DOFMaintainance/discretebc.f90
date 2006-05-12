!##############################################################################
!# ****************************************************************************
!# <name> DiscreteBC </name>
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
!# - Dirichlet boundary conditions are typically represented as a
!#   list of DOF's where a value must be prescribed, plus the value that
!#   must be described in that special DOF.
!#   
!# </purpose>
!##############################################################################

MODULE discretebc

  USE fsystem
  USE scalarbc
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

!</constantblock>

!<constantblock description="Type identifiers for the callback routine during discretisation of BC's">
  
  ! Calculate the function value in a point on the boundary
  INTEGER, PARAMETER :: DISCBC_NEEDFUNC       = 0
  
  ! Calculate the x- and y-derivative in a point on the boundary
  INTEGER, PARAMETER :: DISCBC_NEEDDERIV      = 1

  ! Calculate the integral mean value over an edge on the boundary
  INTEGER, PARAMETER :: DISCBC_NEEDINTMEAN    = 2

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
    
    ! Number of Dirichlet nodes; may be different from the length of the array!
    INTEGER(PREC_DOFIDX)               :: nDOF              = 0
    
    ! Handle to array with all DOF's that refer to Dirichlet nodes
    !   array [1..*] of integer
    INTEGER :: h_IdirichletDOFs   = ST_NOHANDLE
    
    ! Handle to array with the Dirichlet value that should be imposed in these nodes
    !   array [1..*] of double
    INTEGER :: h_IdirichletValues = ST_NOHANDLE
    
  END TYPE
  
  !</typeblock>
  
  !<typeblock>
  
  ! This describes the basic structure for discrete boundary conditions.
  ! A type identifier decides on which boundary conditions this structure
  ! describes. Depending on the type, one of the information blocks
  ! is filled with data about the discrete BC's.
  
  TYPE t_discreteBC
    
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
    
  END TYPE
  
  !</typeblock>

!</types>
  
END MODULE

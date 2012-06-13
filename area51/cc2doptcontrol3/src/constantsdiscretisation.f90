!##############################################################################
!# ****************************************************************************
!# <name> constantsdiscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains program wide constants.
!# </purpose>
!##############################################################################

module constantsdiscretisation

  use fsystem
  
  implicit none
  
  private
  
!<constants>

!<constantblock description="General constant specifying the primal and/or dual space.">

  ! undefined
  integer, parameter, public :: CCSPACE_UNDEFINED = 0

  ! primal space
  integer, parameter, public :: CCSPACE_PRIMAL = 1

  ! dual space
  integer, parameter, public :: CCSPACE_DUAL = 2
  
  ! Control space
  integer, parameter, public :: CCSPACE_CONTROL = 4

!</constantblock>

!<constantblock description = "Supported equations">
  
  ! Navier-Stokes equations, 2D
  integer, parameter, public :: CCEQ_NAVIERSTOKES2D = 0

  ! Stokes equations, 2D
  integer, parameter, public :: CCEQ_STOKES2D = 1

  ! Heat equation, 2D
  integer, parameter, public :: CCEQ_HEAT2D = 2
  
!</constantblock>

!<constantblock description="Identifiers for the IUPWIND stabilisation parameter.">

  ! Streamline diffusion; configured by dupsam
  integer, parameter, public :: CCSTAB_STREAMLINEDIFF    = 0

  ! 1st-order upwind; configured by dupsam
  integer, parameter, public :: CCSTAB_UPWIND            = 1
  
  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'
  integer, parameter, public :: CCSTAB_EDGEORIENTED      = 2

  ! Streamline diffusion; configured by dupsam, new implementation
  integer, parameter, public :: CCSTAB_STREAMLINEDIFF2   = 3

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation
  integer, parameter, public :: CCSTAB_EDGEORIENTED2     = 4

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation.
  ! Precomputed matrix.
  integer, parameter, public :: CCSTAB_EDGEORIENTED3     = 5
!</constantblock>

!</constants>

end module

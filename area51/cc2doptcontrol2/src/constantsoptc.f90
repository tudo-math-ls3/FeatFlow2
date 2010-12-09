!##############################################################################
!# ****************************************************************************
!# <name> constantsoptc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains program wide constants.
!# </purpose>
!##############################################################################

module constantsoptc

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
  
  ! primal and dual space
  integer, parameter, public :: CCSPACE_PRIMALDUAL = 3

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

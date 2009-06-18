!##############################################################################
!# ****************************************************************************
!# <name> transport_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global vairables
!# which are required to solve a conservation law for a scalar variable
!#
!# </purpose>
!##############################################################################

module transport_basic

  use fparser
  use fsystem
  use statistics

  implicit none

  private 

!<constants>

!<constantblock description="Global type of mass matrix">

  ! no mass matrix, i.e., no time derivative
  integer, parameter, public :: MASS_ZERO       = 0
  
  ! consistent mass matrix
  integer, parameter, public :: MASS_CONSISTENT = 1

  ! lumped mass matrix
  integer, parameter, public :: MASS_LUMPED     = 2

!</constantblock>


!<constantblock description="Global type of flow velocities">

  ! zero velocity profile v=0
  integer, parameter, public :: VELOCITY_ZERO              = 0

  ! constant velocity profile v=v(x)
  integer, parameter, public :: VELOCITY_CONSTANT          = 1

  ! linear time-dependent velocity profile v=v(x,t)
  integer, parameter, public :: VELOCITY_TIMEDEP           = 2

  ! nonlinear Burgers' equation in space-time
  integer, parameter, public :: VELOCITY_BURGERS_SPACETIME = 3

  ! nonlinear Buckley-Leverett equation in space-time
  integer, parameter, public :: VELOCITY_BUCKLEV_SPACETIME = 4

  ! nonlinear Burgers' equation in 1D
  integer, parameter, public :: VELOCITY_BURGERS1D         = 5

  ! nonlinear Burgers' equation in 2D
  integer, parameter, public :: VELOCITY_BURGERS2D         = 6

  ! nonlinear Burgers' equation in 3D
  integer, parameter, public :: VELOCITY_BURGERS3D         = 7

  ! nonlinear Buckley-Leverett equation in 1D
  integer, parameter, public :: VELOCITY_BUCKLEV1D         = 8

!</constantblock>


!<constantblock description="Global type of diffusion">

  ! zero diffusion
  integer, parameter, public :: DIFFUSION_ZERO        = 0

  ! isotropic diffusion
  integer, parameter, public :: DIFFUSION_ISOTROPIC   = 1

  ! anisotropic diffusion
  integer, parameter, public :: DIFFUSION_ANISOTROPIC = 2

  ! variable diffusion
  integer, parameter, public :: DIFFUSION_VARIABLE    = 3

!</constantblock>


!<constantblock description="Global type of reaction">

  ! zero reaction
  integer, parameter, public :: REACTION_ZERO = 0

!</constantblock>


!<constantblock description="Global type of right-hand side">

  ! zero right-hand side
  integer, parameter, public :: RHS_ZERO     = 0

  ! analytical right-hand side
  integer, parameter, public :: RHS_ANALYTIC = 1

!</constantblock>

!<constantblock description="Global type of initial solution">

  ! zero initial solution
  integer, parameter, public :: SOLUTION_ZERO     = 0

  ! analytical initial solution
  integer, parameter, public :: SOLUTION_ANALYTIC = 1

  ! graymap profile for initial solution
  integer, parameter, public :: SOLUTION_GRAYMAP  = 2

!</constantblock>


!<constantblock description="Global type of target functional">

  ! zero target functional
  integer, parameter, public :: TFUNC_ZERO     = 0

  ! volume integral target functional
  integer, parameter, public :: TFUNC_VOLINTG  = 1

  ! surface integral target functional
  integer, parameter, public :: TFUNC_SURFINTG = 2

  ! mixed volume + surface integral target functional
  integer, parameter, public :: TFUNC_MIXINTG  = 3

  ! analytical value for target functional
  integer, parameter, public :: TFUNC_ANALYTIC = 4

!</constantblock>


!<constantblock description="Global type of recovery-based error estimation">

  ! L2-projection
  integer, parameter, public :: ERREST_L2PROJECTION  = 1

  ! Superconvergent patch recovery (vertex-based)
  integer, parameter, public :: ERREST_SPR_VERTEX    = 2

  ! Superconvergent patch recovery (element-based)
  integer, parameter, public :: ERREST_SPR_ELEMENT   = 3

  ! Superconvergent patch recovery (face-based)
  integer, parameter, public :: ERREST_SPR_FACE      = 4
 
  ! Limited averaging gradient recovery
  integer, parameter, public :: ERREST_LIMAVR        = 5

  ! Second-difference indicator (by Loehner)
  integer, parameter, public :: ERREST_SECONDDIFF    = 6

!</constantblock>


!<constantblock description="Global constants for error redistribution">

  ! Use error 'as is'
  integer, parameter, public :: ERREST_ASIS          = 0

  ! Equidistribution of error
  integer, parameter, public :: ERREST_EQUIDIST      = 1

  ! Logarithmic equidistribution of error
  integer, parameter, public :: ERREST_LOGEQUIDIST   = 2

  ! Automatic treshold based on RMS
  integer, parameter, public :: ERREST_AUTORMS       = 3

!</constantblock>


!<constantblock description="Global types of perturbation parameters">

  ! Perturbation parameter is chosen as in the NITSOL package
  integer, parameter, public :: PERTURB_NITSOL  = 1

  ! Perturbation parameter is chosen as SQRT(machine precision)
  integer, parameter, public :: PERTURB_SQRTEPS = 2

!</constantblock>
  
!</constants>

end module transport_basic

!##############################################################################
!# ****************************************************************************
!# <name> timestepaux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures for time-stepping algorithms.
!#
!# </purpose>
!##############################################################################

module timestepaux

  use fsystem
  use linearsystemblock

!<constants>

!<constantblock description="Global time-stepping types">

  ! Two-level theta-scheme
  integer, parameter, public :: TSTEP_THETA_SCHEME = 1

  ! Two-level Runge-Kutta scheme
  integer, parameter, public :: TSTEP_RK_SCHEME    = 2
!</constantblock>


!<constantblock description="Adaptive time-stepping types">

  ! No adaptive time-stepping
  integer, parameter, public :: TSTEP_NOADAPT   = 0

  ! Adaptive time-stepping by PID controller
  integer, parameter, public :: TSTEP_PIDADAPT  = 1

  ! Adaptive time-stepping by truncation error analysis
  integer, parameter, public :: TSTEP_AUTOADAPT = 2

  ! Adaptive time-stepping by switched evolution relaxation
  integer, parameter, public :: TSTEP_SERADAPT  = 3
!</constantblock>


!<constantblock description="Global flags for information output">

  ! Silent run, no output at all
  integer, parameter, public :: TSTEP_IOLEVEL_SILENT  = 0

  ! Output only errors
  integer, parameter, public :: TSTEP_IOLEVEL_ERROR   = 1

  ! Output only errors and warnings
  integer, parameter, public :: TSTEP_IOLEVEL_WARNING = 2

  ! Output errors, warnings and information
  integer, parameter, public :: TSTEP_IOLEVEL_INFO    = 3

  ! Output errors, warnings and verbose information
  integer, parameter, public :: TSTEP_IOLEVEL_VERBOSE = 4

!</constantblock>

!</constants>

  ! *****************************************************************************

!<types>

!<typeblock>

  ! This data structure contains settings/parameters for time-stepping. Before
  ! calling the time-stepping algorithm the basic settings need to be initialised.

  type t_timestep

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Name of time-stepping object
    character(LEN=SYS_STRLEN) :: sName = ''

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Identifies the type of the time-stepping algorithm
    integer :: ctimestepType = 0

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! This determines the output level of the time-stepping algorithm
    integer :: ioutputLevel = 0

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Norm in which relative changes of solution should be measured
    integer :: isolNorm = 0

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Initial simulation time
    real(DP) :: dinitialTime = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Final simulation time
    real(DP) :: dfinalTime = 0.0_DP

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Current time instant
    real(DP) :: dTime = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Lower bound for the admissible time step
    real(DP) :: dminStep = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Upper bound for the admissible time step
    real(DP) :: dmaxStep = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Initial value for time step
    real(DP) :: dinitialStep = 0.0_DP

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Current value for time step
    real(DP) :: dStep = 0.0_DP

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Value for previous time step
    real(DP) :: dStep1 = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Factor by which the current time step should be scaled
    ! if the current time step failed; must be < 1.0
    real(DP) :: dstepReductionFactor = 1.0_DP

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Total number of performed time steps
    integer :: nSteps = 0

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Total number of rejected time steps
    integer :: nrejectedSteps = 0

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Implicitness parameter for two-level theta-scheme
    real(DP) :: theta = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Number of multi-steps for Runge-Kutta method
    integer :: multisteps = 0

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Algorithm for adaptive time step control
    integer :: iadaptTimestep = 0

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Initial time for performing adaptive time stepping
    real(DP) :: dadaptTime = 0.0_DP

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Relative changes of the solution
    real(DP) :: drelChange = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Relative stopping criterion. Stop time-stepping algorithm if
    ! !! U^{N+1} - U^{N} !! <= EPSSTAG * dStep
    ! =0: ignore, do not stop time-stepping until final time is reached
    real(DP) :: depsSteady = 0.0_DP

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Weights for multistage methods
    real(DP), dimension(:), pointer :: DmultistepWeights => null()

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Temporal vectors for old solutions
    type(t_vectorBlock), dimension(:), pointer :: RtempVectors => null()

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Substructure for PID controller
    type(t_pidController), pointer :: p_rpidController => null()

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Substructure for automatic controller
    type(t_autoController), pointer :: p_rautoController => null()

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Substructure for switched evolution relaxation controller
    type(t_serController), pointer :: p_rserController => null()

  end type t_timestep

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all settings/parameters for
  ! proportional-integral-derivative (PID) controlling.

  type t_pidController

    ! INPUT PARAMETER FOR PID CONTROLLER:
    ! Exponent for proportional term.
    real(DP) :: dProportionalExponent = 0.075_DP

    ! INPUT PARAMETER FOR PID CONTROLLER:
    ! Exponent for integral term.
    real(DP) :: dIntegralExponent =  0.175_DP

    ! INPUT PARAMETER FOR PID CONTROLLER:
    ! Exponent for derivative term.
    real(DP) :: dDerivativeExponent = 0.01_DP

    ! INPUT PARAMETER FOR PID CONTROLLER:
    ! Factor by which the time step is allowed to increase.
    real(DP) :: dIncreaseFactor = 2.0_DP

    ! INPUT PARAMETER FOR PID CONTROLLER:
    ! Factor by which the time step is allowed to decrease.
    real(DP) :: dDecreaseFactor = 0.5_DP

    ! INPUT PARAMETER FOR PID CONTROLLER:
    ! Tolerance of relative changes.
    real(DP) :: depsRel = 0.0_DP

    ! INPUT PARAMETER FOR PID CONTROLLER:
    ! Maximum tolerance for relative changes
    real(DP) :: dmaxRel = 0.0_DP

    ! OUTPUT PARAMETER FOR PID CONTROLLER:
    ! Measure of the control variable from one time step before.
    real(DP) :: dcontrolValue1 = 0.0_DP

    ! OUTPUT PARAMETER FOR PID CONTROLLER:
    ! Measure of the solution change from two time steps before.
    real(DP) :: dcontrolValue2 = 0.0_DP

  end type t_pidController

!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all settings/parameters
  ! for automatic time step controlling.

  type t_autoController

    ! INPUT PARAMETER FOR AUTOMATIC CONTROLLER:
    ! Factor by which the time step is allowed to decrease.
    real(DP) :: dDecreaseFactor = 0.0_DP

    ! INPUT PARAMETER FOR AUTOMATIC CONTROLLER:
    ! Tolerance of relative changes.
    real(DP) :: depsRel = 0.0_DP

    ! INTERNAL PARAMETER FOR AUTOMATIC CONTROLLER
    ! Temporal vectors for old solutions
    type(t_vectorBlock), dimension(:), pointer :: RtempVectors => null()

  end type t_autoController

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all settings/parameters
  ! for switched evolution relaxation time step controlling.

  type t_serController

    ! INPUT PARAMETER FOR SER CONTROLLER:
    ! Factor by which the time step is allowed to increase.
    real(DP) :: dIncreaseFactor = 0.0_DP

    ! INPUT PARAMETER FOR SER CONTROLLER:
    ! Factor by which the time step is allowed to decrease.
    real(DP) :: dDecreaseFactor = 0.0_DP

    ! OUTPUT PARAMETER FOR SER CONTROLLER
    ! Norm of the stationary defect vector
    real(DP) :: dsteadyDefect = 0.0_DP

    ! OUTPUT PARAMETER FOR SER CONTROLLER
    ! Norm of the stationary defect vector from previous time step
    real(DP) :: dsteadyDefect1 = 0.0_DP

  end type t_serController

!</typeblock>

!</type>

end module timestepaux

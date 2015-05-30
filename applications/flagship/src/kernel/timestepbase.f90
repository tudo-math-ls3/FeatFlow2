!##############################################################################
!# ****************************************************************************
!# <name> timestepbase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic data structures for time-stepping algorithms.
!#
!# The following routines are available:
!#
!# 1.) tstep_igetID
!#     -> Converts a string with a time stepping identifier to the
!#        corresponding time stepping id.
!#
!# 2.) tstep_getType
!#     -> Returns the type of a time stepping scheme.
!#
!# 3.) tstep_getFamily
!#     -> Returns the family of a time stepping scheme.
!#
!# 4.) tstep_getTypeName
!#     -> Returns the name of a time stepping scheme as a string.
!#
!# 5.) tstep_getFamilyName
!#     -> Returns the name of a time stepping family as a string.
!#
!# 6.) tstep_genButcherTableau
!#     -> Generates the Butcher tableau for a given Runge-Kutta scheme
!#
!# 7.) tstep_genAlphaBetaTableau
!#     -> Generates the alpha-beta tableau for a given Runge-Kutta scheme
!#
!#  FAQ - Some explainations
!# -------------------------
!#
!#
!# 1.) What is the difference between 'family' and 'type' of the time
!#     stepping scheme and how is it realised?
!#
!#     The 'family' of the time stepping scheme refers to the whole
!#     class of schemes, e.g., theta-type schemes or Runge-Kutta
!#     schemes. The 'type' is a concrete time stepping scheme. For
!#     example, TSTEP_SSPERK_2_2 represents the explicit strongly
!#     positivity preserving 2-stage Runge-Kutta scheme of order 2.
!#
!#     The timestep constant follows a special pattern oriented on a
!#     bitfield to encode as much information into an integer as
!#     possible. Every timestep identifier consists of a 32 bit
!#     integer, whcih is coded as follows:
!#
!#
!# <verb>
!# Bit | 31 ... 16 | 15           | 14 .......... 0 |
!# -------------------------------------------------------------------------
!#     |    ****   | =1: implicit | stepping scheme |
!#
!#   Bits 0..14  specifies the time stepping scheme/family
!#   Bit  15     specifies if the scheme is explicit (=0) or implicit (=1)
!#   Bits 16..31 specifies details of the time stepping scheme
!3
!#               THETA time stepping schemes:
!#               Bit 16 (=1): Fractional-step theta scheme
!#               Bits 17..31: unused
!#
!#               Runge-Kutta time stepping schemes:
!#               Bit 16 (=1): Strongly positivity preserving Runge-Kutta scheme
!#               Bit 17 (=1): Semi-explicit Runge-Kutta time stepping family of scheme
!#               Bit 18 (=1): Semi-implicit Runge-Kutta time stepping family of scheme
!#               Bit 19 (=1): Diagonally-implicit Runge-Kutta time stepping family of scheme
!#               Bit 20 (=1): Singly-diagonally implicit Runge-Kutta time stepping family of scheme
!#               Bits 21..28: unused
!#               Bit 29 (=1): use embedded Butcher tableau
!#               Bit 30 (=1): use Butcher tableau for representing the coefficients
!#                      (=0): use alpha-beta form for representing the coefficients
!#               
!# </verb>
!# </purpose>
!##############################################################################

module timestepbase

#include "flagship.h"

!$ use omp_lib
  use fsystem
  use genoutput
  use linearsystemblock

  implicit none

  private
  public :: t_timestep
  public :: t_thetaScheme
  public :: t_RungeKuttaScheme
  public :: t_pidController
  public :: t_autoController
  public :: t_serController
  public :: tstep_igetID
  public :: tstep_getType
  public :: tstep_getFamily
  public :: tstep_getTypeName
  public :: tstep_getFamilyName
  public :: tstep_genButcherTableau
  public :: tstep_genAlphaBetaTableau
  
!<constants>

!<constantblock description="Identifiers for time stepping families.">

  ! Unspecified time-stepping scheme
  integer(I32), parameter, public :: TSTEP_UNDEFINED             = -1_I32
  
  ! Two-step theta-time stepping family of schemes
  integer(I32), parameter, public :: TSTEP_THETA                 = 2_I32**4
 
  ! Runge-Kutta time stepping family of schemes
  integer(I32), parameter, public :: TSTEP_RUNGE_KUTTA           = 2_I32**5
  
  ! Bitmask specifying an implicit time stepping scheme
  integer(I32), parameter, public :: TSTEP_IMPLICIT              = 2_I32**15

  ! Bitmask for time stepping type (Bits 0..3)
  integer(I32), parameter, public :: TSTEP_TYPE                  = 15_I32
  
  ! Bitmask for time stepping family (Bits 4..14)
  integer(I32), parameter, public :: TSTEP_FAMILY                = 32752_I32
!</constantblock>
  
!<constantblock description="Identifiers for theta family of time stepping schemes">

  ! Fractional-step theta-scheme
  integer(I32), parameter, public :: TSTEP_FRACTIONAL_STEP       = 2_I32**16

  ! Bitmask for theta-type time stepping schemes
  integer(I32), parameter, public :: TSTEP_THETA_TYPE            = TSTEP_THETA + TSTEP_TYPE&
                                                              + TSTEP_FRACTIONAL_STEP
  
  ! ID for general two-step theta scheme
    integer(I32), parameter, public :: TSTEP_THETA_GEN           = TSTEP_THETA + 0_I32
  
  ! ID for forward Euler scheme (=fully explicit two-step theta-scheme)
  integer(I32), parameter, public :: TSTEP_FORWARD_EULER         = TSTEP_THETA + 1_I32
  
  ! ID for backward Euler scheme (=fully implicit two-step theta-scheme)
  integer(I32), parameter, public :: TSTEP_BACKWARD_EULER        = TSTEP_THETA + 2_I32
  
  ! ID for Crank-Nicholson scheme (=semi-implicit second-order two-step theta-scheme,
  !                                 also termed trapezoidal rule)
  integer(I32), parameter, public :: TSTEP_CRANK_NICHOLSON       = TSTEP_THETA + 3_I32
  integer(I32), parameter, public :: TSTEP_TRAPEZOIDAL           = TSTEP_CRANK_NICHOLSON
  integer(I32), parameter, public :: TSTEP_THETA_2               = TSTEP_CRANK_NICHOLSON
  
  ! ID for general fractional-step theta-scheme
  integer(I32), parameter, public :: TSTEP_FSTHETA_GEN           = TSTEP_THETA + TSTEP_FRACTIONAL_STEP + 4_I32

  ! ID for second-order fractional-step theta scheme
  integer(I32), parameter, public :: TSTEP_FSTHETA_2             = TSTEP_THETA + TSTEP_FRACTIONAL_STEP + 5_I32
!</constantblock>

!<constantblock description="Identifiers for Runge-Kutta family of time stepping schemes">

  ! Strongly stability preserving Runge-Kutta time stepping family of schemes
  integer(I32), parameter, public :: TSTEP_SSPRK                 = 2_I32**16
  
  ! Semi-explicit Runge-Kutta time stepping family of schemes
  integer(I32), parameter, public :: TSTEP_SERP                  = 2_I32**17

  ! Semi-implicit Runge-Kutta time stepping family of schemes
  integer(I32), parameter, public :: TSTEP_SIRP                  = 2_I32**18

  ! Diagonally-implicit Runge-Kutta time stepping family of schemes
  integer(I32), parameter, public :: TSTEP_DIRK                  = 2_I32**19
  
  ! Singly-diagonally implicit Runge-Kutta time stepping family of schemes
  integer(I32), parameter, public :: TSTEP_SDIRK                 = 2_I32**20
                                                              
  ! Bitmask specifying embedded Runge-Kutta scheme
  integer(I32), parameter, public :: TSTEP_EMBEDDED              = 2_I32**29
  
  ! Bitmask specifying Butcher tableau form of Runge-Kutta scheme
  integer(I32), parameter, public :: TSTEP_BUTCHER_TABLEAU       = 2_I32**30

  ! Bitmask for Runge-Kutta-type time stepping schemes
  integer(I32), parameter, public :: TSTEP_RUNGE_KUTTA_TYPE      = TSTEP_RUNGE_KUTTA + TSTEP_TYPE&
                                                              + TSTEP_SSPRK&
                                                              + TSTEP_SERP&
                                                              + TSTEP_SIRP&
                                                              + TSTEP_DIRK&
                                                              + TSTEP_SDIRK&
                                                              + TSTEP_EMBEDDED
    
  ! ID for forward Euler method (implemented as Runge-Kutta method)
  integer(I32), parameter, public :: TSTEP_ERK_FORWARD_EULER     = TSTEP_RUNGE_KUTTA + 0_I32

  ! ID for explicit midpoint rule
  integer(I32), parameter, public :: TSTEP_ERK_MIDPOINT_RULE     = TSTEP_RUNGE_KUTTA + 1_I32

  ! ID for Heun`s method
  integer(I32), parameter, public :: TSTEP_ERK_HEUN              = TSTEP_RUNGE_KUTTA + 2_I32

  ! ID for Ralston`s method
  integer(I32), parameter, public :: TSTEP_ERK_RALSTON           = TSTEP_RUNGE_KUTTA + 3_I32

  ! ID for classical RK4 method
  integer(I32), parameter, public :: TSTEP_ERK_RK4               = TSTEP_RUNGE_KUTTA + 4_I32

  ! ID for 3/8-rule RK method
  integer(I32), parameter, public :: TSTEP_ERK_38_RULE           = TSTEP_RUNGE_KUTTA + 5_I32

  
  ! ID for embedded Heun-Euler method
  integer(I32), parameter, public :: TSTEP_EERK_HEUN_EULER       = TSTEP_RUNGE_KUTTA + TSTEP_EMBEDDED + 0_I32

  ! ID for embedded Bogacki–Shampine method
  integer(I32), parameter, public :: TSTEP_EERK_BOGACKI_SHAMPINE = TSTEP_RUNGE_KUTTA + TSTEP_EMBEDDED + 1_I32

  ! ID for embedded Fehlberg method
  integer(I32), parameter, public :: TSTEP_EERK_FEHLBERG         = TSTEP_RUNGE_KUTTA + TSTEP_EMBEDDED + 2_I32

  ! ID for embedded Cash-Karp method
  integer(I32), parameter, public :: TSTEP_EERK_CASH_KARP        = TSTEP_RUNGE_KUTTA + TSTEP_EMBEDDED + 3_I32

  ! ID for embedded Dormand–Prince method
  integer(I32), parameter, public :: TSTEP_EERK_DORMAND_PRINCE   = TSTEP_RUNGE_KUTTA + TSTEP_EMBEDDED + 4_I32

  ! ID for explicit SSP-(1,1)-RK scheme
  integer(I32), parameter, public :: TSTEP_SSPERK_1_1            = TSTEP_RUNGE_KUTTA + TSTEP_SSPRK + 0_I32

  ! ID for explicit SSP-(2,2)-RK scheme
  integer(I32), parameter, public :: TSTEP_SSPERK_2_2            = TSTEP_RUNGE_KUTTA + TSTEP_SSPRK + 1_I32

  ! ID for explicit SSP-(3,3)-RK scheme
  integer(I32), parameter, public :: TSTEP_SSPERK_3_3            = TSTEP_RUNGE_KUTTA + TSTEP_SSPRK + 2_I32
  
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

!<constantblock description="Strongly stability-preserving explicit Runge-Kutta Schemes">
  
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
    ! Identifier of the time-stepping algorithm
    integer(I32) :: ctimestep = TSTEP_UNDEFINED
    
    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! This determines the output level of the time-stepping algorithm
    integer :: ioutputLevel = 0

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    integer(I32) :: coutputModeError   = 0_I32
    integer(I32) :: coutputModeWarning = 0_I32
    integer(I32) :: coutputModeInfo    = 0_I32
    integer(I32) :: coutputModeVerbose = 0_I32

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Norm in which relative changes of solution should be measured
    integer :: isolNorm = 0

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Scaling parameter for explicit parts
    real(DP) :: dscaleExplicit = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Scaling parameter for implicit parts
    real(DP) :: dscaleImplicit = 0.0_DP
    
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
    real(DP) :: dStepPrevious = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Factor by which the current time step should be scaled
    ! if the current time step failed; must be < 1.0
    real(DP) :: dStepReductionFactor = 1.0_DP

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Total number of performed time steps
    integer :: nSteps = 0

    ! OUTPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Total number of rejected time steps
    integer :: nrejectedSteps = 0

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
    !   || U^{N+1} - U^{N} || <= EPSSTAG * dStep
    ! =0: ignore, do not stop time-stepping until final time is reached
    real(DP) :: depsSteady = 0.0_DP

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Two-step theta scheme
    type(t_thetaScheme), pointer :: p_rthetaScheme => null()

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Runge-Kutta scheme
    type(t_RungeKuttaScheme), pointer :: p_rRungeKuttaScheme => null()
    
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

  ! This data structure contains all settings/parameters for the
  ! two-step theta scheme.
  type t_thetaScheme

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Implicitness parameter for two-step theta-scheme
    real(DP) :: theta = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Implicitness parameter for fractional-step theta-scheme
    real(DP) :: alpha = 0.0_DP

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Implicitness parameter for fractional-step theta-scheme
    real(DP) :: gamma = 0.0_DP
    
    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Temporal vectors for old solutions
    type(t_vectorBlock), dimension(:), pointer :: RtempVectors => null()
    
  end type t_thetaScheme
  
!</typeblock>

  ! *****************************************************************************

!<typeblock>

  ! This data structure contains all settings/parameters for 
  ! multi-step Runge Kutta schemes
  type t_RungeKuttaScheme

    ! INPUT PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Number of stages for Runge-Kutta method
    integer :: nstages = 0
    
    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Coefficients of the Butcher tableau
    real(DP), dimension(:,:), pointer :: Da => null()
    real(DP), dimension(:), pointer   :: Db => null()
    real(DP), dimension(:), pointer   :: Dc => null()

    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Coefficients of the alpha-beta form
    real(DP), dimension(:,:), pointer :: Dalpha => null()
    real(DP), dimension(:,:), pointer :: Dbeta  => null()
    
    ! INTERNAL PARAMETER FOR THE TIME-STEPPING ALGORITHM
    ! Temporal vectors for old solutions
    type(t_vectorBlock), dimension(:), pointer :: RtempVectors => null()
    
  end type t_RungeKuttaScheme
  
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

contains

  ! *****************************************************************************

!<function>

  integer(I32) function tstep_igetID(stimestepName, bcheck)

!<description>   
  ! This routine returns the time stepping id to a given time stepping
  ! name. It is case-insensitive.
!</description>

!<result>
  ! id of the time stepping scheme
!</result>

!<input>
  ! Element name - one of the TSTEP_xxxx constants.
  character (LEN=*) :: stimestepName

  ! Check the time steping type. If set to TRUE and the time stepping
  ! name is invalid, the program is not stopped, but 0 is returned.
  logical, intent(in), optional :: bcheck
!</input>

!</function>

    character(len=len(stimestepName)+1) :: stimestep
    logical :: bchk

    ! SELECT CASE is not allowed for strings (although supported by a majority
    ! of compilers), therefore we have to use a couple of if-commands :(
    ! select case(trim(sys_upcase(stimestepName)))

    stimestep = trim(sys_upcase(stimestepName))

    ! -= Two-step theta scheme =-
    if (stimestep .eq. "TSTEP_THETA_GEN") then
      tstep_igetID = TSTEP_THETA_GEN
    elseif(stimestep .eq. "TSTEP_FORWARD_EULER") then
      tstep_igetID = TSTEP_FORWARD_EULER
    elseif(stimestep .eq. "TSTEP_BACKWARD_EULER") then
      tstep_igetID = TSTEP_BACKWARD_EULER
    elseif(stimestep .eq. "TSTEP_CRANK_NICHOLSON" .or.&
           stimestep .eq. "TSTEP_TRAPEZOIDAL" .or.&
           stimestep .eq. "TSTEP_TSTEP2_THETA") then
      tstep_igetID = TSTEP_CRANK_NICHOLSON
    elseif(stimestep .eq. "TSTEP_FSTHETA_GEN") then
      tstep_igetID = TSTEP_FSTHETA_GEN
    elseif(stimestep .eq. "TSTEP_FSTHETA_2") then
      tstep_igetID = TSTEP_FSTHETA_2

    ! -= explicit SSP-Runge-Kutta schemes =-  
    elseif (stimestep .eq. "TSTEP_SSPERK_1_1") then
      tstep_igetID = TSTEP_SSPERK_1_1
    elseif (stimestep .eq. "TSTEP_SSPERK_2_2") then
      tstep_igetID = TSTEP_SSPERK_2_2
    elseif (stimestep .eq. "TSTEP_SSPERK_3_3") then
      tstep_igetID = TSTEP_SSPERK_3_3

    ! -= explicit Runge-Kutta schemes =-

    ! -= diagonally-implicit Runge-Kutta schemes =-

    ! -= implicit Runge-Kutta schemes =-

    ! Unknown time stepping scheme
    else
      bchk = .false.
      if (present(bcheck)) bchk = bcheck
      if (.not. bchk) then
        call output_line('Error: Unknown time stepping scheme: ' // stimestepName, &
                        OU_CLASS_ERROR,OU_MODE_STD,'tstep_igetID')
        call sys_halt()
      else
        tstep_igetID = 0
      end if
    end if
    
  end function tstep_igetID

  !************************************************************************

!<function>

  elemental integer(I32) function tstep_getType (ctimestep) result (iresult)

!<description>
  ! Determines the 'type' of the time stepping scheme.
!</description>

!<result>
  ! The identifier of the 'type' of the time stepping scheme, ctimestep refers to.
!</result>

!<input>
  ! Time stepping type qualifier
  integer(I32), intent(in) :: ctimestep
!</input>

!</function>

    select case(iand(ctimestep,TSTEP_FAMILY))
    case (TSTEP_THETA)
      iresult = iand(ctimestep,TSTEP_THETA_TYPE)
    case (TSTEP_RUNGE_KUTTA)
      iresult = iand(ctimestep,TSTEP_RUNGE_KUTTA_TYPE)
    case default
      iresult = TSTEP_UNDEFINED
    end select
    
  end function tstep_getType
  
  !************************************************************************

!<function>

  elemental integer(I32) function tstep_getFamily (ctimestep) result (iresult)

!<description>
  ! Determines the 'family' of the time stepping scheme.
!</description>

!<result>
  ! The identifier of the 'family' of the time stepping scheme, ctimestep refers to.
!</result>

!<input>
  ! Time stepping type qualifier
  integer(I32), intent(in) :: ctimestep
!</input>

!</function>

    iresult = iand(ctimestep,TSTEP_FAMILY)
  
  end function tstep_getFamily

  ! ***************************************************************************

!<function>

  character(len=SYS_STRLEN) function tstep_getTypeName(ctimestep) result(sname)

!<description>
  ! This function returns a string which represents the name of the
  ! time stepping scheme 
!</description>

!<input>
  ! The time stepping scheme whose name is to be returned.
  integer(I32), intent(in) :: ctimestep
!</input>

!<result>
  ! A string representing the name of the time stepping scheme.
!</result>

!</function>
  
    select case(tstep_getType(ctimestep))
    ! -= Two-step theta scheme =-
    case(TSTEP_THETA_GEN)
      sname = 'TSTEP_THETA_GEN'
    case(TSTEP_FORWARD_EULER)
      sname = 'TSTEP_FORWARD_EULER'
    case(TSTEP_BACKWARD_EULER)
      sname = 'TSTEP_BACKWARD_EULER'
    case (TSTEP_CRANK_NICHOLSON)
      sname = 'TSTEP_CRANK_NICHOLSON'
    case (TSTEP_FSTHETA_GEN)
      sname = 'TSTEP_FSTHETA_GEN'
    case (TSTEP_FSTHETA_2)
      sname = 'TSTEP_FSTHETA_2'

    ! -= explicit SSP-Runge-Kutta schemes =-
    case (TSTEP_SSPERK_1_1)
      sname = 'TSTEP_SSPERK_1_1'
    case (TSTEP_SSPERK_2_2)
      sname = 'TSTEP_SSPERK_2_2'
    case (TSTEP_SSPERK_3_3)
      sname = 'TSTEP_SSPERK_3_3'

    ! -= explicit Runge-Kutta schemes =-

    ! -= diagonally-implicit Runge-Kutta schemes =-

    ! -= implicit Runge-Kutta schemes =-
      
    case default
      ! unknown element id
      sname = 'unknown'

    end select
    
  end function tstep_getTypeName

  ! ***************************************************************************

!<function>

  character(len=SYS_STRLEN) function tstep_getFamilyName(ctimestep) result(sname)

!<description>
  ! This function returns a string which represents the name of the
  ! time stepping family 
!</description>

!<input>
  ! The time stepping scheme whose name is to be returned.
  integer(I32), intent(in) :: ctimestep
!</input>

!<result>
  ! A string representing the name of the time stepping family.
!</result>

!</function>
  
    select case(tstep_getFamily(ctimestep))
    case (TSTEP_THETA)
      sname = 'TSTEP_THETA'
      
    case (TSTEP_RUNGE_KUTTA)
      sname = 'TSTEP_RUNGE_KUTTA'
      
    case default
      ! unknown element id
      sname = 'unknown'
      
    end select
  end function tstep_getFamilyName

  ! ***************************************************************************

!<subroutine>

  subroutine tstep_genButcherTableau(ctimestep,p_Da,p_Db,p_Dc,&
      ballocate,bcheck)

!<description>
    ! This subroutine generates the Butcher tableau for a given
    ! Runge-Kutta time stepping scheme   
!</description>

!<input>
  ! Switch to indicate whether the pointers should be allocated
  logical, intent(in) :: ballocate
!</input>

!<inputoutput>
  ! Identifier of the time stepping scheme 
  integer(I32), intent(inout) :: ctimestep
!</inputoutput>

!<output>
  ! Pointer to the arrays of the Butcher tableau
  real(DP), dimension(:,:), intent(out), pointer :: p_Da
  real(DP), dimension(:), intent(out), pointer   :: p_Db,p_Dc

  ! Check existence of Butcher tableau. If a Butcher tableau is
  ! available for the specified time stepping scheme, then bcheck is
  ! set to TRUE, otherwise to FALSE. If bcheck is not present and no
  ! Butcher tableau is available, then an error is thrown.
  logical, intent(out), optional :: bcheck
!</output>

!</subroutine>

    if (present(bcheck)) bcheck = .true.
    select case(tstep_getType(ctimestep))
    
    case default
      if (present(bcheck)) then
        bcheck = .false.
        return
      end  if
      
      call output_line('No Butcher tableau available for specified scheme!',&
          OU_CLASS_ERROR,OU_MODE_STD,'tstep_genButcherTableau')
      call sys_halt()
    end select
    
  end subroutine tstep_genButcherTableau

  ! ***************************************************************************

!<subroutine>

  subroutine tstep_genAlphaBetaTableau(ctimestep,p_Dalpha,p_Dbeta,p_Dc,&
      ballocate,bcheck)

!<description>
    ! This subroutine generates the alpha-beta tableau for a given
    ! Runge-Kutta time stepping scheme   
!</description>

!<input>
  ! Switch to indicate whether the pointers should be allocated
  logical, intent(in) :: ballocate
!</input>

!<inputoutput>
  ! Identifier of the time stepping scheme 
  integer(I32), intent(inout) :: ctimestep
!</inputoutput>
  
!<output>
  ! Pointer to the arrays of the alpha-beta tableau
  real(DP), dimension(:,:), intent(out), pointer :: p_Dalpha,p_Dbeta
  real(DP), dimension(:), intent(out), pointer   :: p_Dc

  ! Check existence of alpha-beta tableau. If an alpha-beta tableau is
  ! available for the specified time stepping scheme, then bcheck is
  ! set to TRUE, otherwise to FALSE. If bcheck is not present and no
  ! alpha-beta tableau is available, then an error is thrown.
  logical, intent(out), optional :: bcheck
!</output>

!</subroutine>

    if (present(bcheck)) bcheck = .true.
    select case(tstep_getType(ctimestep))
    case (TSTEP_SSPERK_1_1)

      if (ballocate) allocate(p_Dalpha(1,1),p_Dbeta(1,1),p_Dc(1))
      p_Dalpha = reshape((/ 1.0_DP /), shape(p_Dalpha))
      p_Dbeta  = reshape((/ 1.0_DP /), shape(p_Dbeta))
      p_Dc     = (/ 1.0_DP /)

    case (TSTEP_SSPERK_2_2)

      if (ballocate) allocate(p_Dalpha(2,2),p_Dbeta(2,2),p_Dc(2))
      p_Dalpha = reshape((/ 1.0_DP, 0.5_DP,&
                            0.0_DP, 0.5_DP /), shape(p_Dalpha))
      p_Dbeta  = reshape((/ 1.0_DP, 0.0_DP,&
                            0.0_DP, 0.5_DP /), shape(p_Dbeta))
      p_Dc     = (/ 1.0_DP, 1.0_DP /)
      
    case (TSTEP_SSPERK_3_3)

      if (ballocate) allocate(p_Dalpha(3,3),p_Dbeta(3,3),p_Dc(3))
      p_Dalpha = reshape((/ 1.0_DP, 0.75_DP, 1.0_DP/3.0_DP,&
                            0.0_DP, 0.25_DP, 0.0_DP,&
                            0.0_DP, 0.0_DP , 2.0_DP/3.0_DP /), shape(p_Dalpha))
      p_Dbeta  = reshape((/ 1.0_DP, 0.0_DP , 0.0_DP,&
                            0.0_DP, 0.25_DP, 0.0_DP,&
                            0.0_DP, 0.0_DP , 2.0_DP/3.0_DP /), shape(p_Dbeta))
      p_Dc     = (/ 1.0_DP, 1.0_DP, 1.0_DP /)
      
    case default
      if (present(bcheck)) then
        bcheck = .false.
        return
      end  if
      
      call output_line('No alpha-beta tableau available for specified scheme!',&
          OU_CLASS_ERROR,OU_MODE_STD,'tstep_genButcherTableau')
      call sys_halt()
    end select

    
    
  end subroutine tstep_genAlphaBetaTableau
  
end module timestepbase

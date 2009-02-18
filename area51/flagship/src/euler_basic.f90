!##############################################################################
!# ****************************************************************************
!# <name> euler_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global variables
!# which are required to solve the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!#  1.) euler_getNVAR
!#      Returns the number of variables depending on the spatial dimension
!#
!# </purpose>
!##############################################################################

module euler_basic

  use fparser
  use fsystem
  use problem
  use statistics
  use thermodynamics
  use triangulation

  implicit none

  private
  public :: t_euler
  public :: euler_getNVAR

  interface euler_getNVAR
    module procedure euler_getNvarProblemLevel
    module procedure euler_getNvarDescriptor
  end interface

!<constants>

!<constantblock description="Global type of mass matrix">

  ! no mass matrix, i.e., no time derivative
  integer, parameter, public :: MASS_ZERO       = 0
  
  ! consistent mass matrix
  integer, parameter, public :: MASS_CONSISTENT = 1

  ! lumped mass matrix
  integer, parameter, public :: MASS_LUMPED     = 2

!</constantblock>


!<constantblock description="Global type of system coupling approach">

  ! time-dependent flow
  integer, parameter, public :: SYSTEM_SEGREGATED = 0

  ! steady-state flow
  integer, parameter, public :: SYSTEM_ALLCOUPLED = 1

!</constantblock>

!<constantblock description="Global type of dissipation">

  ! Employ scalar dissipation (default)
  integer, parameter, public :: DISSIPATION_SCALAR = 1

  ! Employ tensorial dissipation
  integer, parameter, public :: DISSIPATION_TENSOR = 2

!</constantblock>


!<constantblock description="Global type of right-hand side">

  ! zero right-hand side
  integer, parameter, public :: RHS_ZERO     = 0

  ! analytical right-hand side
  integer, parameter, public :: RHS_ANALYTIC = 1

!</constantblock>


!<constantblock description="Global types of perturbation parameters">

  ! Perturbation parameter is chosen as in the NITSOL package
  integer, parameter, public :: PERTURB_NITSOL  = 1

  ! Perturbation parameter is chosen as SQRT(machine precision)
  integer, parameter, public :: PERTURB_SQRTEPS = 2

!</constantblock>


!<constantblock description="Global constants from Gas Dynamic">

  ! ratio of specific heats
  real(DP), parameter, public :: GAMMA = GAMMA_AIR
  real(DP), parameter, public :: RGAS  = R_AIR

  ! auxiliary parameters related to gamma
  real(DP), parameter, public :: G1  = GAMMA-1.0
  real(DP), parameter, public :: G2  = (GAMMA-1.0)/2.0
  real(DP), parameter, public :: G3  = 1.0/(GAMMA-1.0)
  real(DP), parameter, public :: G4  = 1.0/GAMMA
  real(DP), parameter, public :: G5  = GAMMA/(GAMMA-1.0)
  real(DP), parameter, public :: G6  = GAMMA-2.0
  real(DP), parameter, public :: G7  = (GAMMA-1.0)/(2.0*GAMMA)
  real(DP), parameter, public :: G8  = (GAMMA+1.0)/(2.0*GAMMA)
  real(DP), parameter, public :: G9  = 2.0/(GAMMA-1.0)
  real(DP), parameter, public :: G10 = 2.0/(GAMMA+1.0)
  real(DP), parameter, public :: G11 = 2.0*GAMMA/(GAMMA-1.0)
  real(DP), parameter, public :: G12 = (GAMMA-1.0)/(GAMMA+1.0)

!</constantblock>


!<constantblock description="Global constants for number of variables">

  ! number of solution components in 1D
  integer, parameter, public :: NVAR1D = 3

  ! number of solution components in 2D
  integer, parameter, public :: NVAR2D = 4

  ! number of solution components in 3D
  integer, parameter, public :: NVAR3D = 5

!</constantblock>

!</constants>

  !*****************************************************************************

!<types>
  
!<typeblock>

  ! This structure contains all required data to describe an instance
  ! of the compressible Euler flow benchmark application
  type t_euler

    ! Dimension of the problem
    integer :: ndimension

    ! Type of the mass matrix
    integer :: imasstype

    ! Type of mass antidiffusion
    integer :: imassantidiffusiontype

    ! Type of dissipation
    integer :: idissipationtype

    ! Type of coupling
    integer :: icoupled

    ! Type of preconditioner
    integer :: iprecond

    ! Type of the right-hand side
    integer :: irhstype

    ! Type of the element(s)
    integer :: ieltype

    ! Matrix format
    integer :: imatrixFormat

    ! System format
    integer :: isystemFormat
    
    ! Function parser for the right-hand side vector
    type(t_fparser) :: rfparserRHS

    ! Timer for the solution process
    type(t_timer) :: rtimerSolution

    ! Timer for the adaptation process
    type(t_timer) :: rtimerAdaptation

    ! Timer for the error estimation process
    type(t_timer) :: rtimerErrorEstimation

    ! Timer for the triangulation process
    type(t_timer) :: rtimerTriangulation

    ! Timer for the assembly of constant coefficient matrices
    type(t_timer) :: rtimerAssemblyCoeff
    
    ! Timer for the assembly of system matrices
    type(t_timer) :: rtimerAssemblyMatrix

    ! Timer for the assembly of residual/right-hand side vectors
    type(t_timer) :: rtimerAssemblyVector

    ! Timer for pre- and post-processing
    type(t_timer) :: rtimerPrePostprocess
    
  end type t_euler

!</typeblock>

contains

  !*****************************************************************************

!<function>

  pure function euler_getNvarProblemLevel(rproblemLevel) result(NVAR)

!<description>
    ! This function returns the number of flow variables
    ! using the information from the problem level structure
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel
!</input>

!<result>
    ! number of variables
    integer :: NVAR
!</result>
!</function>

    NVAR = rproblemLevel%rtriangulation%ndim + 2
    
  end function euler_getNvarProblemLevel

  !*****************************************************************************

!<function>

  pure function euler_getNvarDescriptor(rappDescriptor) result(NVAR)

!<description>
    ! This function returns the number of flow variables
    ! using the information from the application descriptor
!</description>

!<input>
    ! application descriptor
    type(t_euler), intent(IN) :: rappDescriptor
!</input>

!<result>
    ! number of variables
    integer :: NVAR
!</result>
!</function>

    NVAR = rappDescriptor%ndimension + 2
    
  end function euler_getNvarDescriptor
end module euler_basic

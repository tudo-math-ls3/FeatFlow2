!##############################################################################
!# ****************************************************************************
!# <name> structuresoperatorasm </name>
!# ****************************************************************************
!#
!# <purpose>
!# Encapsules structures used for the assembly of space-time operators
!# on one level of a discretisation hierarchy.
!# </purpose>
!##############################################################################

module structuresoperatorasm

  use fsystem
  
  use spatialdiscretisation
  use timediscretisation
  
  use analyticsolution
  
  use structuresgeneral
  use structuresdiscretisation
  use structuresoptcontrol
  
  use assemblytemplates

  implicit none

  private
  
!<constants>

!<constantblock description = "A list of operator tyes for the assembly">

  ! Navier-Stokes operator, primal equation
  integer, parameter, public :: OPTP_PRIMAL = 0

  ! Simple linearised Navier-Stokes, primal equation
  integer, parameter, public :: OPTP_PRIMALLIN_SIMPLE = 1
  
  ! Full linearised Navier-Stokes (Newton), primal equation
  integer, parameter, public :: OPTP_PRIMALLIN = 2
  
  ! Navier-Stokes operator, dual equation
  integer, parameter, public :: OPTP_DUAL = 3

  ! Navier-Stokes operator, linearised dual equation
  integer, parameter, public :: OPTP_DUALLIN_SIMPLE = 4

  ! Navier-Stokes operator, linearised dual equation (Newton)
  integer, parameter, public :: OPTP_DUALLIN = 5
  
  ! Full linearised Navier-Stokes (Newton), dual equation,
  ! nonlinear part in the RHS.
  integer, parameter, public :: OPTP_DUALLIN_RHS = 6

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This type encapsules structures necessary for the
  ! assembly of space-time operators.
  type t_spacetimeOperatorAsm
  
    ! Space discretisation, primal space
    type(t_blockDiscretisation), pointer :: p_rspaceDiscrPrimal => null()

    ! Space discretisation, dual space
    type(t_blockDiscretisation), pointer :: p_rspaceDiscrDual => null()

    ! Space discretisation, control space
    type(t_blockDiscretisation), pointer :: p_rspaceDiscrControl => null()
    
    ! Time discretisation. Defines the time stepping scheme. Primal space
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrPrimal => null()

    ! Time discretisation. Defines the time stepping scheme. Dual space
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrDual => null()

    ! Time discretisation. Defines the time stepping scheme. Control space
    type(t_timeDiscretisation), pointer :: p_rtimeDiscrControl => null()
    
    ! Assembly templates corresponding to the above space discretisation
    type(t_staticSpaceAsmTemplates), pointer :: p_rasmTemplates => null()

    ! Physics of the problem
    type(t_settings_physics), pointer :: p_rphysics => null()
    
    ! Parameters for the discretisation in space
    type(t_settings_discr), pointer :: p_rsettingsDiscr => null()

    ! Optimal-control parameters
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl => null()
    
    ! Reference to the analytic solution defining the RHS of the primal equation
    type(t_anSolution), pointer :: p_rrhsPrimal => null()

    ! Reference to the analytic solution defining the RHS of the dual equation
    type(t_anSolution), pointer :: p_rrhsDual => null()

    ! Reference to the analytic solution defining the target
    type(t_anSolution), pointer :: p_rtargetFlow => null()

    ! Reference to global data
    type(t_globalData), pointer :: p_rglobalData => null()
    
    ! Reference to the debug flags
    type(t_optcDebugFlags), pointer :: p_rdebugFlags => null()
  end type

!</typeblock>

  public :: t_spacetimeOperatorAsm

!<typeblock>

  ! A hierarchy of t_spacetimeOperatorAsm structures.
  type t_spacetimeOpAsmHierarchy

    ! Number of levels in the hierarchy.
    integer :: nlevels = 0
  
    ! The level info structures on all levels.
    type(t_spacetimeOperatorAsm), dimension(:), pointer :: p_RopAsmList => null()
  
  end type

!</typeblock>

  public :: t_spacetimeOpAsmHierarchy

!</types>

end module

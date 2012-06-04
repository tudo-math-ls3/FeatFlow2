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
  use structuresboundaryconditions
  
  use assemblytemplates
  use fespacehierarchybase
  use timescalehierarchy
  use spacetimehierarchy

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

  ! Encapsules analytic data for setting up discrete operators
  ! on the space-time cylinder
  type t_spacetimeOpAsmAnalyticData
  
    ! Definition of the boundary conditions
    type(t_optcBDC), pointer :: p_roptcBDC => null()

    ! Physics of the problem
    type(t_settings_physics), pointer :: p_rphysics => null()
    
    ! Parameters for the discretisation in space
    type(t_settings_spacediscr), pointer :: p_rsettingsSpaceDiscr => null()

    ! Optimal-control parameters
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl => null()
    
    ! Reference to the analytic solution defining the RHS of the primal equation
    type(t_anSolution), pointer :: p_rrhsPrimal => null()

    ! Reference to the analytic solution defining the RHS of the dual equation
    type(t_anSolution), pointer :: p_rrhsDual => null()

    ! Reference to global data
    type(t_globalData), pointer :: p_rglobalData => null()
    
    ! Reference to the debug flags
    type(t_optcDebugFlags), pointer :: p_rdebugFlags => null()
  end type

!</typeblock>

  public :: t_spacetimeOpAsmAnalyticData

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
    
    ! Analytic data necessary to set up operators.
    type(t_spacetimeOpAsmAnalyticData), pointer :: p_ranalyticData => null()
  end type

!</typeblock>

  public :: t_spacetimeOperatorAsm

!<typeblock>

  ! A hierarchy of t_spacetimeOperatorAsm structures.
  type t_spacetimeOpAsmHierarchy

    ! Analytic data necessary to set up operators.
    type(t_spacetimeOpAsmAnalyticData) :: ranalyticData
    
    ! A hierarchy of time levels, primal space
    type(t_timescaleHierarchy), pointer :: p_rtimeHierarchyPrimal => null()

    ! A hierarchy of time levels, dual space
    type(t_timescaleHierarchy), pointer :: p_rtimeHierarchyDual => null()

    ! A hierarchy of time levels, control space
    type(t_timescaleHierarchy), pointer :: p_rtimeHierarchyControl => null()
    
    ! Hierarchy of FEM spaces, primal space.
    type(t_feHierarchy), pointer :: p_rfeHierarchyPrimal => null()

    ! Hierarchy of FEM spaces, dual space.
    type(t_feHierarchy), pointer :: p_rfeHierarchyDual => null()

    ! Hierarchy of FEM spaces, control space.
    type(t_feHierarchy), pointer :: p_rfeHierarchyControl => null()
    
    ! Hierarchy of space assembly structures
    type(t_staticSpaceAsmHierarchy), pointer :: p_rstaticSpaceAsmHier => null()
  
    ! Boundary condition hierarchy for all space levels, primal and dual space
    type(t_optcBDCSpaceHierarchy), pointer :: p_roptcBDCSpaceHierarchy => null()

  end type

!</typeblock>

  public :: t_spacetimeOpAsmHierarchy

!</types>

  ! Get a space-time operator assembly structure based on a space- and time-level
  public :: stoh_getOpAsm_slvtlv

  ! Get a space-time operator assembly strcuture based on a space-time level
  public :: stoh_getOpAsm
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine stoh_getOpAsm_slvtlv (roperatorAsm,roperatorAsmHier,ispacelevel,itimelevel)

!<description>
  ! Creates a space-time operator assembly structure based on a hierarchy
  ! of operator assembly structures, a space-level and a time-level.
!</description>
 
!<input>
  ! A space-time operator assembly hierarchy.
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
  
  ! Space level.
  integer, intent(in) :: ispacelevel

  ! Time level.
  integer, intent(in) :: itimelevel
!</input>

!<outout>
  ! Space-time operator assembly structure to be created.
  type(t_spacetimeOperatorAsm), intent(out) :: roperatorAsm
!</outout>
  
!</subroutine>

    ! Fetch pointers to the structures from the hierarchy.
    
    roperatorAsm%p_rspaceDiscrPrimal => &
        roperatorAsmHier%p_rfeHierarchyPrimal%p_rfeSpaces(ispacelevel)%p_rdiscretisation
    
    roperatorAsm%p_rspaceDiscrDual => &
        roperatorAsmHier%p_rfeHierarchyDual%p_rfeSpaces(ispacelevel)%p_rdiscretisation
    
    roperatorAsm%p_rspaceDiscrControl => &
        roperatorAsmHier%p_rfeHierarchyControl%p_rfeSpaces(ispacelevel)%p_rdiscretisation
    
    roperatorAsm%p_rtimeDiscrPrimal => &
        roperatorAsmHier%p_rtimeHierarchyPrimal%p_RtimeLevels(itimelevel)
    
    roperatorAsm%p_rtimeDiscrDual => &
        roperatorAsmHier%p_rtimeHierarchyDual%p_RtimeLevels(itimelevel)
    
    roperatorAsm%p_rtimeDiscrControl => &
        roperatorAsmHier%p_rtimeHierarchyControl%p_RtimeLevels(itimelevel)
    
    roperatorAsm%p_rasmTemplates => &
        roperatorAsmHier%p_rstaticSpaceAsmHier%p_RasmTemplList(ispacelevel)
    
    roperatorAsm%p_ranalyticData => roperatorAsmHier%ranalyticData

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stoh_getOpAsm (roperatorAsm,rsthierarchy,ilevel,roperatorAsmHier)

!<description>
  ! Creates a space-time operator assembly structure based on a hierarchy
  ! of operator assembly structures, a space-level and a time-level.
!</description>
 
!<input>
  ! A space-time hierarchy
  type(t_spacetimeHierarchy), intent(in) :: rsthierarchy
  
  ! Level in the space-time hierarchy
  integer, intent(in) :: ilevel

  ! A space-time operator assembly hierarchy.
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
!</input>

!<outout>
  ! Space-time operator assembly structure to be created.
  type(t_spacetimeOperatorAsm), intent(out) :: roperatorAsm
!</outout>
  
!</subroutine>
  
    integer :: ispacelevel, itimelevel
    
    ! Get the space- and time-level from the space-time hierarchy
    call sth_getLevel (rsthierarchy,ilevel,&
        ispaceLevel=ispaceLevel,itimeLevel=itimeLevel)
        
    ! Set up the structure
    call stoh_getOpAsm_slvtlv (roperatorAsm,roperatorAsmHier,ispacelevel,itimelevel)

  end subroutine

end module

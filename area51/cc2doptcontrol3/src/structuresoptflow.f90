!##############################################################################
!# ****************************************************************************
!# <name> structuresoptc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying structures of the space-time optimal control solver.
!# </purpose>
!##############################################################################

module structuresoptflow

  use fsystem
  use storage
  use boundary
  use triangulation
  use cubature
  use paramlist
  use discretebc
  use discretefbc
  use fparser
  use linearsystemscalar
  use multilevelprojection
  
  use collection

  use spatialdiscretisation
  use timediscretisation
  use timescalehierarchy

  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimehierarchy

  use spacetimevectors
  use analyticsolution
  use spacetimeinterlevelprojection
  
  use assemblytemplates

  use constantsdiscretisation
  use structuresgeneral
  use structuresdiscretisation
  use structuresboundaryconditions
  use structuresoptcontrol
  use structuresoperatorasm
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! Structure collecting all settings of the space-time optflow solver.
  type t_settings_optflow

    !<!-- ------------------------- -->
    !<!-- MESH, DOMAIN, HIERARCHIES -->
    !<!-- ------------------------- -->

    ! An object for saving the domain
    type(t_boundary) :: rboundary

    ! Physics of the problem
    type(t_settings_physics) :: rphysicsPrimal

    ! Settings controlling the spatial discretisation (element, 
    ! cubature, stabilisation)
    type(t_settings_discr) :: rsettingsSpaceDiscr

    ! Parameters of the optimal control problem
    type(t_settings_optcontrol) :: rsettingsOptControl

    ! Space coarse mesh without refinement
    type(t_triangulation) :: rtriaCoarse
    
    ! Time coarse mesh without refinement
    type(t_timeDiscretisation) :: rtimeCoarse

    ! Settings that define the refinement in space
    type(t_settings_refinement) :: rrefinementSpace

    ! Settings that define the refinement in time
    type(t_settings_refinement) :: rrefinementTime
    
    ! A mesh hierarchy with all available space meshes.
    type(t_meshHierarchy) :: rmeshHierarchy
    
    ! A hierarchy of time levels
    type(t_timescaleHierarchy) :: rtimeHierarchy
    
    ! A level info hierarchy for the assembly of stuff on all levels.
    type(t_staticSpaceAsmHierarchy) :: rspaceAsmHierarchy
    
    ! A hierarchy of operator assembly structures for all levels.
    type(t_spacetimeOpAsmHierarchy) :: roperatorAsmHier
    
    ! A hierarchy of space levels for the primal space
    type(t_feHierarchy) :: rfeHierarchyPrimal

    ! A hierarchy of space levels for the dual space
    type(t_feHierarchy) :: rfeHierarchyDual

    ! A hierarchy of space levels for the control space
    type(t_feHierarchy) :: rfeHierarchyControl
    
    ! Projection hierarchy for the interlevel projection in the primal space.
    type(t_interlevelProjectionHier) :: rprjHierSpacePrimal

    ! Projection hierarchy for the interlevel projection in the dual space.
    type(t_interlevelProjectionHier) :: rprjHierSpaceDual

    ! Projection hierarchy for the interlevel projection in the control space.
    type(t_interlevelProjectionHier) :: rprjHierSpaceControl
    
    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierPrimal
    
    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierDual

    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierControl

    ! Projection hierarchy for the interlevel projection in space/time, primal space.
    type(t_sptiProjHierarchy) :: rprjHierSpaceTimePrimal

    ! Projection hierarchy for the interlevel projection in space/time, dual space.
    type(t_sptiProjHierarchy) :: rprjHierSpaceTimeDual

    ! Projection hierarchy for the interlevel projection in space/time, control space.
    type(t_sptiProjHierarchy) :: rprjHierSpaceTimeControl
    
    ! All debug flags used by the application
    type(t_optcDebugFlags) :: rdebugFlags

    !<!-- ------------------------------------------ -->
    !<!-- INITIAL CONDITION, BOUNDARY CONDITION, RHS -->
    !<!-- ------------------------------------------ -->

    ! Analytic solution defining the initial condition.
    type(t_anSolution) :: rinitialCondition

    ! Analytic solution defining the right-hand side, primal solution
    type(t_anSolution) :: rrhsPrimal

    ! Analytic solution defining the right-hand side, dual solution
    type(t_anSolution) :: rrhsDual

    ! The boundary conditions, analytic definition
    type(t_optcBDC) :: roptcBDC

    !<!-- ------------------------------- -->
    !<!-- APPLICATION SPECIFIC PARAMETERS -->
    !<!-- ------------------------------- -->

    ! An application specific parameter list that allows to pass parameters
    ! to callback routines.
    type(t_parlist), pointer :: p_rparlist => null()
    
    ! User defined global data.
    type(t_globalData) :: rglobalData
    
  end type
  
!</typeblock>

  public :: t_settings_optflow
  
  ! Realises a pointer to the problem structure
  type p_settings_optflow
    
    ! Pointer to the problem structure
    type(t_settings_optflow), pointer :: p_rsettings => null()
    
  end type
  
  public :: p_settings_optflow
  
!</types>

end module

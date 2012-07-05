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
  use spacetimeinterlevelprj
  
  use assemblytemplates

  use constantsdiscretisation
  use structuresgeneral
  use structuresdiscretisation
  use structuresboundaryconditions
  use structuresoptcontrol
  use structuresoperatorasm
  use structurespostproc
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! Structure collecting all settings of the space-time optflow solver.
  type t_settings_optflow

    !<!-- ------------------------- -->
    !<!-- GENERAL SETTINGS          -->
    !<!-- ------------------------- -->

    ! Physics of the problem
    type(t_settings_physics) :: rphysics

    ! Settings controlling the spatial discretisation (element, 
    ! cubature, stabilisation)
    type(t_settings_spacediscr) :: rsettingsSpaceDiscr

    ! Parameters of the optimal control problem
    type(t_settings_optcontrol) :: rsettingsOptControl

    ! Settings that define the refinement in space
    type(t_settings_refinement) :: rrefinementSpace

    ! Settings that define the refinement in time
    type(t_settings_refinement) :: rrefinementTime
    
    ! All debug flags used by the application
    type(t_optcDebugFlags) :: rdebugFlags
    
    ! Postprocessing settings
    type(t_optcPostprocessing) :: rpostproc

    !<!-- --------------------------------- -->
    !<!-- MESH, DOMAIN, TIME-DISCRETISATION -->
    !<!-- --------------------------------- -->

    ! An object for saving the domain
    type(t_boundary) :: rboundary

    ! Space coarse mesh without refinement
    type(t_triangulation) :: rtriaCoarse
    
    ! Time coarse mesh without refinement
    type(t_timeDiscretisation) :: rtimeCoarse

    !<!-- ------------------------- -->
    !<!-- HIERARCHIES IN SPACE      -->
    !<!-- ------------------------- -->

    ! A mesh hierarchy with all available space meshes.
    type(t_meshHierarchy) :: rmeshHierarchy
    
    ! A hierarchy of space levels for the primal space
    type(t_feHierarchy) :: rfeHierarchyPrimal

    ! A hierarchy of space levels for the dual space
    type(t_feHierarchy) :: rfeHierarchyDual

    ! A hierarchy of space levels for the control space
    type(t_feHierarchy) :: rfeHierarchyControl
    
    ! A level info hierarchy for the assembly of stuff on all levels.
    type(t_staticSpaceAsmHierarchy) :: rspaceAsmHierarchy
    
    ! Projection hierarchy for the interlevel projection in the primal space.
    type(t_interlevelProjectionHier) :: rprjHierSpacePrimal

    ! Projection hierarchy for the interlevel projection in the dual space.
    type(t_interlevelProjectionHier) :: rprjHierSpaceDual

    ! Projection hierarchy for the interlevel projection in the control space.
    type(t_interlevelProjectionHier) :: rprjHierSpaceControl
    
    !<!-- ------------------------- -->
    !<!-- HIERARCHIES IN TIME       -->
    !<!-- ------------------------- -->

    ! A hierarchy of time levels
    type(t_timescaleHierarchy) :: rtimeHierarchy

    !<!-- ------------------------- -->
    !<!-- HIERARCHIES IN SPACE-TIME -->
    !<!-- ------------------------- -->
    
    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierPrimal
    
    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierDual

    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy) :: rspaceTimeHierControl

    ! Projection hierarchy for the interlevel projection in space/time, primal space.
    type(t_sptiProjHierarchyBlock) :: rprjHierSpaceTimePrimal

    ! Projection hierarchy for the interlevel projection in space/time, dual space.
    type(t_sptiProjHierarchyBlock) :: rprjHierSpaceTimeDual

    ! Projection hierarchy for the interlevel projection in space/time, control space.
    type(t_sptiProjHierarchyBlock) :: rprjHierSpaceTimeControl
    
    ! A hierarchy of operator assembly structures for all levels.
    ! This is a central discretisation structure passed to all assembly routines.
    type(t_spacetimeOpAsmHierarchy) :: roperatorAsmHier
    
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
    
    ! Discrete initial condition on the maximum space level
    type(t_discreteInitCond) :: rdiscreteInitCond

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

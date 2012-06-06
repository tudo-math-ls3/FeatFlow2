!##############################################################################
!# ****************************************************************************
!# <name> kktsystemhierarchy </name>
!# ****************************************************************************
!#
!# <purpose>
!# Realises a hierarchy of solutions of a KKT system.
!# </purpose>
!##############################################################################

module kktsystemhierarchy

  use fsystem
  use genoutput
  use spatialdiscretisation
  use timediscretisation

  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral

  use fespacehierarchybase  
  use spacetimehierarchy

  use kktsystemspaces
  use kktsystem
  
  use structuresoperatorasm
  use spacematvecassembly

  implicit none
  
  private

!<types>

!<typeblock>
  
  ! A hierarchy of solutions of a KKT system
  type t_kktsystemHierarchy
  
    ! A space-time hierarchy based on the primal space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierPrimal => null()
    
    ! A space-time hierarchy based on the dual space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierDual => null()

    ! A space-time hierarchy based on the control space
    type(t_spaceTimeHierarchy), pointer :: p_rspaceTimeHierControl => null()

    ! Hierarchy of space-time operator assembly structures for the discretisation
    type(t_spacetimeOpAsmHierarchy), pointer :: p_roperatorAsmHier => null()

    ! Number of levels in the hierarchy
    integer :: nlevels = 0

    ! Array of KKT system (solutions) for all the levels
    type(t_kktsystem), dimension(:), pointer :: p_RkktSystems => null()
    
    ! Alternative solution on the maximum level or NULL if there is none.
    type(t_kktsystem), pointer :: p_rkktSystemMaxLevel => null()
  
  end type

!</typeblock>

  public :: t_kktsystemHierarchy

!<typeblock>
  
  ! A hierarchy of directional derivatives of a KKT system
  type t_kktsystemDirDerivHierarchy
  
    ! A connected KKT system hierarchy specifying the evaluation
    ! point of the derivative.
    type(t_kktsystemHierarchy), pointer :: p_rkktsystemHierarchy => null()

    ! Number of levels in the hierarchy
    integer :: nlevels = 0

    ! Array of KKT system (solutions) for all the levels
    type(t_kktsystemDirDeriv), dimension(:), pointer :: p_RkktSysDirDeriv => null()
  
  end type

!</typeblock>

  public :: t_kktsystemDirDerivHierarchy

!</types>

  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
  public :: kkth_initPrimalSol
  
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
  public :: kkth_initDualSol

  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
  public :: kkth_initControl
  
  ! Initialises a KKT system structure for level ilevel of the hierarchy
  public :: kkth_initKKTSystem
  
  ! Initialises a hierarchy of KKT systems.
  public :: kkth_initHierarchy

  ! Cleans up a hierarchy of KKT systems.
  public :: kkth_doneHierarchy
  
  ! Get the KKT system at a level
  public :: kkth_getKKTsystem

  ! Initialises a KKT system derivative structure for level ilevel of the hierarchy
  public :: kkth_initKKTSystemDirDeriv

  ! Initialises a hierarchy of KKT system directional derivatives
  public :: kkth_initHierarchyDirDeriv

  ! Cleans up a hierarchy of KKT systems directional derivatives
  public :: kkth_doneHierarchyDirDeriv
  
  ! Get the KKT system directional derivative at a level
  public :: kkth_getKKTsystemDirDeriv
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initPrimalSol (rvector,rkktsystemHierarchy,ilevel)
  
!<description>
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy
  
  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Vector to be created
  type(t_primalSpace), intent(out) :: rvector
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev
    type(t_feSpaceLevel), pointer :: p_rspaceDiscr
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierPrimal,ilev,&
        p_rspaceDiscr,p_rtimeDiscr)

    ! Create the vector
    call kktsp_initPrimalVector (rvector,&
        p_rspaceDiscr%p_rdiscretisation,p_rtimeDiscr)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initDualSol (rvector,rkktsystemHierarchy,ilevel)
  
!<description>
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy
  
  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Vector to be created
  type(t_dualSpace), intent(out) :: rvector
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev
    type(t_feSpaceLevel), pointer :: p_rspaceDiscr
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierDual,ilev,&
        p_rspaceDiscr,p_rtimeDiscr)

    ! Create the vector
    call kktsp_initDualVector (rvector,&
        p_rspaceDiscr%p_rdiscretisation,p_rtimeDiscr)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initControl (rvector,rkktsystemHierarchy,ilevel)
  
!<description>
  ! Initialises a primal solution according to level ilevel in the
  ! KKT hierarchy.
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy
  
  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Vector to be created
  type(t_controlSpace), intent(out) :: rvector
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev
    type(t_feSpaceLevel), pointer :: p_rspaceDiscr
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierControl,ilev,&
        p_rspaceDiscr,p_rtimeDiscr)

    ! Create the vector
    call kktsp_initControlVector (rvector,&
        p_rspaceDiscr%p_rdiscretisation,p_rtimeDiscr)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initKKTSystem (rkktsystem,rkktsystemHierarchy,ilevel,roperatorAsmHier)
  
!<description>
  ! Initialises a KKT system structure for level ilevel of the hierarchy
!</description>

!<input>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(in) :: rkktsystemHierarchy

  ! Level, the KKT system should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
  
  ! Hierarchy of all possible space-time operators
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
!</input>

!<output>
  ! KKT sytem to be created
  type(t_kktsystem), intent(out) :: rkktsystem
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev, ispacelevel, itimelevel
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemHierarchy%nlevels-ilevel

    ! Get the spatial and time discretisations of that level
    call sth_getLevel (rkktsystemHierarchy%p_rspaceTimeHierPrimal,ilev,&
        ispacelevel=ispacelevel, itimelevel=itimelevel)
    
    ! Initialise the KKT system on level i
    call kkt_initKKTsystem (rkktsystem,roperatorAsmHier,&
        ispacelevel, itimelevel)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initHierarchy (rkktsystemHierarchy,roperatorAsmHier,&
      rspaceTimeHierPrimal,rspaceTimeHierDual,rspaceTimeHierControl,ballocate,&
      rkktSystemMaxLevel)
  
!<description>
  ! Initialises a hierarchy of KKT systems.
!</description>

!<input>
  ! Hierarchy of space-time operator assembly structures.
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
  
  ! Space-time hierarchy in the primal space
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierPrimal

  ! Space-time hierarchy in the dual space
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierDual

  ! Space-time hierarchy in the control space
  type(t_spacetimeHierarchy), intent(in), target :: rspaceTimeHierControl
  
  ! Whether or not to allocate memory for the solutions
  logical, intent(in) :: ballocate
  
  ! OPTIONAL: A KKT system solution on the maximum level.
  ! If specified, this solution is used on the maximum level, and
  ! no memory is allocated in the hierarchy for the maximum level.
  type(t_kktSystem), intent(in), optional, target :: rkktSystemMaxLevel
!</input>

!<output>
  ! Hierarchy of KKT system solutions
  type(t_kktsystemHierarchy), intent(out) :: rkktsystemHierarchy
!</output>

!</subroutine>
    
    ! local variables
    integer :: i,nlevelallocate

    ! Remember the structures
    rkktsystemHierarchy%p_rspaceTimeHierPrimal => rspaceTimeHierPrimal
    rkktsystemHierarchy%p_rspaceTimeHierDual => rspaceTimeHierDual
    rkktsystemHierarchy%p_rspaceTimeHierControl => rspaceTimeHierControl
    rkktsystemHierarchy%p_roperatorAsmHier => roperatorAsmHier
    
    ! Numer of levels in the hierarchy?
    rkktsystemHierarchy%nlevels = rspaceTimeHierPrimal%nlevels
    
    ! Alternative solution on the max. level?
    nlevelallocate = rkktsystemHierarchy%nlevels
    if (present(rkktSystemMaxLevel)) then
      rkktsystemHierarchy%p_rkktSystemMaxLevel => rkktSystemMaxLevel
      nlevelallocate = nlevelallocate - 1
    end if

    nullify(rkktsystemHierarchy%p_RkktSystems)
    
    if (ballocate) then
    
      ! Create the KKT system solutions
      allocate(rkktsystemHierarchy%p_RkktSystems(rkktsystemHierarchy%nlevels))
      do i=1,nlevelallocate
        ! Initialise the KKT system on level i
        call kkth_initKKTSystem (&
            rkktsystemHierarchy%p_RkktSystems(i),rkktsystemHierarchy,i,&
            roperatorAsmHier)
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_doneHierarchy (rkktsystemHierarchy)
  
!<description>
  ! Cleans up a hierarchy of KKT systems.
!</description>

!<inputoutput>
  ! Hierarchy of KKT system solutions to be cleaned up
  type(t_kktsystemHierarchy), intent(inout) :: rkktsystemHierarchy
!</inputoutput>

!</subroutine>
    
    ! local variables
    integer :: i

    ! Remember the structures
    nullify(rkktsystemHierarchy%p_rspaceTimeHierPrimal)
    nullify(rkktsystemHierarchy%p_rspaceTimeHierDual)
    nullify(rkktsystemHierarchy%p_rspaceTimeHierControl)
    nullify(rkktsystemHierarchy%p_roperatorAsmHier)
    nullify(rkktsystemHierarchy%p_rkktSystemMaxLevel)
    
    if (associated(rkktsystemHierarchy%p_RkktSystems)) then
      ! Clean up the solutions on every level
      do i=size(rkktsystemHierarchy%p_RkktSystems),1,-1
        call kkt_doneKKTsystem (rkktsystemHierarchy%p_RkktSystems(i))
      end do

      deallocate(rkktsystemHierarchy%p_RkktSystems)
    end if
    
    rkktsystemHierarchy%nlevels = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_getKKTsystem (rkktsystemHierarchy,ilevel,p_rkktSystem)
  
!<description>
  ! Returns a pointer to the KKT system on level ilevel.
!</description>

!<input>
  ! Hierarchy of KKT system solutions.
  type(t_kktsystemHierarchy), intent(inout) :: rkktsystemHierarchy
  
  ! Level to which a pointer should be set up.
  ! If the value is <= 0, the KKT system at level NLMAX+ilevel is returned.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Pointer to the KKT system on level ilevel
  type(t_kktsystem), pointer :: p_rkktSystem
!</inputoutput>

!</subroutine>

    integer :: ilev
    
    ilev = ilevel
    if (ilev .eq. 0) ilev = rkktsystemHierarchy%nlevels+ilev
    
    ! Get the pointer
    if ((ilev .eq. rkktsystemHierarchy%nlevels) .and. &
        associated(rkktsystemHierarchy%p_rkktSystemMaxLevel)) then
      ! Take the alternative solution on the maximum level
      p_rkktSystem => rkktsystemHierarchy%p_rkktSystemMaxLevel
    else
      p_rkktSystem => rkktsystemHierarchy%p_RkktSystems(ilev)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initKKTSystemDirDeriv (&
      rkktsystemDirDeriv,rkktsystemDirDerivHier,ilevel)
  
!<description>
  ! Initialises a KKT system derivative for level ilevel of the hierarchy
!</description>

!<input>
  ! Hierarchy structure for the directional derivative.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsystemDirDerivHier

  ! Level, the KKT system derivative should be created corresponding to.
  ! If this is set to <= 0, rkktsystem will be created at level NLMAX-ilevel
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! KKT sytem derivative structure to be created
  type(t_kktsystemDirDeriv), intent(out) :: rkktsystemDirDeriv
!</output>

!</subroutine>
    
    ! local variables
    integer :: ilev
    type(t_kktSystem), pointer :: p_rkktsystem
    
    ilev = ilevel
    if (ilevel .le. 0) ilev = rkktsystemDirDerivHier%nlevels-ilevel

    ! Get the KKT system solution on that level
    call kkth_getKKTsystem (&
        rkktsystemDirDerivHier%p_rkktsystemHierarchy,ilevel,p_rkktsystem)

    ! Initialise the KKT system on level i
    call kkt_initKKTsystemDirDeriv (rkktsystemDirDeriv,p_rkktsystem)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_initHierarchyDirDeriv (rkktsystemDirDerivHier,rkktsystemHierarchy)
  
!<description>
  ! Initialises a hierarchy of KKT system directional derivatives
!</description>

!<input>
  ! Evaluation point of the derivative
  type(t_kktsystemHierarchy), intent(in), target :: rkktsystemHierarchy
!</input>

!<output>
  ! Structure for the directional derivative in the point rkktsystemHierarchy.
  type(t_kktsystemDirDerivHierarchy), intent(out) :: rkktsystemDirDerivHier
!</output>

!</subroutine>
    
    ! local variables
    integer :: i
    
    rkktsystemDirDerivHier%p_rkktsystemHierarchy => rkktsystemHierarchy
    rkktsystemDirDerivHier%nlevels = rkktsystemHierarchy%nlevels
    
    ! Allocate memory for the solutions
    allocate(rkktsystemDirDerivHier%p_RkktSysDirDeriv(rkktsystemDirDerivHier%nlevels))
    do i=1,rkktsystemDirDerivHier%nlevels
      ! Initialise the derivative on level i
      call kkth_initKKTSystemDirDeriv (&
          rkktsystemDirDerivHier%p_RkktSysDirDeriv(i),rkktsystemDirDerivHier,i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_doneHierarchyDirDeriv (rkktsystemDirDerivHier)
  
!<description>
  ! Cleans up a hierarchy of KKT systems.
!</description>

!<inputoutput>
  ! Hierarchy structure for the directional derivative.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsystemDirDerivHier
!</inputoutput>

!</subroutine>
    
    ! local variables
    integer :: i

    ! Release memory
    if (associated(rkktsystemDirDerivHier%p_RkktSysDirDeriv)) then
      ! Clean up the solutions on every level
      do i=size(rkktsystemDirDerivHier%p_RkktSysDirDeriv),1,-1
        call kkt_doneKKTsystemDirDeriv (rkktsystemDirDerivHier%p_RkktSysDirDeriv(i))
      end do

      deallocate(rkktsystemDirDerivHier%p_RkktSysDirDeriv)
    end if
    
    ! Final cleanup
    nullify(rkktsystemDirDerivHier%p_rkktsystemHierarchy)
    rkktsystemDirDerivHier%nlevels = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkth_getKKTsystemDirDeriv (rkktsystemDirDerivHier,ilevel,p_rkktSystemDirDeriv)
  
!<description>
  ! Returns a pointer to the KKT system on level ilevel.
!</description>

!<input>
  ! Hierarchy structure for the directional derivative.
  type(t_kktsystemDirDerivHierarchy), intent(inout) :: rkktsystemDirDerivHier
  
  ! Level to which a pointer should be set up.
  ! If the value is <= 0, the KKT system at level NLMAX+ilevel is returned.
  integer, intent(in) :: ilevel
!</input>

!<inputoutput>
  ! Pointer to the KKT system on level ilevel
  type(t_kktsystemDirDeriv), pointer :: p_rkktSystemDirDeriv
!</inputoutput>

!</subroutine>

    integer :: ilev
    
    ilev = ilevel
    if (ilev .eq. 0) ilev = rkktsystemDirDerivHier%nlevels+ilev
    
    ! Get the pointer
    p_rkktSystemDirDeriv => rkktsystemDirDerivHier%p_RkktSysDirDeriv(ilev)

  end subroutine

end module

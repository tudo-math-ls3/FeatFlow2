!##############################################################################
!# ****************************************************************************
!# <name> initoptflow </name>
!# ****************************************************************************
!#
!# <purpose>
!# Initialisation of the OptFlow solver.
!# </purpose>
!##############################################################################

module initoptflow

  use fsystem
  use storage
  use genoutput
  use paramlist
  use basicgeometry
  use boundary
  use triangulation
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  
  use spatialdiscretisation
  use timediscretisation
  
  use cubature
  use filtersupport
  use collection
  
  use analyticsolution
  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimehierarchy
  use timescalehierarchy
  
  use constantsdiscretisation
  use structuresdiscretisation
  use structuresboundaryconditions
  use structuresgeneral
  use structuresmain
  use structuresoptflow
  use structuresoptcontrol
  use hierarchies
  
  use spacematvecassembly
  use initmatrices
  use spacetimeinterlevelprojection
  
  use user_callback
  
  implicit none
  
  private
  
  ! Initialise the OptFlow solver
  public :: init_initOptFlow
  
  ! Clean up the OptFlow solver
  public :: init_doneOptFlow
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine init_initOptFlow (rparlist,rsettings,rsettingsSolver,ioutputLevel)
  
!<description>
  ! Reads standard information from the DAT file and initialises rsettingsSolver.
!</description>
  
!<input>
  ! Parameter list
  type(t_parlist), intent(in), target :: rparlist
  
  ! Basic program settings
  type(t_settings_main), intent(in) :: rsettings

  ! Amount of output during the initialisation.
  ! =0: no output. =1: basic output.
  integer, intent(in) :: ioutputLevel
!</input>
  
!<output>
  ! Settings structure for the optimal control solver.
  type(t_settings_optflow), intent(out), target :: rsettingsSolver
!</output>

!</subroutine>

    character(len=SYS_STRLEN) :: sstr,sstr2,sstr3
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    type(t_feSpaceLevel), pointer :: p_rfeSpacePrimal,p_rfeSpaceDual,p_rfeSpaceControl
    integer :: isuccess
    
    ! Put refrences to the parameter list to the settings structure.
    rsettingsSolver%p_rparlist => rparlist

    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Reading debug flags.")
    end if

    ! At first, read all debug flags.
    call struc_getDebugFlags (rparlist,rsettings%ssectionDebugFlags,&
        rsettingsSolver%rdebugFlags)

    ! Basic space discretisation settings
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising discretisation rules.")
    end if
    call struc_getSpaceDiscrSettings (rparlist,rsettings%ssectionDiscrSpace,&
        rsettingsSolver%rsettingsSpaceDiscr)

    ! Now, read the mesh and the domain
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising coarse mesh (space / time).")
    end if
    call init_initParamTria (rparlist,rsettingsSolver%rboundary,&
        rsettingsSolver%rtriaCoarse,rsettingsSolver%rtimeCoarse,&
        rsettings%ssectionParamTriang,rsettings%ssectionTimeDiscretisation)

    ! Get information about the refinement
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising refinement.")
    end if
    call struc_getRefinementParams (rparlist,&
        rsettingsSolver%rrefinementSpace,rsettingsSolver%rrefinementTime,&
        rsettings%ssectionDiscrSpace,rsettings%ssectionTimeDiscretisation)

    ! Generate hierarchies, in space...
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Refining meshes in space.")
    end if
    call init_initSpaceHierarchy (rsettingsSolver%rboundary,rsettingsSolver%rrefinementSpace,&
        rsettingsSolver%rtriaCoarse,rsettingsSolver%rmeshHierarchy,ioutputlevel)

    ! and time.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Refining meshes in time.")
    end if
    call init_initTimeHierarchy (rsettingsSolver%rrefinementTime,&
        rsettingsSolver%rtimeCoarse,rsettingsSolver%rtimeHierarchy,ioutputlevel)

    ! We now have all space and time meshes and their hierarchies.
    ! Now it is time to generate the discretisation structures that define
    ! which FEM spaces to use for the discretisation on each level.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising space discretisation hierarchy.")
    end if
    call init_initSpaceDiscrHier (&
        rsettingsSolver%rfeHierarchyPrimal,rsettingsSolver%rfeHierarchyDual,&
        rsettingsSolver%rfeHierarchyControl,&
        rsettingsSolver%rphysics,rsettingsSolver%rsettingsOptControl,&
        rsettingsSolver%rsettingsSpaceDiscr,rsettingsSolver%rmeshHierarchy,&
        rsettingsSolver%rboundary,ioutputLevel)
        
    ! Create interlevel projection structures for prlongation/restriction in space.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising projection hierarchy.")
    end if
    call init_initSpacePrjHierarchy (rsettingsSolver%rprjHierSpacePrimal,&
        rsettingsSolver%rfeHierarchyPrimal,rparlist,rsettings%ssectionProlRestSpacePrimal)
    
    call init_initSpacePrjHierarchy (rsettingsSolver%rprjHierSpaceDual,&
        rsettingsSolver%rfeHierarchyDual,rparlist,rsettings%ssectionProlRestSpaceDual)
    
    call init_initSpacePrjHierarchy (rsettingsSolver%rprjHierSpaceControl,&
        rsettingsSolver%rfeHierarchyControl,rparlist,rsettings%ssectionProlRestSpaceControl)

    ! Initialise the underlying space-time hierarchies of the solver.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising space-time hierarchy.")
    end if
    call init_initSpaceTimeHierarchy (rparlist,rsettings%ssectionSpaceTimeRefinement,&
        rsettingsSolver%rrefinementSpace,rsettingsSolver%rrefinementTime,&
        rsettingsSolver%rfeHierarchyPrimal,rsettingsSolver%rfeHierarchyDual,&
        rsettingsSolver%rfeHierarchyControl,rsettingsSolver%rtimeHierarchy,&
        rsettingsSolver%rspaceTimeHierPrimal,rsettingsSolver%rspaceTimeHierDual,&
        rsettingsSolver%rspaceTimeHierControl)
    if (ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line("Space-time hierarchy statistics:")
      call output_line("--------------------------------")
      call output_lbrk()
      call output_line("Primal space:")
      call output_lbrk()
      call sth_printHierStatistics(rsettingsSolver%rspaceTimeHierPrimal)
      call output_lbrk()
      call output_line("Dual space:")
      call output_lbrk()
      call sth_printHierStatistics(rsettingsSolver%rspaceTimeHierDual)
      call output_lbrk()
      call output_line("Control space:")
      call output_lbrk()
      call sth_printHierStatistics(rsettingsSolver%rspaceTimeHierControl)
    end if
        
    ! Initialise the space-time interlevel projection structures
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising space-time projection operators.")
    end if
    
    call init_initSpaceTimePrjHierarchy (rsettingsSolver%rprjHierSpaceTimePrimal,&
        rsettingsSolver%rspaceTimeHierPrimal,&
        rsettingsSolver%rprjHierSpacePrimal,rsettingsSolver%rphysics,&
        rparlist,rsettings%ssectionTimeDiscretisation)
        
    call init_initSpaceTimePrjHierarchy (rsettingsSolver%rprjHierSpaceTimeDual,&
        rsettingsSolver%rspaceTimeHierDual,&
        rsettingsSolver%rprjHierSpaceDual,rsettingsSolver%rphysics,&
        rparlist,rsettings%ssectionTimeDiscretisation)

    call init_initSpaceTimePrjHierarchy (rsettingsSolver%rprjHierSpaceTimeControl,&
        rsettingsSolver%rspaceTimeHierControl,&
        rsettingsSolver%rprjHierSpaceControl,rsettingsSolver%rphysics,&
        rparlist,rsettings%ssectionTimeDiscretisation)
        
    ! Init+Allocate memory for the matrices on all levels and create all
    ! static matrices.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising template matrices.")
    end if
    call inmat_initStaticAsmTemplHier (&
        rsettingsSolver%rspaceAsmHierarchy,&
        rsettingsSolver%rsettingsSpaceDiscr,&
        rsettingsSolver%rfeHierarchyPrimal)
        
    ! Set up operator assembly structures for the assembly of KKT stuff
    call smva_createOpAsmHier(rsettingsSolver%roperatorAsmHier,rsettingsSolver)

    if (ioutputLevel .ge. 2) then
      call output_lbrk()
      call output_line ("Matrix assembly template statistics:")
      call output_lbrk()

      call inmat_printTmplHierStatistic(rsettingsSolver%rspaceAsmHierarchy)
    end if

    ! Calculate template matrices, independent of the physics
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Calculating template matrices.")
      call output_line ("Level: [",bnoLineBreak=.true.)
    end if
    
    call inmat_calcStaticLevelAsmHier (&
        rsettingsSolver%rspaceAsmHierarchy,&
        rsettingsSolver%rsettingsSpaceDiscr,&
        ioutputlevel .ge. 1)

    if (ioutputLevel .ge. 1) then
      call output_line ("]",cdateTimeLogPolicy = OU_DTP_NONE)
    end if

    ! Get the boudary conditions.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising boundary conditions.")
    end if
    
    ! There is one section for the expressions,
    ! one for the primal and one for the dual boundary conditions
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "ssectionBoundaryExpressions",sstr,bdequote=.true.)
        
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "ssectionBoundaryCondPrimal",sstr2,bdequote=.true.)
        
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "ssectionBoundaryCondDual",sstr3,bdequote=.true.)
    
    call struc_initBDC (rsettingsSolver%roptcBDC,rparlist,&
        rsettingsSolver%rphysics,sstr,sstr2,sstr3)

    ! Initialise the physics parameter that describe the
    ! physical equation
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising physics.")
    end if
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "ssectionPhysics",sstr,bdequote=.true.)
    call struc_initPhysics (rparlist,rsettingsSolver%rphysics,sstr)

    ! Ok, now we can initialise the Optimal-Control settings, read in the
    ! target function etc.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising optimal control parameters.")
    end if
    call init_initOptControl (rparlist,rsettings%ssectionOptControl,&
        rsettingsSolver%rsettingsOptControl)

    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising target function.")
    end if
    
    select case (rsettingsSolver%rphysics%cequation)
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      ! Stokes, Navier-Stokes, 2D
      call init_initOptControlTargetFunc2D (rparlist,rsettings%ssectionOptControl,&
          rsettingsSolver%rphysics,rsettingsSolver%rsettingsSpaceDiscr,&
          rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
          rsettingsSolver%rfeHierarchyPrimal,rsettingsSolver%rtimeHierarchy,&
          rsettingsSolver%rboundary,rsettingsSolver%rsettingsOptControl)
    end select
    
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising constraints.")
    end if
    
    select case (rsettingsSolver%rphysics%cequation)
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      ! Stokes, Navier-Stokes, 2D
      call init_initOptControlConstraints (rparlist,rsettings%ssectionOptControl,&
          rsettingsSolver%rsettingsSpaceDiscr,&
          rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
          rsettingsSolver%rfeHierarchyPrimal,rsettingsSolver%rtimeHierarchy,&
          rsettingsSolver%rboundary,rsettingsSolver%rsettingsOptControl)
    end select

    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising RHS.")
    end if
    
    ! Initialise the RHS.
    
    ! Primal RHS
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "srhsPrimal",sstr,bdequote=.true.)
    select case (rsettingsSolver%rphysics%cequation)
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      ! Stokes, Navier-Stokes, 2D
      call init_initFunction (rparlist,sstr,rsettingsSolver%rrhsPrimal,&
          rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
          rsettingsSolver%rsettingsSpaceDiscr,rsettingsSolver%rfeHierarchyPrimal,&
          rsettingsSolver%rboundary,rsettingsSolver%rphysics,isuccess)
    end select
    if (isuccess .eq. 1) then
      call output_line ('Functions created by simulation not yet supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initStandardSolver')
      call sys_halt()
    end if

    ! Dual RHS
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "srhsDual",sstr,bdequote=.true.)
    select case (rsettingsSolver%rphysics%cequation)
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      ! Stokes, Navier-Stokes, 2D
      call init_initFunction (rparlist,sstr,rsettingsSolver%rrhsDual,&
          rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
          rsettingsSolver%rsettingsSpaceDiscr,rsettingsSolver%rfeHierarchyPrimal,&
          rsettingsSolver%rboundary,rsettingsSolver%rphysics,isuccess)
    end select
    if (isuccess .eq. 1) then
      call output_line ('Functions created by simulation not yet supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initStandardSolver')
      call sys_halt()
    end if
       
    ! Get the discretisation of the maximum level.
    p_rtimeDiscr => &
        rsettingsSolver%rtimeHierarchy%p_rtimeLevels(rsettingsSolver%rtimeHierarchy%nlevels)
        
    p_rfeSpacePrimal => rsettingsSolver%rfeHierarchyPrimal% &
        p_rfeSpaces(rsettingsSolver%rfeHierarchyPrimal%nlevels)

    p_rfeSpaceDual => rsettingsSolver%rfeHierarchyDual% &
        p_rfeSpaces(rsettingsSolver%rfeHierarchyDual%nlevels)

    p_rfeSpaceControl => rsettingsSolver%rfeHierarchyControl% &
        p_rfeSpaces(rsettingsSolver%rfeHierarchyControl%nlevels)

    ! Initialise the postprocessing
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising postprocessing.")
    end if
!    call init_initPostprocessing (rparlist,rsettings%ssectionSpaceTimePostprocessing,&
!        rsettingsSpaceDiscr,rsettingsSolver%roptcBDC,&
!        p_rtimeDiscr,p_rfeSpacePrimalDual%p_rdiscretisation,p_rfeSpacePrimal%p_rdiscretisation,&
!        rpostproc,rsettingsSolver)

    ! Call the first user defined init routine to fetch parameters from
    ! the settings structure to the global data structure.
    ! From that point on, we can call assembly routines that use callback
    ! routines.
    call user_initGlobalData (rsettingsSolver,rsettingsSolver%rglobalData)

    ! Initialise the initial condition.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising initial condition.")
    end if
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "sinitialCondition",sstr,bdequote=.true.)
    select case (rsettingsSolver%rphysics%cequation)
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
      ! Stokes, Navier-Stokes, 2D
      call init_initFunction (rparlist,sstr,rsettingsSolver%rinitialCondition,&
          rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
          rsettingsSolver%rsettingsSpaceDiscr,rsettingsSolver%rfeHierarchyPrimal,&
          rsettingsSolver%rboundary,rsettingsSolver%rphysics,isuccess)
    end select
    if (isuccess .eq. 1) then
      call output_line ('Functions created by simulation not yet supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initStandardSolver')
      call sys_halt()
    end if
    
    ! Invoke a routine to discretise the initial condition.
    call smva_createDiscreteInitCond (rsettingsSolver%rdiscreteInitCond,&
        rsettingsSolver%rinitialCondition,&
        rsettingsSolver%rfeHierarchyPrimal%nlevels,rsettingsSolver%rtimeHierarchy%nlevels,&
        rsettingsSolver%roperatorAsmHier,ioutputlevel)

!    ! Now we calculate all assembly template data.
!    !
!    if (ioutputLevel .ge. 1) then
!      call output_lbrk()
!      call output_line ("Calculating template optimal control matrices.")
!      call output_line ("Level: [",bnoLineBreak=.true.)
!    end if
!    
!    call inmat_initStaticTemplHierOptC (rsettingsSolver%rspaceAsmHierarchyOptC,&
!        rsettingsSolver%rfeHierPrimal)
!    
!    call inmat_calcStaticLvlAsmHierOptC (rsettingsSolver%rspaceAsmHierarchy,&
!        rsettingsSolver%rspaceAsmHierarchyOptC,rsettingsSolver%rphysics,&
!        rsettingsSolver%rstabilPrimal,rsettingsSolver%rstabilDual,&
!        rsettingsSolver,ioutputlevel .ge. 1)
!
!    if (ioutputLevel .ge. 1) then
!      call output_line ("]",cdateTimeLogPolicy = OU_DTP_NONE)
!    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_doneOptFlow (rsettings)
  
!<description>
  ! Cleans up the settings structure rsettings.
!</description>
  
!<inputoutput>
  ! Settings structure to clean up.
  type(t_settings_optflow), intent(inout) :: rsettings
!</inputoutput>

!</subroutine>
    ! Release all data from the global data structure.
    call user_doneGlobalData (rsettings%rglobalData)

    ! Release projection and FE hierarchies in space and space/time
    call sptipr_doneProjection(rsettings%rprjHierSpaceTimeControl)
    call sptipr_doneProjection(rsettings%rprjHierSpaceTimeDual)
    call sptipr_doneProjection(rsettings%rprjHierSpaceTimePrimal)
        
    call mlprj_releasePrjHierarchy(rsettings%rprjHierSpaceControl)
    call mlprj_releasePrjHierarchy(rsettings%rprjHierSpaceDual)
    call mlprj_releasePrjHierarchy(rsettings%rprjHierSpacePrimal)

    ! Release the space-time hierarchies.
    call sth_doneHierarchy(rsettings%rspaceTimeHierControl)
    call sth_doneHierarchy(rsettings%rspaceTimeHierDual)
    call sth_doneHierarchy(rsettings%rspaceTimeHierPrimal)

    ! Release all matrices.
    call inmat_doneStaticAsmTemplHier(rsettings%rspaceAsmHierarchy)

    ! Release the initial condition / RHS
    call ansol_done(rsettings%rinitialCondition)

    ! Release optimal control parameters: Target function,...
    call init_doneOptControl (rsettings%rsettingsOptControl)
    
    ! Release all discretisation hierarchies.
    call init_doneSpaceDiscrHier (&
        rsettings%rfeHierarchyPrimal,&
        rsettings%rfeHierarchyDual,&
        rsettings%rfeHierarchyControl)

    ! Release the time meshes...
    call tmsh_releaseHierarchy(rsettings%rtimeHierarchy)

    ! and the space meshes.
    call mshh_releaseHierarchy(rsettings%rmeshHierarchy)

    ! Finally release the coarse mesh, bondary definition, etc.
    call init_doneParamTria (&
        rsettings%rboundary,rsettings%rtriaCoarse,rsettings%rtimeCoarse)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initOptControl (rparlist,ssectionOptC,roptcontrol)
  
!<description>
  ! Initialises the optimal control data structure based on parameters
  ! in a parameter list
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the optimal control can be found.
  character(len=*), intent(in) :: ssectionOptC
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>
    
    ! Initialise the parameters from the parameter list.
    call soptc_initParOptControl (rparlist,ssectionOptC,roptcontrol)

    ! Nothing else to do.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initOptControlTargetFunc2D (rparlist,ssectionOptC,&
      rphysics,rsettingsSpaceDiscr,&
      rtriaCoarse,rrefinementSpace,rfeHierPrimal,rtimeHierarchy,rboundary,roptcontrol)
  
!<description>
  ! Initialises the target function for the optimal control problem
  ! based on the parameters in the DAT file.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the optimal control can be found.
  character(len=*), intent(in) :: ssectionOptC

  ! Physics of the problem
  type(t_settings_physics), intent(in), target :: rphysics

  ! Structure with space discretisation settings
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Underlying spatial coarse mesh of the problem.
  type(t_triangulation), intent(in) :: rtriaCoarse
  
  ! Definition of the underlying domain.
  type(t_boundary), intent(in) :: rboundary

  ! Description of the refinement of rtriaCoarse.
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target function matches the one of the primary
  ! function, this hierarchy is used to save memory.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal
  
  ! Underlying hierarchy of the time discretisation
  type(t_timescaleHierarchy), intent(in) :: rtimeHierarchy
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: stargetFunction
    integer :: isuccess
    
    ! Initialise the target function.
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        'stargetFunction',stargetFunction,bdequote=.true.)
        
    call init_initFunction (rparlist,stargetFunction,roptcontrol%rtargetFunction,&
        rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,&
        rphysics,isuccess)
    
    if (isuccess .eq. 1) then
      call output_line ('Function created by simulation not yet supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initOptControlTargetFunc2D')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initOptControlConstraints (rparlist,ssectionOptC,rsettingsSpaceDiscr,&
      rtriaCoarse,rrefinementSpace,rfeHierPrimal,rtimeHierarchy,rboundary,roptcontrol)
  
!<description>
  ! Initialises the constraint functions for the optimal control problem
  ! based on the parameters in the DAT file.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the optimal control can be found.
  character(len=*), intent(in) :: ssectionOptC

  ! Structure with space discretisation settings
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Underlying spatial coarse mesh of the problem.
  type(t_triangulation), intent(in) :: rtriaCoarse
  
  ! Definition of the underlying domain.
  type(t_boundary), intent(in) :: rboundary

  ! Description of the refinement of rtriaCoarse.
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target function matches the one of the primary
  ! function, this hierarchy is used to save memory.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal
  
  ! Underlying hierarchy of the time discretisation
  type(t_timescaleHierarchy), intent(in) :: rtimeHierarchy
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sfunction
    integer :: isuccess
    
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        'cdistVelConstraints',roptcontrol%rconstraints%cdistVelConstraints,&
        roptcontrol%rconstraints%cdistVelConstraints)

    ! DEPRECATED: ccontrolConstraintsType = cconstraintsType
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        'cdistVelConstType',roptcontrol%rconstraints%cdistVelConstType,&
        roptcontrol%rconstraints%cdistVelConstType)

!    if (roptcontrol%rconstraints%cdistVelConstType .eq. 1) then
!      allocate (roptcontrol%rconstraints%p_rumin1)
!      allocate (roptcontrol%rconstraints%p_rumax1)
!      allocate (roptcontrol%rconstraints%p_rumin2)
!      allocate (roptcontrol%rconstraints%p_rumax2)
!      
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionumin1',sfunction,"",bdequote=.true.)
!          
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rumin1,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptControlConstraints')
!        call sys_halt()
!      end if
!
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionumax1',sfunction,"",bdequote=.true.)
!
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rumax1,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptControlConstraints')
!        call sys_halt()
!      end if
!
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionumin2',sfunction,"",bdequote=.true.)
!
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rumin2,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptControlConstraints')
!        call sys_halt()
!      end if
!
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionumax2',sfunction,"",bdequote=.true.)
!
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rumax2,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptControlConstraints')
!        call sys_halt()
!      end if
!
!    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initOptStateConstraints (rparlist,ssectionOptC,rsettingsSpaceDiscr,&
      rtriaCoarse,rrefinementSpace,rfeHierPrimal,rtimeHierarchy,rboundary,roptcontrol)
  
!<description>
  ! Initialises the state constraint functions for the optimal control problem
  ! based on the parameters in the DAT file.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the optimal control can be found.
  character(len=*), intent(in) :: ssectionOptC

  ! Structure with space discretisation settings
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Underlying spatial coarse mesh of the problem.
  type(t_triangulation), intent(in) :: rtriaCoarse
  
  ! Definition of the underlying domain.
  type(t_boundary), intent(in) :: rboundary

  ! Description of the refinement of rtriaCoarse.
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target function matches the one of the primary
  ! function, this hierarchy is used to save memory.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal
  
  ! Underlying hierarchy of the time discretisation
  type(t_timescaleHierarchy), intent(in) :: rtimeHierarchy
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sfunction
    integer :: isuccess
    
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        'cstateConstraints',roptcontrol%rconstraints%cstateConstraints,&
        roptcontrol%rconstraints%cstateConstraints)

    call parlst_getvalue_int (rparlist,ssectionOptC,&
        'cstateConstraintsType',roptcontrol%rconstraints%cstateConstraintsType,0)

!    if (roptcontrol%rconstraints%ccontrolConstraintsType .eq. 1) then
!      allocate (roptcontrol%rconstraints%p_rymin1)
!      allocate (roptcontrol%rconstraints%p_rymax1)
!      allocate (roptcontrol%rconstraints%p_rymin2)
!      allocate (roptcontrol%rconstraints%p_rymax2)
!      
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionymin1',sfunction,"",bdequote=.true.)
!          
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rymin1,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptStateConstraints')
!        call sys_halt()
!      end if
!
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionymax1',sfunction,"",bdequote=.true.)
!
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rymax1,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptStateConstraints')
!        call sys_halt()
!      end if
!
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionymin2',sfunction,"",bdequote=.true.)
!
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rymin2,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptStateConstraints')
!        call sys_halt()
!      end if
!
!      call parlst_getvalue_string (rparlist,ssectionOptC,&
!          'ssectionymax2',sfunction,"",bdequote=.true.)
!
!      call init_initFunction (rparlist,sfunction,roptcontrol%rconstraints%p_rymax2,&
!          rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
!      if (isuccess .eq. 1) then
!        call output_line ('Cannot set up constraints! Invalid section: '//trim(sfunction), &
!            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptStateConstraints')
!        call sys_halt()
!      end if
!
!    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_doneOptControl (roptcontrol)
  
!<description>
  ! Cleans up information in the optimal control structure.
!</description>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Release control constraints
!    if (roptcontrol%rconstraints%ccontrolConstraintsType .eq. 1) then
!      call ansol_done(roptcontrol%rconstraints%p_rumin1)
!      call ansol_done(roptcontrol%rconstraints%p_rumax1)
!      call ansol_done(roptcontrol%rconstraints%p_rumin2)
!      call ansol_done(roptcontrol%rconstraints%p_rumax2)
!      
!      deallocate (roptcontrol%rconstraints%p_rumin1)
!      deallocate (roptcontrol%rconstraints%p_rumax1)
!      deallocate (roptcontrol%rconstraints%p_rumin2)
!      deallocate (roptcontrol%rconstraints%p_rumax2)
!    end if

    ! Release state constraints
!    if (roptcontrol%rconstraints%cstateConstraintsType .eq. 1) then
!      call ansol_done(roptcontrol%rconstraints%p_rymin1)
!      call ansol_done(roptcontrol%rconstraints%p_rymax1)
!      call ansol_done(roptcontrol%rconstraints%p_rymin2)
!      call ansol_done(roptcontrol%rconstraints%p_rymax2)
!      
!      deallocate (roptcontrol%rconstraints%p_rymin1)
!      deallocate (roptcontrol%rconstraints%p_rymax1)
!      deallocate (roptcontrol%rconstraints%p_rymin2)
!      deallocate (roptcontrol%rconstraints%p_rymax2)
!    end if

    ! Release the target function.
    call ansol_done(roptcontrol%rtargetFunction)
    
    ! Clean up the parameters
    call soptc_doneParOptControl (roptcontrol)
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  recursive subroutine init_initDiscreteAnalytFunc2D (rfunction,&
      ielementType,rboundary,smesh,rtriaCoarse,rrefinementSpace,&
      rsettingsSpaceDiscr,rfeHierarchy,ilevel,rphysics)
  
!<description>
  ! Creates a analytical-function structure that resembles a discrete function.
!</description>

!<input>
  ! Physics of the problem
  type(t_settings_physics), intent(in), target :: rphysics
  
  ! Element to use.
  ! =-1: use the default element in rsettingsSpaceDiscr.
  integer, intent(in) :: ielementType

  ! Definition of the domain
  type(t_boundary), intent(in), target :: rboundary

  ! Underlying mesh. May be ="", in this case the default mesh in rtriaCoarse is used.
  character(len=*), intent(in) :: smesh
  
  ! Default mesh.
  type(t_triangulation), intent(in) :: rtriaCoarse

  ! Description of the refinement in space, how rfeHierarchy was created
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! Space discretisation settings
  type(t_settings_spacediscr), intent(in), target :: rsettingsSpaceDiscr

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target function matches the one of the primary
  ! function, this hierarchy is used to save memory.
  type(t_feHierarchy), intent(in) :: rfeHierarchy
  
  ! Refinement level, the analytical function should resemble.
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Function structure to create.
  type(t_anSolution), intent(out) :: rfunction
!</output>

!</subroutine>

    ! local variables
    integer :: ieltype,iavaillevel
    type(t_collection) :: rcollection
    type(t_spaceDiscrParams) :: rdiscrParams
    type(t_settings_spacediscr), target :: rsettingsDiscrLocal

    ! Initialise the local discretisation structure
    rsettingsDiscrLocal = rsettingsSpaceDiscr
    if (ielementType .ne. -1) then
      ! Change the element type
      rsettingsDiscrLocal%ielementType = ielementType
    end if

    ! Initialise the structure encapsuling the parameters of the discretisation.
    ! The structure is placed in the collection
    
    rdiscrParams%p_rphysics => rphysics
    
    ! Create the hierarchy for the primal space
    rdiscrParams%cspace = CCSPACE_PRIMAL
    rdiscrParams%p_rfeHierarchyPrimal => null()

    ! Can we reuse our hierarchy?
    if (smesh .eq. "") then
    
      ! Is the refinement level smaller than NLMIN?
      if (ilevel .lt. rrefinementSpace%npreref+1) then
        ! We have to create that level.

        rdiscrParams%p_rsettingsSpaceDiscr => rsettingsDiscrLocal
        rcollection%DquickAccess(:) = transfer(rdiscrParams,rcollection%DquickAccess(:))
        
        ! Create the basic analytic solution: Mesh, discretisation structure,...
        call ansol_init (rfunction,ilevel,&
            rtriaCoarse,1,rcollection%IquickAccess(1),&
            fgetDist1LvDiscr,rcollection,rboundary)
      else
      
        ! Ok, here we can hope to reuse existing data structures.

        ! Get the element type from the discretisation
        ieltype = rsettingsSpaceDiscr%ielementType
        
        if ((ielementType .eq. -1) .or. &
            (ielementType .eq. ieltype)) then
          
          ! Very nice, the element type matches.
          rcollection%IquickAccess(1) = ieltype

          ! Get the maximum available level in rfeHierarchy
          iavaillevel = min(rfeHierarchy%nlevels,ilevel-rrefinementSpace%npreref)
          
          ! Use the default discrezisation
          rdiscrParams%p_rsettingsSpaceDiscr => rsettingsSpaceDiscr
          rcollection%DquickAccess(:) = transfer(rdiscrParams,rcollection%DquickAccess(:))
          
          ! And now create the basic function. Only in case ilevel>NLMAX,
          ! new levels are generated.
          call ansol_init (rfunction,ilevel-rrefinementSpace%npreref,&
              rfeHierarchy%p_rfeSpaces(iavaillevel)%p_rdiscretisation,iavaillevel,&
              rcollection%IquickAccess(1),fgetDist1LvDiscr,rcollection)
          
        else
        
          ! Ok, a little bit different, the element type is different.
          ! That means, we can reuse the triangulation but not the discretisation.
          !
          ! Get the maximum available level in rfeHierarchy
          iavaillevel = min(rfeHierarchy%nlevels,ilevel-rrefinementSpace%npreref)

          ! Use the local discretisation withthe new element
          rdiscrParams%p_rsettingsSpaceDiscr => rsettingsDiscrLocal
          rcollection%DquickAccess(:) = transfer(rdiscrParams,rcollection%DquickAccess(:))
          
          ! And now create the basic function. Only in case ilevel>NLMAX,
          ! new levels are generated.
          call ansol_init (rfunction,ilevel-rrefinementSpace%npreref,&
              rfeHierarchy%rmeshHierarchy%p_Rtriangulations(iavaillevel),iavaillevel,&
              rcollection%IquickAccess(1),fgetDist1LvDiscr,rcollection)
          
        end if
        
      end if
    
    else
    
      ! Mesh is different. Then we have to do everything by hand...

      rdiscrParams%p_rsettingsSpaceDiscr => rsettingsDiscrLocal
      rcollection%DquickAccess(:) = transfer(rdiscrParams,rcollection%DquickAccess(:))
      
      ! Create the basic analytic solution: Mesh, discretisation structure,...
      ! As we pass the mesh file here, we also pass the original ilevel
      ! because the mesh has to be refined up to that.
      call ansol_init (rfunction,ilevel,&
          NDIM2D,smesh,rcollection%IquickAccess(1),&
          fgetDist1LvDiscr,rcollection,rboundary)
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initFunction (rparlist,ssection,rfunction,&
      rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierarchy,rboundary,&
      rphysics,ierror)
  
!<description>
  ! Reads in and sets up a function. ssection is the name of a section
  ! containing all parameters that configure the function. The routine creates
  ! a structure rfunction based on these parameters.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found that specify the target function.
  character(len=*), intent(in) :: ssection

  ! Space discretisation settings
  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr

  ! Coarse mesh, corresponding to rfeHierarchy.
  type(t_triangulation), intent(in) :: rtriaCoarse

  ! Description of the refinement in space, how rfeHierarchy was created
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target function matches the one of the primary
  ! function, this hierarchy is used to save memory.
  type(t_feHierarchy), intent(in) :: rfeHierarchy

  ! Definition of the domain.
  type(t_boundary), intent(in), target :: rboundary

  ! Physics of the problem
  type(t_settings_physics), intent(in), target :: rphysics
!</input>

!<output>
  ! Function structure to create.
  type(t_anSolution), intent(out) :: rfunction
  
  ! Success-flag. Returns whether this routine successfully created the function.
  ! =0: Function was successfully created.
  ! =1: Function is defined as forward simulation and was not possible to create.
  !     The caller must invoke init_initFunctionBySimulation to create the function.
  integer, intent(out) :: ierror
  
!</output>

!</subroutine>

    ! local variables
    integer :: ctype,iid,i
    real(dp) :: dstartTime,dtimeMax
    character(len=PARLST_LENLINEBUF), dimension(:), allocatable :: Sexpressions
    character(len=SYS_STRLEN) :: smesh,sfunctionFile
    integer :: ilevel,ielementType,idelta
    integer :: ntimesteps,ncomponents,ifileformat
    type(t_parlstSection), pointer :: p_rsection
    
    ! Is the section available?
    call parlst_querysection(rparlist, ssection, p_rsection)
    if (.not. associated(p_rsection)) then
      call output_line ("Section does not exist: "//trim(ssection), &
          OU_CLASS_ERROR,OU_MODE_STD,"init_initFunction")
      call sys_halt()
    end if
    
    ! Type of the function?
    call parlst_getvalue_int (rparlist,ssection,"ctype",ctype,-1)
        
    ierror = 0
        
    if (ctype .eq. 4) then
      ! Forget it, there are not enough parameters here to create the function.
      ! This situation is more complicated!
      ierror = 1
      return
    end if
    
    ! Get the number of components.
    call parlst_getvalue_int (rparlist,ssection,"ncomponents",ncomponents,0)
    
    if (ncomponents .le. 0) then
      call output_line ("Invalid number of components!", &
          OU_CLASS_ERROR,OU_MODE_STD,"init_initFunction")
      call sys_halt()
    end if

    if (ncomponents .gt. 0) then
      ! Read the expressions
      allocate (Sexpressions(ncomponents))
      
      do i=1,ncomponents
        call parlst_getvalue_string (rparlist,ssection,&
            "sexpression",Sexpressions(i),"",isubstring=i,bdequote=.true.)
      end do
      
    end if
    
    ! Id, level, element type

    call parlst_getvalue_int (rparlist,ssection,&
        "iid",iid,0)

    call parlst_getvalue_int (rparlist,ssection,&
        "ilevel",ilevel,0)

    call parlst_getvalue_int (rparlist,ssection,&
        "ielementType",ielementType,-1)

    ! File format. Formatted or unformatted.

    call parlst_getvalue_int (rparlist,ssection,&
        "ifileformat",ifileformat,1)

    ! Mesh, the function defining files

    call parlst_getvalue_string (rparlist,ssection,&
        "smesh",smesh,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfunctionFile",sfunctionFile,"",bdequote=.true.)

    ! Time discretisation

    call parlst_getvalue_int (rparlist,ssection,&
        "idelta",idelta,1)

    ! Set it up...
    select case (ctype)
    case (-1)
      ! Analytic function given by a callback routine
      call ansol_init(rfunction,ncomponents,0)
      rfunction%iid = iid
      return
    
    case (0)
      ! Zero function
      call ansol_init(rfunction,ncomponents)
      rfunction%iid = iid
      return

    case (3)
      ! Analytical expression as function.
      call ansol_init (rfunction,ncomponents)
      call ansol_configExpressions (rfunction,Sexpressions)
      rfunction%iid = iid
      
      ! Release the text expressions again
      if (ncomponents .gt. 0) then
        deallocate(Sexpressions)
      end if
      
      return

    end select
    
    ! Create an analytical function that resembles a discrete function.
    call init_initDiscreteAnalytFunc2D (rfunction,ielementType,rboundary,&
        smesh,rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierarchy,&
        ilevel,rphysics)
        
    rfunction%iid = iid
    
    ! Finally, set up the function.
    select case (ctype)
      
    case (1)
    
      ! Read stationary function from hard disc
      call ansol_configStationaryFile (rfunction,sfunctionFile,ifileformat .eq. 1)
      
    case (2)
    
      ! Read nonstationary function from hard disc
      
      call parlst_getvalue_double (rparlist,ssection,&
          "dstartTime",dstartTime)

      call parlst_getvalue_double (rparlist,ssection,&
          "dtimeMax",dtimeMax)
      
      call parlst_getvalue_int (rparlist,ssection,&
          "ntimesteps",ntimesteps)

      call ansol_configNonstationaryFile (rfunction, &
          dstartTime,dtimeMax,ntimesteps,&
          "("""//trim(sfunctionFile)//"."",I5.5)",&
          0,idelta,ifileformat .eq. 1)
          
    end select

  end subroutine

!  ! ***************************************************************************
!  
!!<subroutine>
!
!  subroutine init_initPostprocessing (rparlist,ssection,rsettingsSpaceDiscr,&
!      rboundaryConditions,rtimeDiscr,rspaceDiscr,rspaceDiscrPrimal,rpostproc,rsettings)
!
!!<description>
!  ! Initialises rdiscrData with standard values from rsettings.
!!</description>
!
!!<input>
!  ! Parameter list
!  type(t_parlist), intent(in) :: rparlist
!  
!  ! Section where the parameters can be found.
!  character(len=*), intent(in) :: ssection
!  
!  ! Structure with space discretisation settings
!  type(t_settings_spacediscr), intent(in) :: rsettingsSpaceDiscr
!
!  ! Boundary conditions to use.
!  type(t_optcBDC), intent(in), target  :: rboundaryConditions
!
!  ! Underlying space discretisation
!  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr
!
!  ! Underlying space discretisation in the primal space
!  type(t_blockDiscretisation), intent(in), target :: rspaceDiscrPrimal
!
!  ! Underlying time discretisation
!  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
!  
!  ! Global settings structure.
!  type(t_settings_optflow), intent(in), target :: rsettings
!!</input>
!
!!<output>
!  ! Postprocessing structure, to set up with data.
!  type(t_optcPostprocessing), intent(out) :: rpostproc
!!</output>
!
!!</subroutine>
!
!    ! local variables
!    character(len=SYS_STRLEN) :: sstr, sparam
!    integer :: isuccess, npoints, i
!
!    ! Basic initialisation of the postprocessing structure
!    call optcpp_initpostprocessing (rpostproc,rsettings%rphysics,CCSPACE_PRIMALDUAL,&
!        rboundaryConditions,rtimeDiscr,rspaceDiscr,rspaceDiscrPrimal)
!        
!    ! Read remaining parameters from the DAT file.
!    !
!    ! UCD export
!    call parlst_getvalue_int (rparlist,ssection,&
!        "ioutputUCD",rpostproc%ioutputUCD,0)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        "sfilenameUCD",rpostproc%sfilenameUCD,"",bdequote=.true.)
!
!    ! Export of solution and control.
!    call parlst_getvalue_int (rparlist,ssection,&
!        "cwriteFinalSolution",rpostproc%cwriteFinalSolution,1)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        "sfinalSolutionFileName",rpostproc%sfinalSolutionFileName,&
!        "",bdequote=.true.)
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "cwriteFinalControl",rpostproc%cwriteFinalControl,1)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        "sfinalControlFileName",rpostproc%sfinalControlFileName,&
!        "",bdequote=.true.)
!
!    ! function value calculation
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "icalcFunctionalValues",rpostproc%icalcFunctionalValues,0)
!
!    ! Body forces
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "icalcForces",rpostproc%icalcForces,0)
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "ibodyForcesBdComponent",rpostproc%ibodyForcesBdComponent,2)
!
!    call parlst_getvalue_double (rparlist,ssection,&
!        "dbdForcesCoeff1",rpostproc%dbdForcesCoeff1,rsettings%rphysics%dnuConst)
!
!    call parlst_getvalue_double (rparlist,ssection,&
!        "dbdForcesCoeff2",rpostproc%dbdForcesCoeff2,0.1_DP * 0.2_DP**2)
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "iwriteBodyForces",rpostproc%iwriteBodyForces,0)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        "sfilenameBodyForces",rpostproc%sfilenameBodyForces,&
!        "",bdequote=.true.)
!
!    ! Flux
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "icalcFlux",rpostproc%icalcFlux,0)
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "iwriteFlux",rpostproc%iwriteFlux,0)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        'sfilenameFlux',rpostproc%sfilenameFlux,"",bdequote=.true.)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        'dfluxline',sstr,"",bdequote=.true.)
!        
!    call parlst_getvalue_string (rparlist,ssection,&
!        'dfluxline',sstr,"",bdequote=.true.)
!    if (sstr .ne. "") then
!      ! Read the start/end coordinates
!      read(sstr,*) rpostproc%Dfluxline(1),rpostproc%Dfluxline(2),&
!          rpostproc%Dfluxline(3),rpostproc%Dfluxline(4)
!    end if
!
!    ! internal Energy
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "icalcKineticEnergy",rpostproc%icalcKineticEnergy,1)
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "iwriteKineticEnergy",rpostproc%iwriteKineticEnergy,0)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        'sfilenameKineticEnergy',rpostproc%sfilenameKineticEnergy,"",bdequote=.true.)
!
!    ! Error calculation
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "icalcError",rpostproc%icalcError,0)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        "ssectionReferenceFunction",sstr,"",bdequote=.true.)
!    
!    if (sstr .eq. "") &
!        rpostproc%icalcError = 0
!    
!    if (rpostproc%icalcError .ne. 0) then
!    
!      isuccess = 1
!      select case (rsettings%rphysics%cequation)
!      case (0,1)
!        ! Stokes, Navier-Stokes, 2D
!    
!        ! Initialise the reference solution
!        call init_initFunction (rparlist,sstr,rpostproc%ranalyticRefFunction,&
!            rsettings%rtriaCoarse,rsettings%rrefinementSpace,&
!            rsettingsSpaceDiscr,rsettings%rfeHierPrimal,rsettings%rboundary,isuccess)
!      end select
!      
!      if (isuccess .eq. 1) then
!        call output_line ("Reference solution could not be calculated!", &
!            OU_CLASS_ERROR,OU_MODE_STD,"init_initPostprocessing")
!        call sys_halt()
!      end if
!    end if
!    
!    ! Init the points to evaluate
!    npoints = parlst_querysubstrings (rparlist, ssection, &
!        "CEVALUATEPOINTVALUES")
! 
!    if (npoints .gt. 0) then
!      allocate (rpostproc%p_DcoordsPointEval(NDIM2D,npoints))
!      allocate (rpostproc%p_ItypePointEval(NDIM2D,npoints))
!    
!      ! Read the points
!      do i=1,npoints
!        call parlst_getvalue_string (rparlist, ssection, &
!            "CEVALUATEPOINTVALUES", sparam, "", i)
!        read (sparam,*) rpostproc%p_DcoordsPointEval(1,i),rpostproc%p_DcoordsPointEval(2,i),&
!            rpostproc%p_ItypePointEval(1,i),rpostproc%p_ItypePointEval(2,i)
!      end do
!    end if
!
!    call parlst_getvalue_int (rparlist,ssection,&
!        "iwritePointValues",rpostproc%iwritePointValues,0)
!
!    call parlst_getvalue_string (rparlist,ssection,&
!        'sfilenamePointValues',rpostproc%sfilenamePointValues,"",bdequote=.true.)
!
!  end subroutine
!
!  ! ***************************************************************************
!  
!!<subroutine>
!
!  subroutine init_donePostprocessing (rpostproc)
!
!!<description>
!  ! Creans up the postprocessing structure.
!!</description>
!
!!<inputoutput>
!  ! Postprocessing structure, to set up with data.
!  type(t_optcPostprocessing), intent(inout) :: rpostproc
!!</inputoutput>
!
!!</subroutine>
!
!    ! Release the reference solution
!    if (rpostproc%icalcError .ne. 0) then
!      call ansol_done(rpostproc%ranalyticRefFunction)
!    end if
!    
!    ! Release point coordinates for point evaluation
!    if (associated(rpostproc%p_DcoordsPointEval)) then
!      deallocate (rpostproc%p_DcoordsPointEval)
!      deallocate (rpostproc%p_ItypePointEval)
!    end if
!    
!    call optcpp_donepostprocessing(rpostproc)
!
!  end subroutine

end module

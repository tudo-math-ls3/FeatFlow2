!##############################################################################
!# ****************************************************************************
!# <name> initsolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains initialisation routines that initialise settings
!# based on parameters in the DAT files.
!#
!# Routines in this module:
!#
!# 1.) init_initStandardSolver
!#     -> Read in all parameters, create mesh and FEM hierarchies, create a 
!#        nonlinear solver structure ready to use.
!#
!# 2.) init_doneStandardSolver
!#     -> Release a nonlinear solver structure
!#
!# 3.) init_initStartVector
!#     -> Create a start vector to be used in the nonlinear solver.
!#
!# 4.) init_implementInitCondRHS
!#     -> Implement the initial condition corresponding to a space-time
!#        solution into a discretised space-time RHS vector.
!#
!# 5.) init_discretiseRHS
!#     -> Discretise an analytically given RHS.
!#
!# 6.) init_getSpaceDiscrSettings
!#     -> Reads space discretisation settings from the DAT file.
!#
!# Auxiliary routines:
!#
!# </purpose>
!##############################################################################

module initsolver

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
  use fespacehierarchy
  use spacetimehierarchy
  
  use constantsoptc
  use structuresoptc
  use structuresspacetimelinsol
  use structuresoptcspacetimenlsol
  use structuresoptflow
  
  use structuresmain
  
  use spacediscretisation
  
  use spatialbcdef
  use initmatrices
  use postprocessing
  use user_callback
  
  use timescalehierarchy
  use spacematvecassembly
  use spacetimevectors
  use spacetimelinearsystem
  use forwardbackwardsimulation
  use spacetimeinterlevelprojection
  use nonlinearoneshotspacetimesolver
  use timeboundaryconditions
  use timerhsevaluation
  
  implicit none
  
  private
  
  public :: init_initStandardSolver
  public :: init_doneStandardSolver
  public :: init_initStartVector
  public :: init_implementInitCondRHS
  public :: init_discretiseRHS
  public :: init_getSpaceDiscrSettings
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine init_initStandardSolver (rparlist,rsettings,&
      rsettingsSolver,rnlstsolver,rpostproc,rrhs,ioutputLevel)
  
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
  
  ! Nonlinear space-time solver structure configuring the solver.
  type(t_nlstsolver), intent(out) :: rnlstsolver
  
  ! Postprocessing data.    
  type(t_optcPostprocessing), intent(out) :: rpostproc

  ! Analytic solution defining the RHS of the equation.
  type(t_anSolution), intent(out) :: rrhs
!</output>

!</subroutine>

    character(len=SYS_STRLEN) :: sstr,sstr2
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    type(t_feSpaceLevel), pointer :: p_rfeSpacePrimalDual
    type(t_settings_discr) :: rsettingsSpaceDiscr
    integer :: isuccess
    
    ! Put refrences to the parameter list to the settings structure.
    rsettingsSolver%p_rparlist => rparlist

    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Reading debug flags.")
    end if

    ! At first, read all debug flags.
    call init_getDebugFlags (rparlist,rsettings%ssectionDebugFlags,&
        rsettingsSolver%rdebugFlags)

    ! Read the stabilisation that is necessary for the discretisation.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising stabilisation.")
    end if
    call init_initStabil (rparlist,&
        rsettingsSolver%rstabilPrimal,rsettingsSolver%rstabilDual,&
        rsettings%ssectionDiscrSpace)
    
    ! Now, read the mesh and the domain
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising coarse mesh.")
    end if
    call init_initParamTria (rparlist,rsettingsSolver%rboundary,&
        rsettingsSolver%rtriaCoarse,rsettingsSolver%rtimeCoarse,&
        rsettings%ssectionParamTriang,rsettings%ssectionTimeDiscretisation)

    ! Get information about the refinement
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising refinement.")
    end if
    call init_getRefinementParams (rparlist,&
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

    ! Fetch the space discretisation parameters.
    call init_getSpaceDiscrSettings (rparlist,rsettings%ssectionDiscrSpace,&
        rsettingsSpaceDiscr)

    ! We now have all space and time meshes and their hierarchies.
    ! Now it is time to generate the discretisation structures that define
    ! which FEM spaces to use for the discretisation on each level.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising space discretisation hierarchy.")
    end if
    call init_initSpaceDiscrHier (rparlist,rsettingsSpaceDiscr,&
        rsettingsSolver%rboundary,rsettingsSolver%rmeshHierarchy,&
        rsettingsSolver,ioutputLevel)
        
    ! Create interlevel projection structures for prlongation/restriction in space.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising projection hierarchy.")
    end if
    call init_initSpacePrjHierarchy (rsettingsSolver%rprjHierSpacePrimal,&
        rsettingsSolver%rfeHierPrimal,rparlist,rsettings%ssectionProlRestSpace)
    call init_initSpacePrjHierarchy (rsettingsSolver%rprjHierSpacePrimalDual,&
        rsettingsSolver%rfeHierPrimalDual,rparlist,rsettings%ssectionProlRestSpace)

    ! Initialise the underlying space-time hierarchies of the solver.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising space-time hierarchy.")
    end if
    call init_initSpaceTimeHierarchy (rparlist,"SPACETIME-REFINEMENT",&
        rsettingsSolver%rrefinementSpace,rsettingsSolver%rrefinementTime,&
        rsettingsSolver%rfeHierPrimal,rsettingsSolver%rfeHierPrimalDual,&
        rsettingsSolver%rtimeHierarchy,&
        rsettingsSolver%rspaceTimeHierPrimal,rsettingsSolver%rspaceTimeHierPrimalDual)    
    if (ioutputLevel .ge. 2) then
      call output_line("Space-time hierarchy statistics:")
      call output_lbrk()
      call sth_printHierStatistics(rsettingsSolver%rspaceTimeHierPrimalDual)
    end if
        
    ! Initialise the space-time interlevel projection structures
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising space-time projection operators.")
    end if
    call init_initSpaceTimePrjHierarchy (rsettingsSolver%rprjHierSpaceTimePrimal,&
        rsettingsSolver%rtimeCoarse,rsettingsSolver%rspaceTimeHierPrimal,&
        rsettingsSolver%rprjHierSpacePrimal,rparlist,rsettings%ssectionTimeDiscretisation)
    call init_initSpaceTimePrjHierarchy (rsettingsSolver%rprjHierSpaceTimePrimalDual,&
        rsettingsSolver%rtimeCoarse,rsettingsSolver%rspaceTimeHierPrimalDual,&
        rsettingsSolver%rprjHierSpacePrimalDual,rparlist,rsettings%ssectionTimeDiscretisation)
        
    ! Init+Allocate memory for the matrices on all levels and create all 
    ! static matrices.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising template matrices.")
    end if
    call inmat_initStaticAsmTemplHier (rsettingsSolver%rspaceAsmHierarchy,&
        rsettingsSolver%rfeHierPrimal,rsettingsSolver%rfeHierMass,&
        rsettingsSolver%rstabilPrimal,rsettingsSolver%rstabilDual)

    ! Calculate template matrices, independent of the physics
    if (ioutputLevel .ge. 1) then
      call output_line ("Calculating template matrices.")
      call output_line ("Level: [",bnoLineBreak=.true.)
    end if
    
    call inmat_calcStaticLevelAsmHier (rsettingsSolver%rspaceAsmHierarchy,&
        rsettingsSolver,ioutputlevel .ge. 1)

    if (ioutputLevel .ge. 1) then
      call output_line ("]")
    end if

    ! Initialise the nonlinear solver
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising nonlinear space-time solver.")
    end if
    call init_initNlSpaceTimeSolver (rsettingsSolver,rparlist,&
        rsettings%ssectionSpaceTimeSolver,rsettingsSolver%rtimeCoarse,&
        rnlstsolver)

    ! Get the boudary conditions.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising boundary conditions.")
    end if
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "ssectionBoundaryExpressions",sstr,bdequote=.true.)
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "ssectionBoundaryConditions",sstr2,bdequote=.true.)
    call init_initBoundaryConditions (rparlist,sstr,sstr2,&
        rsettingsSolver%roptcBDC)

    ! Initialise the physics parameter that describe the 
    ! physical equation
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising physics.")
    end if
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "ssectionPhysics",sstr,bdequote=.true.)
    call init_initPhysics (rparlist,rsettingsSolver%rphysicsPrimal,sstr)

    ! Ok, now we can initialise the Optimal-Control settings, read in the
    ! target flow etc.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising optimal control parameters.")
    end if
    call init_initOptControl (rparlist,rsettings%ssectionOptControl,&
        rsettingsSolver%rsettingsOptControl)

    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising target flow.")
    end if
    call init_initOptControlTargetFlow (rparlist,rsettings%ssectionOptControl,&
        rsettingsSpaceDiscr,&
        rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
        rsettingsSolver%rfeHierPrimal,rsettingsSolver%rtimeHierarchy,&
        rsettingsSolver%rboundary,rsettingsSolver%rsettingsOptControl)

    ! Initialise the initial condition and RHS.
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising initial condition.")
    end if
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "sinitialCondition",sstr,bdequote=.true.)
    call init_initFlow (rparlist,sstr,rsettingsSolver%rinitialCondition,&
        rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
        rsettingsSpaceDiscr,rsettingsSolver%rfeHierPrimal,rsettingsSolver%rboundary,&
        isuccess)
    if (isuccess .eq. 1) then
      call output_line ('Flow created by simulation not yet supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initStandardSolver')
      call sys_halt()
    end if

    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising RHS.")
    end if
    call parlst_getvalue_string (rparlist,rsettings%ssectionOptControl,&
        "srhs",sstr,bdequote=.true.)
    call init_initFlow (rparlist,sstr,rrhs,&
        rsettingsSolver%rtriaCoarse,rsettingsSolver%rrefinementSpace,&
        rsettingsSpaceDiscr,rsettingsSolver%rfeHierPrimal,rsettingsSolver%rboundary,&
        isuccess)
    if (isuccess .eq. 1) then
      call output_line ('Flow created by simulation not yet supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initStandardSolver')
      call sys_halt()
    end if
       
    ! Get the discretisation of the maximum level.
    p_rtimeDiscr => &
        rsettingsSolver%rtimeHierarchy%p_rtimeLevels(rsettingsSolver%rtimeHierarchy%nlevels)
    p_rfeSpacePrimalDual => rsettingsSolver%rfeHierPrimalDual% &
        p_rfeSpaces(rsettingsSolver%rfeHierPrimalDual%nlevels)

    ! Initialise the postprocessing
    if (ioutputLevel .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising postprocessing.")
    end if
    call init_initPostprocessing (rparlist,rsettings%ssectionSpaceTimePostprocessing,&
        rsettingsSpaceDiscr,rsettingsSolver%roptcBDC,&
        p_rtimeDiscr,p_rfeSpacePrimalDual%p_rdiscretisation,&
        rpostproc,rsettingsSolver)

    ! Call the first user defined init routine to fetch parameters from
    ! the settings structure to the global data structure.
    ! From that point on, we can call assembly routines that use callback
    ! routines.
    call user_initGlobalData (rsettingsSolver,rrhs,rsettingsSolver%rglobalData)

    ! Now we calculate all assembly template data.
    !
    if (ioutputLevel .ge. 1) then
      call output_line ("Calculating template optimal control matrices.")
      call output_line ("Level: [",bnoLineBreak=.true.)
    end if
    
    call inmat_initStaticTemplHierOptC (rsettingsSolver%rspaceAsmHierarchyOptC,&
        rsettingsSolver%rfeHierPrimal)
    
    call inmat_calcStaticLvlAsmHierOptC (rsettingsSolver%rspaceAsmHierarchy,&
        rsettingsSolver%rspaceAsmHierarchyOptC,rsettingsSolver%rphysicsPrimal,&
        rsettingsSolver%rstabilPrimal,rsettingsSolver%rstabilDual,&
        rsettingsSolver,ioutputlevel .ge. 1)

    if (ioutputLevel .ge. 1) then
      call output_line ("]")
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_doneStandardSolver (rsettings,rnlstsolver,rpostproc)
  
!<description>
  ! Cleans up the settings structure rsettings.
!</description>
  
!<inputoutput>
  ! Settings structure to clean up.
  type(t_settings_optflow), intent(inout) :: rsettings

  ! Nonlinear space-time solver structure configuring the solver,
  ! to be cleaned up.
  type(t_nlstsolver), intent(inout) :: rnlstsolver

  ! Postprocessing data.    
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!</subroutine>
    ! Release all data from the global data structure.
    call user_doneGlobalData (rsettings%rglobalData)

    ! Clean up the solver.
    call stnlsinit_doneSolver (rnlstsolver)
 
    ! Release postprocessing structures
    call init_donePostprocessing (rpostproc)

    ! Release projection and FE hierarchies in space and space/time
    call sptipr_doneProjection(rsettings%rprjHierSpaceTimePrimalDual)
    call sptipr_doneProjection(rsettings%rprjHierSpaceTimePrimal)
        
    call mlprj_releasePrjHierarchy(rsettings%rprjHierSpacePrimalDual)
    call mlprj_releasePrjHierarchy(rsettings%rprjHierSpacePrimal)

    ! Release the space-time hierarchies.
    call sth_doneHierarchy(rsettings%rspaceTimeHierPrimalDual)
    call sth_doneHierarchy(rsettings%rspaceTimeHierPrimal)

    ! Release all matrices.
    call inmat_doneStaticTemplHierOptC(rsettings%rspaceAsmHierarchyOptC)
    call inmat_doneStaticAsmTemplHier(rsettings%rspaceAsmHierarchy)

    ! Release the initial condition / RHS
    call ansol_done(rsettings%rinitialCondition)

    ! Release optimal control parameters: Target flow,...
    call init_doneOptControl (rsettings%rsettingsOptControl)
    
    ! Release all discretisation hierarchies.
    call init_doneSpaceDiscrHier (rsettings)

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

  subroutine init_initPhysics (rparlist,rphysics,ssection)
  
!<description>
  ! Reads information about physics from the DAT file and saves them to rphysics.
!</description>
  
!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! A structure receiving the physics of the equation
  type(t_settings_physics), intent(out) :: rphysics
!</output>

!</subroutine>

    ! Get the viscosity parameter.
    ! Note that the parameter in the DAT file is 1/nu !
    call parlst_getvalue_double (rparlist,ssection,&
        'RE',rphysics%dnu,1000.0_DP)
    rphysics%dnu = 1E0_DP/rphysics%dnu

    ! Which type of problem to discretise? (Stokes, Navier-Stokes,...)
    call parlst_getvalue_int (rparlist,ssection,&
        'iequation',rphysics%iequation,0)

    ! Type of subproblem (gradient tensor, deformation tensor,...)
    call parlst_getvalue_int (rparlist,ssection,&
        'isubEquation',rphysics%isubEquation,0)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initBoundaryConditions (rparlist,ssectionExpr,ssectionBDC,roptcBDC)
  
!<description>
  ! Basic initialisation of the boundary condition structure.
!</description>
  
!<input>
  ! Parameter list
  type(t_parlist), intent(in), target :: rparlist
  
  ! Section in the parameter list defining the boundary expressions
  character(len=*), intent(in) :: ssectionExpr
  
  ! Section in the parameter list defining the boundary conditions
  character(len=*), intent(in) :: ssectionBDC
!</input>
  
!<output>
  ! Boundary condition structure.
  type(t_optcBDC), intent(out) :: roptcBDC
!</output>

!</subroutine>

    ! Save a pointer to the parameter list and the name of the sections
    ! containing the BDC's.
    roptcBDC%p_rparamListBDC => rparlist
    roptcBDC%ssectionBdExpressions = ssectionExpr 
    roptcBDC%ssectionBdConditions = ssectionBDC   

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initStabil (rparlist,rstabilPrimal,rstabilDual,ssection)
  
!<description>
  ! Reads information about stabilisation from the DAT file and saves them 
  ! to rstabilPrimal/rstabilDual.
!</description>
  
!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! Stabilisation parameters for the primal equation
  type(t_settings_stabil), intent(out) :: rstabilPrimal

  ! Stabilisation parameters for the dual equation
  type(t_settings_stabil), intent(out) :: rstabilDual
!</output>

!</subroutine>

    ! Get the parameters
    call parlst_getvalue_int (rparlist, ssection, &
        'IUPWIND1', rstabilPrimal%cupwind)
    call parlst_getvalue_int (rparlist, ssection, &
        'IUPWIND2', rstabilDual%cupwind)
    call parlst_getvalue_double (rparlist, ssection, &
        'DUPSAM1', rstabilPrimal%dupsam)
    call parlst_getvalue_double (rparlist, ssection, &
        'DUPSAM2', rstabilDual%dupsam)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initForwardSimPostproc (rparlist,ssection,rpostprocessing)
  
!<description>
  ! Reads information about postprocessing in a forward simulation from the DAT 
  ! file and saves them to rpostprocessing.
!</description>
  
!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  type(t_fbsimPostprocessing), intent(out) :: rpostprocessing
!</output>

!</subroutine>

    ! Get the parameters
    call parlst_getvalue_int (rparlist, ssection, &
        'ioutputUCD', rpostprocessing%ioutputUCD)
    call parlst_getvalue_string (rparlist, ssection, &
        'sfilenameUCD', rpostprocessing%sfilenameUCD,"",bdequote=.true.)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initParamTria (rparlist,rboundary,rtriangulation,rdiscrTime,&
      ssectionSpace,ssectionTime)
  
!<description>
  ! Reads in / create the basic domain and coarse mesh in space and time.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the spatial mesh/domain can be found.
  character(len=*), intent(in) :: ssectionSpace

  ! Section where the parameters of the time mesh can be found.
  character(len=*), intent(in) :: ssectionTime
!</input>

!<output>
  ! Description of the boundary
  type(t_boundary), intent(out) :: rboundary
  
  ! Spatial coarse mesh. The result is a raw mesh.
  type(t_triangulation), intent(out) :: rtriangulation
  
  ! Time coarse mesh.
  type(t_timeDiscretisation) :: rdiscrTime
!</output>

!</subroutine>

  ! local variables
  integer :: imeshType,ncellsX
  
    ! Variable for a filename:  
    character(LEN=60) :: sPRMFile, sTRIFile
    
    integer :: niterations
    real(dp) :: dtimeInit,dtimeMax,dtimeStepTheta
    integer :: ctimeStepScheme
    
    ! -----------------
    ! Space mesh/domain

    ! Get the .prm and the .tri file from the parameter list.
    ! note that parlst_getvalue_string returns us exactly what stands
    ! in the parameter file, so we have to apply READ to get rid of
    ! probable ""!
    call parlst_getvalue_string (rparlist,ssectionSpace,'sParametrisation',&
        sPRMFile,bdequote=.true.)
                              
    call parlst_getvalue_string (rparlist,ssectionSpace,'sMesh',&
        sTRIFile,bdequote=.true.)
    
    ! Read in the parametrisation of the boundary and save it to rboundary.
    call boundary_read_prm(rboundary, sPrmFile)
        
    ! Now set up the basic triangulation. Which type of triangulation should
    ! be set up?
    call parlst_getvalue_int (rparlist,ssectionSpace,'imeshType',imeshType,0)
    select case (imeshType)
    case (0)
      ! Standard mesh, specified by a TRI file.
      call tria_readTriFile2D (rtriangulation, sTRIFile, rboundary)
    case (1)
      ! Sliced QUAD mesh with ncellsX cells on the coarse grid.
      call parlst_getvalue_int (rparlist,ssectionSpace,'ncellsX',ncellsX,1)
      call cc_generateSlicedQuadMesh(rtriangulation,ncellsX)
    case default
      call output_line ('Unknown mesh type!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initParamTria')
      call sys_halt()
    end select
    
    ! Create a standard mesh.
    call tria_initStandardMeshFromRaw(rtriangulation, rboundary)
    
    ! --------------------------
    ! Time mesh / discretisation
    
    ! Get the parameters
    call parlst_getvalue_int (rparlist,ssectionTime,'niterations', &
        niterations, 1000)
    call parlst_getvalue_double (rparlist,ssectionTime,'dtimestart', &
        dtimeInit, 0.0_DP)
    call parlst_getvalue_double (rparlist,ssectionTime,'dtimemax', &
        dtimeMax, 20.0_DP)
    call parlst_getvalue_int (rparlist, ssectionTime,'ctimeStepScheme', &
        ctimeStepScheme, 0)
    call parlst_getvalue_double (rparlist,ssectionTime,'dtimeStepTheta', &
        dtimeStepTheta, 1.0_DP)
    
    ! Initialise the coarse time discretisation
    select case (ctimeStepScheme)
    case (0)
      ! One-step scheme.
      call tdiscr_initOneStepTheta (rdiscrTime, &
          dtimeInit, dtimeMax, niterations, dtimeStepTheta)
    case (1)
      ! FS-Theta scheme.
      call tdiscr_initFSTheta (rdiscrTime, dtimeInit, dtimeMax, niterations)

    case (2)
      ! dG(0)
      call tdiscr_initdG0 (rdiscrTime, dtimeInit, dtimeMax, niterations)
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_doneParamTria (rboundary,rtriangulation,rdiscrTime)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! Description of the boundary
  type(t_boundary), intent(inout) :: rboundary
  
  ! Coarse mesh. The result is a raw mesh.
  type(t_triangulation), intent(inout) :: rtriangulation

  ! Time coarse mesh.
  type(t_timeDiscretisation) :: rdiscrTime
!</inputoutput>

!</subroutine>

    ! Release the data.
    call tria_done (rtriangulation)
    call boundary_release (rboundary)
    call tdiscr_done(rdiscrTime)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_getRefinementParams (rparlist,rrefinementSpace,rrefinementTime,&
      ssectionSpace,ssectionTime)
  
!<description>
  ! Reads in parameters that define the refinement.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the spatial mesh/domain can be found.
  character(len=*), intent(in) :: ssectionSpace

  ! Section where the parameters of the time mesh can be found.
  character(len=*), intent(in) :: ssectionTime
!</input>

!<output>
  ! Description of the refinement in space
  type(t_settings_refinement), intent(out) :: rrefinementSpace

  ! Description of the refinement in time
  type(t_settings_refinement), intent(out) :: rrefinementTime
!</output>

!</subroutine>

    ! local variables
    integer :: nlmin,nlmax

    ! Get the level numbers.
    call parlst_getvalue_int (rparlist,ssectionSpace,'NLMIN',nlmin,1)
    call parlst_getvalue_int (rparlist,ssectionSpace,'NLMAX',nlmax,2)

    ! Do some correction to the level numbers
    if (nlmin .le. 0) nlmin = nlmax+nlmin
    nlmin = min(nlmin,nlmax)
    nlmax = max(nlmin,nlmax)
    
    ! Calculate the total number of levels
    rrefinementSpace%crefType = 0
    rrefinementSpace%nlevels = nlmax-nlmin+1
    rrefinementSpace%npreref = nlmin-1

    ! The same for the time mesh.
    call parlst_getvalue_int (rparlist,ssectionTime,'TIMENLMIN',nlmin,1)
    call parlst_getvalue_int (rparlist,ssectionTime,'TIMENLMAX',nlmax,2)

    ! Do some correction to the level numbers
    if (nlmin .le. 0) nlmin = nlmax+nlmin
    nlmin = min(nlmin,nlmax)
    nlmax = max(nlmin,nlmax)
    
    ! Calculate the total number of levels
    rrefinementTime%crefType = 0
    rrefinementTime%nlevels = nlmax-nlmin+1
    rrefinementTime%npreref = nlmin-1

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceHierarchy (rboundary,rrefinement,&
      rtriaCoarse,rmeshHierarchy,ioutputLevel)
  
!<description>
  ! Initialises a space hierarchy based on a coarse mesh and a refinement
  ! strategy specified in rrefinement.
!</description>

!<input>
  ! Type of refinement
  type(t_settings_refinement), intent(in) :: rrefinement

  ! Description of the boundary
  type(t_boundary), intent(in) :: rboundary
  
  ! Spatial coarse mesh. The result is a raw mesh.
  type(t_triangulation), intent(in), target :: rtriaCoarse
  
  ! Output level during initialisation
  integer, intent(in) :: ioutputLevel
!</input>

!<output>
  ! Refinement specification.
  type(t_meshHierarchy), intent(out) :: rmeshHierarchy
!</output>

!</subroutine>

    ! Create the hierarchy.
    if (ioutputLevel .ge. 1) then
      call output_line ("Pre-refinement to level "//&
          trim(sys_siL(rrefinement%npreref+1,10))//".",bnolinebreak=.true.)
    end if

    call mshh_initHierarchy (rmeshHierarchy,rtriaCoarse,&
        rrefinement%npreref,rrefinement%nlevels,rboundary)
        
    if (ioutputLevel .ge. 1) then
      call output_line (" Done. Creating Hierarchy-Level: [1",bnolinebreak=.true.)
    end if

    ! Refine the coarse mesh.
    call mshh_refineHierarchy2lv (rmeshHierarchy,rrefinement%nlevels,&
        rboundary=rboundary,bprint=ioutputLevel .ge. 1)

    if (ioutputLevel .ge. 1) then
      call output_line ("]",bnolinebreak=.true.)
    end if

    if (ioutputLevel .ge. 2) then
      call output_lbrk ()
      call output_line ('Mesh hierarchy statistics:')
      call output_lbrk ()
      call mshh_printHierStatistics (rmeshHierarchy)
    end if
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initTimeHierarchy (rrefinement,rtimeCoarse,rtimeHierarchy,ioutputLevel)
  
!<description>
  ! Initialises a space hierarchy based on a coarse mesh and a refinement
  ! strategy specified in rrefinement.
!</description>

!<input>
  ! Type of refinement
  type(t_settings_refinement), intent(in) :: rrefinement

  ! Spatial coarse mesh. The result is a raw mesh.
  type(t_timeDiscretisation), intent(in), target :: rtimeCoarse

  ! Output level during initialisation
  integer, intent(in) :: ioutputLevel
!</input>

!<output>
  ! Hierarchy of meshes in space.
  type(t_timescaleHierarchy), intent(out) :: rtimeHierarchy
!</output>

!</subroutine>

    if (ioutputLevel .ge. 1) then
      call output_line ("Creating Level: [1",bnolinebreak=.true.)
    end if

    ! Create the hierarchy.
    call tmsh_createHierarchy (rtimeCoarse,rtimeHierarchy,&
        rrefinement%npreref,rrefinement%nlevels,ioutputLevel .ge. 1)

    if (ioutputLevel .ge. 1) then
      call output_line ("]",bnolinebreak=.true.)
    end if

    if (ioutputLevel .ge. 2) then
      call output_lbrk ()
      call output_line ('Time hierarchy statistics:')
      call output_lbrk ()
      call tmsh_printHierStatistics (rtimeHierarchy)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceDiscrHier (rparlist,rsettingsSpaceDiscr,rboundary,&
      rmeshHierarchy,rsettings,ioutputLevel)

!<description>
  ! Initialises the hierarchies for the space discretisation on all levels
  ! based on the parameters in the parameter list.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Structure with space discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr

  ! Description of the boundary
  type(t_boundary), intent(in) :: rboundary

  ! A mesh hierarchy with all available space meshes.
  type(t_meshHierarchy), intent(in) :: rmeshHierarchy

  ! Output level during initialisation
  integer, intent(in) :: ioutputLevel
!</input>

!<output>
  ! Settings structure for the optimal control solver.
  ! The hierarchy substructures are created.
  type(t_settings_optflow), intent(inout) :: rsettings
!</output>

!</subroutine>
  
    ! local variables
    type(t_collection) :: rcollection
    integer :: icubM
    integer :: i,j,k
    integer(I32) :: celement
    
    ! Read the parameters that define the underlying discretisation.
    ! We use fget1LevelDiscretisation to create the basic spaces.
    ! This routines expects the rcollection%IquickAccess array to be initialised
    ! as follows:
    !   rcollection%IquickAccess(1) = ieltype
    !   rcollection%IquickAccess(2) = nequations (=3:primal space. =6: primal+dual space)
    !   rcollection%IquickAccess(3) = SPDISC_CUB_AUTOMATIC or icubA
    !   rcollection%IquickAccess(4) = SPDISC_CUB_AUTOMATIC or icubB
    !   rcollection%IquickAccess(5) = SPDISC_CUB_AUTOMATIC or icubF
    !
    ! Get parameters about the discretisation from the parameter list
    rcollection%IquickAccess(3) = rsettingsSpaceDiscr%icubStokes
    rcollection%IquickAccess(4) = rsettingsSpaceDiscr%icubB
    rcollection%IquickAccess(5) = rsettingsSpaceDiscr%icubF

    ! We want primal+dual space.
    rcollection%IquickAccess(2) = 6
    
    ! Element type
    rcollection%IquickAccess(1) = rsettingsSpaceDiscr%ielementType

    ! Create an FE space hierarchy based on the existing mesh hierarchy.
    call fesph_createHierarchy (rsettings%rfeHierPrimalDual,&
        rmeshHierarchy%nlevels,rmeshHierarchy,&
        fget1LevelDiscretisation,rcollection,rboundary)
        
    ! Extract the data of the primal space (Block 1..3) and create a new FE 
    ! space hierarchy based on that.
    call fesph_deriveFeHierarchy (rsettings%rfeHierPrimalDual,&
        rsettings%rfeHierPrimal,1,3)
        
    ! Extract the data of the primal velocity space (Block 1..2) and
    ! create a new FE space hierarchy for mass matrices
    call fesph_deriveFeHierarchy (rsettings%rfeHierPrimalDual,&
        rsettings%rfeHierMass,1,3)

    ! Which cubature rule to use for mass matrices?
    icubM = rsettingsSpaceDiscr%icubMass
    
    ! Set the cubature rule in the mass matrix discretisation structure.
    if (icubM .ne. SPDISC_CUB_AUTOMATIC) then
      do i=1,size(rsettings%rfeHierMass%p_rfeSpaces)
        do j=1,size(rsettings%rfeHierMass%p_rfeSpaces(i)%p_rdiscretisation%RspatialDiscr)
          rsettings%rfeHierMass%p_rfeSpaces(i)%p_rdiscretisation%RspatialDiscr(j)%&
              RelementDistr(:)%ccubTypeBilForm = icubM
        end do
      end do
    else
      do i=1,size(rsettings%rfeHierMass%p_rfeSpaces)
        do j=1,size(rsettings%rfeHierMass%p_rfeSpaces(i)%p_rdiscretisation%RspatialDiscr)
          do k=1,size(rsettings%rfeHierMass%p_rfeSpaces(i)%p_rdiscretisation%&
              RspatialDiscr(j)%RelementDistr)
            celement = rsettings%rfeHierMass%p_rfeSpaces(i)%p_rdiscretisation%&
                RspatialDiscr(j)%RelementDistr(k)%celement
          
            rsettings%rfeHierMass%p_rfeSpaces(i)%p_rdiscretisation%RspatialDiscr(j)%&
                RelementDistr(:)%ccubTypeBilForm = spdiscr_getStdCubature(celement)
          end do
        end do
      end do
    end if
    
    if (ioutputLevel .ge. 2) then
      ! Print statistics about the discretisation
      call output_line ("Space discretisation hierarchy statistics:")
      call output_lbrk ()
      call fesph_printHierStatistics (rsettings%rfeHierPrimalDual)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_doneSpaceDiscrHier (rsettings)
  
!<description>
  ! Cleans up the discretisation hierarchies in rsettings.
!</description>

!<output>
  ! Settings structure for the optimal control solver.
  ! The hierarchy substructures are created.
  type(t_settings_optflow), intent(inout) :: rsettings
!</output>

!</subroutine>

    ! Release all discretisation hierarchies.
    call fesph_releaseHierarchy(rsettings%rfeHierMass)
    call fesph_releaseHierarchy(rsettings%rfeHierPrimal)
    call fesph_releaseHierarchy(rsettings%rfeHierPrimalDual)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpacePrjHierarchy (rprjHierarchy,rfehierarchy,rparlist,ssection)
  
!<description>
  ! Creates projection hierarchies for the interlevel projection in space.
!</description>
  
!<input>
  ! Underlying FEM hierarchy.
  type(t_feHierarchy), intent(in) :: rfehierarchy
  
  ! Parameter list with parameters about the projection.
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where parameters about the projection can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! The projection hierarchy to create.
  type(t_interlevelProjectionHier), intent(out) :: rprjHierarchy
!</output>

!</subroutine>

    ! local variables
    integer :: ilev

    ! Initialise a standard interlevel projection structure for every level
    call mlprj_initPrjHierarchy(rprjHierarchy,1,rfehierarchy%nlevels)
        
    do ilev=1,rfehierarchy%nlevels
    
      call mlprj_initPrjHierarchyLevel(rprjHierarchy,ilev,&
          rfeHierarchy%p_rfeSpaces(ilev)%p_rdiscretisation)

      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      if (ilev .gt. 1) then
        call getProlRest (rprjHierarchy%p_Rprojection(ilev), &
            rparlist, ssection)
      end if
          
    end do
    
    call mlprj_commitPrjHierarchy(rprjHierarchy)

  contains
  
    subroutine getProlRest (rprojection, rparamList, sname)
    
    ! Initialises an existing interlevel projection structure rprojection
    ! with parameters from the INI/DAT files. sname is the section in the
    ! parameter list containing parameters about prolongation restriction.

    ! Parameter list that contains the parameters from the INI/DAT file(s).
    type(t_parlist), intent(IN) :: rparamList
    
    ! Name of the section in the parameter list containing the parameters
    ! of the prolongation/restriction.
    character(LEN=*), intent(IN) :: sname

    ! An interlevel projection block structure containing an initial
    ! configuration of prolongation/restriction. The structure is modified
    ! according to the parameters in the INI/DAT file(s).
    type(t_interlevelProjectionBlock), intent(INOUT) :: rprojection

      ! local variables
      type(t_parlstSection), pointer :: p_rsection
      integer :: i1
      real(DP) :: d1

      ! Check that there is a section called sname - otherwise we
      ! cannot create anything!
      
      call parlst_querysection(rparamList, sname, p_rsection) 

      if (.not. associated(p_rsection)) then
        ! We use the default configuration; stop here.
        return
      end if
      
      ! Now take a look which parameters appear in that section.

      ! Prolongation/restriction order for velocity components
      call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
      
      if (i1 .ne. -1) then
        ! Initialise order of prolongation/restriction for velocity components
        rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
        rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
        rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
      end if

      ! Prolongation/restriction order for pressure
      call parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
      
      if (i1 .ne. -1) then
        ! Initialise order of prolongation/restriction for pressure components
        rprojection%RscalarProjection(:,NDIM2D+1)%iprolongationOrder  = i1
        rprojection%RscalarProjection(:,NDIM2D+1)%irestrictionOrder   = i1
        rprojection%RscalarProjection(:,NDIM2D+1)%iinterpolationOrder = i1
      end if
      
      ! Prolongation/restriction variant for velocity components
      ! in case of Q1~ discretisation
      call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
      
      if (i1 .ne. -1) then
        rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
        rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1
      end if
      
      ! Aspect-ratio indicator in case of Q1~ discretisation
      ! with extended prolongation/restriction
      call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
      
      if (i1 .ne. 1) then
        rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
        rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
      end if

      ! Aspect-ratio bound for switching to constant prolongation/restriction
      ! in case of Q1~ discretisation with extended prolongation/restriction
      call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
      
      if (d1 .ne. 20.0_DP) then
        rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
        rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
      end if

    end subroutine

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceTimePrjHierarchy (rprjHierarchy,rtimeDisr,rhierarchy,&
      rprojHierarchySpace,rparlist,ssection)
  
!<description>
  ! Creates projection hierarchies for the interlevel projection in space.
!</description>
  
!<input>
  ! Underlying time (coarse) grid.
  type(t_timeDiscretisation), intent(in) :: rtimeDisr

  ! Underlying space-time hierarchy.
  type(t_spaceTimeHierarchy), intent(in) :: rhierarchy
  
  ! Projection hierarchy in space
  type(t_interlevelProjectionHier), intent(in), target :: rprojHierarchySpace
  
  ! Parameter list with parameters about the projection.
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where parameters about the projection can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! The projection hierarchy to create.
  type(t_sptiProjHierarchy), intent(out) :: rprjHierarchy
!</output>

!</subroutine>

    ! local variables
    integer :: iorderTimeProlRest

    ! Type of prolongation/restriction in time
    call parlst_getvalue_int (rparlist, ssection, &
        'iorderTimeProlRest', iorderTimeProlRest, -1)
        
    if (iorderTimeProlRest .eq. -1) then
      ! Automatic mode. Check the time stepping theta to determine
      ! if we have 1st or 2nd order in time.
      iorderTimeProlRest = 1
      if (rtimeDisr%dtheta .eq. 0.5_DP) then
        ! 2nd order Crank Nicolson
        iorderTimeProlRest = 2
      end if
      if (rtimeDisr%ctype .eq. TDISCR_FSTHETA) then
        ! 2nd order FS-Theta
        iorderTimeProlRest = 2
      end if
    end if
    
    call sptipr_initProjection (rprjHierarchy,rhierarchy,&
        rprojHierarchySpace,iorderTimeProlRest)
    
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

    ! Alpha/Gamma parameters
    call parlst_getvalue_double (rparlist,ssectionOptC,&
        'dalphaC',roptcontrol%dalphaC,1.0_DP)
        
    call parlst_getvalue_double (rparlist,ssectionOptC,&
        'dgammaC',roptcontrol%dgammaC,0.0_DP)
        
    ! Type of the formulation
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        'ispaceTimeFormulation',roptcontrol%ispaceTimeFormulation,0)
    
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        'iconvectionExplicit',roptcontrol%iconvectionExplicit,0)
    
    ! Parameters defining the constraints
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        'ccontrolConstraints',roptcontrol%rconstraints%ccontrolConstraints,0)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        'dumin1',roptcontrol%rconstraints%dumin1,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        'dumax1',roptcontrol%rconstraints%dumax1,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        'dumin2',roptcontrol%rconstraints%dumin2,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        'dumax2',roptcontrol%rconstraints%dumax2,1.0E10_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initOptControlTargetFlow (rparlist,ssectionOptC,rsettingsSpaceDiscr,&
      rtriaCoarse,rrefinementSpace,rfeHierPrimal,rtimeHierarchy,rboundary,roptcontrol)
  
!<description>
  ! Initialises the target flow for the optimal control problem
  ! based on the parameters in the DAT file.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the optimal control can be found.
  character(len=*), intent(in) :: ssectionOptC

  ! Structure with space discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr

  ! Underlying spatial coarse mesh of the problem.
  type(t_triangulation), intent(in) :: rtriaCoarse
  
  ! Definition of the underlying domain.
  type(t_boundary), intent(in) :: rboundary

  ! Description of the refinement of rtriaCoarse.
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target flow matches the one of the primary
  ! flow, this hierarchy is used to save memory.
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
    character(len=SYS_STRLEN) :: stargetFlow
    integer :: isuccess
    
    ! Initialise the target flow.
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        'stargetFlow',stargetFlow,bdequote=.true.)
    call init_initFlow (rparlist,stargetFlow,roptcontrol%rtargetFlow,&
        rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,isuccess)
    if (isuccess .eq. 1) then
      call output_line ('Flow created by simulation not yet supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initOptControlTargetFlow')
      call sys_halt()
    end if


  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_doneOptControl (roptcontrol)
  
!<description>
  ! Cleans up information about constraints
!</description>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Release the target flow.
    call ansol_done(roptcontrol%rtargetFlow)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  recursive subroutine init_initDiscreteAnalytFlow (ielementType,rboundary,&
      smesh,rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,&
      ilevel,rflow)
  
!<description>
  ! Creates a analytical-flow structure that resembles a discrete flow.
!</description>

!<input>
  ! Element to use.
  ! =-1: use the default element in rsettingsSpaceDiscr.
  integer, intent(in) :: ielementType

  ! Definition of the domain
  type(t_boundary), intent(in), target :: rboundary

  ! Underlying mesh. May be ="", in this case the default mesh in rtriaCoarse is used.
  character(len=*), intent(in) :: smesh
  
  ! Default mesh.
  type(t_triangulation), intent(in) :: rtriaCoarse

  ! Description of the refinement in space, how rfeHierPrimal was created
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! Space discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target flow matches the one of the primary
  ! flow, this hierarchy is used to save memory.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal
  
  ! Refinement level, the analytical flow should resemble.
  integer, intent(in) :: ilevel
!</input>

!<output>
  ! Flow structure to create.
  type(t_anSolution), intent(out) :: rflow
!</output>

!</subroutine>

    ! local variables
    integer :: ieltype,iavaillevel
    type(t_collection) :: rcollection

    ! Can we reuse our hierarchy?
    if (smesh .eq. "") then
    
      ! Is the refinement level smaller than NLMIN?
      if (ilevel .lt. rrefinementSpace%npreref+1) then
        ! We have to create that level.
        !
        ! Set up the collection for the creation of appropriate discretisation structures
        rcollection%IquickAccess(1) = ielementType
        rcollection%IquickAccess(2) = NDIM2D+1
        rcollection%IquickAccess(3) = SPDISC_CUB_AUTOMATIC
        rcollection%IquickAccess(4) = SPDISC_CUB_AUTOMATIC
        rcollection%IquickAccess(5) = SPDISC_CUB_AUTOMATIC
        
        if (ielementType .eq. -1) then
          ! Get the element type from the discretisation section
          rcollection%IquickAccess(1) = rsettingsSpaceDiscr%ielementType
        end if
        
        ! Create the basic analytic solution: Mesh, discretisation structure,...
        call ansol_init (rflow,ilevel,&
            rtriaCoarse,1,rcollection%IquickAccess(1),&
            fget1LevelDiscretisation,rcollection,rboundary)
      else
      
        ! Ok, here we can hope to reuse existing data structures.
        !
        ! Set up the collection for the creation of appropriate discretisation structures
        rcollection%IquickAccess(2) = NDIM2D+1
        rcollection%IquickAccess(3) = SPDISC_CUB_AUTOMATIC
        rcollection%IquickAccess(4) = SPDISC_CUB_AUTOMATIC
        rcollection%IquickAccess(5) = SPDISC_CUB_AUTOMATIC

        ! Get the element type from the discretisation
        ieltype = rsettingsSpaceDiscr%ielementType
        
        if ((ielementType .eq. -1) .or. &
            (ielementType .eq. ieltype)) then
          
          ! Very nice, the element type matches.
          rcollection%IquickAccess(1) = ieltype

          ! Get the maximum available level in rfeHierPrimal
          iavaillevel = min(rfeHierPrimal%nlevels,ilevel-rrefinementSpace%npreref)
          
          ! And now create the basic flow. Only in case ilevel>NLMAX,
          ! new levels are generated.
          call ansol_init (rflow,ilevel-rrefinementSpace%npreref,&
              rfeHierPrimal%p_rfeSpaces(iavaillevel)%p_rdiscretisation,iavaillevel,&
              rcollection%IquickAccess(1),fget1LevelDiscretisation,rcollection)
          
        else
        
          ! Ok, a little bit different, the element type is different.
          rcollection%IquickAccess(1) = ielementType
          
          ! That means, we can reuse the triangulation but not the discretisation.
          !
          ! Get the maximum available level in rfeHierPrimal
          iavaillevel = min(rfeHierPrimal%nlevels,ilevel-rrefinementSpace%npreref)
          
          ! And now create the basic flow. Only in case ilevel>NLMAX,
          ! new levels are generated.
          call ansol_init (rflow,ilevel-rrefinementSpace%npreref,&
              rfeHierPrimal%rmeshHierarchy%p_Rtriangulations(iavaillevel),iavaillevel,&
              rcollection%IquickAccess(1),fget1LevelDiscretisation,rcollection)
          
        end if
        
      end if
    
    else
    
      ! Mesh is different. Then we have to do everything by hand...
      !
      ! Set up the collection for the creation of appropriate discretisation structures
      rcollection%IquickAccess(1) = ielementType
      rcollection%IquickAccess(2) = NDIM2D+1
      rcollection%IquickAccess(3) = SPDISC_CUB_AUTOMATIC
      rcollection%IquickAccess(4) = SPDISC_CUB_AUTOMATIC
      rcollection%IquickAccess(5) = SPDISC_CUB_AUTOMATIC
      
      if (ielementType .eq. -1) then
        ! Get the element type from the discretisation
        rcollection%IquickAccess(1) = rsettingsSpaceDiscr%ielementType
      end if
      
      ! Create the basic analytic solution: Mesh, discretisation structure,...
      ! As we pass the mesh file here, we also pass the original ilevel
      ! because the mesh has to be refined up to that.
      call ansol_init (rflow,ilevel,&
          NDIM2D,smesh,rcollection%IquickAccess(1),&
          fget1LevelDiscretisation,rcollection,rboundary)
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initFlow (rparlist,ssectionFlow,rflow,&
      rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,rboundary,&
      isuccess)
  
!<description>
  ! Reads in and sets up a flow. ssectionFlow is the name of a section
  ! containing all parameters that configure the flow. The routine creates
  ! a structure rflow based on these parameters.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found that specify the target flow.
  character(len=*), intent(in) :: ssectionFlow

  ! Space discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr

  ! Coarse mesh, corresponding to rfeHierPrimal.
  type(t_triangulation), intent(in) :: rtriaCoarse

  ! Description of the refinement in space, how rfeHierPrimal was created
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! A hierarchy of space levels for velocity+pressure (primal/dual space).
  ! If the element of the target flow matches the one of the primary
  ! flow, this hierarchy is used to save memory.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal

  ! Definition of the domain.
  type(t_boundary), intent(in), target :: rboundary
!</input>

!<output>
  ! Flow structure to create.
  type(t_anSolution), intent(out) :: rflow
  
  ! Success-flag. Returns whether this routine successfully created the flow.
  ! =0: Flow was successfully created.
  ! =1: Flow is defined as forward simulation and was not possible to create.
  !     The caller must invoke init_initFlowBySimulation to create the flow.
  integer, intent(out) :: isuccess
  
!</output>

!</subroutine>

    ! local variables
    integer :: ctype,iid
    real(dp) :: dstartTime,dtimeMax
    character(len=PARLST_LENLINEBUF), dimension(6) :: Sexpressions
    character(len=SYS_STRLEN) :: smesh,sflowFile
    integer :: ilevel,ielementType,idelta
    integer :: ntimesteps
    type(t_parlstSection), pointer :: p_rsection
    
    ! Is the section available?
    call parlst_querysection(rparlist, ssectionFlow, p_rsection) 
    if (.not. associated(p_rsection)) then
      call output_line ('Flow section does not exist: '//trim(ssectionFlow), &
          OU_CLASS_ERROR,OU_MODE_STD,'init_initFlow')
      call sys_halt()
    end if
    
    ! Type of the flow?    
    call parlst_getvalue_int (rparlist,ssectionFlow,'ctype',ctype,-1)
        
    if (ctype .eq. 4) then
      ! Forget it, there are not enough parameters here to create the flow.
      ! This situation is more complicated!
      isuccess = 1
      return
    end if

    ! Read some more parameters, we may need them.
    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'sflowExpressionY1',Sexpressions(1),"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'sflowExpressionY2',Sexpressions(2),"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'sflowExpressionP',Sexpressions(3),"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'sflowExpressionL1',Sexpressions(4),"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'sflowExpressionL2',Sexpressions(5),"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'sflowExpressionXi',Sexpressions(6),"",bdequote=.true.)
    
    ! Id, level, element type

    call parlst_getvalue_int (rparlist,ssectionFlow,&
        'iid',iid,0)

    call parlst_getvalue_int (rparlist,ssectionFlow,&
        'ilevel',ilevel,0)

    call parlst_getvalue_int (rparlist,ssectionFlow,&
        'ielementType',ielementType,-1)

    ! Mesh, the flow diles

    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'smesh',smesh,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssectionFlow,&
        'sflowFile',sflowFile,"",bdequote=.true.)

    ! Time discretisation

    call parlst_getvalue_int (rparlist,ssectionFlow,&
        'idelta',idelta,1)

    ! Set it up...
    select case (ctype)
    case (-1)
      ! Analytic flow given by a callback routine
      call ansol_init(rflow,6,0)
      rflow%iid = iid
      return
    
    case (0)
      ! Zero flow
      call ansol_init(rflow,6)
      rflow%iid = iid
      return

    case (3)
      ! Analytical expression as flow.
      call ansol_init (rflow,6)
      call ansol_configExpressions (rflow,Sexpressions)
      rflow%iid = iid
      return

    end select
    
    ! Create an analytical flow that resembles a discrete flow.
    call init_initDiscreteAnalytFlow (ielementType,rboundary,&
        smesh,rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rfeHierPrimal,&
        ilevel,rflow)  
        
    rflow%iid = iid  
    
    ! Finally, set up the flow.
    select case (ctype)
      
    case (1)
    
      ! Read stationary flow from hard disc
      call ansol_configStationaryFile (rflow,sflowfile,.true.)
      
    case (2)
    
      ! Read nonstationary flow from hard disc
      
      call parlst_getvalue_double (rparlist,ssectionFlow,&
          'dstartTime',dstartTime)

      call parlst_getvalue_double (rparlist,ssectionFlow,&
          'dtimeMax',dtimeMax)
      
      call parlst_getvalue_int (rparlist,ssectionFlow,&
          'ntimesteps',ntimesteps)

      call ansol_configNonstationaryFile (rflow, &
          dstartTime,dtimeMax,ntimesteps,&
          "("""//trim(sflowfile)//"."",I5.5)",&
          0,idelta,.true.)
          
    end select

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  recursive subroutine init_initFlowBySimulation (rparlist,ssectionFlow,rflow,&
!      rsettingsSpaceDiscr,rsettingsSolver)
!  
!!<description>
!  ! Sets up a flow by performing a forward simulation.
!!</description>
!
!!<input>
!  ! Parameter list
!  type(t_parlist), intent(in) :: rparlist
!  
!  ! Section where the parameters can be found that specify the target flow.
!  character(len=*), intent(in) :: ssectionFlow
!
!  ! Space discretisation settings
!  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr
!
!  ! Settings structure for the optimal control solver.
!  type(t_settings_optflow), intent(out), target :: rsettingsSolver
!!</input>
!
!!<output>
!  ! Flow structure to create.
!  type(t_anSolution), intent(out) :: rflow
!!</output>
!
!!</subroutine>
!
!    ! local variables
!    integer :: ctype,iavaillevel,iid,isuccess
!    type(t_collection) :: rcollection
!    real(dp) :: dstartTime,dtimeMax,dtimeStepTheta
!    character(len=SYS_STRLEN) :: ssectionInitCond
!    character(len=SYS_STRLEN) :: smesh,sflowFile
!    integer :: ilevel,ielementType,idelta,ieltype
!    integer :: ntimesteps,ctimeStepScheme
!    type(t_anSolution) :: rinitcondflow
!    
!    ! Read some parameters, we may need them.
!    
!    ! Id, level, element type
!
!    call parlst_getvalue_int (rparlist,ssectionFlow,&
!        'iid',iid,0)
!
!    call parlst_getvalue_int (rparlist,ssectionFlow,&
!        'ilevel',ilevel,0)
!
!    call parlst_getvalue_int (rparlist,ssectionFlow,&
!        'ielementType',ielementType,-1)
!
!    ! Mesh, initial condition
!
!    call parlst_getvalue_string (rparlist,ssectionFlow,&
!        'smesh',smesh,"",bdequote=.true.)
!
!    call parlst_getvalue_string (rparlist,ssectionFlow,&
!        'ssectionInitSol',ssectionInitSol,bdequote=.true.)
!
!    ! Time discretisation
!
!    call parlst_getvalue_double (rparlist,ssectionFlow,&
!        'dstartTime',dstartTime)
!
!    call parlst_getvalue_double (rparlist,ssectionFlow,&
!        'dtimeMax',dtimeMax)
!    
!    call parlst_getvalue_int (rparlist,ssectionFlow,&
!        'ntimesteps',ntimesteps)
!
!    call parlst_getvalue_int (rparlist,ssectionFlow,&
!        'ctimeStepScheme',ctimeStepScheme)
!
!    call parlst_getvalue_int (rparlist,ssectionFlow,&
!        'dtimeStepTheta',dtimeStepTheta)
!
!    ! Prepare the flow
!    call init_initDiscreteAnalytFlow (ielementType,rboundary,&
!        smesh,rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,&
!        rsettingsSolver%rfeHierPrimal,ilevel,rflow)
!
!    rflow%iid = iid
!    
!    ! What is with the initial condition?  
!    call init_initFlow (rparlist,ssectionInitSol,rinitcondflow,&
!        rtriaCoarse,rrefinementSpace,rsettingsSpaceDiscr,rsettingsSolver%rfeHierPrimal,
!        rboundary,isuccess)
!    
!    if (isuccess .eq. 1) then
!      ! Oops, that one must be calculated by simulation.
!      call init_initFlowBySimulation (rparlist,ssectionInitSol,rinitcondflow,&
!          rsettingsSpaceDiscr,rsettingsSolver)
!    end if
!    
!    ! Solution at the start time in rinitcondflow is our initial condition.
!    
!    
!    
!    
!    
!    ! Release the initial condition
!    call ansol_done (rinitcondflow)
!
!  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initSpaceTimeHierarchy (rparlist,ssection,&
      rrefinementSpace,rrefinementTime,&
      rfeHierPrimal,rfeHierPrimalDual,rtimeHierarchy,&
      rspaceTimeHierPrimal,rspaceTimeHierPrimalDual)
  
!<description>
  ! Creates a space-time hierarchy based on the parameters in the parameter list.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the parameters how to set up the
  ! space-time hierarchies.
  character(len=*), intent(in) :: ssection

  ! Settings that define the refinement in space
  type(t_settings_refinement), intent(in) :: rrefinementSpace

  ! Settings that define the refinement in time
  type(t_settings_refinement), intent(in) :: rrefinementTime

  ! A mesh hierarchy with all available space meshes, only primal space.
  type(t_feHierarchy), intent(in) :: rfeHierPrimal
  
  ! A mesh hierarchy with all available space meshes, primal + dual space.
  type(t_feHierarchy), intent(in) :: rfeHierPrimalDual

  ! A hierarchy of time levels
  type(t_timescaleHierarchy), intent(in) :: rtimeHierarchy 
!</input>

!<inputoutput>
  ! A space-time hierarchy based on the primal/dual space
  type(t_spaceTimeHierarchy), intent(out) :: rspaceTimeHierPrimal
  
  ! A space-time hierarchy based on the primal+dual space
  type(t_spaceTimeHierarchy), intent(out) :: rspaceTimeHierPrimalDual
!</inputoutput>

!</subroutine>

    integer :: ispacelevelcoupledtotimelevel,nmaxSimulRefLevel
    real(DP) :: dspacetimeRefFactor
    
    ! Get the parameters from the parameter list.
    !
    ! Should we couple space and time coarsening/refinement?
    call parlst_getvalue_int (rparlist,ssection,&
        'ispacelevelcoupledtotimelevel',ispacelevelcoupledtotimelevel,1)

    ! Parameters of the refinement?
    call parlst_getvalue_int (rparlist,ssection,&
        'nmaxSimulRefLevel',nmaxSimulRefLevel,0)
    
    if (nmaxSimulRefLevel .le. 0) then
      nmaxSimulRefLevel = rrefinementTime%nlevels+&
          nmaxSimulRefLevel-rrefinementTime%npreref
    end if
    nmaxSimulRefLevel = min(rrefinementTime%nlevels,max(1,nmaxSimulRefLevel))
        
    call parlst_getvalue_double (rparlist,ssection,&
        'dspacetimeRefFactor',dspacetimeRefFactor,1.0_DP)

    ! Create the hierarchies.
    call sth_initHierarchy (rspaceTimeHierPrimal,&
        rfeHierPrimal,rtimeHierarchy)

    call sth_initHierarchy (rspaceTimeHierPrimalDual,&
        rfeHierPrimalDual,rtimeHierarchy)
        
    select case (ispacelevelcoupledtotimelevel)
    case (0)
      ! Only in time, space level stays at max.

      call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimal,&
          rrefinementSpace%nlevels,rrefinementSpace%nlevels,&
          1,rrefinementTime%nlevels,dspacetimeRefFactor)

      call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimalDual,&
          rrefinementSpace%nlevels,rrefinementSpace%nlevels,&
          1,rrefinementTime%nlevels,dspacetimeRefFactor)
          
    case (1)
      ! Simultaneous refinement in space+time

      call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimal,&
          1,rrefinementSpace%nlevels,&
          1,rrefinementTime%nlevels,dspacetimeRefFactor)

      call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimalDual,&
          1,rrefinementSpace%nlevels,&
          1,rrefinementTime%nlevels,dspacetimeRefFactor)

    case (2)
      ! Only in space, time level stays at max.

      call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimal,&
          1,rrefinementSpace%nlevels,&
          rrefinementTime%nlevels,rrefinementTime%nlevels,dspacetimeRefFactor)

      call sth_defineHierarchyByCoarsening (rspaceTimeHierPrimalDual,&
          1,rrefinementSpace%nlevels,&
          rrefinementTime%nlevels,rrefinementTime%nlevels,dspacetimeRefFactor)

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateSlicedQuadMesh (rtriangulation,ncellsX)
  
!<description>
  ! This routine generates a standard [0,1]^2 QUAD mesh with ncellsX cells in
  ! X-direction. By nature, these cells show an anisotropy of ncellsX:1.
!</description>
  
!<input>
  ! Number of cells in X-direction.
  integer, intent(IN) :: ncellsX
!</input>
    
!<output>
  ! Triangulation structure that receives the triangulation.
  type(t_triangulation), intent(OUT) :: rtriangulation
!</output>
    
    ! local variables
    real(DP), dimension(:,:), pointer :: p_Ddata2D
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer(I32), dimension(:,:), pointer :: p_Idata2D
    integer(I32), dimension(:), pointer :: p_Idata,p_IverticesAtBoundary,p_IboundaryCpIdx
    integer(I32) :: ivt, iel
    integer(I32), dimension(2) :: Isize
    
    ! Initialise the basic mesh
    rtriangulation%ndim = NDIM2D
    rtriangulation%NEL = ncellsX
    rtriangulation%NVT = (ncellsX+1)*2
    rtriangulation%NMT = 0
    rtriangulation%NNVE = 4
    rtriangulation%NNEE = 4
    rtriangulation%NBCT = 1
    rtriangulation%InelOfType(:) = 0
    rtriangulation%InelOfType(TRIA_NVEQUAD2D) = rtriangulation%NEL
    
    ! Allocate memory for the basic arrays on the heap
    ! 2d array of size(NDIM2D, NVT)
    Isize = (/NDIM2D,int(rtriangulation%NVT,I32)/)
    call storage_new ('tria_read_tri2D', 'DCORVG', Isize, ST_DOUBLE, &
        rtriangulation%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointers to the coordinate array
    ! p_Ddata2D is the pointer to the coordinate array
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
    
    ! Initialise the point coordinates.
    ! Odd vertices on the bottom, even vertices on top of the QUAD mesh,
    ! numbered from left to right.
    do ivt=0,ncellsX
      p_Ddata2D(1,2*ivt+1) = real(ivt,DP)/real(ncellsX,DP)
      p_Ddata2D(2,2*ivt+1) = 0.0_DP

      p_Ddata2D(1,2*ivt+2) = real(ivt,DP)/real(ncellsX,DP)
      p_Ddata2D(2,2*ivt+2) = 1.0_DP
    end do
    
    ! Allocate memory for IverticesAtElement
    ! build the old KVERT...
    ! 2d array of size(NVE, NEL)
    Isize = (/rtriangulation%NNVE,int(rtriangulation%NEL,I32)/)
    call storage_new ('tria_read_tri2D', 'KVERT', Isize, ST_INT, &
        rtriangulation%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
        
    ! Get the pointer to the IverticesAtElement array and read the array
    call storage_getbase_int2D(&
        rtriangulation%h_IverticesAtElement,p_Idata2D)

    ! Initialise the connectivity for the cells.
    do iel=0,ncellsX-1
      p_Idata2D(1,iel+1) = 2*iel+1
      p_Idata2D(2,iel+1) = 2*iel+3
      p_Idata2D(3,iel+1) = 2*iel+4
      p_Idata2D(4,iel+1) = 2*iel+2
    end do
    
    ! Allocate memory for InodalProperty 
    call storage_new ('tria_read_tri2D', 'KNPR', &
        int(rtriangulation%NVT,I32), ST_INT, &
        rtriangulation%h_InodalProperty, ST_NEWBLOCK_ZERO)
    
    ! Get the pointer to the InodalProperty array
    call storage_getbase_int(&
        rtriangulation%h_InodalProperty,p_Idata)

    ! All vertices are on the boundary
    p_Idata(:) = 1
    
    ! Number of vertices on the boundary -- all of them
    rtriangulation%NVBD = rtriangulation%NVT
    
    ! Allocate memory for IverticesAtBoundary.
    call storage_new ('tria_generateBasicBoundary', &
        'KVBD', int(rtriangulation%NVBD,I32), &
        ST_INT, rtriangulation%h_IverticesAtBoundary, ST_NEWBLOCK_NOINIT)
        
    ! Allocate memory for the boundary component index vector.
    ! Initialise that with zero!
    call storage_new ('tria_generateBasicBoundary', &
        'KBCT', int(rtriangulation%NBCT+1,I32), &
        ST_INT, rtriangulation%h_IboundaryCpIdx, ST_NEWBLOCK_ZERO)
    
    ! Get pointers to the arrays
    call storage_getbase_int (&
        rtriangulation%h_IverticesAtBoundary,p_IverticesAtBoundary)
        
    call storage_getbase_int (&
        rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
    
    ! The first element in p_IboundaryCpIdx is (as the head) always =1.
    p_IboundaryCpIdx(1) = 1
    p_IboundaryCpIdx(2) = 1 + rtriangulation%NVBD

    ! Initialise the numbers of the vertices on the boundary
    do ivt=0,ncellsX
      p_IverticesAtBoundary (ivt+1) = 2*ivt+1
      p_IverticesAtBoundary (2*(ncellsX+1)-ivt) = 2*ivt+2
    end do
    
    ! Allocate memory for  and DvertexParameterValue
    call storage_new ('tria_generateBasicBoundary', &
        'DVBDP', int(rtriangulation%NVBD,I32), &
        ST_DOUBLE, rtriangulation%h_DvertexParameterValue, ST_NEWBLOCK_NOINIT)
    
    ! Get the array where to store boundary parameter values.
    call storage_getbase_double (&
        rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
    call storage_getbase_double2D(&
        rtriangulation%h_DvertexCoords,p_Ddata2D)
        
    ! Initialise the parameter values of the vertices. For the bottommost
    ! edge, they coincide wit the coordinate. For the topmost edge,
    ! that's 3 - x-coordinate.
    do ivt=1,ncellsX+1
      p_DvertexParameterValue (ivt) = p_Ddata2D(1,p_IverticesAtBoundary(ivt))
      p_DvertexParameterValue (ncellsX+1+ivt) = &
          3.0_DP - p_Ddata2D(1,p_IverticesAtBoundary(ncellsX+1+ivt))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_getDebugFlags (rparlist,ssection,rdebugFlags)
  
!<description>
  ! Reads all debug flags from the parameter list.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the parameters.
  character(len=*), intent(in) :: ssection
!</input>

!<output>
  ! Structure with the parameters of the nonlinear soace-time solver.
  type(t_optcDebugFlags), intent(out) :: rdebugFlags
!</output>

!</subroutine>

    call parlst_getvalue_double (rparlist,ssection,&
        'dprimalDualCoupling',rdebugFlags%dprimalDualCoupling,1.0_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        'ddualPrimalCoupling',rdebugFlags%ddualPrimalCoupling,1.0_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        'dterminalCondDecoupled',rdebugFlags%dterminalCondDecoupled,1.0_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        'dtimeCoupling',rdebugFlags%dtimeCoupling,1.0_DP)

    call parlst_getvalue_int (rparlist,ssection,&
        'ipressureFullyImplicit',rdebugFlags%ipressureFullyImplicit,1)

    call parlst_getvalue_int (rparlist,ssection,&
        'cumfpackWriteMatrix',rdebugFlags%cumfpackWriteMatrix,0)

    call parlst_getvalue_string (rparlist,ssection,&
        'sumfpackMatrixFilename',rdebugFlags%sumfpackMatrixFilename,&
        "./matrix.txt",bdequote=.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initNlSpaceTimeSolver (rsettings,rparlist,ssection,rdiscrTime,rnlstsolver)
  
!<description>
  ! Reads the parameters of the nonlinear space-time solver and initialises the
  ! corresponding structure.
!</description>

!<input>
  ! Settings structure for the optimal control solver.
  type(t_settings_optflow), intent(in), target :: rsettings

  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the parameters.
  character(len=*), intent(in) :: ssection

  ! Time coarse mesh.
  type(t_timeDiscretisation) :: rdiscrTime
!</input>

!<output>
  ! Structure with the parameters of the nonlinear soace-time solver.
  type(t_nlstsolver), intent(out) :: rnlstsolver
!</output>

!</subroutine>

    character(len=SYS_STRLEN) :: sstr,sstrSmoother,sstrSmootherPrec
    character(len=SYS_STRLEN) :: sstrSgrSolver,sstrSgrPrecond
    type(t_settings_nlstprec) :: rsettingsPrecond
    integer :: cpreconditioner

    ! Figure out which type of solver we should use to solve
    ! the problem. 
    call parlst_getvalue_int (rparlist,ssection,&
        'cpreconditioner',cpreconditioner,0)
        
    call parlst_getvalue_int (rparlist, ssection, &
        'ioutputLevel', rnlstsolver%ioutputLevel, 2)
                              
    ! Depending on the type of the solver, get the parameters of the subsolver(s).
    select case (cpreconditioner)
    case (0)
      ! Not implemented
      
    case (1)
      ! Simple single-grid solver.
      rsettingsPrecond%ctypeSolver = 0
      
      ! Get the section with the parameters
      call parlst_getvalue_string (rparlist,ssection,&
          'ssectionSingleGridSolver',sstrSgrSolver,&
          "TIME-SINGLEGRIDSOLVER",bdequote=.true.)
          
      call parlst_getvalue_string (rparlist,ssection,&
          'ssectionSinglePreconditioner',sstrSgrPrecond,&
          "TIME-SINGLEGRIDPRECOND",bdequote=.true.)
          
      ! Get the parameters from that subsection. This configures the solver.
      call init_getParamsPrec_sgrsolver (rparlist,sstrSgrSolver,&
          rsettingsPrecond%ronegridSolver)
      
      ! Get parameters for a possible preconditioner
      call init_getParamsPrec_sgrprec (rparlist,sstrSgrPrecond,&
          rsettingsPrecond%ronegridPrecond)

    case (2)
      ! Multigrid solver
      rsettingsPrecond%ctypeSolver = 1
    
      ! Get the section with the parameters
      call parlst_getvalue_string (rparlist,ssection,&
          'ssectionMultigrid',sstr,"TIME-MULTIGRID",bdequote=.true.)
          
      call parlst_getvalue_string (rparlist,ssection,&
          'ssectionSmoother',sstrSmoother,"TIME-SMOOTHER",bdequote=.true.)

      call parlst_getvalue_string (rparlist,ssection,&
          'ssectionSmootherPrecond',sstrSmootherPrec,&
          "TIME-SMOOTHERPRECOND",bdequote=.true.)
          
      call parlst_getvalue_string (rparlist,ssection,&
          'ssectionCgrSolver',sstrSgrSolver,&
          "TIME-SINGLEGRIDSOLVER",bdequote=.true.)
          
      call parlst_getvalue_string (rparlist,ssection,&
          'ssectionCgrPreconditioner',sstrSgrPrecond,&
          "TIME-SINGLEGRIDPRECOND",bdequote=.true.)
          
      ! Get the parameters from that subsections.
      call cc_getParamsPrec_mgsolver (rparlist,&
          rsettingsPrecond%rmgSolver,rdiscrTime,sstr)
      
      ! Read the parameters of the smoother on every level.
      call init_getParamsPrec_smoother (rparlist,sstrSmoother,&
          rsettingsPrecond%rsmoother)

      ! Read the parameters for the preconditioner of the smoother on every level.
      call init_getParamsPrec_sgrprec (rparlist,sstrSmootherPrec,&
          rsettingsPrecond%rprecsmoother)
      
      ! Go on reading the parameters of the coarse grid solver,
      call init_getParamsPrec_sgrsolver (rparlist,sstrSgrSolver,&
          rsettingsPrecond%ronegridSolver)
      
      ! and its preconditioner (in case it has one).
      call init_getParamsPrec_sgrprec (rparlist,sstrSgrPrecond,&
          rsettingsPrecond%ronegridPrecond)

    end select
    
    ! Using these preconditioner information, create a nonlinear solver structure.
    call stnlsinit_initSolver (rsettings,&
        rsettings%rspaceTimeHierPrimalDual%nlevels,rsettingsPrecond,rnlstsolver)
    
    ! Get parameters that configure the iteration.
    call parlst_getvalue_int (rparlist,ssection,&
        'ctypeNonlinearIteration',rnlstsolver%ctypeNonlinearIteration,1)

    call parlst_getvalue_int (rparlist, ssection, &
        'nminIterations', rnlstsolver%nminIterations, 1)
    
    call parlst_getvalue_int (rparlist, ssection, &
        'nmaxIterations', rnlstsolver%nmaxIterations, 10)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'depsRel', rnlstsolver%depsRel, 1E-5_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'depsAbs', rnlstsolver%depsAbs, 1E-5_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'depsDiff', rnlstsolver%depsDiff, 0.0_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'domega', rnlstsolver%domega, 1.0_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'dinexactNewtonEpsRel', rnlstsolver%dinexactNewtonEpsRel, 1.0E-2_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'dinexactNewtonExponent', rnlstsolver%dinexactNewtonExponent, 2.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getParamsPrec_mgsolver (rparlist,rnlstsubsolver,rdiscrTime,ssection)
  
!<description>
  ! Reads parameters that configure the linear space-time multigrid
  ! preconditioner.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the solver parameters.
  character(len=*), intent(in) :: ssection

  ! Time coarse mesh.
  type(t_timeDiscretisation) :: rdiscrTime
!</input>

!<output>
  ! Structure with the parameters to be initialised.
  type(t_nlstprec_mgrsolver), intent(out) :: rnlstsubsolver
!</output>

!</subroutine>

    ! Get the parameters of the solver.
    call parlst_getvalue_int (rparlist, ssection, &
        'nmaxIterations', rnlstsubsolver%nmaxIterations, 1)
    
    call parlst_getvalue_int (rparlist, ssection, &
        'nminIterations', rnlstsubsolver%nminIterations, 10)
    
    call parlst_getvalue_int (rparlist, ssection, &
        'ioutputLevel', rnlstsubsolver%ioutputLevel, 0)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'depsRel', rnlstsubsolver%depsRel, 1E-5_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'depsAbs', rnlstsubsolver%depsAbs, 1E-5_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'depsDiff', rnlstsubsolver%depsDiff, 0.0_DP)
    
    call parlst_getvalue_int (rparlist, ssection, &
        'icycle', rnlstsubsolver%icycle, 0)
    
    call parlst_getvalue_int (rparlist, ssection, &
        'istoppingCriterion', rnlstsubsolver%istoppingCriterion, 0)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'dalphaMin', rnlstsubsolver%dalphaMin,1.0_DP)
    
    call parlst_getvalue_double (rparlist, ssection, &
        'dalphaMax', rnlstsubsolver%dalphaMax,1.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_getParamsPrec_smoother (rparlist,ssection,rparams)
  
!<description>
  ! Reads parameters that configure the linear space-time smoother.
  !
  ! Note: This method does not read parameters of linear subsolvers
  ! in space. It just reads all the information which is level
  ! independent, so which holds simultaneously for all levels.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the solver parameters.
  character(len=*), intent(in) :: ssection
!</input>

!<output>
  ! Structure with the parameters to be initialised.
  type(t_nlstprec_mgrsmooth), intent(out) :: rparams
!</output>

!</subroutine>
    
    call parlst_getvalue_int (rparlist, ssection, &
        'cspaceTimeSmoother', rparams%cspaceTimeSmoother, 0)

    call parlst_getvalue_double (rparlist, ssection, &
        'domega', rparams%domega, 1.0_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'drelax', rparams%drelax, 1.0_DP)
        
    call parlst_getvalue_int (rparlist, ssection, &
        'nsmpre', rparams%nsmpre, 0)

    call parlst_getvalue_int (rparlist, ssection, &
        'nsmpost', rparams%nsmpost, 1)

    call parlst_getvalue_int (rparlist, ssection, &
        'ioutputLevel', rparams%ioutputLevel, 1)
        
    call parlst_getvalue_double (rparlist, ssection, &
        'depsRel', rparams%depsRel, 0.0_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'depsAbs', rparams%depsAbs, 0.0_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'depsDiff', rparams%depsDiff, 0.0_DP)

    call parlst_getvalue_int (rparlist, ssection, &
        'istoppingCriterion', rparams%istoppingCriterion, 1)

    ! Probably there is a linear subsolver. Get the name with the section
    ! configuring it.
    call parlst_getvalue_string (rparlist, ssection, &
        'slinearSolver', rparams%slinearSpaceSolver, "",bdequote=.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_getParamsPrec_sgrsolver (rparlist,ssection,rparams)
  
!<description>
  ! Reads parameters that configure a linear one-level space-time solver.
  !
  ! Note: This method does not read parameters of linear subsolvers
  ! in space.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section that contains the solver parameters.
  character(len=*), intent(in) :: ssection
!</input>

!<output>
  ! Structure with the parameters to be initialised.
  type(t_nlstprec_sgrsolver), intent(out) :: rparams
!</output>

!</subroutine>

    call parlst_getvalue_int (rparlist, ssection, &
        'ctypeSolver', rparams%ctypeSolver, 0)

    call parlst_getvalue_double (rparlist, ssection, &
        'domega', rparams%domega, 1.0_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'drelax', rparams%drelax, 1.0_DP)
                                
    call parlst_getvalue_int (rparlist, ssection, &
        'nminIterations', rparams%nminIterations, 1)
        
    call parlst_getvalue_int (rparlist, ssection, &
        'nmaxIterations', rparams%nmaxIterations, 100)
        
    call parlst_getvalue_double (rparlist, ssection, &
        'depsRel', rparams%depsRel, 1E-5_DP)
        
    call parlst_getvalue_double (rparlist, ssection, &
        'depsAbs', rparams%depsAbs, 1E-5_DP)
        
    call parlst_getvalue_double (rparlist, ssection, &
        'depsDiff', rparams%depsDiff, 0.0_DP)
        
    call parlst_getvalue_double (rparlist, ssection, &
        'ddivRel', rparams%ddivRel, 1.0_DP)
        
    call parlst_getvalue_int (rparlist, ssection, &
        'istoppingCriterion', rparams%istoppingCriterion, 0)

    call parlst_getvalue_int (rparlist, ssection, &
        'ioutputLevel', rparams%ioutputLevel, 100)

    ! ...

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_getParamsPrec_sgrprec (rparlist,ssection,rparams)
  
!<description>
  ! Reads parameters that configure the linear space-time coarse grid 
  ! preconditioner.
  !
  ! Note: This method does not read parameters of linear subsolvers
  ! in space.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in), target :: rparlist
  
  ! Section that contains the solver parameters.
  character(len=*), intent(in) :: ssection
!</input>

!<output>
  ! Structure with the parameters to be initialised.
  type(t_nlstprec_sgrprec), intent(out) :: rparams
!</output>

!</subroutine>

    call parlst_getvalue_double (rparlist, ssection, &
        'domega', rparams%domega, 1.0_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'drelax', rparams%drelax, 1.0_DP)

    call parlst_getvalue_int (rparlist, ssection, &
        'nminIterations', rparams%nminIterations, 1)

    call parlst_getvalue_int (rparlist, ssection, &
        'nmaxIterations', rparams%nmaxIterations, 100)

    call parlst_getvalue_double (rparlist, ssection, &
        'depsRel', rparams%depsRel, 1E-5_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'depsAbs', rparams%depsAbs, 1E-5_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'depsDiff', rparams%depsDiff, 0.0_DP)

    call parlst_getvalue_double (rparlist, ssection, &
        'ddivRel', rparams%ddivRel, 1.0_DP)

    call parlst_getvalue_int (rparlist, ssection, &
        'istoppingCriterion', rparams%istoppingCriterion, 0)

    call parlst_getvalue_int (rparlist, ssection, &
        'ioutputLevel', rparams%ioutputLevel, 100)

    call parlst_getvalue_int (rparlist, ssection, &
        'ifbSORPartialUpdate', rparams%ifbSORPartialUpdate,0)

    ! Probably there is a linear subsolver. Get the name with the section
    ! configuring it.
    call parlst_getvalue_string (rparlist, ssection, &
        'slinearSolver', rparams%slinearSpaceSolver, "",bdequote=.true.)
        
    ! Remember the parameter list, it contains additional parameters
    ! about the linear solver in space.
    rparams%p_rparlist => rparlist

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_initStartVector (rsettings,rsettingsSpaceDiscr,ilevel,rparlist,ssection,&
      rinitialCondition,rvector,rrhs,ioutputLevel)
  
!<description>
  ! Creates a space-time vector used as start vector.
!</description>

!<input>
  ! Settings structure.
  type(t_settings_optflow), intent(inout) :: rsettings
  
  ! Structure with space discretisation settings. Defines how the space 
  ! discretisation hierarchy in rsettings is set up.
  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr

  ! Level in the global space-time hierarchy where the start vector should
  ! be created.
  integer, intent(in) :: ilevel

  ! Parameter list with parameters configuring subsolvers e.g. in space.
  type(t_parlist), intent(in) :: rparlist
  
  ! Name of the section containing the definition of the start vector.
  character(len=*), intent(in) :: ssection
  
  ! The initial condition.
  type(t_anSolution), intent(inout) :: rinitialCondition
  
  ! Space-time vector containing the RHS of the forward problem.
  type(t_spaceTimeVector), intent(inout) :: rrhs
  
  ! Output level during the initialisation phase
  integer, intent(in) :: ioutputLevel
!</input>

!<inputoutput>
  ! Space-time vector containing the initial solution.
  ! The vector must already have been created. The routine will
  ! initialise the content of the vector.
  type(t_spaceTimeVector), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    integer :: ctypeStartVector,i,ispacelevel,isuccess
    real(dp) :: dtime,dstartVectorWeight
    logical :: bsuccess
    type(t_anSolution) :: rlocalsolution
    character(len=SYS_STRLEN) :: sstartVector,sstartVectorSolver
    character(len=SYS_STRLEN) :: sstartVectorBoundaryConditions
    character(len=SYS_STRLEN) :: sstartVectorPostprocessing
    type(t_vectorBlock) :: rvectorSpace
    type(t_simSolver) :: rsimsolver
    type(t_ccoptSpaceTimeMatrix) :: rspaceTimeMatrix
    type(t_spaceTimeMatrixDiscrData) :: rdiscrData
    type(t_optcBDC) :: rboudaryConditions
    type(t_matrixBlock) :: rmassMatrix

    ! Get the definition of the start vector
    call parlst_getvalue_int (rparlist, ssection, &
        'ctypeStartVector', ctypeStartVector, 0)

    call parlst_getvalue_double (rparlist, ssection, &
        'dstartVectorWeight', dstartVectorWeight, 1.0_DP)

    call parlst_getvalue_string (rparlist, ssection, &
        'sstartVector', sstartVector, "",bdequote=.true.)

    call parlst_getvalue_string (rparlist, ssection, &
        'sstartVectorSolver', sstartVectorSolver, "CC-NONLINEAR",bdequote=.true.)
        
    call parlst_getvalue_string (rparlist, ssection, &
        'sstartVectorBoundaryConditions', sstartVectorBoundaryConditions,&
        "BDCONDITIONS",bdequote=.true.)

    call parlst_getvalue_string (rparlist, ssection, &
        'sstartVectorPostprocessing', sstartVectorPostprocessing, &
        "",bdequote=.true.)
        
    ! Create a temp vector in space in the form of the space-time
    ! vector.
    call sth_getLevel (rsettings%rspaceTimeHierPrimalDual,ilevel,ispaceLevel=ispaceLevel)
    call lsysbl_createVectorBlock(rvector%p_rspaceDiscr,rvectorSpace,.true.)
        
    ! Create a mass matrix for the projection
    call imnat_getL2PrjMatrix(rsettings%rspaceAsmHierarchy%p_RasmTemplList(ispaceLevel),&
        CCSPACE_PRIMAL,rvector%p_rspaceDiscr,rmassMatrix)

    ! What to do?
    select case (ctypeStartVector)
    case (0)
      ! Initialise with zero.
      call sptivec_clearVector(rvector)
      
      ! Put the initial condition to the first timestep.
      call ansol_prjToVector (rinitialCondition,rinitialCondition%rtimeDiscr%dtimeInit,&
          rvectorSpace,1,3,rmassMatrix)
          
      call sptivec_setTimestepData (rvector, 1, rvectorSpace)
      
    case (1)
      ! Get the initial condition
      call ansol_prjToVector (rinitialCondition,rinitialCondition%rtimeDiscr%dtimeInit,&
          rvectorSpace,1,3,rmassMatrix)
          
      ! Propagate to all timesteps
      do i=1,rvector%neqTime
        call sptivec_setTimestepData (rvector, i, rvectorSpace)
      end do
      
    case (2)
      ! Read the analytic solution.
      call init_initFlow (rparlist,sstartVector,rlocalsolution,&
          rsettings%rtriaCoarse,rsettings%rrefinementSpace,&
          rsettingsSpaceDiscr,rsettings%rfeHierPrimal,rsettings%rboundary,isuccess)
      if (isuccess .eq. 1) then
        call output_line ('Flow created by simulation not yet supported!', &
            OU_CLASS_ERROR,OU_MODE_STD,'init_initOptControlTargetFlow')
        call sys_halt()
      end if
      
      ! Project it down to all timesteps.
      do i=1,rvector%neqTime
        call tdiscr_getTimestep(rvector%p_rtimeDiscr,i-1,dtime)
        call ansol_prjToVector (rlocalsolution,dtime,rvectorSpace,1,3,rmassMatrix)
        call sptivec_setTimestepData (rvector, i, rvectorSpace)
      end do
      
      ! Release the analytic flow again.
      call ansol_done(rlocalsolution)
    
    case (3)
      ! Initialise with zero.
      call sptivec_clearVector(rvector)
      
      ! Put the initial solution to the first timestep.
      call ansol_prjToVector (rinitialCondition,rinitialCondition%rtimeDiscr%dtimeInit,&
          rvectorSpace,1,3,rmassMatrix)
          
      call sptivec_setTimestepData (rvector, 1, rvectorSpace)
      
      ! Initialise a flow solver that simulates the flow to calculate the
      ! start vector.
      call fbsim_init (rsettings, rparlist, sstartVectorSolver, &
          1, rsettings%rfeHierPrimalDual%nlevels, FBSIM_SOLVER_NLFORWARD, rsimsolver)
          
      ! Create a nonlinear space-time matrix that resembles the current
      ! forward equation.
      call nlstslv_initStdDiscrData (rsettings,rsettings%rspaceTimeHierPrimalDual%nlevels,&
          rdiscrData)
      call stlin_initSpaceTimeMatrix (&
          rspaceTimeMatrix,MATT_OPTCONTROL,rdiscrData,rvector,&
          rsettings%rglobalData,rsettings%rdebugFlags)
      call fbsim_setMatrix (rsimsolver,rspaceTimeMatrix)
      
      ! Get the boundary conditions
      call init_initBoundaryConditions (rparlist,SEC_SBDEXPRESSIONS,&
          sstartVectorBoundaryConditions,rboudaryConditions)
      
      ! Get settings about postprocessing
      if (sstartVectorPostprocessing .ne. "") then
        call init_initForwardSimPostproc (rparlist,&
            sstartVectorPostprocessing,rsimsolver%rpostprocessing)
      end if
      
      ! Simulate...
      call fbsim_simulate (rsimsolver, rvector, rrhs, &
          rboudaryConditions,1, rvector%p_rtimeDiscr%nintervals,rvector,1.0_DP,bsuccess)
          
      ! Clean up
      call fbsim_done (rsimsolver)
      call stlin_releaseSpaceTimeMatrix(rspaceTimeMatrix)
    
      ! Probably print some statistical data
      if (ioutputLevel .ge. 1) then
        call output_lbrk ()
        call output_line ("Total time for forward simulation:      "//&
            trim(sys_sdL(rsimsolver%dtimeTotal,10)))

        call output_line ("Total time for nonlinear solver:        "//&
          trim(sys_sdL(rsimsolver%dtimeNonlinearSolver,10)))
          
        call output_line ("Total time for defect calculation:      "//&
          trim(sys_sdL(rsimsolver%dtimeDefectCalculation,10)))

        call output_line ("Total time for matrix assembly:         "//&
          trim(sys_sdL(rsimsolver%dtimeMatrixAssembly,10)))
          
        call output_line ("Total time for linear solver:           "//&
          trim(sys_sdL(rsimsolver%dtimeLinearSolver,10)))
          
        call output_line ("Total time for postprocessing:          "//&
          trim(sys_sdL(rsimsolver%dtimePostprocessing,10)))
          
        call output_line ("Total #iterations nonlinear solver:     "//&
          trim(sys_siL(rsimsolver%nnonlinearIterations,10)))
          
        call output_line ("Total #iterations linear solver:        "//&
          trim(sys_siL(rsimsolver%nlinearIterations,10)))
      end if
    
    end select
    
    ! Implement the Dirichlet boundary conditions.
    call tbc_implementBCsolution (rsettings%roptcBDC,rvector,&
        rsettings%rglobalData,rvectorSpace)
        
    ! Release memory
    call lsysbl_releaseMatrix(rmassMatrix)
    call lsysbl_releaseVector (rvectorSpace)

    ! Scale the start vector.
    do i=2,rvector%neqTime
      call sptivec_scaleVector (rvector,dstartVectorWeight,i)
    end do

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine init_initPostprocessing (rparlist,ssection,rsettingsSpaceDiscr,&
      rboundaryConditions,rtimeDiscr,rspaceDiscr,rpostproc,rsettings)

!<description>
  ! Initialises rdiscrData with standard values from rsettings.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection
  
  ! Structure with space discretisation settings
  type(t_settings_discr), intent(in) :: rsettingsSpaceDiscr

  ! Boundary conditions to use.
  type(t_optcBDC), intent(in), target  :: rboundaryConditions

  ! Underlying space discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Underlying time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
  
  ! Global settings structure.
  type(t_settings_optflow), intent(in) :: rsettings
!</input>

!<output>
  ! Postprocessing structure, to set up with data.
  type(t_optcPostprocessing), intent(out) :: rpostproc
!</output>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sstr
    integer :: isuccess

    ! Basic initialisation of the postprocessing structure
    call optcpp_initpostprocessing (rpostproc,CCSPACE_PRIMALDUAL,&
        rboundaryConditions,rtimeDiscr,rspaceDiscr)
        
    ! Read remaining parameters from the DAT file.
    !
    ! UCD export
    call parlst_getvalue_int (rparlist,ssection,&
        "ioutputUCD",rpostproc%ioutputUCD,0)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfilenameUCD",rpostproc%sfilenameUCD,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfinalSolutionFileName",rpostproc%sfinalSolutionFileName,&
        "",bdequote=.true.)

    ! function value calculation

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcFunctionalValues",rpostproc%icalcFunctionalValues,0)

    ! Body forces

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcForces",rpostproc%icalcForces,0)

    call parlst_getvalue_int (rparlist,ssection,&
        "ibodyForcesBdComponent",rpostproc%ibodyForcesBdComponent,2)

    call parlst_getvalue_double (rparlist,ssection,&
        "dbdForcesCoeff1",rpostproc%dbdForcesCoeff1,rsettings%rphysicsPrimal%dnu)

    call parlst_getvalue_double (rparlist,ssection,&
        "dbdForcesCoeff2",rpostproc%dbdForcesCoeff2,0.1_DP * 0.2_DP**2)

    call parlst_getvalue_int (rparlist,ssection,&
        "iwriteBodyForces",rpostproc%iwriteBodyForces,0)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfilenameBodyForces",rpostproc%sfilenameBodyForces,&
        "",bdequote=.true.)

    ! Flux

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcFlux",rpostproc%icalcFlux,0)

    call parlst_getvalue_int (rparlist,ssection,&
        "iwriteFlux",rpostproc%iwriteFlux,0)

    call parlst_getvalue_string (rparlist,ssection,&
        'sfilenameFlux',rpostproc%sfilenameFlux,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        'dfluxline',sstr,"",bdequote=.true.)
        
    call parlst_getvalue_string (rparlist,ssection,&
        'dfluxline',sstr,"",bdequote=.true.)
    if (sstr .ne. "") then
      ! Read the start/end coordinates
      read(sstr,*) rpostproc%Dfluxline(1),rpostproc%Dfluxline(2),&
          rpostproc%Dfluxline(3),rpostproc%Dfluxline(4)
    end if

    ! internal Energy

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcKineticEnergy",rpostproc%icalcKineticEnergy,1)

    call parlst_getvalue_int (rparlist,ssection,&
        "iwriteKineticEnergy",rpostproc%iwriteKineticEnergy,0)

    call parlst_getvalue_string (rparlist,ssection,&
        'sfilenameKineticEnergy',rpostproc%sfilenameKineticEnergy,"",bdequote=.true.)

    ! Error calculation

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcError",rpostproc%icalcError,0)

    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionReferenceFunction",sstr,"",bdequote=.true.)
    
    if (sstr .eq. "") &
        rpostproc%icalcError = 0
    
    if (rpostproc%icalcError .ne. 0) then
      ! Initialise the reference solution
      call init_initFlow (rparlist,sstr,rpostproc%ranalyticRefFunction,&
          rsettings%rtriaCoarse,rsettings%rrefinementSpace,&
          rsettingsSpaceDiscr,rsettings%rfeHierPrimal,rsettings%rboundary,isuccess)
      if (isuccess .eq. 1) then
        call output_line ("Flow created by simulation not yet supported!", &
            OU_CLASS_ERROR,OU_MODE_STD,"init_initOptControlTargetFlow")
        call sys_halt()
      end if
    end if

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine init_donePostprocessing (rpostproc)

!<description>
  ! Creans up the postprocessing structure.
!</description>

!<inputoutput>
  ! Postprocessing structure, to set up with data.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!</subroutine>

    ! Release the reference solution    
    if (rpostproc%icalcError .ne. 0) then
      call ansol_done(rpostproc%ranalyticRefFunction)
    end if
    
    call optcpp_donepostprocessing(rpostproc)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine init_discretiseRHS (rsettings,rrhs,rrhsDiscrete)
  
!<description>
  ! Discretises the RHS according to the space/time discretisation scheme
  ! in rrhsDiscrete. The result is summed up to rrhsDiscrete.
  ! Note: The boundary conditions are not implemented!
!</description>
  
!<input>
  ! All settings of the optimal control flow solver.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! Analytic solution defining the RHS of the equation.
  type(t_anSolution), intent(inout) :: rrhs
!</input>

!<inputoutput>
  ! Discrete RHS. The calculated RHS is added here.
  type(t_spacetimevector), intent(inout) :: rrhsDiscrete
!</inputoutput>

!</subroutine>

    ! CAll the assembly routine to do that task.
    call trhsevl_assembleRHS (rsettings%rglobalData, rrhs, rrhsDiscrete, &
        rsettings%rsettingsOptControl)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine init_generateInitCondRHS(rsettings,rtimediscr,rx,rb)

!<description>
  ! Generates the RHS vector used for the initial condition.
!</description>

!<input>
  ! All settings of the optimal control flow solver.
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Structure defining the time discretisation
  type(t_timeDiscretisation), intent(in) :: rtimediscr
  
  ! Space solution at the start time.
  type(t_vectorBlock), intent(in), target :: rx
!</input>

!<inputoutput>
  ! RHS at the initial time.
  type(t_vectorBlock), intent(inout) :: rb
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    logical :: bconvectionExplicit
    real(DP) :: dequationType,dtstep
    type(t_spatialMatrixDiscrData) :: rmatrixDiscr
    type(t_spatialMatrixNonlinearData), target :: rnonlinearity
    real(dp), dimension(:), pointer :: p_Ddata
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr
    
!    ! If the following constant is set from 1.0 to 0.0, the primal system is
!    ! decoupled from the dual system!
!    real(DP), parameter :: dprimalDualCoupling = 1.0_DP
!    
!    ! If the following constant is set from 1.0 to 0.0, the dual system is
!    ! decoupled from the primal system!
!    real(DP), parameter :: ddualPrimalCoupling = 1.0_DP
!    
!    ! If the following parameter is set from 1.0 to 0.0, the time coupling
!    ! is disabled, resulting in a stationary simulation in every timestep.
!    real(DP), parameter :: dtimeCoupling = 1.0_DP

    real(DP) :: dprimalDualCoupling,ddualPrimalCoupling,dtimeCoupling
    
    dprimalDualCoupling = rsettings%rdebugFlags%dprimalDualCoupling
    ddualPrimalCoupling = rsettings%rdebugFlags%ddualPrimalCoupling
    dtimeCoupling = rsettings%rdebugFlags%dtimeCoupling
    
    ! Create two temp vectors representing the first timestep.
    p_rspaceDiscr => rx%p_rblockDiscr
        
    ! DEBUG!!!
    call lsysbl_getbase_double (rb,p_Ddata)

    ! Form a t_spatialMatrixNonlinearData structure that encapsules the nonlinearity
    ! of the spatial matrix.
    rnonlinearity%p_rvector1 => rx
    rnonlinearity%p_rvector2 => rx
    rnonlinearity%p_rvector3 => rx

    ! The initial condition is implemented as:
    !
    !   (M/dt + A) y_0  =  b_0 := (M/dt + A) y^0
    !
    ! i.e. we take the solution vector of the 0th timestep. multiply it by the
    ! (Navier--)Stokes equation and what we receive is the RHS for the
    ! terminal condition. We only have to take care of bondary conditions.

    call smva_getDiscrData (rsettings,rsettings%rfeHierPrimalDual%nlevels,rmatrixDiscr)
    call smva_initNonlinMatrix (rnonlinearSpatialMatrix,rmatrixDiscr,rnonlinearity)
    
    ! Disable the submatices for the dual solution and the coupling.
    ! We only want to generate the RHS for the primal solution.
    call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
    call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
    call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,2)

    ! Set up the matrix weights the matrix in the 0th timestep.
    ! We set up only the stuff for the primal equation; for setting up
    ! the RHS, there is no dual equation and also no coupling between the
    ! primal and dual solutions.
    bconvectionExplicit = rsettings%rsettingsOptControl%iconvectionExplicit .ne. 0

    dequationType = 1.0_DP
    if (rsettings%rsettingsOptControl%ispaceTimeFormulation .ne. 0) &
      dequationType = -1.0_DP
      
    call tdiscr_getTimestep(rtimediscr,1,dtstep=dtstep)
    rnonlinearSpatialMatrix%Dalpha(1,1) = dtimeCoupling * 1.0_DP/dtstep
    rnonlinearSpatialMatrix%Dtheta(1,1) = 1.0_DP
    
    if (.not. bconvectionExplicit) then
      rnonlinearSpatialMatrix%Dgamma(1,1) = real(1-rsettings%rphysicsPrimal%iequation,DP)
    end if

    rnonlinearSpatialMatrix%Deta(1,1) = 1.0_DP
    rnonlinearSpatialMatrix%Dtau(1,1) = 1.0_DP
        
    ! Create by substraction: rd = 0*rd - (- A11 x1) = A11 x1.
    ! Clear the primal RHS for that purpose.
    call lsyssc_clearVector(rb%RvectorBlock(1))
    call lsyssc_clearVector(rb%RvectorBlock(2))
    call lsyssc_clearVector(rb%RvectorBlock(3))
    call smva_assembleDefect (rnonlinearSpatialMatrix,rx,rb,-1.0_DP)
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine init_implementInitCondRHS (rsettings,rsolution,rrhs)

!<description>
  ! Generates the RHS vector used for the initial condition and implements
  ! it into the space-time vector rrhs.
!</description>

!<input>
  ! All settings of the optimal control flow solver.
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Space time solution vector.
  type(t_spaceTimeVector), intent(in), target :: rsolution
!</input>

!<inputoutput>
  ! Vector where the first entry receives the RHS for the initial condition.
  ! The vector must have been initialised. The vector data is overwritten
  ! by the RHS data.
  type(t_spaceTimeVector), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr
    type(t_vectorBlock), target :: rx,rb
    
    ! Create two temp vectors representing the first timestep.
    p_rtimeDiscr => rsolution%p_rtimeDiscr
    p_rspaceDiscr => rsolution%p_rspaceDiscr
        
    call lsysbl_createVectorBlock(p_rspaceDiscr,rx)
    call lsysbl_createVectorBlock(p_rspaceDiscr,rb)
    call sptivec_getTimestepData (rsolution, 1, rx)
    call sptivec_getTimestepData (rrhs, 1, rb)

    ! DEBUG!!!
    call lsysbl_getbase_double (rb,p_Ddata)

    call init_generateInitCondRHS(rsettings,p_rtimediscr,rx,rb)

    ! Save the new RHS.
    call sptivec_setTimestepData (rrhs, 1, rb)
    
    ! Release the temp vectors.
    call lsysbl_releaseVector(rb)
    call lsysbl_releaseVector(rx)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine init_getSpaceDiscrSettings (rparlist,ssection,rsettingsSpaceDiscr)

!<description>
  ! Extracts main discretisation settings in space from the DAT file.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section containing the parameters
  character(len=*), intent(in) :: ssection
!</input>

!<output>
   ! Structure receiving the main discretisation settings
   type(t_settings_discr), intent(out) :: rsettingsSpaceDiscr
!</output>

!</subroutine>

    character(len=SYS_STRLEN) :: sstr

    call parlst_getvalue_int (rparlist,ssection,&
        'ielementType',rsettingsSpaceDiscr%ielementType,3)
                              
    call parlst_getvalue_string (rparlist,ssection,'scubStokes',sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          'icubStokes',rsettingsSpaceDiscr%icubStokes,int(SPDISC_CUB_AUTOMATIC))
    else
      rsettingsSpaceDiscr%icubStokes = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rparlist,ssection,'scubB',sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          'icubB',rsettingsSpaceDiscr%icubB,int(SPDISC_CUB_AUTOMATIC))
    else
      rsettingsSpaceDiscr%icubB = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rparlist,ssection,'scubF',sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          'icubF',rsettingsSpaceDiscr%icubF,int(SPDISC_CUB_AUTOMATIC))
    else
      rsettingsSpaceDiscr%icubF = cub_igetID(sstr)
    end if
    
    ! Which cubature rule to use for mass matrices?
    call parlst_getvalue_string (rparlist,ssection,'scubMass',sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          'icubM',rsettingsSpaceDiscr%icubMass,int(SPDISC_CUB_AUTOMATIC))
    else
      rsettingsSpaceDiscr%icubMass = cub_igetID(sstr)
    end if    

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_initSpacePrec_mg (rparlist,ssection,nequations,&
!      nlevels,RfilterChain,p_rsolverNode)
!  
!!<description>
!  ! Creates a node for the space multigrid solver and initialises it with
!  ! parameters from a DAT file.
!!</description>
!
!!<input>
!  ! Parameter list
!  type(t_parlist), intent(in) :: rparlist
!  
!  ! Section that contains the solver parameters.
!  character(len=*), intent(in) :: ssection
!  
!  ! Type of equation. =3 for primal or dual equation. =6 for primal+dual equation
!  
!  ! Number of levels
!  integer, intent(in) :: nlevels
!  
!  ! Filter chain to be applied during the application
!  type(t_filterChain), dimension(:), intent(in), target :: RfilterChain
!!</input>
!
!!<output>
!  ! Pointer to a space solver node.
!  type(t_linsolNode), pointer :: p_rsolverNode
!!</output>
!
!!</subroutine>
!
!    ! Create the solver
!    call linsol_initMultigrid2 (p_rsolverNode,nlevels,RfilterChain)
!    
!    ! Manually trim the coarse grid correction in Multigrid to multiply the 
!    ! pressure equation with -1. This (un)symmetrises the operator and gives
!    ! much better convergence rates.
!    call cgcor_release(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection)
!    call cgcor_init(p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection,nequations)
!    p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(3) &
!        = -1.0_DP
!    if (nequations .gt. 3) then
!      p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%p_DequationWeights(6) &
!          = -1.0_DP
!    end if
!
!    ! Init standard solver parameters and extended multigrid parameters
!    ! from the DAT file.
!    call linsolinit_initParams (p_rsolverNode,rparlist,ssection,&
!        LINSOL_ALG_UNDEFINED)
!    call linsolinit_initParams (p_rsolverNode,rparlist,ssection,&
!        LINSOL_ALG_MULTIGRID2)
!
!  end subroutine

!Nach der 1. Iteration ist der Def. nicht mehr 0 tritz Stokes.
!Das Startresiduum passt nicht

end module

!##############################################################################
!# ****************************************************************************
!# <name> kktsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the realisation of the discrete KKT system
!# which stems from the discretisation of the optimisation problem.
!# </purpose>
!##############################################################################

module kktsystem

  use fsystem
  use storage
  use genoutput
  
  use mprimitives
  use spatialdiscretisation
  use timediscretisation
  use element
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  
  use scalarpde
  use dofmapping
  use bcassemblybase
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  use linearsolver
  use feevaluation
  use statistics
  use numbersets
  
  use spacetimevectors
  use analyticsolution
  
  use constantsdiscretisation
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresboundaryconditions
  use structuresgeneral
  use assemblytemplates
  
  use structuresoperatorasm
  use spacematvecassembly
  use spatialbc
  
  use kktsystemspaces
  use spacesolver
  
  use newtonderivative
  
  implicit none
  
  private

!<types>

!<typeblock>

  ! This structure encapsules the discrete KKT system. It holds
  ! a solution of the primal equation, the dual equation and the
  ! corresponding control (if not directly calculated from the dual
  ! solution on-demand).
  type t_kktsystem
  
    ! Solution of the primal equation of the KKT system.
    type(t_primalSpace), pointer :: p_rprimalSol => null()
    
    ! Solution of the dual equation of the KKT system.
    type(t_dualSpace), pointer :: p_rdualSol => null()
    
    ! Solution of the control equation of the KKT system.
    type(t_controlSpace), pointer :: p_rcontrol => null()

    ! "Intermediate" control. This is a control-like vector
    !    u~ = -1/alpha lambda
    ! computed from the dual solution lambda.
    type(t_controlSpace), pointer :: p_rintermedControl => null()

    ! Underlying space-time operator assembly hierarchy
    ! specifying all possible space and time discretisations / levels.
    type(t_spacetimeOpAsmHierarchy), pointer :: p_roperatorAsmHier => null()
    
    ! Boundary conditions of the system.
    type(t_optcBDC), pointer :: p_roptcBDC => null()

    ! Space-level in the global space-time hierarchy, the solver should be applied to.
    integer :: ispacelevel = 0

    ! Time-level in the global space-time hierarchy, the solver should be applied to.
    integer :: itimelevel = 0
    
  end type

!</typeblock>

  public :: t_kktsystem

!<typeblock>

  ! Encapsules a directional derivative of a KKT system.
  type t_kktsystemDirDeriv
  
    ! Reference to the underlying KKT system. Defines the evaluation
    ! point where the KKT system is linearised.
    type(t_kktsystem), pointer :: p_rkktsystem => null()
  
    ! Solution of the linearised primal equation of the KKT system.
    ! Specifies the directional derivative of the primal equation
    ! into a given direction p_rprimalDirection.
    type(t_primalSpace), pointer :: p_rprimalSolLin => null()
    
    ! Solution of the linearised dual equation of the KKT system.
    ! Specifies the directional derivative of the dual equation
    ! into the direction specified by p_rdualDirection.
    type(t_dualSpace), pointer :: p_rdualSolLin => null()

    ! Solution of the linearised control equation of the KKT system.
    ! Specifies the directional derivative of the control equation
    ! into the direction specified by p_rcontrolDirection.
    type(t_controlSpace), pointer :: p_rcontrolLin => null()
    
  end type

!</typeblock>

  public :: t_kktsystemDirDeriv

!</types>

  ! Initialises a KKT system.
  public :: kkt_initKKTsystem

  ! Cleans up a KKT system.
  public :: kkt_doneKKTsystem

  ! Initialises the structure for a directional derivative of the
  ! solutions of a KKT system.
  public :: kkt_initKKTsystemDirDeriv

  ! Cleans up the structure for a directional derivative of the
  ! solutions of a KKT system.
  public :: kkt_doneKKTsystemDirDeriv

  ! Solve the primal equation
  public :: kkt_solvePrimal

  ! Solve the dual equation
  public :: kkt_solveDual

  ! Calculate the control from the solution of the primal/dual equation
  public :: kkt_dualToControl
  
  ! Calculate the residual of the control equation(s)
  public :: kkt_calcControlRes

  ! Solve the primal equation of the linearised KKT system
  public :: kkt_solvePrimalDirDeriv

  ! Solve the dual equation of the linearised KKT system
  public :: kkt_solveDualDirDeriv

  ! Calculate the control of the linearised KKT system 
  ! from the solution of the primal/dual equation
  public :: kkt_dualToControlDirDeriv

  ! Calculate the residual of the control equation(s) in the linearised KKT system
  public :: kkt_calcControlResDirDeriv
  
  ! Applies the operator of the control equation in the linearised KKT system
  public :: kkt_applyControlDirDeriv

  ! Clear a KKT structure
  public :: kkt_clear

  ! Clear a KKT derivative structure
  public :: kkt_clearDirDeriv

  ! Calculates the norm of a residual in the control space,
  ! weighted by the corresponding weighting factors.
  public :: kkt_controlResidualNorm

  ! Calculates the control at a given point in time.
  public :: kkt_getControlAtTime

  ! Saves the solutions of a KKT system to file sequences
  public :: kkt_saveToFiles

  ! Saves the solutions of a KKT system to file sequences
  public :: kkt_saveDirDerivToFiles

contains

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_initKKTsystem (rkktsystem,&
      roperatorAsmHier,ispacelevel,itimelevel,roptcBDC)
  
!<description>
  ! Initialises a KKT system structure.
!</description>

!<input>  
  ! Parameters for the assembly of space-time operators
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier
  
  ! Boudary conditions of the KKT system
  type(t_optcBDC), intent(in), target :: roptcBDC

  ! Space-level in the global space-time hierarchy, the solver should be applied to.
  integer, intent(in) :: ispacelevel

  ! Time-level in the global space-time hierarchy, the solver should be applied to.
  integer, intent(in) :: itimelevel
!</input>

!<output>
  ! Structure defining the KKT system.
  type(t_kktsystem), intent(out) :: rkktsystem
!</output>

!</subroutine>

    ! local variables
    type(t_spacetimeOperatorAsm) :: roperatorAsm

    ! Remember the structures
    rkktsystem%p_roperatorAsmHier => roperatorAsmHier
    rkktsystem%p_roptcBDC => roptcBDC
    rkktsystem%ispacelevel = ispacelevel
    rkktsystem%itimelevel = itimelevel
    
    ! Get the underlying space and time discretisation structures.
    call stoh_getOpAsm_slvtlv (roperatorAsm,roperatorAsmHier,ispacelevel,itimelevel)
    
    ! Allocate memory for the solutions of the KKT system.
    allocate (rkktsystem%p_rprimalSol)
    call kktsp_initPrimalVector (rkktsystem%p_rprimalSol,&
        roperatorAsm%p_rspaceDiscrPrimal,roperatorAsm%p_rtimeDiscrPrimal)

    allocate (rkktsystem%p_rdualSol)
    call kktsp_initDualVector (rkktsystem%p_rdualSol,&
        roperatorAsm%p_rspaceDiscrDual,roperatorAsm%p_rtimeDiscrDual)

    allocate (rkktsystem%p_rcontrol)
    call kktsp_initControlVector (rkktsystem%p_rcontrol,&
        roperatorAsm%p_rspaceDiscrControl,roperatorAsm%p_rtimeDiscrControl)
  
    allocate (rkktsystem%p_rintermedControl)
    call kktsp_initControlVector (rkktsystem%p_rintermedControl,&
        roperatorAsm%p_rspaceDiscrControl,roperatorAsm%p_rtimeDiscrControl)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_doneKKTsystem (rkktsystem)
  
!<description>
  ! Cleans up a KKT system structure.
!</description>

!<inputoutput>
  ! Structure defining the KKT system.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>

    ! Clean up the structures
    nullify(rkktsystem%p_roperatorAsmHier)
    
    ! Release memory
    call kktsp_donePrimalVector (rkktsystem%p_rprimalSol)
    deallocate (rkktsystem%p_rprimalSol)

    call kktsp_doneDualVector (rkktsystem%p_rdualSol)
    deallocate (rkktsystem%p_rdualSol)

    call kktsp_doneControlVector (rkktsystem%p_rintermedControl)
    deallocate (rkktsystem%p_rintermedControl)

    call kktsp_doneControlVector (rkktsystem%p_rcontrol)
    deallocate (rkktsystem%p_rcontrol)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_initKKTsystemDirDeriv (rkktsystemDirDeriv,rkktsystem)
  
!<description>
  ! Initialises the structure for a directional derivative of the
  ! solutions of a KKT system.
!</description>

!<input>  
  ! Structure defining the KKT system.
  type(t_kktsystem), intent(in), target :: rkktsystem
!</input>

!<output>
  ! Structure defining the directional derivative.
  type(t_kktsystemDirDeriv), intent(out), target :: rkktsystemDirDeriv
!</output>

!</subroutine>
  
    ! Remember the structures
    rkktsystemDirDeriv%p_rkktsystem => rkktsystem
    
    ! Allocate memory for the solutions of the KKT system.
    allocate (rkktsystemDirDeriv%p_rprimalSolLin)
    call kktsp_initPrimalVector (rkktsystemDirDeriv%p_rprimalSolLin,&
        rkktsystem%p_rprimalSol%p_rvector%p_rspaceDiscr,&
        rkktsystem%p_rprimalSol%p_rvector%p_rtimeDiscr)

    allocate (rkktsystemDirDeriv%p_rdualSolLin)
    call kktsp_initDualVector (rkktsystemDirDeriv%p_rdualSolLin,&
        rkktsystem%p_rdualSol%p_rvector%p_rspaceDiscr,&
        rkktsystem%p_rdualSol%p_rvector%p_rtimeDiscr)

    allocate (rkktsystemDirDeriv%p_rcontrolLin)
    call kktsp_initControlVector (rkktsystemDirDeriv%p_rcontrolLin,&
        rkktsystem%p_rcontrol%p_rvector%p_rspaceDiscr,&
        rkktsystem%p_rcontrol%p_rvector%p_rtimeDiscr)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_doneKKTsystemDirDeriv (rkktsystemDirDeriv)
  
!<description>
  ! Cleans up the structure for a directional derivative of the
  ! solutions of a KKT system.
!</description>

!<inputoutput>
  ! Structure defining the directional derivative.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
!</inputoutput>

!</subroutine>

    ! Clean up the structures
    nullify(rkktsystemDirDeriv%p_rkktsystem)
    
    ! Release memory
    call kktsp_donePrimalVector (rkktsystemDirDeriv%p_rprimalSolLin)
    deallocate (rkktsystemDirDeriv%p_rprimalSolLin)

    call kktsp_doneDualVector (rkktsystemDirDeriv%p_rdualSolLin)
    deallocate (rkktsystemDirDeriv%p_rdualSolLin)

    call kktsp_doneControlVector (rkktsystemDirDeriv%p_rcontrolLin)
    deallocate (rkktsystemDirDeriv%p_rcontrolLin)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solvePrimal (rkktsystem,rspaceSolver,cspatialInitCondPolicy,rstatistics)
  
!<description>
  ! Solves the primal equation in the KKT system based on the control in the
  ! rkktsystem structure.
!</description>
  
!<input>
  ! Defines a policy how to generate the initial condition of a timestep.
  ! =0: Always take zero
  ! =1: Propagate the solution of the previous/next timestep to the
  !     current one. (Default)
  ! =2: Take the solution of the last space-time iteration
  integer, intent(in) :: cspatialInitCondPolicy
!</input>
  
!<inputoutput>
  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics
!<output>

!</subroutine>

    ! local variables
    integer :: idofTime, ierror
    type(t_spaceslSolverStat) :: rstatLocal
    
    call stat_startTimer (rstatistics%rtotalTime)

    ! -------------------------------------------------------------------------
    ! Basic description
    !
    ! What we have to do here is a loop through all timesteps.
    ! For a nonlinear state equation, the problems to solve in each timestep
    ! are nonlinear. The control is given in rkktsystem and acts as right-hand
    ! side, on the boundary or whereever.
    !
    ! For the nonlinear Navier-Stokes equations with distributed control
    ! for example, the equation to solve here reads
    !
    !    y_t - Laplace(y) + (y grad) y + grad(p) = u
    !                                     -div y = 0
    !
    ! with the control u being given in rkktsystem.
    ! -------------------------------------------------------------------------
    ! All the timestepping weights and so on are realised in the 
    ! matrix-vector assembly routines. Here, we only have to
    ! apply a loop over all unknowns in time.
    ! -------------------------------------------------------------------------
    
    ! Initialise basic solver structures
    call spaceslh_initStructure (rspaceSolver, &
        rkktsystem%ispacelevel, rkktsystem%itimelevel, &
        rkktsystem%p_roperatorAsmHier,rstatLocal,ierror)

    ! Sum up statistics
    call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)

    if (ierror .ne. 0) then
      call output_line("Error initialising the solver structures.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_solvePrimal")
      call sys_halt()
    end if
    
    ! -----------------------
    ! Initial condition
    ! -----------------------
    ! Take the initial condition from the structure with the discrete
    ! initial condition
    call smva_implementInitCond (rkktsystem%p_rprimalSol,&
        rkktSystem%p_roperatorAsmHier%p_rdiscreteInitCond)

    ! -----------------------
    ! Loop over all timesteps
    ! -----------------------
    do idofTime = 1,rkktsystem%p_rprimalSol%p_rvector%NEQtime
    
      ! Apply the solver to update the solution in timestep idofTime.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call spaceslh_solve (rspaceSolver,idofTime,cspatialInitCondPolicy,SPACESLH_EQNF_DEFAULT,&
          rstatLocal,rkktsystem%ispacelevel,rkktsystem%p_rprimalSol,rcontrol=rkktsystem%p_rcontrol)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      ! Sum up statistics
      call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)
      
    end do ! step
   
    call spaceslh_doneStructure (rspaceSolver)
    
    call stat_stopTimer (rstatistics%rtotalTime)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solveDual (rkktsystem,rspaceSolver,cspatialInitCondPolicy,rstatistics)
  
!<description>
  ! Solves the dual equation in the KKT system based on the control and the
  ! primal solution in the rkktsystem structure.
!</description>
  
!<input>
  ! Defines a policy how to generate the initial condition of a timestep.
  ! =0: Always take zero
  ! =1: Propagate the solution of the previous/next timestep to the
  !     current one. (Default)
  ! =2: Take the solution of the last space-time iteration
  integer, intent(in) :: cspatialInitCondPolicy
!</input>
  
!<inputoutput>
  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem

  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics
!<output>

!</subroutine>

    ! local variables
    integer :: idofTime, ierror
    type(t_spaceslSolverStat) :: rstatLocal

    call stat_startTimer (rstatistics%rtotalTime)

    ! -------------------------------------------------------------------------
    ! The solution of the dual equation is rather similar to the primal
    ! equation, but even simpler. The timeloop loops backward through the
    ! timesteps, and in every timestep, a linear problem has to be solved.
    ! -------------------------------------------------------------------------
    
    ! Initialise basic solver structures
    call spaceslh_initStructure (rspaceSolver, &
        rkktsystem%ispacelevel, rkktsystem%itimelevel, &
        rkktsystem%p_roperatorAsmHier,rstatLocal,ierror)

    ! Sum up statistics
    call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)

    if (ierror .ne. 0) then
      call output_line("Error initialising the solver structures.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_solveDual")
      call sys_halt()
    end if
    
    ! ----------------------------------
    ! Loop over all timesteps, backwards
    ! ----------------------------------
    do idofTime = rkktsystem%p_rprimalSol%p_rvector%NEQtime,1,-1
    
      ! Apply the solver to update the solution in timestep idofTime.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call spaceslh_solve (rspaceSolver,idofTime,cspatialInitCondPolicy,SPACESLH_EQNF_DEFAULT,&
          rstatLocal,rkktsystem%ispacelevel,rkktsystem%p_rprimalSol,rdualSol=rkktsystem%p_rdualSol)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      ! Sum up statistics
      call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)

    end do ! step
   
    call spaceslh_doneStructure (rspaceSolver)
    
    call stat_stopTimer (rstatistics%rtotalTime)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_dualToControl (rkktsystem,rcontrol)
  
!<description>
  ! Calculates the intermediate control
  !    u_intermed = u - (lambda + alpha u)
  ! from the dual solution lambda and the old control.
  !
  ! Applies the projection and saves the intermediate control
  ! in rcontrol:
  !
  !   rcontrol = P(u_intermed)
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system.
  ! The control, primal and dual variable in this structure are used to
  ! calculate the residual.
  type(t_kktsystem), intent(inout), target :: rkktsystem

  ! Control vector that receives the new control.
  type(t_controlSpace), intent(inout) :: rcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: icomp,istep
    real(DP) :: dtheta,dwmin,dwmax,dtime
    type(t_vectorBlock), pointer :: p_rdualSpace, p_rcontrolSpace, p_rintermedControl
    type(t_vectorBlock), pointer :: p_rcontrolSpaceOutput
    type(t_spaceTimeVector), pointer :: p_rdualSol
    type(t_optcBDCSpace) :: roptcBDCspace

    type(t_settings_physics), pointer :: p_rphysics
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl

    type(t_spacetimeOperatorAsm) :: roperatorAsm

    ! Fetch some structures
    p_rphysics => &
        rkktsystem%p_roperatorAsmHier%ranalyticData%p_rphysics
    p_rsettingsOptControl => &
        rkktsystem%p_roperatorAsmHier%ranalyticData%p_rsettingsOptControl

    ! Get the underlying space and time discretisation structures.
    call stoh_getOpAsm_slvtlv (roperatorAsm,&
        rkktsystem%p_roperatorAsmHier,rkktsystem%ispacelevel,rkktsystem%itimelevel)

    ! This is strongly equation and problem dependent
    ! and may imply a projection to the admissible set.
    !
    ! We apply a loop over all steps and construct the
    ! control depending on the timestep scheme.
    !
    ! Which timestep scheme do we have?
    
    p_rdualSol => rkktsystem%p_rdualSol%p_rvector
    
    ! Timestepping technique?
    select case (p_rdualSol%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = p_rdualSol%p_rtimeDiscr%dtheta
      
      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (p_rdualSol%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)
      
        ! Loop over all timesteps.
        do istep = 1,p_rdualSol%p_rtimeDiscr%nintervals+1
        
          ! Fetch the dual and control vectors.
          call sptivec_getVectorFromPool (&
              rkktsystem%p_rdualSol%p_rvectorAccess,istep,p_rdualSpace)

          call sptivec_getVectorFromPool (&
              rkktsystem%p_rcontrol%p_rvectorAccess,istep,p_rcontrolSpace)

          call sptivec_getVectorFromPool (&
              rkktsystem%p_rintermedControl%p_rvectorAccess,istep,p_rintermedControl)

          call sptivec_getVectorFromPool (&
              rcontrol%p_rvectorAccess,istep,p_rcontrolSpaceOutput)

          ! icomp counts the component in the control
          icomp = 0
          
          ! Which equation do we have?
          select case (p_rphysics%cequation)
          
          ! -------------------------------------------------------------
          ! Stokes/Navier Stokes.
          ! -------------------------------------------------------------
          case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
            
            ! Which type of control is applied?
            
            ! -----------------------------------------------------------
            ! Distributed control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

              ! The first two components of the control read
              !
              !    u_intermed = u - (lambda + alpha u)
              !               = (1-alpha) u - lambda
              !
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpace%RvectorBlock(icomp),p_rdualSpace%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaDistC,-1.0_DP,&
                  p_rintermedControl%RvectorBlock(icomp))

              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpace%RvectorBlock(icomp),p_rdualSpace%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaDistC,-1.0_DP,&
                  p_rintermedControl%RvectorBlock(icomp))
                  
              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
                  call output_line("Alpha=0 not possible without contraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
                  call sys_halt()
                end if
                
                icomp = icomp - 2                

                icomp = icomp + 1
                call lsyssc_copyVector (&
                    p_rintermedControl%RvectorBlock(icomp),p_rcontrolSpaceOutput%RvectorBlock(icomp))

                icomp = icomp + 1
                call lsyssc_copyVector (&
                    p_rintermedControl%RvectorBlock(icomp),p_rcontrolSpaceOutput%RvectorBlock(icomp))

              ! ----------------------------------------------------------
              ! Box constraints, implemented by DOF
              ! ----------------------------------------------------------
              case (1)
              
                ! Applying the projection to the intermediate control gives the control:
                !
                !   u = P(u_intermed)
                
                icomp = icomp - 2
                
                dwmin = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmin1
                dwmax = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmax1
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax)

                dwmin = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmin2
                dwmax = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmax2
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax)

              case default          
                call output_line("Unknown constraints",&
                    OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
                call sys_halt()

              end select ! constraints

            end if ! alphaDistC
          
            ! -----------------------------------------------------------
            ! L2 Boundary control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaL2BdC .ge. 0.0_DP) then

              ! No control in the initial solution
              if (istep .gt. 1) then

                ! Characteristics of the current timestep.
                call tdiscr_getTimestep(roperatorasm%p_rtimeDiscrPrimal,istep-1,dtime)

                ! Calculate the region where boundary control is applied
                call sbc_assembleBDconditions (rkktSystem%p_roptcBDC,roptcBDCSpace,dtime,&
                    p_rphysics%cequation,OPTP_PRIMAL,SBC_DIRICHLETBCC,&
                    p_rintermedControl%p_rblockDiscr)

                ! The first two components of the control read
                !
                !    u_intermed  =  u - [ alpha u - (nu dn lambda - xi n) ]
                !                =  (1-alpha) u + nu dn lambda - xi n
                !
                ! Calculate "nu dn lambda - xi n"
                call kkt_calcL2BdCNavSt (roperatorAsm%p_rasmTemplates,p_rphysics,&
                    p_rdualSpace,p_rintermedControl,icomp+1,roptcBDCSpace)
                
                ! Calculate
                !    u_intermed = (1-alpha) u + u_intermed
                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rcontrolSpace%RvectorBlock(icomp),p_rintermedControl%RvectorBlock(icomp),&
                    1.0_DP-p_rsettingsOptControl%dalphaL2BdC,1.0_DP,&
                    p_rintermedControl%RvectorBlock(icomp))

                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rcontrolSpace%RvectorBlock(icomp),p_rintermedControl%RvectorBlock(icomp),&
                    1.0_DP-p_rsettingsOptControl%dalphaL2BdC,1.0_DP,&
                    p_rintermedControl%RvectorBlock(icomp))
                    
                ! Do we have constraints?
                select case (p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%cconstraints)

                ! ----------------------------------------------------------
                ! No constraints
                ! ----------------------------------------------------------
                case (0)

                  if (p_rsettingsOptControl%dalphaL2BdC .eq. 0.0_DP) then
                    call output_line("Alpha=0 not possible without contraints",&
                        OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
                    call sys_halt()
                  end if
                  
                  icomp = icomp - 2                

                  icomp = icomp + 1
                  call lsyssc_copyVector (&
                      p_rintermedControl%RvectorBlock(icomp),p_rcontrolSpaceOutput%RvectorBlock(icomp))

                  icomp = icomp + 1
                  call lsyssc_copyVector (&
                      p_rintermedControl%RvectorBlock(icomp),p_rcontrolSpaceOutput%RvectorBlock(icomp))

                ! ----------------------------------------------------------
                ! Box constraints, implemented by DOF
                ! ----------------------------------------------------------
                case (1)
                
                  ! Applying the projection to the intermediate control gives the control:
                  !
                  !   u = P(u_intermed)
                  
                  icomp = icomp - 2
                  
                  dwmin = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmin1
                  dwmax = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmax1
                  icomp = icomp + 1
                  call nwder_applyMinMaxProjByDof (&
                      p_rcontrolSpaceOutput%RvectorBlock(icomp),1.0_DP,&
                      1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                      1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax)

                  dwmin = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmin2
                  dwmax = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmax2
                  icomp = icomp + 1
                  call nwder_applyMinMaxProjByDof (&
                      p_rcontrolSpaceOutput%RvectorBlock(icomp),1.0_DP,&
                      1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                      1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax)

                case default          
                  call output_line("Unknown constraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
                  call sys_halt()

                end select ! constraints
                
              end if

            end if ! alphaL2BdC

          ! -------------------------------------------------------------
          ! Heat equation
          ! -------------------------------------------------------------
          case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
            
            ! Which type of control is applied?
            
            ! -----------------------------------------------------------
            ! Distributed control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

              ! The first two components of the control read
              !
              !    u_intermed = u - (lambda + alpha u)
              !               = (1-alpha) u - lambda
              !
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpace%RvectorBlock(icomp),p_rdualSpace%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaDistC,-1.0_DP,&
                  p_rintermedControl%RvectorBlock(icomp))

              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
                  call output_line("Alpha=0 not possible without contraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
                  call sys_halt()
                end if
                
                icomp = icomp - 1

                icomp = icomp + 1
                call lsyssc_copyVector (&
                    p_rintermedControl%RvectorBlock(icomp),p_rcontrolSpaceOutput%RvectorBlock(icomp))

              ! ----------------------------------------------------------
              ! Box constraints, implemented by DOF
              ! ----------------------------------------------------------
              case (1)
              
                ! rcontrol contains the intermediate control as well.
                ! Applying the projection gives the control:
                !
                !   u = P(u_intermed)

                icomp = icomp - 1
              
                dwmin = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmin1
                dwmax = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmax1
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax)

              case default          
                call output_line("Unknown constraints",&
                    OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
                call sys_halt()

              end select ! constraints

            end if ! alpha

            ! -----------------------------------------------------------
            ! L2 Boundary control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaL2BdC .ge. 0.0_DP) then

              call output_line("Boundary control not available.",&
                  OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
              call sys_halt()

            end if
            
          end select ! equation
          
          ! Save the new control
          call sptivec_commitVecInPool (rkktsystem%p_rintermedControl%p_rvectorAccess,istep)
          call sptivec_commitVecInPool (rcontrol%p_rvectorAccess,istep)
        
        end do ! istep

      end select
    
    end select    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_calcL2BdCNavSt (rasmTemplates,rphysics,rdualSol,rcontrol,icomp,roptcBDCspace)
  
!<description>
  ! Calculates the L2 boundary control term
  !      nu dn lambda - xi n
  ! from the dual solution.
!</description>
  
!<input>
  ! Assembly template structure
  type (t_staticSpaceAsmTemplates), intent(in) :: rasmTemplates
  
  ! The physics of the problem.
  type(t_settings_physics), intent(in), target :: rphysics
  
  ! Dual solution
  type(t_vectorBlock), intent(in) :: rdualSol
  
  ! Component in the control from which on the boundary control term
  ! should be saved to.
  integer, intent(in) :: icomp
  
  ! Structure defining boundary conditions
  type(t_optcBDCSpace), intent(in) :: roptcBDCspace
!</input>

!<inputoutput>
  ! Control vector. The boundary control term is saved to component
  ! icomp, icomp+1,...
  type(t_vectorBlock), intent(inout) :: rcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    real(DP), dimension(:), allocatable :: p_DparValues
    real(DP), dimension(:), allocatable :: p_DnormalX,p_DnormalY
    real(DP), dimension(:), allocatable :: p_DlambdaX,p_DlambdaY,p_Dxi
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer :: ibct,iseg,iidx1,iidx2,i,NEQ,nvbd,NEQlocal
    real(DP) :: dnu
    
    p_rdiscr => rdualSol%p_rblockDiscr%RspatialDiscr(1)
    
    select case (rphysics%cequation)
    case (CCEQ_NAVIERSTOKES2D,CCEQ_STOKES2D)
      select case (rphysics%cviscoModel)
        case (0)
          ! This is ok.
          dnu = rphysics%dnuConst
        case default
          call output_line("Unsupported equation.",&
              OU_CLASS_ERROR,OU_MODE_STD,"kkt_calcL2BdCNavSt")
          call sys_halt()
      end select
    case default
      call output_line("Unsupported equation.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_calcL2BdCNavSt")
      call sys_halt()
    end select

    ! Unfortunately, this is element dependent...
    
    if (p_rdiscr%inumFESpaces .ne. 1) then
      call output_line("Unsupported discretisation.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_calcL2BdCNavSt")
      call sys_halt()
    end if

    select case (p_rdiscr%RelementDistr(1)%celement)
    case (EL_Q2)
    
      ! Set up a list of points where to evaluate on the boundary.
      ! All vertices and edge midpoints on the boundary.
      
      call storage_getbase_double (&
          p_rdiscr%p_rtriangulation%h_DvertexParameterValue,p_DvertexParameterValue)
      call storage_getbase_double (&
          p_rdiscr%p_rtriangulation%h_DedgeParameterValue,p_DedgeParameterValue)
      call storage_getbase_int (&
          p_rdiscr%p_rtriangulation%h_IboundaryCpIdx,p_IboundaryCpIdx)
          
      ! Length of the destination vector
      NEQ = rcontrol%Rvectorblock(icomp)%NEQ
      
      ! Number of vertices/edges on the boundary
      nvbd = p_rdiscr%p_rtriangulation%NVBD

      ! Temp memory          
      allocate (p_DparValues(NEQ))
      allocate (p_DnormalX(NEQ))
      allocate (p_DnormalY(NEQ))
          
      ! Loop over the boundary components
      do ibct = 1,p_rdiscr%p_rtriangulation%nbct
      
        ! Loop over the segments
        iidx1 = p_IboundaryCpIdx(ibct)
        iidx2 = p_IboundaryCpIdx(ibct+1)-1
        do iseg = iidx1,iidx2
        
          ! Get the parameter values
          p_DparValues(iseg) = p_DvertexParameterValue(iseg)
          p_DparValues(iseg+nvbd) = p_DedgeParameterValue(iseg)
          
        end do
      
        ! Get he normal vectors in these points.
        call boundary_getNormalVec2D (p_rdiscr%p_rboundary, ibct, &
            p_DparValues(iidx1:iidx2),p_DnormalX(iidx1:iidx2), p_DnormalY(iidx1:iidx2))

        call boundary_getNormalVec2D (p_rdiscr%p_rboundary, ibct, &
            p_DparValues(iidx1+nvbd:iidx2+nvbd),p_DnormalX(iidx1+nvbd:iidx2+nvbd), &
            p_DnormalY(iidx1+nvbd:iidx2+nvbd))

      end do
      
      ! Temp memory for intermediate calculations
      allocate (p_DlambdaX(NEQ))
      allocate (p_DlambdaY(NEQ))
      allocate (p_Dxi(NEQ))
      
      ! Loop over the boundary components
      do ibct = 1,p_rdiscr%p_rtriangulation%nbct
      
        ! Get the position of the boundary segment
        iidx1 = p_IboundaryCpIdx(ibct)
        iidx2 = p_IboundaryCpIdx(ibct+1)-1
        NEQlocal = iidx2-iidx1

        ! Evaluare XI
        call fevl_evaluateBdr2d (DER_FUNC, p_Dxi(iidx1:iidx2), rdualSol%RvectorBlock(3), &
            p_DparValues(iidx1:iidx2),ibct,BDR_PAR_01)

        call fevl_evaluateBdr2d (DER_FUNC, p_Dxi(iidx1+nvbd:iidx2+nvbd), rdualSol%RvectorBlock(3), &
            p_DparValues(iidx1+nvbd:iidx2+nvbd),ibct,BDR_PAR_01)

        ! Evaluare D(lambda1)
        call fevl_evaluateBdr2d (DER_DERIV2D_X, p_DlambdaX(iidx1:iidx2), rdualSol%RvectorBlock(1), &
            p_DparValues(iidx1:iidx2),ibct,BDR_PAR_01)

        call fevl_evaluateBdr2d (DER_DERIV2D_Y, p_DlambdaY(iidx1:iidx2), rdualSol%RvectorBlock(1), &
            p_DparValues(iidx1:iidx2),ibct,BDR_PAR_01)

        call fevl_evaluateBdr2d (DER_DERIV2D_X, p_DlambdaX(iidx1+nvbd:iidx2+nvbd), rdualSol%RvectorBlock(1), &
            p_DparValues(iidx1+nvbd:iidx2+nvbd),ibct,BDR_PAR_01)

        call fevl_evaluateBdr2d (DER_DERIV2D_Y, p_DlambdaY(iidx1+nvbd:iidx2+nvbd), rdualSol%RvectorBlock(1), &
            p_DparValues(iidx1+nvbd:iidx2+nvbd),ibct,BDR_PAR_01)
            
        call lsyssc_getbase_double (rcontrol%RvectorBlock(icomp),p_Ddata)
            
        ! Calculate the control in X-direction
        do i=1,NEQlocal
          ! vertices
          p_Ddata(iidx1-1+i) = &
              dnu * p_DlambdaX(i)*p_DnormalX(i)  +  dnu * p_DlambdaY(i)*p_DnormalY(i) - &
              p_Dxi(i)*p_DnormalX(i)
        end do

        do i=nvbd+1,nvbd+NEQlocal
          ! edges
          p_Ddata(iidx1-1+i) = &
              dnu * p_DlambdaX(i)*p_DnormalX(i)  +  dnu * p_DlambdaY(i)*p_DnormalY(i) - &
              p_Dxi(i)*p_DnormalX(i)
        end do
            
        ! Evaluare D(lambda2)
        call fevl_evaluateBdr2d (DER_DERIV2D_X, p_DlambdaX(iidx1:iidx2), rdualSol%RvectorBlock(2), &
            p_DparValues(iidx1:iidx2),ibct,BDR_PAR_01)

        call fevl_evaluateBdr2d (DER_DERIV2D_Y, p_DlambdaY(iidx1:iidx2), rdualSol%RvectorBlock(2), &
            p_DparValues(iidx1:iidx2),ibct,BDR_PAR_01)

        call fevl_evaluateBdr2d (DER_DERIV2D_X, p_DlambdaX(iidx1+nvbd:iidx2+nvbd), rdualSol%RvectorBlock(2), &
            p_DparValues(iidx1+nvbd:iidx2+nvbd),ibct,BDR_PAR_01)

        call fevl_evaluateBdr2d (DER_DERIV2D_Y, p_DlambdaY(iidx1+nvbd:iidx2+nvbd), rdualSol%RvectorBlock(2), &
            p_DparValues(iidx1+nvbd:iidx2+nvbd),ibct,BDR_PAR_01)
            
        call lsyssc_getbase_double (rcontrol%RvectorBlock(icomp+1),p_Ddata)
            
        ! Calculate the control in Y-direction
        do i=1,NEQlocal
          ! vertices
          p_Ddata(iidx1-1+i) = &
              dnu * p_DlambdaX(i)*p_DnormalX(i)  +  dnu * p_DlambdaY(i)*p_DnormalY(i) - &
              p_Dxi(i)*p_DnormalY(i)
        end do

        do i=nvbd+1,nvbd+NEQlocal
          ! edges
          p_Ddata(iidx1-1+i) = &
              dnu * p_DlambdaX(i)*p_DnormalX(i)  +  dnu * p_DlambdaY(i)*p_DnormalY(i) - &
              p_Dxi(i)*p_DnormalY(i)
        end do
      
      end do
        
      deallocate (p_DparValues)
      deallocate (p_DnormalX)
      deallocate (p_DnormalY)
      deallocate (p_DlambdaX)
      deallocate (p_DlambdaY)
      deallocate (p_Dxi)
      
      ! Restrict the computed control to the L2 boundary control region
      call kkt_restrictControlToBDRegion (&
          rcontrol,icomp,p_rdiscr,roptcBDCspace%rdirichletControlBoundary)
      call kkt_restrictControlToBDRegion (&
          rcontrol,icomp+1,p_rdiscr,roptcBDCspace%rdirichletControlBoundary)
        
    case default
    
      call output_line("Unsupported discretisation.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_solvePrimal")
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_restrictControlToBDRegion (rcontrol,icomp,rspatialDiscr,rbdRegionList)
  
!<description>
  ! Restricts a control vector to a set of boundary regions defined
  ! by rbdRegionList. DOFs not in this region are set to zero.
!</description>
  
!<input>
  ! Boundary regino list, the control is to be restricted to.
  type(t_boundaryRegionList), intent(in) :: rbdRegionList
  
  ! Spatial discretisation structure of the 2D space.
  type(t_spatialDiscretisation), intent(in) :: rspatialDiscr
  
  ! Number of the component in the control which is to be restricted.
  integer, intent(in) :: icomp
!</input>

!<inputoutput>
  ! Control vector.
  type(t_vectorBlock), intent(inout) :: rcontrol
!</inputoutput>

!</subroutine>

    type(t_bdRegionEntry), pointer :: p_rbdEntry
    type(t_directAccessIntSet) :: rset
    integer :: h_Idofs
    integer :: i,nvbd
    integer, dimension(:), pointer :: p_Idofs
    integer, dimension(:), pointer :: p_IverticesOnBoundary
    integer, dimension(:), pointer :: p_IedgesOnBoundary
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Allocate an integer set for canceling out entries
    call nsets_initDASet (rset,dof_igetNDofGlob(rspatialDiscr))

    ! Loop over the boundary regions.
    p_rbdEntry => rbdRegionList%p_rbdHead
    do while (associated(p_rbdEntry))
      
      ! Figure out the DOFs there.
      h_Idofs = ST_NOHANDLE
      call bcasm_getDOFsInBDRegion (rspatialDiscr,p_rbdEntry%rboundaryRegion, h_Idofs)
      call storage_getbase_int (h_Idofs,p_Idofs)
      
      ! Mark these DOFs.
      call nsets_putElements (rset,p_Idofs)
      
      ! Release the DOF-list
      call storage_free (h_Idofs)
      
      ! Next
      p_rbdEntry => p_rbdEntry%p_nextBdRegion
      
    end do
    
    ! Quick and dirty implementation...
    
    ! Element type?
    select case (rcontrol%p_rblockDiscr%RspatialDiscr(icomp)%RelementDistr(1)%celement)
    
    case (EL_P2_1D)
      ! Next step: Loop through all vertices/edges on the boundary = DOFs in the
      ! control vector (in that order).

      nvbd = rspatialDiscr%p_rtriangulation%nvbd
      call lsyssc_getbase_double (rcontrol%RvectorBlock(icomp),p_Ddata)

      call storage_getbase_int (rspatialDiscr%p_rtriangulation%h_IverticesAtBoundary,&
          p_IverticesOnBoundary)
      call storage_getbase_int (rspatialDiscr%p_rtriangulation%h_IedgesAtBoundary,&
          p_IedgesOnBoundary)
          
      do i=1,nvbd
      
        ! If this vertex is not in the set, set the corresponding vector entry to zero.
        if (nsets_DASetContains (rset,p_IverticesOnBoundary(i)) .eq. 0) then
          p_Ddata(i) = 0.0_DP
        end if

        ! The DOF of the edge is shifted by nvbd. Check the edge
        if (nsets_DASetContains (rset,&
            rspatialDiscr%p_rtriangulation%NVT+p_IedgesOnBoundary(i)) .eq. 0) then
          p_Ddata(nvbd+i) = 0.0_DP
        end if
      
      end do
    
    
    case default
      call output_line ("Unsupported element.", &
          OU_CLASS_ERROR,OU_MODE_STD,"cc_getDirBCNavSt2D")
      call sys_halt()
    end select
    
    ! Release the set, finish
    call nsets_doneDASet (rset)

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine kkt_projectVector (rmassMatrix,rsource,rdest,rcubatureInfo)
!  
!!<description>
!  ! Projects a scalar vector into a different FEM space.
!!</description>
!  
!!<input>
!  ! Mass matrix in the target space.
!  type(t_matrixScalar), intent(in) :: rmassMatrix
!  
!  ! Source vector
!  type(t_vectorScalar), intent(in) :: rsource
!
!  ! Cubature information structure for applying cubature.
!  type(t_scalarCubatureInfo), intent(in) :: rcubatureInfo
!!</input>
!
!!<inputoutput>
!  ! Destination vector
!  type(t_vectorScalar), intent(inout) :: rdest
!!</inputoutput>
!
!!</subroutine>
!
!    type(t_vectorBlock) :: rsourceTemp,rdestTemp,rrhsTemp,rtemp
!    type(t_matrixBlock) :: rmassTemp
!    type(t_fev2Vectors) :: revalVectors
!    type(t_linsolNode), pointer :: p_rsolverNode, p_rpreconditioner
!    
!    ! Get block vectors/matrices
!    call lsysbl_createMatFromScalar (rmassMatrix,rmassTemp)
!    call lsysbl_createVecFromScalar (rsource,rsourceTemp)
!    call lsysbl_createVecFromScalar (rdest,rdestTemp)
!    
!    ! Create an appropriate RHS.
!    call lsysbl_copyVector (rdestTemp,rrhsTemp)
!    call lsysbl_copyVector (rdestTemp,rtemp)
!    call lsysbl_clearVector (rrhsTemp)
!    
!    call fev2_addVectorToEvalList(revalVectors,rrhsTemp%RvectorBlock(1))
!    
!    call bma_buildVector (rrhsTemp,BMA_CALC_STANDARD,&
!        bma_fcalc_rhsFE, revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
!    
!    ! Apply an L2 projection.
!    call linsol_initJacobi (p_rpreconditioner)
!    call linsol_initDefCorr (p_rsolverNode,p_rpreconditioner)
!    
!    p_rsolverNode%domega = 0.7_DP
!    p_rsolverNode%nmaxIterations = 1000
!    p_rsolverNode%depsRel = 1E-14
!    p_rsolverNode%depsAbs = 1E-12
!    
!    call linsol_setMatrix (p_rsolverNode, rmassTemp)
!    call linsol_initStructure (p_rsolverNode)
!    call linsol_initData (p_rsolverNode)
!    
!    call lsysbl_clearVector (rdest)
!    call linsol_solveAdaptively (p_rsolverNode,rdest,rrhsTemp,rtemp)
!    
!    call linsol_doneData (p_rsolverNode)
!    call linsol_doneStructure (p_rsolverNode)
!    call linsol_done (p_rsolverNode)
!    
!    ! Release
!    call fev2_releaseVectorList (revalVectors)
!    call lsysbl_releaseMatrix (rmassTemp)
!    call lsysbl_releaseVector (rtemp)
!    call lsysbl_releaseVector (rrhsTemp)
!    call lsysbl_releaseVector (rsourceTemp)
!    call lsysbl_releaseVector (rdestTemp)
!
!  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_calcControlRes (rkktsystem,rresidual,dres,iresnorm)
  
!<description>
  ! Calculates the residual of the control equation
!</description>
  
!<input>
  ! Structure defining the KKT system.
  ! The control, primal and dual variable in this structure are used to
  ! calculate the residual.
  type(t_kktsystem), intent(inout), target :: rkktsystem
  
  ! type of norm. A LINALG_NORMxxxx constant.
  integer, intent(in) :: iresnorm
!</input>

!<inputoutput>
  ! Receives the residual in the control space.
  type(t_controlSpace), intent(inout) :: rresidual
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>

    ! The control equation reads
    !
    !   J'(u)  =  u - P ( u - ( lambda + alpha u) )  =  0
    !
    ! with P(.) the projection operator to the admissible space
    ! of controls. In the case of no constraints, this reduces to
    !
    !   J'(u)  =  lambda + alpha u
    !
    ! The residual of the control equation is its negative.
    !
    !   d  =  -J'(u)  =  P ( u - ( lambda + alpha u) )  -  u
    !
    ! The "intermediate" control, i.e. the term in the projection,
    !
    !      u_intermed = u - ( lambda + alpha u )
    !
    ! is already present as in rkktsystem.
    ! Transfer it to rresidual and apply the projection.
    call kkt_dualToControl (rkktsystem,rresidual)
    
    ! Add -u:   rresidual = rresidual - u
    ! Calculate the norm of the residual.
    call kktsp_controlLinearComb (&
        rkktsystem%p_rcontrol,-1.0_DP,rresidual,1.0_DP)
   
    call kkt_controlResidualNorm (&
        rkktsystem%p_roperatorAsmHier%ranalyticData,&
        rresidual,dres,iresnorm)
   
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,rspaceSolver,&
      cspatialInitCondPolicy,ceqnflags,rstatistics)
  
!<description>
  ! Solves the linearised primal equation in the KKT system.
!</description>
  
!<input>
  ! Defines a policy how to generate the initial condition of a timestep.
  ! =0: Always take zero
  ! =1: Propagate the solution of the previous/next timestep to the
  !     current one. (Default)
  ! =2: Take the solution of the last space-time iteration
  integer, intent(in) :: cspatialInitCondPolicy
  
  ! Equation flags that specify modifications to the equation to solve.
  ! One of the SPACESLH_EQNF_xxxx constants.
  integer, intent(in) :: ceqnflags
!</input>

!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv

  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics
!<output>

!</subroutine>

    ! local variables
    type(t_kktSystem), pointer :: p_rkktSystem
    integer :: ierror, idoftime
    type(t_spaceslSolverStat) :: rstatLocal
    type(t_timer) :: rtimer
    
    call stat_startTimer (rstatistics%rtotalTime)
    
    p_rkktSystem => rkktsystemDirDeriv%p_rkktsystem
   
    ! Initialise basic solver structures
    call stat_startTimer (rtimer)
    call spaceslh_initStructure (rspaceSolver, &
        p_rkktsystem%ispacelevel, &
        p_rkktsystem%itimelevel, &
        p_rkktsystem%p_roperatorAsmHier,rstatLocal,ierror)
    call stat_stopTimer (rtimer)

    ! Sum up statistics
    call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)

    if (ierror .ne. 0) then
      call output_line("Error initialising the solver structures.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_solvePrimal")
      call sys_halt()
    end if
    
    ! -----------------------
    ! Loop over all timesteps
    ! -----------------------
    do idofTime = 1,rkktsystemDirDeriv%p_rprimalSolLin%p_rvector%NEQtime
    
      ! Apply the solver to update the solution in timestep idofTime.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call spaceslh_solve (rspaceSolver,idofTime,cspatialInitCondPolicy,&
          ceqnflags,rstatLocal,p_rkktsystem%ispacelevel,&
          p_rkktsystem%p_rprimalSol,&
          rprimalSolLin=rkktsystemDirDeriv%p_rprimalSolLin,&
          rcontrolLin=rkktsystemDirDeriv%p_rcontrolLin)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      ! Sum up statistics
      call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)
      
    end do ! step
   
    call spaceslh_doneStructure (rspaceSolver)

    call stat_stopTimer (rstatistics%rtotalTime)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solveDualDirDeriv (rkktsystemDirDeriv,rspaceSolver,&
      cspatialInitCondPolicy,ceqnflags,rstatistics)
  
!<description>
  ! Solves the linearised dual equation in the KKT system.
!</description>
  
!<input>
  ! Defines a policy how to generate the initial condition of a timestep.
  ! =0: Always take zero
  ! =1: Propagate the solution of the previous/next timestep to the
  !     current one. (Default)
  ! =2: Take the solution of the last space-time iteration
  integer, intent(in) :: cspatialInitCondPolicy

  ! Equation flags that specify modifications to the equation to solve.
  ! One of the SPACESLH_EQNF_xxxx constants.
  integer, intent(in) :: ceqnflags
!</input>

!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv

  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!<output>
  ! Statistics structure
  type(t_spaceslSolverStat), intent(out) :: rstatistics
!<output>

!</subroutine>
   
    ! local variables
    type(t_kktSystem), pointer :: p_rkktSystem
    integer :: ierror, idoftime
    type(t_spaceslSolverStat) :: rstatLocal
    
    call stat_startTimer (rstatistics%rtotalTime)

    p_rkktSystem => rkktsystemDirDeriv%p_rkktsystem
   
    ! Initialise basic solver structures
    call spaceslh_initStructure (rspaceSolver, &
        p_rkktsystem%ispacelevel, &
        p_rkktsystem%itimelevel, &
        p_rkktsystem%p_roperatorAsmHier,rstatLocal,ierror)

    ! Sum up statistics
    call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)

    if (ierror .ne. 0) then
      call output_line("Error initialising the solver structures.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_solvePrimal")
      call sys_halt()
    end if
    
    ! ----------------------------------
    ! Loop over all timesteps, backwards
    ! ----------------------------------
    do idofTime = rkktsystemDirDeriv%p_rprimalSolLin%p_rvector%NEQtime,1,-1
    
      ! Apply the solver to update the solution in timestep idofTime.
      output_iautoOutputIndent = output_iautoOutputIndent + 2
      call spaceslh_solve (rspaceSolver,idofTime,cspatialInitCondPolicy,&
          ceqnflags,rstatLocal,p_rkktsystem%ispacelevel,&
          p_rkktsystem%p_rprimalSol,&
          rdualSol=p_rkktsystem%p_rdualSol,&
          rprimalSolLin=rkktsystemDirDeriv%p_rprimalSolLin,&
          rdualSolLin=rkktsystemDirDeriv%p_rdualSolLin)
      output_iautoOutputIndent = output_iautoOutputIndent - 2
      
      ! Sum up statistics
      call spacesl_sumStatistics(rstatLocal,rstatistics,.false.)
      
    end do ! step
   
    call spaceslh_doneStructure (rspaceSolver)

    call stat_stopTimer (rstatistics%rtotalTime)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_dualToControlDirDeriv (rkktsystemDirDeriv,rcontrolLin)
  
!<description>
  ! From the solution of the linearised primal and dual problem, this routine
  ! calculates the corresponding linearised control
  ! "DP(u_intermed) ( u~ - (lambda~ + alpha u~) )".
!</description>
  
!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv
!</inputoutput>

!<inputoutput>
  ! This receives the control of the linearised control equation.
  type(t_controlSpace), intent(inout) :: rcontrolLin
!</inputoutput>

!</subroutine>
   
    ! local variables
    integer :: icomp,istep
    real(DP) :: dtheta,dwmin,dwmax,dtime
    type(t_vectorBlock), pointer :: p_rdualSpaceLin, p_rcontrolSpaceLin
    type(t_vectorBlock), pointer :: p_rcontrolSpaceLinOutput, p_rintermedControl
    type(t_spaceTimeVector), pointer :: p_rdualSolLin
    type(t_optcBDCSpace) :: roptcBDCspace

    type(t_settings_physics), pointer :: p_rphysics
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl

    type(t_spacetimeOperatorAsm) :: roperatorAsm
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddual,p_Dcontrol,p_DcontrolOut,p_DintermedC

    ! Fetch some structures
    p_rphysics => &
        rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData%p_rphysics
    p_rsettingsOptControl => &
        rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData%p_rsettingsOptControl

    ! Get the underlying space and time discretisation structures.
    call stoh_getOpAsm_slvtlv (roperatorAsm,&
        rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier,&
        rkktsystemDirDeriv%p_rkktsystem%ispacelevel,&
        rkktsystemDirDeriv%p_rkktsystem%itimelevel)

    ! This is strongly equation and problem dependent
    ! and may imply a projection to the admissible set.
    !
    ! We apply a loop over all steps and construct the
    ! control depending on the timestep scheme.
    !
    ! Which timestep scheme do we have?
    
    p_rdualSolLin => rkktsystemDirDeriv%p_rdualSolLin%p_rvector
    
    ! Timestepping technique?
    select case (p_rdualSolLin%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = p_rdualSolLin%p_rtimeDiscr%dtheta
      
      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (p_rdualSolLin%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)
      
        ! Loop over all timesteps.
        do istep = 1,p_rdualSolLin%p_rtimeDiscr%nintervals+1
        
          ! Fetch the dual and control vectors.
          call sptivec_getVectorFromPool (&
              rkktsystemDirDeriv%p_rdualSolLin%p_rvectorAccess,istep,p_rdualSpaceLin)

          call sptivec_getVectorFromPool (&
              rkktsystemDirDeriv%p_rcontrolLin%p_rvectorAccess,istep,p_rcontrolSpaceLin)

          call sptivec_getVectorFromPool (&
              rcontrolLin%p_rvectorAccess,istep,p_rcontrolSpaceLinOutput)

          call sptivec_getVectorFromPool (&
              rkktsystemDirDeriv%p_rkktSystem%p_rintermedControl%p_rvectorAccess,istep,p_rintermedControl)

          ! DEBUG!!!
          call lsysbl_getbase_double (p_rdualSpaceLin,p_Ddual)
          call lsysbl_getbase_double (p_rcontrolSpaceLin,p_Dcontrol)
          call lsysbl_getbase_double (p_rcontrolSpaceLinOutput,p_DcontrolOut)
          call lsysbl_getbase_double (p_rintermedControl,p_DintermedC)

          ! icomp counts the component in the control
          icomp = 0
          
          ! Which equation do we have?
          select case (p_rphysics%cequation)
          
          ! -------------------------------------------------------------
          ! Stokes/Navier Stokes.
          ! -------------------------------------------------------------
          case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
            
            ! Which type of control is applied?
            
            ! -----------------------------------------------------------
            ! Distributed control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

              ! The linearised control equation reads
              !
              !    u~ - DP(u_intermed) ( u~ - (lambda~ + alpha u~) )  =  0
              !
              ! The "linearised intermediate control" reads
              !
              !    u_intermed~ = u~ - (lambda~ + alpha u~)
              !                = (1-alpha) u~ - lambda~
              !
              ! Calculate that.
              !
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpaceLin%RvectorBlock(icomp),p_rdualSpaceLin%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaDistC,-1.0_DP,&
                  p_rcontrolSpaceLinOutput%RvectorBlock(icomp))

              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpaceLin%RvectorBlock(icomp),p_rdualSpaceLin%RvectorBlock(icomp),&
                 1.0_DP-p_rsettingsOptControl%dalphaDistC,-1.0_DP,&
                 p_rcontrolSpaceLinOutput%RvectorBlock(icomp))

              ! The actual linearised control is calculated by
              ! applying an appropriate projection:
              !
              !    u~ = DP(u_intermed) (u_intermed~) 
              !       = DP(u_intermed) ( u~ - (lambda~ + alpha u~) )

              ! ----------------------------------------------------------
              ! Constraints in the distributed control
              ! ----------------------------------------------------------

              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
                  call output_line("Alpha=0 not possible without contraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                  call sys_halt()
                end if
              
              ! ----------------------------------------------------------
              ! Box constraints, implemented by DOF.
              ! ----------------------------------------------------------
              case (1)
              
                icomp = icomp - 2

                ! Create the "restricted" control.
                dwmin = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmin1
                dwmax = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmax1
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceLinOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp),0.0_DP,0.0_DP)

                dwmin = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmin2
                dwmax = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmax2
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceLinOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp),0.0_DP,0.0_DP)

              case default          
                call output_line("Unknown constraints",&
                    OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                call sys_halt()

              end select ! constraints

            end if ! alpha
          
            ! -----------------------------------------------------------
            ! L2 Boundary control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaL2BdC .ge. 0.0_DP) then

              ! No control in the initial solution
              if (istep .gt. 1) then

                ! Characteristics of the current timestep.
                call tdiscr_getTimestep(roperatorasm%p_rtimeDiscrPrimal,istep-1,dtime)

                ! Calculate the region where boundary control is applied
                call sbc_assembleBDconditions (rkktsystemDirDeriv%p_rkktSystem%p_roptcBDC,&
                    roptcBDCSpace,dtime,p_rphysics%cequation,OPTP_PRIMAL,SBC_DIRICHLETBCC,&
                    p_rintermedControl%p_rblockDiscr)

                ! The first two components of the control read
                !
                !    u_intermed~  =  u~ - [ alpha u~ - (nu dn lambda~ - xi~ n) ]
                !                 =  (1-alpha) u~ + nu dn lambda~ - xi~ n
                !
                ! Calculate "nu dn lambda - xi n"
                call kkt_calcL2BdCNavSt (roperatorAsm%p_rasmTemplates,p_rphysics,&
                    p_rdualSpaceLin,p_rcontrolSpaceLinOutput,icomp+1,roptcBDCSpace)
                    
                ! Calculate
                !    u_intermed~ = (1-alpha) u~ + u_intermed
                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rcontrolSpaceLin%RvectorBlock(icomp),p_rcontrolSpaceLinOutput%RvectorBlock(icomp),&
                    1.0_DP-p_rsettingsOptControl%dalphaL2BdC,1.0_DP,&
                    p_rcontrolSpaceLinOutput%RvectorBlock(icomp))

                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rcontrolSpaceLin%RvectorBlock(icomp),p_rcontrolSpaceLinOutput%RvectorBlock(icomp),&
                    1.0_DP-p_rsettingsOptControl%dalphaL2BdC,1.0_DP,&
                    p_rcontrolSpaceLinOutput%RvectorBlock(icomp))
                    
                ! Do we have constraints?
                select case (p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%cconstraints)

                ! ----------------------------------------------------------
                ! No constraints
                ! ----------------------------------------------------------
                case (0)

                  if (p_rsettingsOptControl%dalphaL2BdC .eq. 0.0_DP) then
                    call output_line("Alpha=0 not possible without contraints",&
                        OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                    call sys_halt()
                  end if
                
                ! ----------------------------------------------------------
                ! Box constraints, implemented by DOF.
                ! ----------------------------------------------------------
                case (1)
                
                  icomp = icomp - 2

                  ! Create the "restricted" control.
                  dwmin = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmin1
                  dwmax = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmax1
                  icomp = icomp + 1
                  call nwder_applyMinMaxProjByDof (&
                      p_rcontrolSpaceLinOutput%RvectorBlock(icomp),1.0_DP,&
                      1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                      1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp),0.0_DP,0.0_DP)

                  dwmin = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmin2
                  dwmax = p_rsettingsOptControl%rconstraints%rconstraintsL2BdC%dmax2
                  icomp = icomp + 1
                  call nwder_applyMinMaxProjByDof (&
                      p_rcontrolSpaceLinOutput%RvectorBlock(icomp),1.0_DP,&
                      1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                      1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp),0.0_DP,0.0_DP)

                case default          
                  call output_line("Unknown constraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                  call sys_halt()

                end select ! constraints
                
              end if

            end if ! alphaL2BdC

          ! -------------------------------------------------------------
          ! Heat equation
          ! -------------------------------------------------------------
          case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
            
            ! Which type of control is applied?
            
            ! -----------------------------------------------------------
            ! Distributed control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

              ! The linearised control equation reads
              !
              !    u~ - DP(u_intermed) ( u~ - (lambda~ + alpha u~) )  =  0
              !
              ! The "linearised intermediate control" reads
              !
              !    u_intermed~ = u~ - (lambda~ + alpha u~)
              !                = (1-alpha) u~ - lambda~
              !
              ! Calculate that.
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpaceLin%RvectorBlock(icomp),p_rdualSpaceLin%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaDistC,-1.0_DP,&
                  p_rcontrolSpaceLinOutput%RvectorBlock(icomp))

              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
                  call output_line("Alpha=0 not possible without contraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                  call sys_halt()
                end if
              
              ! ----------------------------------------------------------
              ! Box constraints
              ! ----------------------------------------------------------
              case (1)
              
                icomp = icomp - 1
              
                dwmin = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmin1
                dwmax = p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%dmax1
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceLinOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp),0.0_DP,0.0_DP)

              case default          
                call output_line("Unknown constraints",&
                    OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                call sys_halt()

              end select ! constraints

            end if ! alpha

            ! -----------------------------------------------------------
            ! L2 Boundary control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaL2BdC .ge. 0.0_DP) then

              call output_line("Boundary control not available.",&
                  OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
              call sys_halt()

            end if

          end select ! equation
          
          ! Save the new linearised control
          call sptivec_commitVecInPool (rcontrolLin%p_rvectorAccess,istep)
        
        end do ! istep

      end select
    
    end select    
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_calcControlResDirDeriv (rkktsystemDirDeriv,rrhs,rresidual,&
      dres,iresnorm)
  
!<description>
  ! Calculates the residual of the control equation of the linearised
  ! KKT system  "J''(u) u~ = d".
!</description>
  
!<input>
  ! Structure defining the KKT system.
  ! The control / primal / dual variables in this structure
  ! shall define the value of the functional "J''(u) u~".
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv

  ! The right-hand side "d" of the control equation in the linearised
  ! KKT system "J''(u) u~ = d".
  type(t_controlSpace), intent(inout) :: rrhs

  ! type of norm. A LINALG_NORMxxxx constant.
  integer, intent(in), optional :: iresnorm
!</input>

!<inputoutput>
  ! Receives the residual in the control space.
  type(t_controlSpace), intent(inout) :: rresidual
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out), optional :: dres
!</output>

!</subroutine>

    ! The (linearised) control equation reads:
    !
    !    J''(u) u~  =  u~ - P'(u_intermed) ( u~ - ( lambda~ + alpha u~ ) )  =  rhs
    !
    ! The residual of the control equation is (for the distributed 
    ! control case)
    !
    !   res = rhs - J''(u) u~
    !       = rhs - u~ + P'(u_intermed) ( u~ - ( lambda~ + alpha u~ ) )
    !
    ! The result is written to rresidual, thus, rresidual receives a
    ! fully qualified description of the residual in the control space.
    !
    ! First, add the RHS to the residual.
    ! This is done by creating an appropriate structure.
    !
    ! a) rresidual = P'(u_intermed) ( u~ - ( lambda~ + alpha u~ ) )
    ! We expect rkktsystemDirDeriv to represent u, u_intermed, u~ and lambda~.
    call kkt_dualToControlDirDeriv (rkktsystemDirDeriv,rresidual)

    ! b) rresidual = rresidual + rhs - u~
    call kktsp_controlLinearComb (&
        rrhs,1.0_DP,&
        rkktsystemDirDeriv%p_rcontrolLin,-1.0_DP,&
        rresidual,1.0_DP)
        
    if (present(dres)) then
      call kkt_controlResidualNorm (&
          rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData,&
          rresidual,dres,iresnorm)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_controlResidualNorm (ranalyticData,rcontrolRes,dres,iresnorm)
  
!<description>
  ! Calculates the norm of a residual in the control space.
!</description>
  
!<input>
  ! Structure defining analytic data.
  type(t_spacetimeOpAsmAnalyticData), intent(in), target :: ranalyticData

  ! type of norm. A LINALG_NORMxxxx constant.
  ! Specifies the norm to use for the subcomponents.
  integer, intent(in), optional :: iresnorm
!</input>

!<inputoutput>
  ! Residual in the control space.
  type(t_controlSpace), intent(inout) :: rcontrolRes
!</inputoutput>

!<output>
  ! l2-Norm of the residual. The subcomponents are calculated using the norm
  ! iresnorm.
  real(DP), intent(out), optional :: dres
!</output>

!</subroutine>
   
    ! local variables
    integer :: icomp,istep,itotalcomp
    type(t_vectorBlock), pointer :: p_rcontrolSpace
    real(DP), dimension(:), pointer :: p_Ddata

    type(t_settings_physics), pointer :: p_rphysics
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl

    ! Fetch some structures
    p_rphysics => ranalyticData%p_rphysics
    p_rsettingsOptControl => ranalyticData%p_rsettingsOptControl

    ! This is strongly equation and problem dependent
    ! and may imply a projection to the admissible set.
    !
    ! We apply a loop over all steps and construct the
    ! control depending on the timestep scheme.
    
    dres = 0.0_DP
    itotalcomp = 0
    
    ! Loop over all timesteps.
    do istep = 2,rcontrolRes%p_rvector%p_rtimeDiscr%nintervals+1
    
      ! Get the control vector.
      call sptivec_getVectorFromPool (&
          rcontrolRes%p_rvectorAccess,istep,p_rcontrolSpace)
          
      ! DEBUG!!!
      call lsysbl_getbase_double (p_rcontrolSpace,p_Ddata)
          
      ! icomp counts the component in the control
      icomp = 0
      
      ! Which equation do we have?
      select case (p_rphysics%cequation)
      
      ! -------------------------------------------------------------
      ! Stokes/Navier Stokes.
      ! -------------------------------------------------------------
      case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
        
        ! Which type of control is applied?
        
        ! -----------------------------------------------------------
        ! Distributed control
        ! -----------------------------------------------------------
        if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

          ! Do we have constraints?
          select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

          ! ----------------------------------------------------------
          ! No constraints
          ! ----------------------------------------------------------
          case (0)

            if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
              call output_line("Alpha=0 not possible without contraints",&
                  OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
              call sys_halt()
            end if
          
            ! The first two components of the control read
            !
            !    d = u + 1/alpha lambda
            !
            ! Multiplying with alpha, we get the norm
            !
            !   || alpha d || = || alpha u + lambda ||
            icomp = icomp + 1
            dres = dres + ( & !p_rsettingsOptControl%dalphaDistC * &
                lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

            icomp = icomp + 1
            dres = dres + ( & !p_rsettingsOptControl%dalphaDistC * &
                lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2
                
            itotalcomp = itotalcomp + 2

          ! ----------------------------------------------------------
          ! Box constraints
          ! ----------------------------------------------------------
          case (1)

            ! This scaling can only be done for alpha > 0.
            ! Otherwise, we have to use the non-scaled residual.
            if (p_rsettingsOptControl%dalphaDistC .gt. 0.0_DP) then
              
              ! The first two components of the control read
              !
              !    d = alpha u + lambda
              !
              ! so we get the norm
              !
              !   || alpha d || = || alpha u - alpha P(-1/alpha lambda)) ||
              icomp = icomp + 1
              dres = dres + (& !p_rsettingsOptControl%dalphaDistC * &
                  lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

              icomp = icomp + 1
              dres = dres + (& !p_rsettingsOptControl%dalphaDistC * &
                  lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2
          
            else
              ! The first two components of the control read
              !
              !    d = u + 1/alpha lambda
              !
              ! so we get
              !
              !   || d || = || u + P(1/alpha lambda) ||
              icomp = icomp + 1
              dres = dres + (& !p_rsettingsOptControl%dalphaDistC * &
                  lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

              icomp = icomp + 1
              dres = dres + (& !p_rsettingsOptControl%dalphaDistC * &
                  lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2
            end if
                
            itotalcomp = itotalcomp + 2
                
          case default          
            call output_line("Unknown constraints",&
                OU_CLASS_ERROR,OU_MODE_STD,"kkt_controlResidualNorm")
            call sys_halt()

          end select ! constraints

        end if ! alpha

        ! -----------------------------------------------------------
        ! L2 boundary control
        ! -----------------------------------------------------------
        if (p_rsettingsOptControl%dalphaL2BdC .ge. 0.0_DP) then

          ! ----------------------------------------------------------
          ! No constraints
          ! ----------------------------------------------------------

          if (p_rsettingsOptControl%dalphaL2BdC .eq. 0.0_DP) then
            call output_line("Alpha=0 not possible without contraints",&
                OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
            call sys_halt()
          end if
        
          ! The first two components of the control read
          !
          !    d = alpha u lambda
          !
          ! so we get the norm
          !
          !   || alpha d || = || alpha u + lambda ||
          icomp = icomp + 1
          dres = dres + (& !p_rsettingsOptControl%dalphaL2BdC * &
              lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

          icomp = icomp + 1
          dres = dres + (& !p_rsettingsOptControl%dalphaL2BdC * &
              lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2
              
          itotalcomp = itotalcomp + 2

        end if ! alpha
      
      ! -------------------------------------------------------------
      ! Heat equation
      ! -------------------------------------------------------------
      case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
        
        ! Which type of control is applied?
        
        ! -----------------------------------------------------------
        ! Distributed control
        ! -----------------------------------------------------------
        if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

          ! Do we have constraints?
          select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

          ! ----------------------------------------------------------
          ! No constraints
          ! ----------------------------------------------------------
          case (0)

            if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
              call output_line("Alpha=0 not possible without contraints",&
                  OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
              call sys_halt()
            end if
          
            ! The first two components of the linearised control read
            !
            !    d = u + 1/alpha lambda
            !
            ! Multiuplying with alpha, we get the norm
            !
            !   || alpha d || = || alpha u + lambda ||
            icomp = icomp + 1
            dres = dres + (& !p_rsettingsOptControl%dalphaDistC * &
                lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

            itotalcomp = itotalcomp + 1

          case default          
            call output_line("Unknown constraints",&
                OU_CLASS_ERROR,OU_MODE_STD,"kkt_controlResidualNorm")
            call sys_halt()

          end select ! constraints

        end if ! alpha

      end select ! equation
      
    end do ! istep

    ! Calculate the total l2-norm
    dres = sqrt(dres/real(itotalcomp,DP))
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_applyControlDirDeriv (rkktsystemDirDeriv,rrhs)
  
!<description>
  ! Applies the control equation of the linearised
  ! KKT system:  "d := J''(u) u~"  for a given u~.
!</description>
  
!<input>
  ! Structure defining the KKT system.
  ! The control / primal / dual variables in this structure
  ! shall define the value of the functional "J''(u) u~".
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
!</input>

!<inputoutput>
  ! Receives the result "d" in the control space.
  type(t_controlSpace), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! The (linearised) control equation reads:
    !
    !    J''(u) u~  =  u~ - P'(u_intermed) ( u~ - (lambda~ + alpha u~ ) )
    !
    ! The result is written to rrhs, thus, rrhs receives a
    ! fully qualified description of matrix-vector-product in the control space.
    !
    ! a) rrhs = P'(u_intermed) ( u~ - (lambda~ + alpha u~ ) )
    ! We expect rkktsystemDirDeriv to represent u, u_intermed, u~ and lambda~.
    call kkt_dualToControlDirDeriv (rkktsystemDirDeriv,rrhs)

    ! b) rrhs = -rrhs + u~
    call kktsp_controlLinearComb (&
        rkktsystemDirDeriv%p_rcontrolLin,1.0_DP,&
        rrhs,-1.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_clear (rkktSystem)  
!<description>
  ! Clears a KKT structure.
!</description>

!<inputoutput>
  ! Structure to be cleared.
  type(t_kktsystem), intent(inout), target :: rkktsystem
!</inputoutput>

!</subroutine>

    call kktsp_clearPrimal (rkktsystem%p_rprimalSol)
    call kktsp_clearDual (rkktsystem%p_rdualSol)
    call kktsp_clearControl (rkktsystem%p_rcontrol)
    call kktsp_clearControl (rkktsystem%p_rintermedControl)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_clearDirDeriv (rkktSystemDirDeriv)
!<description>
  ! Clears a KKT structure.
!</description>

!<inputoutput>
  ! Structure to be cleared.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktSystemDirDeriv
!</inputoutput>

!</subroutine>

    call kktsp_clearPrimal (rkktSystemDirDeriv%p_rprimalSolLin)
    call kktsp_clearDual (rkktSystemDirDeriv%p_rdualSolLin)
    call kktsp_clearControl (rkktSystemDirDeriv%p_rcontrolLin)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_getControlAtTime (rkktsystem,dtime,p_rvector)
  
!<description>
  ! Calculates the control at a given point in time.
!</description>

!<input>
  ! Structure defining the KKT system.
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! Point in time where the control should be evaluated
  real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Pointer to a buffer that controls the control at time dtime.
  ! The buffer is taken from the vector buffer in rkktsystem.
  type(t_vectorBlock), pointer :: p_rvector
!</inputoutput>

!</subroutine>

    ! local variables.
    integer :: iindex,icomp
    type(t_vectorBlock), pointer :: p_rvecSource
    real(DP) :: dacttime,dtheta
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    type(t_settings_physics), pointer :: p_rphysics
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl

    ! Fetch some structures
    p_rphysics => &
        rkktsystem%p_roperatorAsmHier%ranalyticData%p_rphysics
    p_rsettingsOptControl => &
        rkktsystem%p_roperatorAsmHier%ranalyticData%p_rsettingsOptControl

    ! Ok, this is a bit tedious. At first allocate and reserve
    ! a vector in the buffer we can use for output.
    iindex = -1
    call sptivec_getFreeBufferFromPool (&
        rkktsystem%p_rcontrol%p_rvectorAccess,iindex,p_rvector)
        
    call sptivec_lockVecInPool (&
        rkktsystem%p_rcontrol%p_rvectorAccess,iindex)
    
    ! Get the time discretisation
    p_rtimeDiscr => rkktsystem%p_rcontrol%p_rvector%p_rtimeDiscr
    
    ! Timestepping technique?
    select case (p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = p_rtimeDiscr%dtheta
      
      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)
      
        ! The time discretisation of the control space is shifted for the
        ! distributed control. At first, depending on dtheta,
        ! calculate the actual time that must be passed to the space-time
        ! vector routines that calculate the control.
        !
        ! The time dtime must be shifted by 1-dtheta to the right.
        call tdiscr_shiftTimestamp(p_rtimeDiscr,dtime,1.0_DP-dtheta,dacttime)
        
        ! Rescale the time into the range [0,1].
        call mprim_linearRescale(dacttime,&
            p_rtimeDiscr%dtimeInit,p_rtimeDiscr%dtimeMax,0.0_DP,1.0_DP,dacttime)
        dacttime = min(max(0.0_DP,dacttime),1.0_DP)
        
        ! Get the control at that time.
        call sptivec_getTimestepDataByTime (&
            rkktsystem%p_rcontrol%p_rvectorAccess, dacttime, p_rvecSource)


        ! icomp counts the component in the control
        icomp = 0
        
        ! Which equation do we have?
        select case (p_rphysics%cequation)
        
        ! -------------------------------------------------------------
        ! Stokes/Navier Stokes.
        ! -------------------------------------------------------------
        case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
          
          ! Which type of control is applied?
          
          ! -----------------------------------------------------------
          ! Distributed control
          ! -----------------------------------------------------------
          if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

            ! Get the source vector with the discributed control.
            

            ! Do we have constraints?
            select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

            ! ----------------------------------------------------------
            ! No constraints, Box constraints
            ! ----------------------------------------------------------
            case (0,1)

              if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
                call output_line("Alpha=0 not possible without contraints",&
                    OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
                call sys_halt()
              end if
            
              ! Copy the distributed control        
              icomp = icomp + 1
              call lsyssc_copyVector (p_rvecSource%RvectorBlock(icomp),p_rvector%RvectorBlock(icomp))

              icomp = icomp + 1
              call lsyssc_copyVector (p_rvecSource%RvectorBlock(icomp),p_rvector%RvectorBlock(icomp))
                  
            case default          
              call output_line("Unknown constraints",&
                  OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
              call sys_halt()
                  
            end select ! constraints

          end if ! alpha

        ! -------------------------------------------------------------
        ! Heat Stokes.
        ! -------------------------------------------------------------
        case (CCEQ_HEAT2D,CCEQ_NL1HEAT2D)
          
          ! Which type of control is applied?
          
          ! -----------------------------------------------------------
          ! Distributed control
          ! -----------------------------------------------------------
          if (p_rsettingsOptControl%dalphaDistC .ge. 0.0_DP) then

            ! Get the source vector with the discributed control.
            

            ! Do we have constraints?
            select case (p_rsettingsOptControl%rconstraints%rconstraintsDistCtrl%cconstraints)

            ! ----------------------------------------------------------
            ! No constraints, Box constraints
            ! ----------------------------------------------------------
            case (0,1)

              if (p_rsettingsOptControl%dalphaDistC .eq. 0.0_DP) then
                call output_line("Alpha=0 not possible without contraints",&
                    OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
                call sys_halt()
              end if
            
              ! Copy the distributed control        
              icomp = icomp + 1
              call lsyssc_copyVector (p_rvecSource%RvectorBlock(icomp),p_rvector%RvectorBlock(icomp))

            case default          
              call output_line("Unknown constraints",&
                  OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
              call sys_halt()
                  
            end select ! constraints

          end if ! alpha

        case default          
          call output_line("Unknown equation",&
              OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
          call sys_halt()

        end select ! equation
        
      case default          
        call output_line("Unknown timestep sub-scheme",&
            OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
        call sys_halt()

      end select ! timestep sub-scheme
      
    case default          
      call output_line("Unknown timestep scheme",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_getControlAtTime")
      call sys_halt()

    end select ! Timstep scheme    

    ! Unlock the vector, finish
    call sptivec_unlockVecInPool (rkktsystem%p_rcontrol%p_rvectorAccess,iindex)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_saveToFiles (rkktSystem,sfilename,cspace)
!<description>
  ! Saves the solutions of a KKT system to file sequences
!</description>

!<input>
  ! The solutions to be saved.
  type(t_kktsystem), intent(inout), target :: rkktsystem
  
  ! Basic filename
  character(len=*), intent(in) :: sfilename

  ! OPTIONAL: Space to be saved. If not specified, all spaces
  ! will be saved.
  integer, intent(in), optional :: cspace
!</input>

!</subroutine>

    logical :: bprimal, bdual, bcontrol, bintermed
    
    bprimal = .false.
    bdual = .false.
    bcontrol = .false.
    bintermed = .false.
    
    if (.not. present(cspace)) then
    
      bprimal = .true.
      bdual = .true.
      bcontrol = .true.

    else
    
      select case (cspace)
      case (CCSPACE_PRIMAL)
        bprimal = .true.
        
      case (CCSPACE_DUAL)
        bdual = .true.
        
      case (CCSPACE_CONTROL)
        bcontrol = .true.

      case (CCSPACE_INTERMEDCONTROL)
        bintermed = .true.
        
      end select
    
    end if

    if (bprimal) &
      call sptivec_saveToFileSequence (&
          rkktsystem%p_rprimalSol%p_rvector,"("""//trim(sfilename)//"_primal."",I5.5)",.true.)
    if (bdual) &
      call sptivec_saveToFileSequence (&
          rkktsystem%p_rdualSol%p_rvector,"("""//trim(sfilename)//"_dual."",I5.5)",.true.)
    if (bcontrol) &
      call sptivec_saveToFileSequence (&
          rkktsystem%p_rcontrol%p_rvector,"("""//trim(sfilename)//"_control."",I5.5)",.true.)
    if (bintermed) &
      call sptivec_saveToFileSequence (&
          rkktsystem%p_rintermedControl%p_rvector,"("""//trim(sfilename)//"_intermedc."",I5.5)",.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_saveDirDerivToFiles (rkktSystemDirDeriv,sfilename,cspace)
!<description>
  ! Saves the directional derivative of a KKT system to file sequences
!</description>

!<input>
  ! The solutions to be saved.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktSystemDirDeriv
  
  ! Basic filename
  character(len=*), intent(in) :: sfilename

  ! OPTIONAL: Space to be saved. If not specified, all spaces
  ! will be saved.
  integer, intent(in), optional :: cspace
!</input>

!</subroutine>

    logical :: bprimal, bdual, bcontrol
    
    bprimal = .false.
    bdual = .false.
    bcontrol = .false.
    
    if (.not. present(cspace)) then
    
      bprimal = .true.
      bdual = .true.
      bcontrol = .true.

    else
    
      select case (cspace)
      case (CCSPACE_PRIMAL)
        bprimal = .true.
        
      case (CCSPACE_DUAL)
        bdual = .true.
        
      case (CCSPACE_CONTROL)
        bcontrol = .true.
        
      end select
    
    end if

    if (bprimal) &
      call sptivec_saveToFileSequence (&
          rkktSystemDirDeriv%p_rprimalSolLin%p_rvector,&
          "("""//trim(sfilename)//"_primallin."",I5.5)",.true.)
    if (bdual) &
      call sptivec_saveToFileSequence (&
          rkktSystemDirDeriv%p_rdualSolLin%p_rvector,&
          "("""//trim(sfilename)//"_duallin."",I5.5)",.true.)
    if (bcontrol) &
      call sptivec_saveToFileSequence (&
          rkktSystemDirDeriv%p_rcontrolLin%p_rvector,&
          "("""//trim(sfilename)//"_controllin."",I5.5)",.true.)

  end subroutine

end module

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
  use genoutput
  
  use mprimitives
  use spatialdiscretisation
  use timediscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  use statistics
  
  use spacetimevectors
  use analyticsolution
  
  use constantsdiscretisation
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use assemblytemplates
  
  use structuresoperatorasm
  use spacematvecassembly
  
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

contains

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_initKKTsystem (rkktsystem,&
      roperatorAsmHier,ispacelevel,itimelevel)
  
!<description>
  ! Initialises a KKT system structure.
!</description>

!<input>  
  ! Parameters for the assembly of space-time operators
  type(t_spacetimeOpAsmHierarchy), intent(in), target :: roperatorAsmHier

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
    real(DP) :: dtheta,dweight,dwmin,dwmax
    type(t_vectorBlock), pointer :: p_rdualSpace, p_rcontrolSpace, p_rintermedControl
    type(t_vectorBlock), pointer :: p_rcontrolSpaceOutput
    type(t_spaceTimeVector), pointer :: p_rdualSol

    type(t_settings_physics), pointer :: p_rphysics
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl

    ! Fetch some structures
    p_rphysics => &
        rkktsystem%p_roperatorAsmHier%ranalyticData%p_rphysics
    p_rsettingsOptControl => &
        rkktsystem%p_roperatorAsmHier%ranalyticData%p_rsettingsOptControl

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
            if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

              ! The first two components of the control read
              !
              !    u_intermed = u - (lambda + alpha u)
              !               = (1-alpha) u - lambda
              !
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpace%RvectorBlock(icomp),p_rdualSpace%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaC,-1.0_DP,&
                  p_rintermedControl%RvectorBlock(icomp))

              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpace%RvectorBlock(icomp),p_rdualSpace%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaC,-1.0_DP,&
                  p_rintermedControl%RvectorBlock(icomp))
                  
              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
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
                
                dwmin = p_rsettingsOptControl%rconstraints%ddistVelUmin1
                dwmax = p_rsettingsOptControl%rconstraints%ddistVelUmax1
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax)

                dwmin = p_rsettingsOptControl%rconstraints%ddistVelUmin2
                dwmax = p_rsettingsOptControl%rconstraints%ddistVelUmax2
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
          
          ! -------------------------------------------------------------
          ! Heat equation
          ! -------------------------------------------------------------
          case (CCEQ_HEAT2D)
            
            ! Which type of control is applied?
            
            ! -----------------------------------------------------------
            ! Distributed control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

              ! The first two components of the control read
              !
              !    u_intermed = u - (lambda + alpha u)
              !               = (1-alpha) u - lambda
              !
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpace%RvectorBlock(icomp),p_rdualSpace%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaC,-1.0_DP,&
                  p_rintermedControl%RvectorBlock(icomp))

              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
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
              
                dwmin = p_rsettingsOptControl%rconstraints%ddistVelUmin1
                dwmax = p_rsettingsOptControl%rconstraints%ddistVelUmax1
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
    ! The 'intermediate' control, i.e. the term in the projection,
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
    real(DP) :: dtheta,dweight,dwmin,dwmax
    type(t_vectorBlock), pointer :: p_rdualSpace
    type(t_vectorBlock), pointer :: p_rdualSpaceLin, p_rcontrolSpaceLin
    type(t_vectorBlock), pointer :: p_rcontrolSpaceLinOutput, p_rintermedControl
    type(t_spaceTimeVector), pointer :: p_rdualSolLin

    type(t_settings_physics), pointer :: p_rphysics
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl

    ! Fetch some structures
    p_rphysics => &
        rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData%p_rphysics
    p_rsettingsOptControl => &
        rkktsystemDirDeriv%p_rkktsystem%p_roperatorAsmHier%ranalyticData%p_rsettingsOptControl

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
            if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

              ! The linearised control equation reads
              !
              !    u~ - DP(u_intermed) ( u~ - (lambda~ + alpha u~) )  =  0
              !
              ! The 'linearised intermediate control' reads
              !
              !    u_intermed~ = u~ - (lambda~ + alpha u~)
              !                = (1-alpha) u~ - lambda~
              !
              ! Calculate that.
              !
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpaceLin%RvectorBlock(icomp),p_rdualSpaceLin%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaC,-1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp))

              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpaceLin%RvectorBlock(icomp),p_rdualSpaceLin%RvectorBlock(icomp),&
                 1.0_DP-p_rsettingsOptControl%dalphaC,-1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp))

              ! The actual linearised control is calculated by
              ! applying an appropriate projection:
              !
              !    u~ = DP(u_intermed) (u_intermed~) 
              !       = DP(u_intermed) ( u~ - (lambda~ + alpha u~) )

              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
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
                dwmin = p_rsettingsOptControl%rconstraints%ddistVelUmin1
                dwmax = p_rsettingsOptControl%rconstraints%ddistVelUmax1
                icomp = icomp + 1
                call nwder_applyMinMaxProjByDof (&
                    p_rcontrolSpaceLinOutput%RvectorBlock(icomp),1.0_DP,&
                    1.0_DP,p_rintermedControl%RvectorBlock(icomp),dwmin,dwmax,&
                    1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp),0.0_DP,0.0_DP)

                dwmin = p_rsettingsOptControl%rconstraints%ddistVelUmin2
                dwmax = p_rsettingsOptControl%rconstraints%ddistVelUmax2
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
          
          ! -------------------------------------------------------------
          ! Heat equation
          ! -------------------------------------------------------------
          case (CCEQ_HEAT2D)
            
            ! Which type of control is applied?
            
            ! -----------------------------------------------------------
            ! Distributed control
            ! -----------------------------------------------------------
            if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

              ! The linearised control equation reads
              !
              !    u~ - DP(u_intermed) ( u~ - (lambda~ + alpha u~) )  =  0
              !
              ! The 'linearised intermediate control' reads
              !
              !    u_intermed~ = u~ - (lambda~ + alpha u~)
              !                = (1-alpha) u~ - lambda~
              !
              ! Calculate that.
              icomp = icomp + 1
              call lsyssc_vectorLinearComb ( &
                  p_rcontrolSpaceLin%RvectorBlock(icomp),p_rdualSpaceLin%RvectorBlock(icomp),&
                  1.0_DP-p_rsettingsOptControl%dalphaC,-1.0_DP,p_rcontrolSpaceLinOutput%RvectorBlock(icomp))

              ! Do we have constraints?
              select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

              ! ----------------------------------------------------------
              ! No constraints
              ! ----------------------------------------------------------
              case (0)

                if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
                  call output_line("Alpha=0 not possible without contraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                  call sys_halt()
                end if
              
              ! ----------------------------------------------------------
              ! Box constraints
              ! ----------------------------------------------------------
              case (1)
              
                icomp = icomp - 1
              
                dwmin = p_rsettingsOptControl%rconstraints%ddistVelUmin1
                dwmax = p_rsettingsOptControl%rconstraints%ddistVelUmax1
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

          end select ! equation
          
          ! Save the new linearised control
          call sptivec_commitVecInPool (rcontrolLin%p_rvectorAccess,istep)
        
        end do ! istep

      end select
    
    end select    
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_calcControlResDirDeriv (rkktsystemDirDeriv,rrhs,rresidual,dres,iresnorm)
  
!<description>
  ! Calculates the residual of the control equation of the linearised
  ! KKT system  "J''(u) g = d".
!</description>
  
!<input>
  ! Structure defining the KKT system.
  ! The control / primal / dual variables in this structure
  ! shall define the value of the functional "J''(u) g".
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv

  ! The right-hand side "d" of the control equation in the linearised
  ! KKT system "J''(u) g = d".
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
    !    J''(u) g  =  g - P'(-1/alpha B lambda) ( -1/alpha [B'(u)]* lambda_g )  =  rhs
    !
    ! The residual of the control equation is (for the distributed 
    ! control case)
    !
    !   res = rhs - J''(u) g
    !       = rhs - g + P'(-1/alpha B lambda) ( -1/alpha [B'(u)]* lambda_g )
    !
    ! The result is written to rresidual, thus, rresidual receives a
    ! fully qualified description of the residual in the control space.
    !
    ! First, add the RHS to the residual.
    ! This is done by creating an appropriate structure.
    !
    ! a) rresidual = P'(-1/alpha B lambda) ( -1/alpha [B'(u)]* lambda_g )
    ! We expect rkktsystemDirDeriv to represent the value "J''(u) g".
    ! To calculate the residual, we need a representation of this value
    ! in the control space.
    call kkt_dualToControlDirDeriv (rkktsystemDirDeriv,rresidual)

    ! b) rresidual = rresidual + rhs - g
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
    do istep = 1,rcontrolRes%p_rvector%p_rtimeDiscr%nintervals+1
    
      ! Get the control vector.
      call sptivec_getVectorFromPool (&
          rcontrolRes%p_rvectorAccess,istep,p_rcontrolSpace)
          
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
        if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

          ! Do we have constraints?
          select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

          ! ----------------------------------------------------------
          ! No constraints
          ! ----------------------------------------------------------
          case (0)

            if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
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
            dres = dres + (p_rsettingsOptControl%dalphaC * &
                lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

            icomp = icomp + 1
            dres = dres + (p_rsettingsOptControl%dalphaC * &
                lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2
                
            itotalcomp = itotalcomp + 2

          ! ----------------------------------------------------------
          ! Box constraints
          ! ----------------------------------------------------------
          case (1)

            ! This scaling can only be done for alpha > 0.
            ! Otherwise, we have to use the non-scaled residual.
            if (p_rsettingsOptControl%dalphaC .gt. 0.0_DP) then
              
              ! The first two components of the control read
              !
              !    d = u + 1/alpha lambda
              !
              ! Multiplying with alpha, we get the norm
              !
              !   || alpha d || = || alpha u - alpha P(-1/alpha lambda)) ||
              icomp = icomp + 1
              dres = dres + (p_rsettingsOptControl%dalphaC * &
                  lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

              icomp = icomp + 1
              dres = dres + (p_rsettingsOptControl%dalphaC * &
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
              dres = dres + (p_rsettingsOptControl%dalphaC * &
                  lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2

              icomp = icomp + 1
              dres = dres + (p_rsettingsOptControl%dalphaC * &
                  lsyssc_vectorNorm(p_rcontrolSpace%RvectorBlock(icomp),iresnorm))**2
            end if
                
            itotalcomp = itotalcomp + 2
                
          case default          
            call output_line("Unknown constraints",&
                OU_CLASS_ERROR,OU_MODE_STD,"kkt_controlResidualNorm")
            call sys_halt()

          end select ! constraints

        end if ! alpha
      
      ! -------------------------------------------------------------
      ! Heat equation
      ! -------------------------------------------------------------
      case (CCEQ_HEAT2D)
        
        ! Which type of control is applied?
        
        ! -----------------------------------------------------------
        ! Distributed control
        ! -----------------------------------------------------------
        if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

          ! Do we have constraints?
          select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

          ! ----------------------------------------------------------
          ! No constraints
          ! ----------------------------------------------------------
          case (0)

            if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
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
            dres = dres + (p_rsettingsOptControl%dalphaC * &
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
  ! KKT system:  "d := J''(u) g"  for a given g
!</description>
  
!<input>
  ! Structure defining the KKT system.
  ! The control / primal / dual variables in this structure
  ! shall define the value of the functional "J''(u) g".
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
!</input>

!<inputoutput>
  ! Receives the result "d" in the control space.
  type(t_controlSpace), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! The (linearised) control equation reads:
    !
    !    J''(u) g  =  g - P'(u) ( -1/alpha [B'(u)]* lambda_g )
    !
    ! The result is written to rrhs, thus, rrhs receives a
    ! fully qualified description of matrix-vector-product in the control space.
    !
    ! a) rrhs = P'(u) ( -1/alpha [B'(u)]* lambda_g )
    ! We expect rkktsystemDirDeriv to represent the value "J''(u) g".
    ! To calculate the residual, we need a representation of this value
    ! in the control space.
    call kkt_dualToControlDirDeriv (rkktsystemDirDeriv,rrhs)

    ! b) rrhs = -rrhs + g
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
          if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

            ! Get the source vector with the discributed control.
            

            ! Do we have constraints?
            select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

            ! ----------------------------------------------------------
            ! No constraints, Box constraints
            ! ----------------------------------------------------------
            case (0,1)

              if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
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
        case (CCEQ_HEAT2D)
          
          ! Which type of control is applied?
          
          ! -----------------------------------------------------------
          ! Distributed control
          ! -----------------------------------------------------------
          if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then

            ! Get the source vector with the discributed control.
            

            ! Do we have constraints?
            select case (p_rsettingsOptControl%rconstraints%cdistVelConstraints)

            ! ----------------------------------------------------------
            ! No constraints, Box constraints
            ! ----------------------------------------------------------
            case (0,1)

              if (p_rsettingsOptControl%dalphaC .eq. 0.0_DP) then
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

end module

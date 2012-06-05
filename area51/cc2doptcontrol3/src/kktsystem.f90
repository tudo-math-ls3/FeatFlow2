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
  
  use spacetimevectors
  use analyticsolution
  
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use assemblytemplates
  
  use structuresoperatorasm
  use spacematvecassembly
  
  use kktsystemspaces
  use spacesolver
  
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
    type(t_dualSpace), pointer :: p_rdualLinSol => null()

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

    allocate (rkktsystemDirDeriv%p_rdualLinSol)
    call kktsp_initDualVector (rkktsystemDirDeriv%p_rdualLinSol,&
        rkktsystem%p_rprimalSol%p_rvector%p_rspaceDiscr,&
        rkktsystem%p_rprimalSol%p_rvector%p_rtimeDiscr)

    allocate (rkktsystemDirDeriv%p_rcontrolLin)
    call kktsp_initControlVector (rkktsystemDirDeriv%p_rcontrolLin,&
        rkktsystem%p_rprimalSol%p_rvector%p_rspaceDiscr,&
        rkktsystem%p_rprimalSol%p_rvector%p_rtimeDiscr)
   
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
  type(t_kktsystemDirDeriv), intent(out), target :: rkktsystemDirDeriv
!</inputoutput>

!</subroutine>

    ! Clean up the structures
    nullify(rkktsystemDirDeriv%p_rkktsystem)
    
    ! Release memory
    call kktsp_donePrimalVector (rkktsystemDirDeriv%p_rprimalSolLin)
    deallocate (rkktsystemDirDeriv%p_rprimalSolLin)

    call kktsp_doneDualVector (rkktsystemDirDeriv%p_rdualLinSol)
    deallocate (rkktsystemDirDeriv%p_rdualLinSol)

    call kktsp_doneControlVector (rkktsystemDirDeriv%p_rcontrolLin)
    deallocate (rkktsystemDirDeriv%p_rcontrolLin)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solvePrimal (rkktsystem,rspaceSolver)
  
!<description>
  ! Solves the primal equation in the KKT system based on the control in the
  ! rkktsystem structure.
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: idofTime, ierror
    
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
        rkktsystem%p_roperatorAsmHier,ierror)

    if (ierror .ne. 0) then
      call output_line("Error initialising the solver structures.",&
          OU_CLASS_ERROR,OU_MODE_STD,"kkt_solvePrimal")
      call sys_halt()
    end if
    
    ! -----------------------
    ! Loop over all timesteps
    ! -----------------------
    do idofTime = 1,rkktsystem%p_rprimalSol%p_rvector%NEQtime
    
      ! Apply the solver to update the solution in timestep idofTime.
      call spaceslh_solve (rspaceSolver,idofTime,rkktsystem%ispacelevel,&
          rkktsystem%p_rprimalSol,rcontrol=rkktsystem%p_rcontrol)
      
    end do ! step
   
    call spaceslh_doneStructure (rspaceSolver)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solveDual (rkktsystem,rspaceSolver)
  
!<description>
  ! Solves the dual equation in the KKT system based on the control and the
  ! primal solution in the rkktsystem structure.
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem

  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: idofTime, ierror

    ! -------------------------------------------------------------------------
    ! The solution of the dual equation is rather similar to the primal
    ! equation, but even simpler. The timeloop loops backward through the
    ! timesteps, and in every timestep, a linear problem has to be solved.
    ! -------------------------------------------------------------------------
    
    ! Initialise basic solver structures
    call spaceslh_initStructure (rspaceSolver, &
        rkktsystem%ispacelevel, rkktsystem%itimelevel, &
        rkktsystem%p_roperatorAsmHier,ierror)

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
      call spaceslh_solve (rspaceSolver,idofTime,rkktsystem%ispacelevel,&
          rkktsystem%p_rprimalSol,rdualSol=rkktsystem%p_rdualSol)
      
    end do ! step
   
    call spaceslh_doneStructure (rspaceSolver)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_dualToControl (rkktsystem,rcontrol)
  
!<description>
  ! From the solution of the primal and dual problem, this routine
  ! calculates the corresponding control.
!</description>
  
!<input>
  ! Structure defining the KKT system.
  type(t_kktsystem), intent(inout), target :: rkktsystem
!</input>

!<inputoutput>
  ! Receives the control.
  type(t_controlSpace), intent(inout), target :: rcontrol
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: icomp,istep
    real(DP) :: dtheta
    type(t_vectorBlock), pointer :: p_rdualSpace, p_rcontrolSpace
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
              rcontrol%p_rvectorAccess,istep,p_rcontrolSpace)
              
          ! icomp counts the component in the control
          icomp = 0
          
          ! Which equation do we have?
          select case (p_rphysics%cequation)
          
          ! -------------------------------------------------------------
          ! Stokes/Navier Stokes.
          ! -------------------------------------------------------------
          case (0,1)
            
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
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControl")
                  call sys_halt()
                end if

                ! The first two components of the control read
                !
                !    u = -1/alpha lambda
                !
                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rdualSpace%RvectorBlock(icomp),p_rcontrolSpace%RvectorBlock(icomp),&
                    -1.0_DP/p_rsettingsOptControl%dalphaC,0.0_DP)

                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rdualSpace%RvectorBlock(icomp),p_rcontrolSpace%RvectorBlock(icomp),&
                    -1.0_DP/p_rsettingsOptControl%dalphaC,0.0_DP)
              
              end select ! constraints

            end if ! alpha
          
          end select ! equation
          
          ! Save the new control
          call sptivec_commitVecInPool (rcontrol%p_rvectorAccess,istep)
        
        end do ! istep

      end select
    
    end select    
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_calcControlRes (rkktsystem,rresidual,dres)
  
!<description>
  ! Calculates the residual of the control equation
!</description>
  
!<input>
  ! Structure defining the KKT system.
  ! The control, primal and dual variable in this structure are used to
  ! calculate the residual.
  type(t_kktsystem), intent(inout), target :: rkktsystem
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
    !   J'(u)  =  u - P ( -1/alpha [B'(u)]* lambda )  =  0
    !
    ! with P(.) the projection operator to the admissible space
    ! of controls.
    !
    ! The residual of the control equation is its negative.
    !
    !   d  =  -J'(u)  =  -u + P ( -1/alpha [B'(u)]* lambda ) 
    !
    ! We call kkt_dualToControl to calculate a new control 
    !
    !      rresidual = P ( -1/alpha [B'(u)]* lambda )
    !
    ! from the primal/dual variables in rkktsystem.
    call kkt_dualToControl (rkktsystem,rresidual)
    
    ! Add -u:   rresidual = rresidual - u
    call kktsp_controlLinearComb (rkktsystem%p_rcontrol,-1.0_DP,rresidual,1.0_DP)
    
    ! Calculate the norm of the residual
    dres = sptivec_vectorNorm (rresidual%p_rvector,LINALG_NORML2)
   
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine kkt_solvePrimalDirDeriv (rkktsystemDirDeriv,rspaceSolver)
  
!<description>
  ! Solves the linearised primal equation in the KKT system.
!</description>
  
!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv

  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!</subroutine>
   
    ! ... to be done
    call sys_halt()

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solveDualDirDeriv (rkktsystemDirDeriv,rspaceSolver)
  
!<description>
  ! Solves the linearised dual equation in the KKT system.
!</description>
  
!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv

  ! Space solver structure used for solving subequations in space
  type(t_spaceSolverHierarchy), intent(inout) :: rspaceSolver
!</inputoutput>

!</subroutine>
   
    ! ... to be done
    call sys_halt()

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_dualToControlDirDeriv (rkktsystemDirDeriv,rcontrolLin)
  
!<description>
  ! From the solution of the linearised primal and dual problem, this routine
  ! calculates the corresponding linearised control.
!</description>
  
!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv
!</inputoutput>

!<inputoutput>
  ! OPTIONAL: If specified, this receives the control.
  ! If not specified, the control in rkktsystem is overwritten.
  type(t_controlSpace), intent(inout) :: rcontrolLin
!</inputoutput>

!</subroutine>
   
    ! local variables
    integer :: icomp,istep
    real(DP) :: dtheta
    type(t_vectorBlock), pointer :: p_rdualSpace, p_rcontrolSpace
    type(t_spaceTimeVector), pointer :: p_rdualLinSol

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
    
    p_rdualLinSol => rkktsystemDirDeriv%p_rdualLinSol%p_rvector
    
    ! Timestepping technique?
    select case (p_rdualLinSol%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = p_rdualLinSol%p_rtimeDiscr%dtheta
      
      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (p_rdualLinSol%p_rtimeDiscr%itag)
      
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
        do istep = 1,p_rdualLinSol%p_rtimeDiscr%nintervals+1
        
          ! Fetch the dual and control vectors.
          call sptivec_getVectorFromPool (&
              rkktsystemDirDeriv%p_rdualLinSol%p_rvectorAccess,istep,p_rdualSpace)

          call sptivec_getVectorFromPool (&
              rcontrolLin%p_rvectorAccess,istep,p_rcontrolSpace)
              
          ! icomp counts the component in the control
          icomp = 0
          
          ! Which equation do we have?
          select case (p_rphysics%cequation)
          
          ! -------------------------------------------------------------
          ! Stokes/Navier Stokes.
          ! -------------------------------------------------------------
          case (0,1)
            
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

                if (p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then
                  call output_line("Alpha=0 not possible without contraints",&
                      OU_CLASS_ERROR,OU_MODE_STD,"kkt_dualToControlDirDeriv")
                  call sys_halt()
                end if
              
                ! The first two components of the linearised control read
                !
                !    u~ = -1/alpha lambda~
                !
                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rdualSpace%RvectorBlock(icomp),p_rcontrolSpace%RvectorBlock(icomp),&
                    -1.0_DP/p_rsettingsOptControl%dalphaC,0.0_DP)

                icomp = icomp + 1
                call lsyssc_vectorLinearComb ( &
                    p_rdualSpace%RvectorBlock(icomp),p_rcontrolSpace%RvectorBlock(icomp),&
                    -1.0_DP/p_rsettingsOptControl%dalphaC,0.0_DP)
              
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

  subroutine kkt_calcControlResDirDeriv (rkktsystemDirDeriv,rrhs,rresidual,dres)
  
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

    ! The (linearised) control equation reads:
    !
    !    J''(u) g  =  g - (-1/alpha) P'(u) ( -1/alpha [B'(u)]* lambda_g )  =  d
    !
    ! The residual of the control equation is (for the distributed 
    ! control case)
    !
    !   res = d - J''(u) g
    !       = d - g - 1/alpha P'(u) ( -1/alpha [B'(u)]* lambda_g )
    !
    ! The result is written to rresidual, thus, rresidual receives a
    ! fully qualified description of the residual in the control space.
    !
    ! First, add the RHS to the residual.
    ! This is done by creating an appropriate structure.
    !
    ! a) rresidual = -1/alpha P'(u) ( -1/alpha [B'(u)]* lambda_g )
    ! We expect rkktsystemDirDeriv to represent the value "J''(u) g".
    ! To calculate the residual, we need a representation of this value
    ! in the control space.
    call kkt_dualToControlDirDeriv (rkktsystemDirDeriv,rresidual)

    ! b) rresidual = rresidual + d - g
    call kktsp_controlLinearComb (&
        rrhs,1.0_DP,&
        rkktsystemDirDeriv%p_rcontrolLin,-1.0_DP,&
        rresidual,1.0_DP,dres)

  end subroutine

end module

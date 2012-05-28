!##############################################################################
!# ****************************************************************************
!# <name> newtoniterationlinear </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises the linear subsolver of the Newton iteration in the 
!# control space of the optimal control problem.
!# </purpose>
!##############################################################################

module newtoniterationlinear

  use fsystem
  use genoutput
  use paramlist
  
  use spatialdiscretisation
  use timediscretisation
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
  
  use spacematvecassembly
  
  use kktsystemspaces
  use kktsystem
  
  private

!<type>

  !<typeblock>
  
  ! Contains the parameters of the Newton iteration.
  type t_linsolParameters
  
    ! OUTPUT: Result
    ! The result of the solution process.
    ! =-1: convergence criterion not reached.
    ! =0: success.
    ! =1: iteration broke down, diverging.
    ! =2: iteration broke down, preconditioner did not work.
    ! =3: error in the parameters.
    integer :: iresult = 0

    ! OUTPUT: Number of performed iterations.
    integer :: iiterations = 0
    
    ! OUTPUT: Initial residual in the control space.
    real(DP) :: dresInit = 0.0_DP

    ! OUTPUT: Final residual in the control space.
    real(DP) :: dresFinal = 0.0_DP

    ! Relative stopping criterion. Stop iteration if
    ! !!defect!! < EPSREL * !!initial defect!!.
    ! =0: ignore, use absolute stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsRel = 1E-5_DP

    ! Absolute stopping criterion. Stop iteration if
    ! !!defect!! < EPSABS.
    ! =0: ignore, use relative stopping criterion; standard = 1E-5
    ! Remark: do not set depsAbs=depsRel=0!
    real(DP) :: depsAbs = 1E-5_DP

    ! Relative divergence criterion.  Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL * !!initial defect!!
    ! A value of SYS_INFINITY_DP disables the relative divergence check.
    ! standard = 1E3
    real(DP) :: ddivRel = 1E6_DP

    ! Absolute divergence criterion. Treat iteration as
    ! diverged if
    !   !!defect!! >= DIVREL
    ! A value of SYS_INFINITY_DP disables the absolute divergence check.
    ! standard = SYS_INFINITY_DP
    real(DP) :: ddivAbs = SYS_INFINITY_DP

    ! Minimum number of iterations top perform
    integer :: nminIterations = 1

    ! Maximum number of iterations top perform
    integer :: nmaxIterations = 50
    
    ! General relaxation parameter
    real(DP) :: domega = 1.0_DP
  
    ! Temporary memory for calculations
    type(t_controlSpace), pointer :: p_rtempVector => null()

  end type
  
  !</typeblock>
  
  public :: t_linsolParameters

!</types>

  ! Apply preconditioning of a defect with the Newton preconditioner
  ! in the control space.
  public :: newtonlin_precondNewton
  
  ! Basic initialisation of the Newton solver
  public :: newtonlin_init
  
  ! Structural initialisation
  public :: newtonlin_initStructure
  
  ! Final initialisation
  public :: newtonlin_initData
  
  ! Cleanup of data
  public :: newtonlin_doneData

  ! Cleanup of structures
  public :: newtonlin_doneStructure

  ! Final cleanup
  public :: newtonlin_done

contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_getResidual (rkktsystemDirDeriv,rrhs,rresidual,dres)
  
!<description>
  ! Calculates the residual in the linearised control equation.
!</description>
  
!<inputoutput>
  ! Structure defining the linearised KKT system.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
  
  ! On output, this structure receives a representation of the search
  ! direction / residual in the Newton iteration
  type(t_controlSpace), intent(inout) :: rresidual
!</inputoutput>

!<output>
  ! L2-Norm of the residual
  real(DP), intent(out) :: dres
!</output>

!</subroutine>

    ! -------------------------------------------------------------
    ! Step 1: Solve the primal and dual system.
    ! -------------------------------------------------------------

    ! Solve the primal equation, update the primal solution.
    call kkt_solvePrimalDirDeriv (rkktsystemDirDeriv)
    
    ! Solve the dual equation, update the dual solution.
    call kkt_solveDualDirDeriv (rkktsystemDirDeriv)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the residual
    ! -------------------------------------------------------------

    ! Take the solution of the linearised primal/dual system and
    ! calculate the residual in the control space.
    call kkt_calcControlResDirDeriv (rkktsystemDirDeriv,rrhs,rresidual,dres)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_richardson (rlinsolParam,rkktsystemDirDeriv,rrhs)
  
!<description>
  ! Applies a Richardson iteration in the control space
  ! for the linearised KKT system.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Structure defining the linearised KKT system.
  ! The linearised control in this structure is used as
  ! initial and target vector.
  type(t_kktsystemDirDeriv), intent(inout), target :: rkktsystemDirDeriv
  
  ! Right-hand side of the linearised control equation
  type(t_controlSpace), intent(inout), target :: rrhs
!</inputoutput>

!</subroutine>

    ! Richardson is a do loop...
    rlinsolParam%iiterations = 0

    do while (.true.)
    
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      ! Compute the residual and its norm.
      call newtonlin_getResidual (rkktsystemDirDeriv,rrhs,&
          rlinsolParam%p_rtempVector,rlinsolParam%dresFinal)

      if (rlinsolParam%iiterations .eq. 1) then
        ! Remember the initial residual
        rlinsolParam%dresInit = rlinsolParam%dresFinal
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (rlinsolParam%iiterations .ge. rlinsolParam%nminIterations) then
        ! Check the residual.
        !
        ! Absolute residual
        if (rlinsolParam%dresFinal .le. rlinsolParam%depsAbs) then
          rlinsolParam%iresult = 0
          exit
        end if
        
        ! Relative residual
        if (rlinsolParam%dresFinal .le. rlinsolParam%depsRel * rlinsolParam%dresInit) then
          rlinsolParam%iresult = 0
          exit
        end if
      end if
      
      ! -------------------------------------------------------------
      ! Check for divergence
      ! -------------------------------------------------------------
      ! Absolute residual
      if (rlinsolParam%dresFinal .ge. rlinsolParam%ddivAbs) then
        rlinsolParam%iresult = 1
        exit
      end if
      
      ! Relative residual
      if (rlinsolParam%dresFinal .ge. rlinsolParam%ddivRel * rlinsolParam%dresInit) then
        rlinsolParam%iresult = 1
        exit
      end if
      
      ! -------------------------------------------------------------
      ! Check other stopping criteria
      ! -------------------------------------------------------------
      if (rlinsolParam%iiterations .ge. rlinsolParam%nmaxIterations) then
        ! Maximum number of iterations reached.
        rlinsolParam%iresult = -1
        exit
      end if
      
      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Add the residual (damped) to the current control.
      ! This updates the current control.
      call kktsp_controlLinearComb (rlinsolParam%p_rtempVector,rlinsolParam%domega,&
          rkktsystemDirDeriv%p_rcontrolLin,1.0_DP)
    
      ! -------------------------------------------------------------
      ! Proceed with the iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rlinsolParam%iiterations = rlinsolParam%iiterations + 1
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_precondNewton (rlinsolParam,rkktsystemDirDeriv,rnewtonDir)
  
!<description>
  ! Calculates the Newton search direction by applying the Newton
  ! preconditioner to the given residual rnewtonDir.
  ! This calculates a correction g with
  !      J''(u) g = d
  ! where d=rnewtonDir. On return of this routine, there is rnewtonDir=g.
!</description>

!<inputoutput>
  ! Parameters for the iteration.
  ! The output parameters are changed according to the iteration.
  type(t_linsolParameters), intent(inout) :: rlinsolParam

  ! Structure defining directional derivatives of the KKT system.
  ! Used for temporary calculations.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv 

  ! This structure contains the search direction to be
  ! preconditioned. Will be replaced by the precoditioned search direction.
  type(t_controlSpace), intent(inout) :: rnewtonDir
!</inputoutput>

!</subroutine>

    ! For the calculation of the Newon search direction, we have to solve
    ! the linear system
    !
    !    J''(u) g = d
    !
    ! On the entry of this routine, there is d=rnewtonDir. When the routine
    ! is left, rnewtonDir receives g.
    !
    ! For this purpose, we need a linear solver based on J''(u) in the
    ! control space. For distributed control, this is multigrid-based,
    ! for boundary control (or general control with a limited number of
    ! control variables) a single-grid solver.
    !
    ! The system to solve is rather similar to the nonlinear system,
    ! but involves a couple of additional terms. It is still an iteration
    ! in u and thus, involves a forward and a backward solve for every 
    ! update -- probably on different refinement levels.
    !
    ! In the following, rnewtonDir is used as right-hand side of the
    ! equation. The solution is calculated in 
    !     rkktsystemDirDeriv%p_rcontrolLin
    ! which starts with zero.
    
    call kktsp_clearPrimal (rkktsystemDirDeriv%p_rprimalSolLin)
    call kktsp_clearDual (rkktsystemDirDeriv%p_rdualSolLin)
    call kktsp_clearControl (rkktsystemDirDeriv%p_rcontrolLin)
    
    ! Call the Richardson iteration to calculate an update
    ! for the control.
    call newtonlin_richardson (rlinsolParam,rkktsystemDirDeriv,rnewtonDir)
    
    ! Overwrite the rnewtonDir with the update.
    call kktsp_controlLinearComb (&
        rkktsystemDirDeriv%p_rcontrolLin,1.0_DP,rnewtonDir,0.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_init (rlinsolParam,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the preconditioner.
!</description>

!<input>
  ! Parameter list that contains the data for the preconditioning.
  type(t_parlist), intent(in) :: rparlist
  
  ! Entry Section in the parameter list containing the data of the 
  ! preconditioner.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(out) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initStructure (rlinsolParam)
  
!<description>
  ! Structural initialisation of the preconditioner
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_initData (rlinsolParam)
  
!<description>
  ! Final preparation of the preconditioner for preconditioning.
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneData (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_doneStructure (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonlin_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonlin_done (rlinsolParam)
  
!<description>
  ! Clean up a preconditioner.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_linsolParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

end module

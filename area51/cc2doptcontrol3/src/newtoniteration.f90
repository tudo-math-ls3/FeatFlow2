!##############################################################################
!# ****************************************************************************
!# <name> newtoniteration </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a (semismooth) Newton iteration in the control 
!# space of the optimal control problem.
!# </purpose>
!##############################################################################

module newtoniteration

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
  
  use newtoniterationlinear
  
  private
  
!<type>

  !<typeblock>
  
  ! Contains the parameters of the Newton iteration.
  type t_newtonParameters
  
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
  
    ! Damping parameter for the Newton iteration
    real(DP) :: domega = 1.0_DP
  
    ! Parameters for the linear subsolver.
    type(t_linsolParameters) :: rlinsolParam

  end type
  
  !</typeblock>
  
  public :: t_newtonParameters

!</types>

  ! Apply a Newton iteration in the control space
  public :: newtonit_solve
  
  ! Basic initialisation of the Newton solver
  public :: newtonit_init
  
  ! Structural initialisation
  public :: newtonit_initStructure
  
  ! Final initialisation
  public :: newtonit_initData
  
  ! Cleanup of data
  public :: newtonit_doneData

  ! Cleanup of structures
  public :: newtonit_doneStructure

  ! Final cleanup
  public :: newtonit_done
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_getResidual (rkktsystem,rresidual,dres)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system. The control in this structure
  ! defines the current 'state' of the Newton algorithm.
  ! On output, the 
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! On output, this structure receives a representation of the search
  ! direction / residual in the Newton iteration.
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
    call kkt_solvePrimal (rkktsystem)
    
    ! Solve the dual equation, update the dual solution.
    call kkt_solveDual (rkktsystem)

    ! -------------------------------------------------------------
    ! Step 2: Calculate the search direction   
    ! -------------------------------------------------------------

    ! The search direction is just the residual in the control equation.
    call kkt_calcControlRes (rkktsystem,rresidual,dres)

  end subroutine

    ! ***************************************************************************

!<subroutine>

  subroutine newtonit_updateControl (rnewtonParam,rkktsystem,rcorrection)
  
!<description>
  ! Calculates the basic (unprecondiotioned) search direction of the 
  ! Newton algorithm.
!</description>  

!<input>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_newtonParameters), intent(in) :: rnewtonParam

  ! The preconditioned search direction
  type(t_controlSpace), intent(inout) :: rcorrection
!</input>
  
!<inputoutput>
  ! Structure defining the KKT system. The control in this structure
  ! is updated according to the search direction.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>

    ! Currectly, this is just a linear combination of the control variables.
    !
    !    u_n+1  =  u_n + g_n
    !
    ! Later, a step length control can be added here.
    
    call kktsp_controlLinearComb (&
        rcorrection,rnewtonParam%domega,rkktsystem%p_rcontrol,1.0_DP)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_solve (rnewtonParam,rkktsystem)
  
!<inputoutput>
  ! Parameters for the Newton iteration.
  ! The output parameters are changed according to the iteration.
  type(t_newtonParameters), intent(inout) :: rnewtonParam

  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>
   
    ! local variables
    type(t_kktsystemDirDeriv) :: rkktsystemDirDeriv
    type(t_controlSpace) :: rdescentDir
    
    ! Prepare a structure that encapsules the directional derivative.
    
    ! Apply the Newton iteration
    rnewtonParam%iiterations = 0
    
    do while (.true.)
    
      ! The Newton iteration reads
      !
      !    u_n+1  =  u_n  -  [J''(u_n)]^-1  J'(u_n)
      !
      ! or in other words,
      !
      !    u_n+1  =  u_n  +  [J''(u_n)]^-1  d_n
      !
      ! with the residual   d_n = -J'(u_n)   specifying a 'search direction'.
      
      ! -------------------------------------------------------------
      ! Get the current residual / search direction
      ! -------------------------------------------------------------

      ! Compute the basic (unpreconditioned) search direction d_n.
      call newtonit_getResidual (rkktsystem,rdescentDir,rnewtonParam%dresFinal)

      if (rnewtonParam%iiterations .eq. 1) then
        ! Remember the initial residual
        rnewtonParam%dresInit = rnewtonParam%dresFinal
      end if

      ! -------------------------------------------------------------
      ! Check for convergence
      ! -------------------------------------------------------------
      if (rnewtonParam%iiterations .ge. rnewtonParam%nminIterations) then
        ! Check the residual.
        !
        ! Absolute residual
        if (rnewtonParam%dresFinal .le. rnewtonParam%depsAbs) then
          rnewtonParam%iresult = 0
          exit
        end if
        
        ! Relative residual
        if (rnewtonParam%dresFinal .le. rnewtonParam%depsRel * rnewtonParam%dresInit) then
          rnewtonParam%iresult = 0
          exit
        end if
      end if
      
      ! -------------------------------------------------------------
      ! Check for divergence
      ! -------------------------------------------------------------
      ! Absolute residual
      if (rnewtonParam%dresFinal .ge. rnewtonParam%ddivAbs) then
        rnewtonParam%iresult = 1
        exit
      end if
      
      ! Relative residual
      if (rnewtonParam%dresFinal .ge. rnewtonParam%ddivRel * rnewtonParam%dresInit) then
        rnewtonParam%iresult = 1
        exit
      end if
      
      ! -------------------------------------------------------------
      ! Check other stopping criteria
      ! -------------------------------------------------------------

      if (rnewtonParam%iiterations .ge. rnewtonParam%nmaxIterations) then
        ! Maximum number of iterations reached.
        rnewtonParam%iresult = -1
        exit
      end if
      
      ! -------------------------------------------------------------
      ! Preconditioning with the Newton matrix
      ! -------------------------------------------------------------
      
      ! Actual Newton iteration. Apply the Newton preconditioner
      ! to get the Newton search direction:
      !
      !    J''(u_n) g_n  =  d_n
      !
      call newtonlin_precondNewton (rnewtonParam%rlinsolParam,&
          rkktsystemDirDeriv,rdescentDir)
      
      ! -------------------------------------------------------------
      ! Update of the solution
      ! -------------------------------------------------------------

      ! Update the control according to the search direction:
      !
      !    u_n+1  =  u_n  +  g_n
      !
      ! or to any configured step-length control rule.
      call newtonit_updateControl (rnewtonParam,rkktsystem,rdescentDir)
      
      ! -------------------------------------------------------------
      ! Proceed with the next iteration
      ! -------------------------------------------------------------
      ! Next iteration
      rnewtonParam%iiterations = rnewtonParam%iiterations + 1
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_init (rlinsolParam,rparlist,ssection)
  
!<description>
  ! Basic initialisation of the Newton solver.
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
  type(t_newtonParameters), intent(out) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initStructure (rlinsolParam)
  
!<description>
  ! Structural initialisation of the Newton solver.
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_initData (rlinsolParam)
  
!<description>
  ! Final preparation of the Newton solver.
!</description>

!<inputoutput>
  ! Structure to be initialised.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneData (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initData.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_doneStructure (rlinsolParam)
  
!<description>
  ! Cleanup of the data initalised in newtonit_initStructure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine newtonit_done (rlinsolParam)
  
!<description>
  ! Clean up the Newton iteration.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_newtonParameters), intent(inout) :: rlinsolParam
!</inputoutput>

!</subroutine>
   
  end subroutine

end module

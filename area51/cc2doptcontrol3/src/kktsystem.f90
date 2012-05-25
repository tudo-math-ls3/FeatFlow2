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
  
  use spacematvecassembly
  
  use kktsystemspaces
  
  implicit none
  
  private

!<types>

!<typeblock>

  ! This structure encapsules the discrete KKT system. It holds
  ! a solution of the primal equation, the dual equation and the
  ! corresponding control (if not directly calculated from the dual
  ! solution on-demand).
  type t_kktsystem
  
    ! Underlying physics
    type(t_settings_physics), pointer :: p_rphysics => null()
    
    ! Structure defining the optimal control problem to calculate
    type(t_settings_optcontrol), pointer :: p_roptControl => null()
  
    ! Solution of the primal equation of the KKT system.
    type(t_primalSpace), pointer :: p_rprimalSol => null()
    
    ! Solution of the dual equation of the KKT system.
    type(t_dualSpace), pointer :: p_rdualSol => null()
    
    ! Solution of the control equation of the KKT system.
    type(t_controlSpace), pointer :: p_rcontrol => null()

    ! Parameters for the assembly of space-time operators
    type(t_spacetimeOperatorAsm), pointer :: p_rspacetimeOperatorAsm => null()
    
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

  ! Solve the primal equation
  public :: kkt_solvePrimal

  ! Solve the dual equation
  public :: kkt_solveDual

  ! Calculate the control from the solution of the primal/dual equation
  public :: kkt_calcControl

  ! Calculate the residual of the control equation(s)
  public :: kkt_calcControlRes

  ! Solve the primal equation of the linearised KKT system
  public :: kkt_solvePrimalDirDeriv

  ! Solve the dual equation of the linearised KKT system
  public :: kkt_solveDualDirDeriv

  ! Calculate the control of the linearised KKT system 
  ! from the solution of the primal/dual equation
  public :: kkt_calcControlDirDeriv

  ! Calculate the residual of the control equation(s) in the linearised KKT system
  public :: kkt_calcControlResDirDeriv

contains

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solvePrimal (rkktsystem)
  
!<description>
  ! Solves the primal equation in the KKT system based on the control in the
  ! rkktsystem structure.
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>

    ! ... to be done
    call sys_halt()
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solveDual (rkktsystem)
  
!<description>
  ! Solves the dual equation in the KKT system based on the control and the
  ! primal solution in the rkktsystem structure.
!</description>
  
!<inputoutput>
  ! Structure defining the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</inputoutput>

!</subroutine>

    ! ... to be done
    call sys_halt()
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_calcControl (rkktsystem,rcontrol)
  
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
  type(t_controlSpace), intent(inout) :: rcontrol
!</inputoutput>

!</subroutine>

    ! This is strongly equation and problem dependent
    ! and may imply a projection to the admissible set.
    ! ... to be done
    call sys_halt()
   
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
    ! We call kkt_calcControl to calculate a new control 
    !
    !      rresidual = P ( -1/alpha [B'(u)]* lambda )
    !
    ! from the primal/dual variables in rkktsystem.
    call kkt_calcControl (rkktsystem,rresidual)
    
    ! Add -u:   rresidual = rresidual - u
    call kktsp_controlLinearComb (rkktsystem%p_rcontrol,-1.0_DP,rresidual,1.0_DP)
    
    ! Calculate the norm of the residual
    dres = sptivec_vectorNorm (rresidual%p_rvector,LINALG_NORML2)
   
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine kkt_solvePrimalDirDeriv (rkktsystemDirDeriv)
  
!<description>
  ! Solves the linearised primal equation in the KKT system.
!</description>
  
!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv
!</inputoutput>

!</subroutine>
   
    ! ... to be done
    call sys_halt()

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_solveDualDirDeriv (rkktsystemDirDeriv)
  
!<description>
  ! Solves the linearised dual equation in the KKT system.
!</description>
  
!<inputoutput>
  ! Structure defining a directional derivative of the KKT system.
  ! The solutions in this structure are taken as initial
  ! values. On exit, the structure contains improved solutions.
  type(t_kktsystemDirDeriv), intent(inout) :: rkktsystemDirDeriv
!</inputoutput>

!</subroutine>
   
    ! ... to be done
    call sys_halt()

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine kkt_calcControlDirDeriv (rkktsystemDirDeriv,rcontrol)
  
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
  type(t_controlSpace), intent(inout) :: rcontrol
!</inputoutput>

!</subroutine>
   
    ! This is strongly equation and problem dependent.
    ! ... to be done
    call sys_halt()
    
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
    call kkt_calcControlDirDeriv (rkktsystemDirDeriv,rresidual)

    ! b) rresidual = rresidual + d - g
    call kktsp_controlLinearComb (&
        rrhs,1.0_DP,&
        rkktsystemDirDeriv%p_rcontrolLin,-1.0_DP,&
        rresidual,1.0_DP,dres)

  end subroutine

end module

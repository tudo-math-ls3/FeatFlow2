!##############################################################################
!# ****************************************************************************
!# <name> afc_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global vairables
!# which are required to solve a conservation law for a scalar variable
!#
!# The following routines are available:
!#
!# 1.) afc_updateSolverMatrix
!#      -> update the matrices in the solver structure
!#
!# </purpose>
!##############################################################################

module afc_basic

  use storage
  use fparser
  use linearsystemscalar
  use paramlist
  use statistics
  use triangulation

  use boundaryfilter
  use problem
  use solver

  implicit none

  private
  public :: afc_updateSolverMatrix
  
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<constants>

!<constantblock description="Global format flags for convection-diffusion-reaction problem">
  
  ! primary boundary description
  integer, parameter, public :: CDEQ_BOUNDARY_PRIMAL  = 1

  ! dual boundary description
  integer, parameter, public :: CDEQ_BOUNDARY_DUAL    = 2


  ! stabilisation for convective part
  integer, parameter, public :: CDEQ_AFCSTAB_CONVECTION       = 1

  ! stabilisation for diffusive part
  integer, parameter, public :: CDEQ_AFCSTAB_DIFFUSION        = 2


  ! primary velocity vector
  integer, parameter, public :: CDEQ_VELOCITY                 = 1

  ! template matrix
  integer, parameter, public :: CDEQ_MATRIX_TEMPLATE          = 1

  ! system matrix
  integer, parameter, public :: CDEQ_MATRIX_A                 = 2

  ! system matrix (extended sparsity pattern)
  integer, parameter, public :: CDEQ_MATRIX_J                 = 3

  ! transport matrix
  integer, parameter, public :: CDEQ_MATRIX_L                 = 4

  ! consistent mass matrix
  integer, parameter, public :: CDEQ_MATRIX_MC                = 5

  ! lumped mass matrix
  integer, parameter, public :: CDEQ_MATRIX_ML                = 6

  ! coefficient matrix Dx(phi_i) Dx(phi_j) + Dy(phi_i) Dy(phi_j)
  integer, parameter, public :: CDEQ_MATRIX_S                 = 7

  ! coefficient matrix phi_i Dx(phi_j)
  integer, parameter, public :: CDEQ_MATRIX_CX                = 8

  ! coefficient matrix phi_i Dy(phi_j)
  integer, parameter, public :: CDEQ_MATRIX_CY                = 9

  ! coefficient matrix phi_i Dz(phi_j)
  integer, parameter, public :: CDEQ_MATRIX_CZ                = 10

!</constantblock>


!<constantblock description="Global type of flow velocities">

  ! zero velocity profile v=0
  integer, parameter, public :: VELOCITY_NONE                 = 0

  ! constant velocity profile v=v(x)
  integer, parameter, public :: VELOCITY_CONSTANT             = 1

  ! linear time-dependent velocity profile v=v(x,t)
  integer, parameter, public :: VELOCITY_TIMEDEP              = 2

  ! nonlinear Burgers' equation in space-time
  integer, parameter, public :: VELOCITY_BURGERS_SPACETIME    = 3

  ! nonlinear Buckley-Leverett equation in space-time
  integer, parameter, public :: VELOCITY_BUCKLEV_SPACETIME    = 4

  ! nonlinear Burgers' equation in 1D
  integer, parameter, public :: VELOCITY_BURGERS1D            = 5

  ! nonlinear Burgers' equation in 2D
  integer, parameter, public :: VELOCITY_BURGERS2D            = 6

  ! nonlinear Burgers' equation in 3D
  integer, parameter, public :: VELOCITY_BURGERS3D            = 7

  ! nonlinear Buckley-Leverett equation in 1D
  integer, parameter, public :: VELOCITY_BUCKLEV1D            = 8

!</constantblock>


!<constantblock description="Global type of flow">

  ! transient flow
  integer, parameter, public :: FLOW_TRANSIENT                = 0

  ! steady-state flow
  integer, parameter, public :: FLOW_STEADYSTATE              = 1

  ! pseudo transient flow
  integer, parameter, public :: FLOW_PSEUDOTRANSIENT          = 2

!</constantblock>


!<constantblock description="Global type of diffusion">

  ! zero diffusion
  integer, parameter, public :: DIFF_NONE                     = 0

  ! isotropic diffusion
  integer, parameter, public :: DIFF_ISOTROPIC                = 1

  ! anisotropic diffusion
  integer, parameter, public :: DIFF_ANISOTROPIC              = 2

!</constantblock>


!<constantblock description="Global type of r.h.s.">

  ! zero right-hand side
  integer, parameter, public :: RHS_ZERO                      = 0

  ! analytical right-hand side
  integer, parameter, public :: RHS_ANALYTIC                  = 1

!</constantblock>


!<constantblock description="Global types of perturbation parameters">

  ! Perturbation parameter is chosen as in the NITSOL package
  integer, parameter, public :: PERTURB_NITSOL                = 1

  ! Perturbation parameter is chosen as SQRT(machine precision)
  integer, parameter, public :: PERTURB_SQRTEPS               = 2

!</constantblock>

!<constantblock description="Flags for matrix update specification bitfield">

  ! Update the matrix only for the current solver
  integer(I32), parameter, public :: UPDMAT_NORECURSIVE            = 2**0

  ! Update the matrix for coarse grid solver
  integer(I32), parameter, public :: UPDMAT_LINEAR_SOLVER          = 2**1

  ! Update the preconditioner for the coarse grid solver
  integer(I32), parameter, public :: UPDMAT_LINEAR_PRECOND         = 2**2

  ! Update the solver of the linear multigrid solver
  integer(I32), parameter, public :: UPDMAT_LINEARMG_SOLVER        = 2**3

  ! Update the smoothers of the linear multigrid solver
  integer(I32), parameter, public :: UPDMAT_LINEARMG_SMOOTHER      = 2**4

  ! Update the coarsegrid solver of the linear multigrid solver
  integer(I32), parameter, public :: UPDMAT_LINEARMG_COARSEGRID    = 2**5

!</constantblock>

!<constantblock description="Predified bitfields for matrix update">

  ! Update everything
  integer(I32), parameter, public :: UPDMAT_ALL = UPDMAT_LINEAR_SOLVER +&
                                                  UPDMAT_LINEAR_PRECOND +&
                                                  UPDMAT_LINEARMG_SOLVER +&
                                                  UPDMAT_LINEARMG_SMOOTHER +&
                                                  UPDMAT_LINEARMG_COARSEGRID
  
  ! Update Jacobian matrix for steady-state flows
  integer(I32), parameter, public :: UPDMAT_JAC_STEADY = UPDMAT_LINEAR_SOLVER +&
                                                         UPDMAT_LINEARMG_SOLVER +&
                                                         UPDMAT_LINEARMG_SMOOTHER 

  ! Update Jacobian matrix for transient flows
  integer(I32), parameter, public :: UPDMAT_JAC_TRANSIENT = UPDMAT_LINEAR_SOLVER +&
                                                            UPDMAT_LINEAR_PRECOND +&
                                                            UPDMAT_LINEARMG_SOLVER +&
                                                            UPDMAT_LINEARMG_SMOOTHER 

!</constantblock>
  
!</constants>

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<globals>

  !*****************************************************************
  ! Parameter list
  type(t_parlist), save, public :: rparlist


  !*****************************************************************
  ! Timers
  type(t_timer), save, public :: rtimer_total
  type(t_timer), save, public :: rtimer_solution
  type(t_timer), save, public :: rtimer_adaptivity
  type(t_timer), save, public :: rtimer_errorestimation
  type(t_timer), save, public :: rtimer_triangulation
  type(t_timer), save, public :: rtimer_assembly_coeff
  type(t_timer), save, public :: rtimer_assembly_matrix
  type(t_timer), save, public :: rtimer_assembly_resrhs
  type(t_timer), save, public :: rtimer_prepostprocess

  
  !*****************************************************************
  ! Benchmark configuration

  ! type of flow
  integer, save, public :: iflowtype = FLOW_TRANSIENT


  !*****************************************************************
  ! R.H.S handling

  ! type of r.h.s
  integer, save, public :: irhstype = RHS_ZERO

  ! Function parser for analytical r.h.s
  type(t_fparser), save, public :: rrhsParser


  !*****************************************************************
  ! Diffusion handling

  ! Diffusion "matrix" in 1D
  real(DP), dimension(1,1), save, public :: DdiffusionMatrix1D

  ! Diffusion matrix in 2D
  real(DP), dimension(2,2), save, public :: DdiffusionMatrix2D

  ! Diffusion matrix in 3D
  real(DP), dimension(3,3), save, public :: DdiffusionMatrix3D

  ! type of diffusion
  integer, save, public :: idiffusiontype = DIFF_ISOTROPIC

  
  !*****************************************************************
  ! Velocity handling
  
  ! type of velocity vector
  integer, save, public :: ivelocitytype = VELOCITY_CONSTANT

  ! Switch to mark velocity vector for update
  logical, save, public :: bvelocityUpdate

  ! Function parser for velocity vector
  type(t_fparser), save, public :: rvelocityParser
   
!</globals>
 
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

contains
  
  !*****************************************************************************

!<subroutine>

  subroutine afc_updateSolverMatrix(rproblemLevel, rsolver, imatrix, iupdflag, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine updates the solver structure by setting the matrices.
    ! If the optional parameters NLMINOPT and NLMAXOPT are not given, then 
    ! only the current level of the given multigrid structure is processed.
!</description>

!<input>
    ! multigrid level to start with
    type(t_problemLevel), intent(IN), target :: rproblemLevel

    ! matrix number
    integer, intent(IN) :: imatrix

    ! flags which is used to specify the matrices to be updated
    integer(I32), intent(IN) :: iupdflag

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: nlmin,nlmax

    ! Set minimal level
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = solver_getNLMIN(rsolver, rproblemLevel%ilev)
    end if

    ! Set maximum level
    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = solver_getNLMAX(rsolver, rproblemLevel%ilev)
    end if
    
    ! Ok, let us update the solver (recursively?)
    call updateMatrix(rproblemLevel, rsolver, imatrix, iupdflag, nlmin, nlmax)

  contains

    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Set the solver matrices recursively. The routine must be called
    ! with a valid solver structure RSOLVER which serves as top-level
    ! solver. Depending on the parameter IUPDFLAG, the matrices
    ! are only set for the given solver structure or recursively for
    ! all subnodes. The parameter IMATRIX denotes the matrix number.
    ! The parameters NLMIN/NLMAX denote the minimum/maximum level to
    ! be considered.
    
    recursive subroutine updateMatrix(rproblemLevel, rsolver, imatrix, iupdflag, nlmin, nlmax)

      type(t_problemLevel), intent(IN), target :: rproblemLevel
      type(t_solver), intent(INOUT) :: rsolver
      integer(I32), intent(IN) :: iupdflag
      integer, intent(IN) :: imatrix,nlmin,nlmax

      ! local variables
      type(t_problemLevel), pointer :: rproblemLevelTmp,rproblemLevelCoarse
      integer :: i


      ! What kind of solver are we?
      select case(rsolver%csolverType)
      case (SV_FMG)
        
        ! There are no matrices for the full multigrid solver.
        if (iand(iupdflag, UPDMAT_NORECURSIVE) .eq. UPDMAT_NORECURSIVE) return
        
        ! Directly proceed to the coarsegrid solver which serves as
        ! nonlinear solver on each level
        call updatematrix(rproblemLevel, rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                          imatrix, iupdflag, nlmin, nlmax)
        

      case (SV_NONLINEARMG)
        
        ! There are no matrices for the nonlinear multigrid solver.
        if (iand(iupdflag, UPDMAT_NORECURSIVE) .eq. UPDMAT_NORECURSIVE) return
        
        ! Proceed to the coarsegrid solver
        call updateMatrix(rproblemLevel, rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                          imatrix, iupdflag, rsolver%p_solverMultigrid%nlmin,&
                          rsolver%p_solverMultigrid%nlmin)

        ! Proceed to the smoothers
        if (associated(rsolver%p_solverMultigrid%p_smoother)) then
          do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                 ubound(rsolver%p_solverMultigrid%p_smoother,1)
            
            call updateMatrix(rproblemLevel, rsolver%p_solverMultigrid%p_smoother(i),&
                              imatrix, iupdflag, rsolver%p_solverMultigrid%nlmin,&
                              rsolver%p_solverMultigrid%nlmin)
          end do
        end if
        
        
      case (SV_NONLINEAR)
        
        ! There are no matrices for the nonlinear single-grid solver.
        if (iand(iupdflag, UPDMAT_NORECURSIVE) .eq. UPDMAT_NORECURSIVE) return

        ! Proceed to the linear solver subnode.
        call updateMatrix(rproblemLevel, rsolver%p_solverSubnode,&
                          imatrix, iupdflag, nlmin, nlmax)
        
        
      case (SV_LINEARMG)
        
        ! Are there multiple levels?
        if (rsolver%p_solverMultigrid%nlmin .eq.&
            rsolver%p_solverMultigrid%nlmax) then

          ! Proceed to single grid solver
          call updateMatrix(rproblemLevel, rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                            imatrix, iupdflag, rsolver%p_solverMultigrid%nlmin,&
                            rsolver%p_solverMultigrid%nlmin)
          
        else

          ! We are on the level of the linear multigrid solver.         
          rproblemLevelTmp    => rproblemLevel
          rproblemLevelCoarse => rproblemLevel
          do while(associated(rproblemLevelTmp))
            
            ! Do we have to set matrices for this level?
            if (rproblemLevelTmp%ilev > nlmax) then
              rproblemLevelTmp => rproblemLevelTmp%p_rproblemLevelCoarse
              cycle
            elseif(rproblemLevelTmp%ilev < nlmin) then
              exit
            end if
            
            if (iand(iupdflag, UPDMAT_LINEARMG_SOLVER) .eq.&
                               UPDMAT_LINEARMG_SOLVER) then
              ! Set the system matrix for the linear solver
              call solver_setSolverMatrix(rsolver, rproblemLevelTmp%Rmatrix(imatrix), rproblemLevelTmp%ilev)
            end if
          
            if (iand(iupdflag, UPDMAT_LINEARMG_SMOOTHER) .eq.&
                               UPDMAT_LINEARMG_SMOOTHER) then
              ! Set the system matrix for the linear smoother
              ! Note that the smoother is not required in the coarsest level
              if (rproblemLevelTmp%ilev .gt. nlmin) then
                call solver_setSmootherMatrix(rsolver, rproblemLevelTmp%Rmatrix(imatrix),&
                                              rproblemLevelTmp%ilev)
              end if
            end if
            
            ! Switch to next coarser level
            rproblemLevelCoarse => rproblemLevelTmp
            rproblemLevelTmp    => rproblemLevelTmp%p_rproblemLevelCoarse
          end do
          
          if (iand(iupdflag, UPDMAT_LINEARMG_COARSEGRID) .eq.&
                             UPDMAT_LINEARMG_COARSEGRID) then
            ! Set the system matrix for the linear coarse grid solver
            call updateMatrix(rproblemLevelCoarse, rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                              imatrix, iupdflag, rsolver%p_solverMultigrid%nlmin,&
                              rsolver%p_solverMultigrid%nlmin)
          end if
        end if


      case (SV_LINEAR)

        ! The solver matrix and preconditioner matrix are only updated if
        ! the current level satisfies nlmin <= ilev <= nlmax
        if (nlmin .eq. rproblemLevel%ilev) then
          
          if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                             UPDMAT_LINEAR_SOLVER) then
            ! Set the system matrix of the single-grid solver
            call solver_setSolverMatrix(rsolver, rproblemLevel%Rmatrix(imatrix))
          end if
          
          if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                             UPDMAT_LINEAR_PRECOND) then
            ! Set the system matrix of the preconditioner
            call solver_setPrecondMatrix(rsolver, rproblemLevel%Rmatrix(imatrix))
          end if

        else
          rproblemLevelTmp => rproblemLevel
          do while(associated(rproblemLevelTmp))
            
            ! Are we on the coarse grid level?
            if (rproblemLevelTmp%ilev > nlmin) then
              rproblemLevelTmp => rproblemLevelTmp%p_rproblemLevelCoarse
              cycle
            elseif(rproblemLevelTmp%ilev .eq. nlmin) then
              exit
            else
              call output_line('Invalid multigrid level!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
              call sys_halt()
            end if
          end do
          
          if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                             UPDMAT_LINEAR_SOLVER) then
            ! Set the system matrix of the single-grid solver
            call solver_setSolverMatrix(rsolver, rproblemLevelTmp%Rmatrix(imatrix))
          end if
          
          if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                             UPDMAT_LINEAR_PRECOND) then
            ! Set the system matrix of the preconditioner
            call solver_setPrecondMatrix(rsolver, rproblemLevelTmp%Rmatrix(imatrix))
          end if
        end if
        
      case DEFAULT
        call output_line('Unsupported solver type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
        call sys_halt()
      end select
    end subroutine updateMatrix
  end subroutine afc_updateSolverMatrix
  
end module afc_basic

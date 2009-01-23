!##############################################################################
!# ****************************************************************************
!# <name> euler_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines and provides all global variables
!# which are required to solve the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!#  1.) euler_getNVAR
!#      Returns the number of variables depending on the spatial dimension
!#
!#  2.) euler_getVariableNodewise
!#      -> Returns individual variable from solution vector in interleave format
!#
!#  3.) euler_getVariableBlockwise
!#      -> Returns individual variable from solution vector in block format
!#
!#  4.) euler_updateSolverMatrix
!#      -> Updates the matrices in the solver structure
!#
!# </purpose>
!##############################################################################

module euler_basic

  use afcstabilisation
  use fsystem
  use linearsystemscalar
  use paramlist
  use statistics
  use triangulation

  use thermodynamics
  use problem
  use solver


  implicit none

  private
  public :: euler_getNVAR
  public :: euler_getVariableNodewise
  public :: euler_getVariableBlockwise
  public :: euler_updateSolverMatrix

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!<constants>

!<constantblock description="Global format flags for boundary treatment">
  
  ! primary boundary description
  integer, parameter, public :: CNSE_BOUNDARY_PRIMAL  = 1

  ! dual boundary description
  integer, parameter, public :: CNSE_BOUNDARY_DUAL    = 2

!</constantblock>


!<constantblock description="Global format flags for spatial stabilization">

  ! stabilisation for inviscid part
  integer, parameter, public :: CNSE_AFCSTAB_INVISCID = 1

  ! stabilisation for viscous part
  integer, parameter, public :: CNSE_AFCSTAB_VISCOUS  = 2

!</constantblock>


!<constantblock description="Global format flags for matrix addressing">

  ! system matrix
  integer, parameter, public :: CNSE_MATRIX_A         = 1

  ! system matrix (extended sparsity pattern)
  integer, parameter, public :: CNSE_MATRIX_J         = 1   ! NOTE: just for compatibility reasons

  ! template matrix
  integer, parameter, public :: CNSE_MATRIX_TEMPLATE  = 2

  ! consistent mass matrix
  integer, parameter, public :: CNSE_MATRIX_MC        = 3

  ! lumped mass matrix
  integer, parameter, public :: CNSE_MATRIX_ML        = 4

  ! coefficient matrix phi_i Dx(phi_j)
  integer, parameter, public :: CNSE_MATRIX_CX        = 5
  
  ! coefficient matrix phi_i Dy(phi_j)
  integer, parameter, public :: CNSE_MATRIX_CY        = 6

  ! coefficient matrix phi_i Dz(phi_j)
  integer, parameter, public :: CNSE_MATRIX_CZ        = 7
  
!</constantblock>


!<constantblock description="Global flags for time measurement">

  ! time for computing the solution
  integer, parameter, public :: CNSE_TIME_SOLUTION             = 1

  ! time for grid adaptivity
  integer, parameter, public :: CNSE_TIME_ADAPTIVITY           = 2

  ! time for error estimation
  integer, parameter, public :: CNSE_TIME_ERRORESTIMATION      = 3

  ! time for generation of triangulation
  integer, parameter, public :: CNSE_TIME_TRIANGULATION        = 4

  ! time for coefficient matrix assembly
  integer, parameter, public :: CNSE_TIME_ASSEMBLY_COEFFICIENT = 5

  ! time for global matrix assembly
  integer, parameter, public :: CNSE_TIME_ASSEMBLY_MATRIX      = 6

  ! time for global defect + rhs assembly
  integer, parameter, public :: CNSE_TIME_ASSEMBLY_RESIDUALRHS = 7

  ! time for pre- and post-processing
  integer, parameter, public :: CNSE_TIME_PREPOSTPROCESS       = 8

!</constantblock>


!<constantblock description="Global type of flow behavior">

  ! time-dependent flow
  integer, parameter, public :: FLOW_TRANSIENT                 = 0

  ! steady-state flow
  integer, parameter, public :: FLOW_STEADYSTATE               = 1

  ! pseudo transient flow
  integer, parameter, public :: FLOW_PSEUDOTRANSIENT           = 2

!</constantblock>


!<constantblock description="Global type of coupling approach">

  ! time-dependent flow
  integer, parameter, public :: FLOW_SEGREGATED                = 0

  ! steady-state flow
  integer, parameter, public :: FLOW_ALLCOUPLED                = 1

!</constantblock>


!<constantblock description="Global type of system structure">

  ! time-dependent flow
  integer, parameter, public :: SYSTEM_INTERLEAVEFORMAT        = 0

  ! steady-state flow
  integer, parameter, public :: SYSTEM_BLOCKFORMAT             = 1

!</constantblock>


!<constantblock description="Global types of perturbation parameters">

  ! Perturbation parameter is chosen as in the NITSOL package
  integer, parameter, public :: PERTURB_NITSOL                 = 1

  ! Perturbation parameter is chosen as SQRT(machine precision)
  integer, parameter, public :: PERTURB_SQRTEPS                = 2

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


!<constantblock description="Global constants from Gas Dynamic">

  ! ratio of specific heats
  real(DP), parameter, public :: GAMMA = GAMMA_AIR
  real(DP), parameter, public :: RGAS  = R_AIR

  ! auxiliary parameters related to gamma
  real(DP), parameter, public :: G1  = GAMMA-1.0
  real(DP), parameter, public :: G2  = (GAMMA-1.0)/2.0
  real(DP), parameter, public :: G3  = 1.0/(GAMMA-1.0)
  real(DP), parameter, public :: G4  = 1.0/GAMMA
  real(DP), parameter, public :: G5  = GAMMA/(GAMMA-1.0)
  real(DP), parameter, public :: G6  = GAMMA-2.0
  real(DP), parameter, public :: G7  = (GAMMA-1.0)/(2.0*GAMMA)
  real(DP), parameter, public :: G8  = (GAMMA+1.0)/(2.0*GAMMA)
  real(DP), parameter, public :: G9  = 2.0/(GAMMA-1.0)
  real(DP), parameter, public :: G10 = 2.0/(GAMMA+1.0)
  real(DP), parameter, public :: G11 = 2.0*GAMMA/(GAMMA-1.0)
  real(DP), parameter, public :: G12 = (GAMMA-1.0)/(GAMMA+1.0)

!</constantblock>


!<constantblock description="Global constants for number of variables">

  ! number of solution components in 1D
  integer, parameter, public :: NVAR1D = 3

  ! number of solution components in 2D
  integer, parameter, public :: NVAR2D = 4

  ! number of solution components in 3D
  integer, parameter, public :: NVAR3D = 5

!</constantblock>


!<constantblock description="Global constants for type of variables">

  ! position of conservative variables in solution vector
  integer, parameter, public :: VAR_DENSITY    = 1
  integer, parameter, public :: VAR_MOMENTUM_X = 2
  integer, parameter, public :: VAR_MOMENTUM_Y = 3
  integer, parameter, public :: VAR_MOMENTUM_Z = 4
!  integer, parameter, public :: VAR_ENERGY     = NVAR2D

  ! auxiliary/derived variables
!  integer, parameter, public :: VAR_PRESSURE   = NVAR2D+1
!  integer, parameter, public :: VAR_MACHNUMBER = NVAR2D+2

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
  ! Possible values are: FLOW_TRANSIENT, FLOW_STEADYSTATE
  integer, save, public :: iflowtype     = FLOW_TRANSIENT

  ! type of coupling
  ! Possible values are: FLOW_SEGREGATED, FLOW_ALLCOUPLED
  integer, save, public :: icoupled      = FLOW_SEGREGATED

  ! type of preconditioner
  ! Possible values are: AFCSTAB_SCALARDISSIPATION,
  ! AFCSTAB_TENSORDISSIPATION, AFCSTAB_GALERKIN
  integer, save, public :: iprecond      = AFCSTAB_SCALARDISSIPATION

  ! type of system block structure
  ! Possible values are: SYSTEM_INTERLEAVEFORMAT, SYSTEM_BLOCKFORMAT
  integer, save, public :: isystemFormat = SYSTEM_INTERLEAVEFORMAT

!</globals>

contains

  !*****************************************************************************

!<function>

  pure function euler_getNVAR(rproblemLevel) result(NVAR)

!<description>
    ! This function returns the number of variables present in the given level
!</description>

!<input>
    ! multigrid level
    type(t_problemLevel), intent(IN) :: rproblemLevel
!</input>

!<result>
    ! number of variables
    integer :: NVAR
!</result>
!</function>

    NVAR = rproblemLevel%rtriangulation%ndim + 2
    
  end function euler_getNVAR

  !*****************************************************************************

!<subroutine>
  
  subroutine euler_getVariableNodewise(NEQ, NVAR, Dsolution, ivariable, Dvalue)

!<description>
    ! This subroutine extracts an individual variable from the
    ! global solution vector store in interleave format.
!</description>

!<input>
    ! Number of equations
    integer, intent(IN) :: NEQ

    ! Number of variables
    integer, intent(IN) :: NVAR
    
    ! Solution vector
    real(DP), dimension(NVAR,NEQ), intent(IN) :: Dsolution

    ! Variable to extract
    integer, intent(IN) :: ivariable
!</input>

!<output>
    ! Value of the extracted variable
    real(DP), dimension(:), intent(OUT) :: Dvalue
!</output>
!</subroutine>

    ! local variables
    real(DP) :: p
    integer :: ieq
    
    
    ! Which variable should be extracted?
    if (ivariable .eq. VAR_DENSITY) then
      ! density
      do ieq = 1, NEQ
        Dvalue(ieq) = Dsolution(VAR_DENSITY, ieq)
      end do
      
    elseif ((ivariable .ge. VAR_MOMENTUM_X) .and.&
            (ivariable .le. NVAR)) then
      ! velocities and energy (located at NVAR)
      do ieq = 1, NEQ
        Dvalue(ieq) = Dsolution(ivariable, ieq)/&
                      Dsolution(VAR_DENSITY, ieq)
      end do
      
    elseif (ivariable .eq. NVAR+1) then
      ! pressure ...
      select case(NVAR)
      case (NVAR1D)
        ! ... in 1D
        do ieq = 1, NEQ
          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
              Dsolution(NVAR1D, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq))
        end do

      case (NVAR2D)
        ! ... in 2D
        do ieq = 1, NEQ
          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
              Dsolution(NVAR2D, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq))
        end do
      
      case (NVAR3D)
        ! ... in 3D
        do ieq = 1, NEQ
          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
              Dsolution(NVAR2D, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Z, ieq)/Dsolution(VAR_DENSITY, ieq))
        end do
      end select

    elseif(ivariable .eq. NVAR+2) then
      ! Mach number ...
      select case(NVAR)
      case (NVAR1D)
        ! ... in 1D
        do ieq = 1, NEQ
          p = thdyn_pressure(GAMMA_AIR,&
              Dsolution(NVAR1D, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq))
          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq))
        end do

      case (NVAR2D)
        ! ... in 2D
        do ieq = 1, NEQ
          p = thdyn_pressure(GAMMA_AIR,&
              Dsolution(NVAR2D, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq))
          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq))
        end do
        
      case (NVAR3D)
        ! ... in 3D
        do ieq = 1, NEQ
          p = thdyn_pressure(GAMMA_AIR,&
              Dsolution(NVAR3D, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Z, ieq)/Dsolution(VAR_DENSITY, ieq))
          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
              Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_X, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Y, ieq)/Dsolution(VAR_DENSITY, ieq),&
              Dsolution(VAR_MOMENTUM_Z, ieq)/Dsolution(VAR_DENSITY, ieq))
        end do
      end select
      
    else
      call output_line('Unsupported variable!',&
                        OU_CLASS_ERROR,OU_MODE_STD,'euler_getVariableNodewise')
      call sys_halt()
    end if
  end subroutine euler_getVariableNodewise
  
  !*****************************************************************************
  
!<subroutine>
  
  subroutine euler_getVariableBlockwise(NEQ, NVAR, Dsolution, ivariable, Dvalue)

!<description>
    ! This subroutine extracts an individual variable from the
    ! global solution vector store in block format.
!</description>

!<input>
    ! Number of equations
    integer, intent(IN) :: NEQ

    ! Number of variables
    integer, intent(IN) :: NVAR
    
    ! Solution vector
    real(DP), dimension(NEQ,NVAR), intent(IN) :: Dsolution

    ! Variable to extract
    integer, intent(IN) :: ivariable
!</input>

!<output>
    ! Value of the extracted variable
    real(DP), dimension(NEQ), intent(OUT) :: Dvalue
!</output>
!</subroutine>     
    
    ! local variables
    real(DP) :: p
    integer :: ieq
    
    ! Which variable should be extracted?
    if (ivariable .eq. VAR_DENSITY) then
      ! density
      do ieq = 1, NEQ
        Dvalue(ieq) = Dsolution(ieq, VAR_DENSITY)
      end do
      
    elseif ((ivariable .ge. VAR_MOMENTUM_X) .and.&
            (ivariable .le. NVAR)) then
      ! velocities and energy (located at NVAR)
      do ieq = 1, NEQ
        Dvalue(ieq) = Dsolution(ieq,ivariable)/&
                      Dsolution(ieq, VAR_DENSITY)
      end do
      
    elseif(ivariable .eq. NVAR+1) then
      ! pressure ...
      select case(NVAR)
      case (NVAR1D)
        ! ... in 1D
        do ieq = 1, NEQ
          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
              Dsolution(ieq, NVAR1D)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY))
        end do

      case (NVAR2D)
        ! ... in 2D
        do ieq = 1, NEQ
          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
              Dsolution(ieq, NVAR2D)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY))
        end do
        
      case (NVAR3D)
        ! ... in 3D
        do ieq = 1, NEQ
          Dvalue(ieq) = thdyn_pressure(GAMMA_AIR,&
              Dsolution(ieq, NVAR3D)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Z)/Dsolution(ieq, VAR_DENSITY))
        end do
      end select
      
    elseif (ivariable .eq. NVAR+2) then
      ! Mach number ...
      select case(NVAR)
      case (NVAR1D)
        ! ... in 1D
        do ieq = 1, NEQ
          p = thdyn_pressure(GAMMA_AIR,&
              Dsolution(ieq, NVAR1D)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY))
          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY))
        end do

      case (NVAR2D)
        ! ... in 2D
        do ieq = 1, NEQ
          p = thdyn_pressure(GAMMA_AIR,&
              Dsolution(ieq, NVAR2D)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY))
          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY))
        end do

      case (NVAR3D)
        ! ... in 3D
        do ieq = 1, NEQ
          p = thdyn_pressure(GAMMA_AIR,&
              Dsolution(ieq, NVAR3D)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Z)/Dsolution(ieq, VAR_DENSITY))
          Dvalue(ieq) = thdyn_Machnumber(GAMMA_AIR, p,&
              Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_X)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Y)/Dsolution(ieq, VAR_DENSITY),&
              Dsolution(ieq, VAR_MOMENTUM_Z)/Dsolution(ieq, VAR_DENSITY))
        end do
      end select
      
    else
      call output_line('Unsupported variable!',&
                        OU_CLASS_ERROR,OU_MODE_STD,'euler_getVariableBlockwise')
      call sys_halt()
    end if
  end subroutine euler_getVariableBlockwise

  !*****************************************************************************

!<subroutine>

  subroutine euler_updateSolverMatrix(rproblemLevel, rsolver, imatrix,&
                                   iupdflag, nlminOpt, nlmaxOpt)

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
    
    recursive subroutine updateMatrix(rproblemLevel, rsolver, imatrix,&
                                      iupdflag, nlmin, nlmax)
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

        ! Directly proceed to the coarsegrid solver
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

        ! Directly proceed to the linear solver subnode.
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
            
            ! What type of matrix format are we
            select case(isystemFormat)
            case (SYSTEM_INTERLEAVEFORMAT)
              
              if (iand(iupdflag, UPDMAT_LINEARMG_SOLVER) .eq.&
                                 UPDMAT_LINEARMG_SOLVER) then
                ! Set the system matrix for the linear solver
                call solver_setSolverMatrix(rsolver, rproblemLevelTmp%Rmatrix(imatrix),&
                                            rproblemLevelTmp%ilev)
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

            case (SYSTEM_BLOCKFORMAT)

              if (iand(iupdflag, UPDMAT_LINEARMG_SOLVER) .eq.&
                                 UPDMAT_LINEARMG_SOLVER) then
                ! Set the system matrix for the linear solver
                call solver_setSolverMatrix(rsolver, rproblemLevelTmp%RmatrixBlock(imatrix),&
                                            rproblemLevelTmp%ilev)
              end if
              
              if (iand(iupdflag, UPDMAT_LINEARMG_SMOOTHER) .eq.&
                                 UPDMAT_LINEARMG_SMOOTHER) then
                ! Set the system matrix for the linear smoother
                ! Note that the smoother is not required in the coarsest level
                if (rproblemLevelTmp%ilev .gt. nlmin) then
                  call solver_setSmootherMatrix(rsolver, rproblemLevelTmp%RmatrixBlock(imatrix),&
                                                rproblemLevelTmp%ilev)
                end if
              end if
              
            case DEFAULT
              call output_line('Unsupported system format!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
              call sys_halt()
            end select
            
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
          
          ! What type of matrix format are we
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)

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

          case (SYSTEM_BLOCKFORMAT)

            if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                               UPDMAT_LINEAR_SOLVER) then
              ! Set the system matrix of the single-grid solver
              call solver_setSolverMatrix(rsolver, rproblemLevel%RmatrixBlock(imatrix))
            end if
          
            if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                               UPDMAT_LINEAR_PRECOND) then
              ! Set the system matrix of the preconditioner
              call solver_setPrecondMatrix(rsolver, rproblemLevel%RmatrixBlock(imatrix))
            end if

          case DEFAULT
            call output_line('Unsupported system format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
            call sys_halt()
          end select

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
          
          ! What type of matrix format are we
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            
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
            
          case (SYSTEM_BLOCKFORMAT)

            if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                               UPDMAT_LINEAR_SOLVER) then
              ! Set the system matrix of the single-grid solver
              call solver_setSolverMatrix(rsolver, rproblemLevelTmp%RmatrixBlock(imatrix))
            end if
          
            if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                               UPDMAT_LINEAR_PRECOND) then
              ! Set the system matrix of the preconditioner
              call solver_setPrecondMatrix(rsolver, rproblemLevelTmp%RmatrixBlock(imatrix))
            end if
            
          case DEFAULT
            call output_line('Unsupported system format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
            call sys_halt()
          end select
        end if
        
      case DEFAULT
        call output_line('Unsupported solver type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
        call sys_halt()
      end select
    end subroutine updateMatrix
  end subroutine euler_updateSolverMatrix

end module euler_basic

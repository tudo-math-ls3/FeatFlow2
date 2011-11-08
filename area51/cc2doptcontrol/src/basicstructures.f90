!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic problem definitions for the cc2dmini_method2
!# solver. The basic structure and content of the different structures
!# are described here.
!# </purpose>
!##############################################################################

module basicstructures

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use timestepping
  use timediscretisation
  use spacetimevectors
  use discretebc
  use discretefbc
  use fparser
  
  use collection
  
  implicit none
  
!<types>

!<typeblock>

  ! A type block specifying all 'static' information which are depending
  ! on a discretisation and a triangulation. Such static information can be
  ! precalculated and is valid until the mesh or the FE spaces change.
  type t_staticLevelInfo
  
    ! A template FEM matrix that defines the structure of Laplace/Stokes/...
    ! matrices. The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateFEM

    ! A template FEM matrix that defines the structure of the pressure
    ! matrices. The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateFEMPressure

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2). The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateGradient

    ! A template FEM matrix that defines the structure of divergence
    ! matrices (D1/D2). The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateDivergence

    ! Precalculated Stokes matrix for that specific level (=nu*Laplace)
    type(t_matrixScalar) :: rmatrixStokes
    
    ! Precalculated B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! Precalculated B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2
    
    ! Precalculated D1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixD1

    ! Precalculated D2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixD2

    ! Precalculated D1^T-matrix for that specific level.
    ! This matrix either coincides with rmatrixB1 (in case D1=B1^T)
    ! or describes an independent D1 matrix.
    type(t_matrixScalar) :: rmatrixD1T

    ! Precalculated D2-matrix for that specific level.
    ! This matrix either coincides with rmatrixB1 (in case D2=B2^T)
    ! or describes an independent D2 matrix.
    type(t_matrixScalar) :: rmatrixD2T

    ! A scalar discretisation structure that specifies how to generate
    ! the mass matrix in the velocity FEM space.
    type(t_spatialDiscretisation) :: rdiscretisationMass

    ! A scalar discretisation structure that specifies how to generate
    ! the mass matrix in the pressure FEM space.
    type(t_spatialDiscretisation) :: rdiscretisationMassPressure

    ! Precalculated mass matrix for the velocity space.
    type(t_matrixScalar) :: rmatrixMass

    ! Precalculated mass matrix for the pressure space.
    type(t_matrixScalar) :: rmatrixMassPressure

    ! An object specifying the discretisation in the velocity space
    ! to be used for (edge) stabilisation. Only used if edge stabilisation
    ! is activated, otherwise this coincides with the default discretisation
    ! of the velocity space.
    type(t_spatialDiscretisation) :: rdiscretisationStabil
    
    ! Matrix with a precomputed EOJ stabilisation operator -- if EOJ is active.
    type(t_matrixScalar) :: rmatrixEOJ1
    type(t_matrixScalar) :: rmatrixEOJ2
    
  end type

!</typeblock>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation

    ! An object specifying the block discretisation to be used for (edge)
    ! stabilisation. Only used if edge stabilisation is activated, otherwise
    ! this coincides with rdiscretisation.
    type(t_blockDiscretisation) :: rdiscretisationStabil

    ! An object specifying the block discretisation structure only for the
    ! primal space. This coincides with the discretisation structure of the dual
    ! space.
    type(t_blockDiscretisation) :: rdiscretisationPrimal

!    ! Three temp vectors for the full system.
!    type(t_vectorBlock) :: rtempVector1
!    type(t_vectorBlock) :: rtempVector2
!    type(t_vectorBlock) :: rtempVector3
!
!    ! Reference ro rtempVector(1:3), which corresponds to the primal solution.
!    type(t_vectorBlock) :: rtempVectorPrimal
!
!    ! Reference ro rtempVector(4:6), which corresponds to the dual solution.
!    type(t_vectorBlock) :: rtempVectorDual
!
    ! This flag signales whether there are Neumann boundary components
    ! visible on the boundary of this level or not. If there are no
    ! Neumann boundary components visible, the equation gets indefinite
    ! for the pressure.
    !logical :: bhasNeumannBoundary
    
    ! A structure containing all static information about this level.
    type(t_staticLevelInfo) :: rstaticInfo

  end type
  
!</typeblock>


!<typeblock>
  
  ! Application-specific type block for the nonstationary Nav.St. problem
  type t_problem_explTimeStepping
  
    ! Number of time steps; former NITNS
    integer :: niterations         = 0
    
    ! Absolute start time of the simulation
    real(DP) :: dtimeInit          = 0.0_DP
    
    ! Maximum time of the simulation
    real(DP) :: dtimeMax           = 0.0_DP
    
    ! Time-stepping scheme (IFRSTP);
    ! 0=one step scheme, 1=fractional step, 2=dG(0)
    integer :: ctimeStepScheme     = 0
    
    ! Parameter for one step scheme (THETA) if ctimeStepScheme=0;
    ! =0:Forward Euler(instable), =1: Backward Euler, =0.5: Crank-Nicolson
    real(DP) :: dtimeStepTheta     = 1.0_DP
    
    ! Parameter for the one-step theta scheme in the pressure.
    ! usually = 1.0, independent of dtimeStepTheta.
    real(dp) :: dtimeStepThetaPressure = 1.0_DP
    
  end type

!</typeblock>


!<typeblock>

  ! This block saves parameters of the optimal control problem
  type t_problem_optcontrol
  
    ! $\alpha$ parameter of the optimal control functional
    real(DP) :: dalphaC = 1.0_DP
    
    ! $\gamma$ parameter of the nonstationary optimal control functional
    real(DP) :: dgammaC = 0.0_DP
  
    ! Type of target flow.
    ! =0: analytically given.
    ! =1: stationary solution, read from file.
    ! =2: nonstationary solution, read from a sequence of files.
    ! =3: stationary solution, read from file. Arbitrary level.
    ! =4: nonstationary solution, read from a sequence of files. Arbitrary level.
    ! =5: analytically given in stargetFlowExpressionX/stargetFlowExpressionY
    integer :: itypeTargetFlow = 0
    
    ! Analytical expression that defines the target flow for the X/Y coordinate
    ! in case itypeTargetFlow = 5.
    character(SYS_STRLEN) :: stargetFlowExpressionX = ''
    character(SYS_STRLEN) :: stargetFlowExpressionY = ''
    
    ! Parser object for the target flow expressions.
    type(t_fparser) :: rparserTargetFlowExpression
    
    ! Formulation of the Space-time problem.
    ! =0: usual formulation as specified in the DFG applicance
    ! =1: Formulation for the generation of reference results from papers
    ! The two formulations differ in a "-"-sign in front of the dual velocity.
    integer :: ispaceTimeFormulation = 0
  
    ! Whether to treat the convection explicitly or implicitly.
    ! =0: Treat the convection implicitely.
    ! =1: Treat the convection explicitely. (This is the formulation of the paper
    !     of Baerwolff and Hinze!)
    integer :: iconvectionExplicit = 0

    ! Name of the file with the target flow
    character(SYS_STRLEN) :: stargetFlow = ''
    
    ! Refinement level of the target flow.
    integer :: ilevelTargetFlow = 0
    
    ! Element type of the target flow
    integer :: ielementTypeTargetFlow = -1
    
    ! Name of the TRI file with the mesh corresponding to target flow.
    ! ='': Use the same mesh as for the computation of the solution
    character(SYS_STRLEN) :: smeshTargetFlow = ''

    ! Underlying mesh for the target flow if itypeTargetFlow<>0.
    type(t_triangulation), pointer :: p_rtriangulationTargetFlow => null()
    
    ! Discretisation structure specifying the discretisation of the
    ! target flow.
    type(t_blockDiscretisation), pointer :: p_rdiscrTargetFlow => null()
    
    ! Solution vector containing a stationary target flow.
    ! Only used if itypeTargetFlow <> 0.
    type(t_vectorBlock) :: rtargetFlow
  
    ! Time discretisation for a nonstationary target flow.
    ! Only used if itypeTargetFlow = 2 or 4
    type(t_timeDiscretisation) :: rtargetTimeDiscr
    
    ! Solution vector containing a nonstationary target flow.
    ! Only used if itypeTargetFlow = 2 or 4
    type(t_spaceTimeVector) :: rtargetFlowNonstat
    
    ! Target flow granilarity for the filenames of a nonstationary target flow.
    ! Specifies how the number in the filename is increased.
    integer :: itargetFlowDelta = 1
    
    ! Number of files that describe the target flow.
    ! -1=As many files as there are solution vectors on the finest time level.
    integer :: itargetFlowTimesteps = -1
  
    ! Type of constraints to apply to the control u.
    ! =0: No constraints.
    ! =1: Constant constraints on u active:
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2
    !     Implemented by cubature point.
    ! =2: Constant constraints on u active:
    !     dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2
    !     Implementation by DOF.
    integer :: ccontrolConstraints = 0

    ! Constraints on u_1
    real(DP) :: dumin1 = -1.0E10
    real(DP) :: dumax1 = 1.0E10

    ! Constraints in u_2
    real(DP) :: dumin2 = -1.0E10
    real(DP) :: dumax2 = 1.0E10
    
  end type

!</typeblock>


!<typeblock>
  
  ! Structure with some pointers to internal variables of a special algorithm.
  ! Used to provide data to the parser so the user can advice the program to do
  ! something with this data.
  ! Only valid during the execution of the command parser.
  type t_problem_oneshot
  
    ! Current global iteration
    integer :: iglobIter
  
    ! Norm of the defect vector
    real(DP) :: ddefNorm
    
    ! Pointer to the current space-time solution vector
    type(t_spaceTimeVector), pointer :: p_rx

    ! Pointer to the current space-time RHS vector
    type(t_spaceTimeVector), pointer :: p_rb
    
  end type

!</typeblock>

!<typeblock>

  ! This type block encapsules all physical constants and configuration
  ! parameters for the primal equation. This includes e.g. the type of the equation,
  ! viscosity parameter etc.
  type t_physicsPrimal
  
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    ! Type of problem.
    ! =0: Navier-Stokes.
    ! =1: Stokes.
    integer :: iequation
    
    ! Type of subproblem of the main problem. Depending on iequationType.
    ! If iequationType=0 or =1:
    ! =0: (Navier-)Stokes with gradient tensor
    ! =1: (Navier-)Stokes with deformation tensor
    integer :: isubEquation
    
  
  end type

!</typeblock>



!<typeblock>

  ! Application-specific type block for the Nav.St. problem
  type t_problem
  
    ! Output level during the initialisation phase.
    integer :: MSHOW_Initialisation
  
    ! Output level of the application.
    integer :: MT_OutputLevel
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: NLMIN
    
    ! Maximum refinement level
    integer :: NLMAX
    
    ! Physics of the problem
    type(t_physicsPrimal) :: rphysicsPrimal
    
    ! Type of right hand side in the primal equation.
    ! =-1: given as analytical expression srhsExpressionX,srhsExpressionY.
    ! >=0: analytically given by callback function
    integer :: irhs
    
    ! Analytical expression that defines the RHS for the X/Y coordinate
    ! in case irhs = -1.
    character(len=SYS_STRLEN) :: srhsExpressionX
    character(len=SYS_STRLEN) :: srhsExpressionY
    
    ! Parser object for the RHS expressions.
    type(t_fparser) :: rrhsParser
        
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer           :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    type(t_problem_lvl), dimension(:), pointer :: RlevelInfo
    
    ! Type of simulation.
    ! =0: stationary simulation.
    ! =1: time-dependent simulation with explicit time stepping configured
    !     by rtimedependence
    integer                               :: itimedependence
    
    ! A parameter block for everything that controls the time dependence.
    ! Only valid if itimedependence=1!
    type(t_problem_explTimeStepping)      :: rtimedependence
    
    ! A parameter block for everything that controls the optimal control
    ! problem.
    type(t_problem_optcontrol)            :: roptcontrol
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection)                    :: rcollection
    
    ! A param list that saves all parameters from the DAT/INI file(s).
    type(t_parlist)                       :: rparamList

    ! A t_problem_oneshot structure that is only valid during the
    ! execution of the command parser.
    type(t_problem_oneshot)               :: rdataOneshot

  end type

!</typeblock>

!</types>

!<globals>

  ! Directory containing the data files.
  character(LEN=SYS_STRLEN), save :: DIR_DATA = "./data";

!</globals>

!******************************************************************************
! Documentation
!******************************************************************************

!******************************************************************************
! The problem structure t_problem
!******************************************************************************
!
! The problem structure t_problem collects all information of importance
! for the main CC2D solver module. It contains:
!  - Analytic definition of the boundary
!  - Analytic boundary conditions
!  - Information from the INI/DAT files
!  - Main solution and RHS vector
!  - Information for preconditioning in nonlinear loop
!  - For every level in the discretisation:
!    - Triangulation
!    - Block discretisation
!    - Matrices and
!    - Temporary vectors
! This problem structure is available only in the top-level routines of
! the CC2D module; it's not available in any callback routines.
! Parameters for callback routines are passed through the collection
! structure rcollection, which is part of the problem structure.
!
!
!******************************************************************************
! The collection structure t_problem%rcollection
!******************************************************************************
!
! The collection structure collects all information that is necessary to be
! passed to callback routines. This e.g. allows to pass matrices/vectors/
! constants from the main problem.
! This structure contains the following informaion, which is added by the
! initialisation routines to the collection:
!
! Global, level independent data:
!
! Name                  | Description
! ----------------------+------------------------------------------------------
! NU                    | Reciprocal 1/RE of parameter RE from the DAT file
! IBOUNDARY             ! =0: stationary Dirichlet boundary conditions
!                       ! =1: nonstationary boundary conditions, possibly with
!                       !     pressure drop
!                       ! =2: nonstationary boundary conditions, possibly with
!                       !     prssure drop and/or Neumann boundary parts
! UPSAM                 | Stabilisation parameter
!
! On every level between NLMIN and NLMAX:
!
! Expressions for boundary conditions are saved in special
! sections in the collection:
!
! Section [BDEXPRESSIONS]:
! ------------------------
! Saves the expressions that are to be evaluated on the boundary.
!
! During the assembly of RHS vectors and boundary conditions, the following
! additional information is valid in the collection:
!
!   collection%IquickAccess(1)   = 0: stationary,
!                                  1: nonstationary with explicit time stepping
!   collection%DquickAccess(1)   = current simulation time
!   collection%DquickAccess(2)   = minimum simulation time
!   collection%DquickAccess(3)   = maximum simulation time
  
end module

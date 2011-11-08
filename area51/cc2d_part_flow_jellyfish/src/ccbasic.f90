!##############################################################################
!# ****************************************************************************
!# <name> ccbasic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic problem definitions for the cc2d solver.
!# The basic structure and content of the different structures
!# are described here.
!# </purpose>
!##############################################################################

module ccbasic

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
  use discretebc
  use discretefbc
  use linearsystemscalar
  use linearsystemblock
  use geometry
  
  use collection
  
  use adaptivetimestep
    
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
    
    ! Precomputed matrix for edge stabilisation. Only active if
    ! iupwind = CCMASM_STAB_FASTEDGEORIENTED.
    type(t_matrixScalar) :: rmatrixStabil
    
    ! A prolongation matrix for the velocity.
    type(t_matrixScalar) :: rmatrixProlVelocity
    
    ! A prolongation matrix for the pressure.
    type(t_matrixScalar) :: rmatrixProlPressure
    
    ! An interpolation matrix for the velocity.
    type(t_matrixScalar) :: rmatrixInterpVelocity
    
    ! An interpolation matrix for the pressure.
    type(t_matrixScalar) :: rmatrixInterpPressure
    
  end type

!</typeblock>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation of the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation

    ! An object specifying the block discretisation to be used for (edge)
    ! stabilisation. Only used if edge stabilisation is activated, otherwise
    ! this coincides with rdiscretisation.
    type(t_blockDiscretisation) :: rdiscretisationStabil

    ! Temporary vector in the size of the RHS/solution vector on that level.
    type(t_vectorBlock) :: rtempVector

    ! A variable describing the discrete boundary conditions fo the velocity.
    ! Points to NULL until the BC`s are discretised for the first time.
    type(t_discreteBC), pointer :: p_rdiscreteBC => null()
  
    ! A structure for discrete fictitious boundary conditions
    ! Points to NULL until the BC`s are discretised for the first time.
    type(t_discreteFBC), pointer :: p_rdiscreteFBC => null()
    
    ! This flag signales whether there are Neumann boundary components
    ! visible on the boundary of this level or not. If there are no
    ! Neumann boundary components visible, the equation gets indefinite
    ! for the pressure.
    logical :: bhasNeumannBoundary
    
    ! A structure containing all static information about this level.
    type(t_staticLevelInfo) :: rstaticInfo
    
  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for the nonstationary Nav.St. problem">

  type t_problem_explTimeStepping
  
    ! Number of current time step; changes during the simulation.
    integer :: itimeStep           = 0
    
    ! Current simulation time; changes during the simulation.
    real(DP) :: dtime              = 0.0_DP
  
    ! Maximum number of time steps; former NITNS
    integer :: niterations         = 0
    
    ! Absolute start time of the simulation
    real(DP) :: dtimeInit          = 0.0_DP
    
    ! Time step size; former TSTEP
    real(DP) :: dtimeStep          = 0.0_DP
    
    ! Maximum time of the simulation
    real(DP) :: dtimeMax           = 0.0_DP
  
    ! Lower limit for the time derivative to be treated as zero. Former EPSNS.
    ! Simulation stops if time derivative drops below this value.
    real(DP) :: dminTimeDerivative = 0.00001_DP
    
    ! Configuration block for the adaptive time stepping.
    type(t_adaptimeTimeStepping) :: radaptiveTimeStepping
    
  end type

!</typeblock>


!<typeblock description="Application-specific type block configuring the stabilisation">

  type t_problem_stabilisation
  
    ! Type of stabilisation. =0: Streamline diffusion. =1: Upwind, 2=edge oriented.
    integer :: iupwind = 0
    
    ! Stabilisation parameter for the nonlinearity.
    ! Standard values: Streamline diffusion: 1.0. Upwind: 0.1. Edge oriented: 0.01.
    real(DP) :: dupsam = 1.0_DP
    
    ! Calculation of local H for streamline diffusion
    integer :: clocalH = 0
    
  end type

!</typeblock>

!<typeblock>

  ! Object that collects statistics gathered during the simulation.
  type t_cc_statistics
    
    ! Total number of nonlinear iterations
    integer :: nnonlinearIterations = 0
    
    ! Total number of iterations in the linear solver
    integer :: nlinearIterations = 0
    
    ! Total number of calculated timesteps
    integer :: ntimesteps = 0

    ! Total time used for nonlinear solver
    real(DP) :: dtimeNonlinearSolver = 0.0_DP
    
    ! Total time used for linear solver
    real(DP) :: dtimeLinearSolver = 0.0_DP

    ! Total time used for symbolical/numerical factorisation
    ! in the preparation of the linear solver
    real(DP) :: dtimeLinearSolverFactorisation = 0.0_DP
    
    ! Total time for matrix assembly
    real(DP) :: dtimeMatrixAssembly = 0.0_DP
    
    ! Total time for defect vector calculation
    real(DP) :: dtimeDefectCalculation = 0.0_DP
    
    ! Total time for grid generation
    real(DP) :: dtimeGridGeneration = 0.0_DP
    
    ! Total time for postprocessing
    real(DP) :: dtimePostprocessing = 0.0_DP
    
    ! Total time for calculating optimal correction factors in the
    ! nonlinear iteration
    real(DP) :: dtimeOptimalCorrection = 0.0_DP
    
    ! Total time for solver (without initial initialisation)
    real(DP) :: dtimeSolver = 0.0_DP
    
    ! Total time
    real(DP) :: dtimeTotal = 0.0_DP
    
  end type

!</typeblock>


!<typeblock description="Application-specific type block for the Nav.St. problem">

  type t_problem
  
    ! Output level during the initialisation phase.
    integer :: MSHOW_Initialisation
  
    ! Output level of the application.
    integer :: MT_OutputLevel
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: NLMIN
    
    ! Maximum refinement level
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
    real(dp) :: dx
    
    real(dp) :: dy
    
    real(dp) :: drad
    
    real(dp) :: drho2
    
    ! Type of problem.
    ! =0: Stokes.
    ! =1: Navier-Stokes.
    integer :: iequation
    
    ! Type of subproblem of the main problem. Depending on iequationType.
    ! If iequationType=0 or =1:
    ! =0: (Navier-)Stokes with gradient tensor
    ! =1: (Navier-)Stokes with deformation tensor
    integer :: isubEquation
    
    ! Model for the viscosity.
    ! =0: Constant viscosity.
    ! =1: Power law nu = nu_0 * z^(dviscoexponent/2 - 1), nu_0 = 1/RE, z=||D(u)||^2+dviscoEps
    integer :: cviscoModel
        
    ! Exponent parameter for the viscosity model
    real(DP) :: dviscoexponent

    ! Epsilon regularisation for the viscosity model
    real(DP) :: dviscoEps

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! iboundary parameter from the DAT file. Specifies whether to update
    ! boundary conditions in nonstationary simulations or not.
    ! =0: stationary Dirichlet boundary conditions
    ! =1: nonstationary boundary conditions, possibly with
    !     pressure drop
    ! =2: nonstationary boundary conditions, possibly with
    !     prssure drop and/or Neumann boundary parts
    integer :: iboundary

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
    
    ! A configuration block for the stabilisation of the convection.
    type(t_problem_stabilisation) :: rstabilisation
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection)                    :: rcollection
    
    ! A param list that saves all parameters from the DAT/INI file(s).
    type(t_parlist)                       :: rparamList

    ! A statistics structure gathering statistical data about the
    ! simulation.
    type(t_cc_statistics)                 :: rstatistics
    
    real(dp),dimension(2) :: DResForceX
    
    real(dp),dimension(2) :: DResForceY

    real(dp),dimension(2) :: DTrqForce
    
    real(dp),dimension(2) :: DVelX

    real(dp),dimension(2) :: DVelY
    
    real(dp),dimension(2) :: DAngVel

    real(dp) :: dCoefficientDrag

    real(dp) :: dCoefficientLift
    
    type(t_particleCollection) :: rparticleCollection

  end type

!</typeblock>

!</types>

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
! the CC2D module; it is not available in any callback routines.
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

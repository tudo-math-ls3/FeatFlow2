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
  
  use collection
  
  use adaptivetimestep
    
  implicit none
  
!<types>

!<typeblock>

  ! A type block specifying all 'template' information which are depending
  ! on a discretisation and a triangulation. Such template information can be
  ! precalculated and is valid until the mesh or the FE spaces change.
  ! It can be used to create vectors and matrices e.g.
  type t_asmTemplates
  
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

    ! Precalculated mass matrix for the velocity space.
    type(t_matrixScalar) :: rmatrixMass

    ! Precalculated penalty matrix for the velocity space.
    type(t_matrixScalar) :: rmatrixPenalty

    ! Precalculated mass matrix for the pressure space.
    type(t_matrixScalar) :: rmatrixMassPressure

    ! An object specifying the discretisation in the velocity space
    ! to be used for (edge) stabilisation. Only used if edge stabilisation
    ! is activated, otherwise this coincides with the default discretisation
    ! of the velocity space.
    type(t_spatialDiscretisation) :: rdiscretisationStabil

    ! An object specifying the discretisation in the velocity space
    ! to be used for penalty matrix.
    type(t_spatialDiscretisation) :: rdiscretisationPenalty
    
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

!<typeblock description="General assembly information">

  type t_matrixAssembly
  
    ! Cubature formula tag for the Stokes matrices
    integer :: icubA = CUB_GEN_AUTO

    ! Cubature formula tag for the B-matrices matrices
    integer :: icubB = CUB_GEN_AUTO

    ! Cubature formula tag for the mass matrices
    integer :: icubM = CUB_GEN_AUTO

    ! Cubature formula tag for the penalty matrices
    integer :: icubMp = CUB_GEN_AUTO

    ! Cubature formula tag for the penalty matrices 
    ! with adaptive cub. form.
    integer :: icubMpa = CUB_GEN_AUTO
    integer :: icubrefin = 1

    ! Use mass lumping for the velocity
    logical :: bmassLumpingVelocity = .false.
    
    ! Type of mass lumping
    integer :: imassLumpTypeVelocity = LSYSSC_LUMP_DIAG
  
  end type

!</typeblock>

!<typeblock>

  ! Type block configuring the RHS.
  type t_rhsAssembly
    
    ! Cubature formula tag for the right hand side vectors
    integer :: icubF = CUB_GEN_AUTO

    ! Type of the RHS.
    ! =0: zero.
    ! =1: stationary RHS, analytically defined by coeff_RHS_x/coeff_RHS_y
    ! =2: nonstationary RHS, analytically defined by coeff_RHS_x/coeff_RHS_y
    ! =3; Stationary RHS prescribed by a file.
    ! =4: Nonstationary RHS prescribed by a sequence of files.
    !     The files are linearly interpolated over the whole current
    !     space-time cylinder.
    integer :: ctype = 0
    
    ! If ctype=3/4: Basic filename of the files configuring the RHS.
    character(len=SYS_STRLEN) :: sfilename = ""
    
    ! If ctype=4: First index in the filename.
    ! The filename sfilename is extended by an index ".00001", ",.00002",...
    integer :: ifirstindex = 0
    
    ! If ctype=4: Number of files in the sequence.
    integer :: inumfiles = 0
    
    ! Whether or not files are formatted.
    integer :: iformatted = 1
    
    ! If ctype=3/4: Block vector representing the stationary RHS from the file.
    type(t_vectorBlock) :: rrhsVector

    ! If ctype=4: A second RHS vector used for interpolating.
    type(t_vectorBlock) :: rrhsVector2
    
    ! If ctype=4: Index of the RHS currently loaded in rrhsVector.
    ! =-1 if this is undetermined.
    integer :: icurrentRhs = -1
    
    ! If ctype=4: Index of the RHS currently loaded in rrhsVector2.
    ! =-1 if this is undetermined.
    integer :: icurrentRhs2 = -1
    
    ! For irhs=4: start time of the RHS.
    real(DP) :: dtimeInit = 0.0_DP

    ! For irhs=4: end time of the RHS.
    real(DP) :: dtimeMax = 0.0_DP
    
    ! Multiplier for the X-component
    real(DP) :: dmultiplyX = 1.0_DP
    
    ! Multiplier for the Y-component
    real(DP) :: dmultiplyY = 1.0_DP

  end type

!</typeblock>

!<typeblock description="Type block defining dynamic information about a level that change in every timestep">

  type t_dynamicLevelInfo
  
    ! A variable describing the discrete boundary conditions fo the velocity.
    type(t_discreteBC) :: rdiscreteBC
  
    ! A structure for discrete fictitious boundary conditions
    type(t_discreteFBC) :: rdiscreteFBC
    
    ! This flag signales whether there are Neumann boundary components
    ! visible on the boundary of this level or not. If there are no
    ! Neumann boundary components visible, the equation gets indefinite
    ! for the pressure.
    logical :: bhasNeumannBoundary = .false.
    
    ! Handle to a list of edges with Dirichlet boundary conditions on one of the
    ! velocity components. =ST_NOHANDLE if there are no Dirichlet boundary segments.
    integer :: hedgesDirichletBC = ST_NOHANDLE
    
    ! Number of edges with Dirichlet boudary conditions.
    integer :: nedgesDirichletBC = 0
    
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

    ! An object specifying the block discretisation to be used for Penalty.
    type(t_blockDiscretisation) :: rdiscretisationPenalty

    ! Temporary vector in the size of the RHS/solution vector on that level.
    type(t_vectorBlock) :: rtempVector

    ! A structure containing all template information about this level.
    type(t_asmTemplates) :: rasmTempl

    ! A structure containing all dynamic information about this level.
    type(t_dynamicLevelInfo) :: rdynamicInfo
    
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

!<typeblock>

  ! This type block encapsules all physical constants and configuration
  ! parameters for the primal equation. This includes e.g. the type of the equation,
  ! viscosity parameter etc.
  type t_problem_physics
  
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu
    
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
    ! =1: Power law: nu = nu_0 * z^(dviscoexponent/2 - 1), nu_0 = 1/RE, z=||D(u)||^2+dviscoEps
    ! =2: Bingham fluid: nu = nu_0 + dviscoyield / sqrt(|D(u)||^2+dviscoEps^2), nu_0 = 1/RE
    integer :: cviscoModel
        
    ! Exponent parameter for the viscosity model
    real(DP) :: dviscoexponent

    ! Epsilon regularisation for the viscosity model
    real(DP) :: dviscoEps
    
    ! Yield stress for Bingham fluid
    real(DP) :: dviscoYield

  end type

!</typeblock>


!<typeblock description="Application-specific type block configuring the stabilisation">

  type t_problem_stabilisation
  
    ! Type of stabilisation. =0: Streamline diffusion. =1: Upwind, 2=edge oriented.
    integer :: iupwind = 0
    
    ! Cubature formula for the EOJ stabilisation.
    integer(I32) :: ccubEOJ = CUB_G4_1D
    
    ! Stabilisation parameter for the nonlinearity.
    ! Standard values: Streamline diffusion: 1.0. Upwind: 0.1. Edge oriented: 0.01.
    real(DP) :: dupsam = 1.0_DP
    
    ! 2nd Relaxation parameter in the jump stabilisation. Standard = 0.0
    real(dp) :: dUpsamStar = 0.0

    ! Exponent for edge length weight in the jump stabilisation. Standard = 2.0
    ! (corresponds to a weight dupsam*h^2).
    real(dp) :: deojEdgeExp = 2.0
    
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

    ! Total time for RHS assembly
    real(DP) :: dtimeRHSAssembly = 0.0_DP
    
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
    
    ! A block containing the physics of the problem.
    type(t_problem_physics) :: rphysics

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
    
    ! Assembly structure for the matrices.
    type(t_matrixAssembly) :: rmatrixAssembly
    
    ! Assembly structure specifying the RHS.
    type(t_rhsAssembly) :: rrhsAssembly

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
    integer :: itimedependence
    
    ! A parameter block for everything that controls the time dependence.
    ! Only valid if itimedependence=1!
    type(t_problem_explTimeStepping) :: rtimedependence
    
    ! A configuration block for the stabilisation of the convection.
    type(t_problem_stabilisation) :: rstabilisation
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection
    
    ! A param list that saves all parameters from the DAT/INI file(s).
    type(t_parlist) :: rparamList

    ! A statistics structure gathering statistical data about the
    ! simulation.
    type(t_cc_statistics) :: rstatistics
    
    ! Penalty parameters
    real(DP) :: dLambda 
    integer :: ipenalty
    
    type(t_particleCollection) :: rparticleCollection
    
    integer  :: iParticles
    
    real(dp) :: dCoefficientDrag

    real(dp) :: dCoefficientLift
    
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

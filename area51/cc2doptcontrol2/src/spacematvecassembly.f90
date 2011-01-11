!##############################################################################
!# ****************************************************************************
!# <name> spacematvecassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the very basic matrix assembly routines for the
!# core equation. It's independent of any nonlinear iteration and provides
!# just one functionality: Assemble a matrix or a vector based on a given 
!# set of parameters.
!#
!# The discretised core equation reads at the moment:
!#
!#  $$           A_1 y   +  \eta_1 B p   +          R lambda      = f_1 $$
!#  $$ \tau_{11} B^T y   +  \kappa_{11} I p                       = f_2 $$
!#
!#  $$ (R_{22}+R2_{22}) y +  A_{22} \lambda  + \eta_{22}   B \xi           = f_3 $$
!#  $$               \tau_{22} B^T \lambda  + \kappa_{22} I \xi           = f_4 $$
!#
!# with
!#
!#   $$ A_1 = \iota_{11} I  +  \alpha_{11} M  +  \theta_{11} L  +  
!#            \gamma_{11} N(y) + dnewton_{11} N* + gammaT_{11} N^t(y) + dnewtonT_{11} N*^t(y)(y)$$
!#   $$ A_2 = \iota_{22} I  +  \alpha_{22} M  +  \theta_{22} L  +  
!#            \gamma_{22} N(y) + dnewton_2 N*(y) + \gammaT_{22} N^t(y) + dnewtonT_{22} N*^t(y)$$
!#   $$ R_1 = \alpha_{12} (P)M(\lambda) $$
!#   $$ R_2 = \alpha_{21} M  +  
!#            \gamma_{21} N(\lambda) + dnewton_{21} N*(\lambda) + 
!#            gammaT_{21} N^t(\lambda) + dnewtonT_{21} N*^t(\lambda) $$
!#   $$ R2_2 =  dnewton_{21} N*(\lambda) + gammaT_{21} N^t(\lambda) $$
!#  
!# and
!#
!#   $I$     = identity matrix,
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity including stabilisation, depending on the 
!#             primal velocity, i.e.
!#                   $$ (y\Delta)\cdot $$
!#   $N*(y)$ = Newton matrix, depending on the primal velocity, i.e.
!#                  $$ (\Delta y)\cdot $$
!#   $N^t(y)$ = Transposed nonlinearity including stabilisation, depending on the 
!#              primal velocity, i.e.
!#                   $$ (y\Delta)^t\cdot $$
!#   $N*^t(y)$ = Transposed Newton matrix, depending on the primal velocity, i.e.
!#                  $$ (\Delta y)^t\cdot $$
!#   
!#   $\iota_{ij}$  = 0/1     - switches the identity matrix on/off,
!#   $\alpha_{ij}$ = 0/1     - switches the mass matrices on/off;
!#                             =0 for stationary problem,
!#   $\theta_{íj}$           - weight for the Laplace matrix,
!#   $\gamma_{ij}$ = 0/1     - Switches the nonlinearity on/off;
!#                             =0 for Stokes system,
!#   $dnewton_{ij} \in R$    - Switches the Newton matrix on/off.
!#   $\eta_{ij}$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau_{ij}$   = 0/1     - Switches the 'B^T'-term on/off,
!#   $\kappa_{ij}$ = 0/1     - Switches of the identity matrix I for the pressure
!#                             in the continuity equation
!#
!# Note that the nonlinear part is always dependent on the primal velocity --
!# for the primal equation as well as for the dual one!
!#
!# To assemble such a matrix, the application has to follow two steps:
!# 
!# 1.) Create a structure of type t_nonlinearSpatialMatrix and set the
!#     parameters in this structure according to the matrix which should
!#     be assembled.
!#
!# 2.) Call the matrix/vector assembly routine with this structure as input.
!#
!# The matrix assembly routine will return an initialised matrix structure
!# of a system matrix based on the provided parameters.
!#
!# The module contains the following routines:
!#
!# 1.) smva_assembleMatrix
!#     -> Assembles a matrix based on a set of input parameters.
!#
!# 2.) smva_assembleDefect
!#     -> Set up a defect vector d:=b-A(x)x
!#
!# 3.) smva_projectControlTimestep
!#     -> Projects the entries of a vector into a range of numbers.
!#
!# 4.) smva_initNonlinMatrix
!#     -> Initialises a nonlinear-matrix structure with basic parameters.
!#
!# 5.) smva_initNonlinearData
!#     -> Initialises a nonlinear-data structure.
!# </purpose>
!##############################################################################

module spacematvecassembly

  use fsystem
  use storage
  use genoutput
  use basicgeometry
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use scalarpde
  use linearsystemscalar
  use linearsystemblock
  use linearsolver
  use bilinearformevaluation
  use linearformevaluation
  use linearsolverautoinitialise
  use matrixrestriction
  use trilinearformevaluation
  use matrixmodification
  use matrixio
  use feevaluation
  use convection
  use collection
  use mprimitives
  
  use constantsoptc
  use assemblytemplates
  use assemblytemplatesoptc
  use structuresoptc
  
  use structuresoptflow

  use optcontrolconvection
    
  implicit none
  
  private
  
!<constants>

!<constantblock description="Identifiers for the 'coperation' input parameter of the matrix assembly routine">

  ! Allocate memory if necessary.
  integer(I32), parameter, public :: CCMASM_ALLOCMEM              = 1
  
  ! Compute all matrix entries.
  integer(I32), parameter, public :: CCMASM_COMPUTE               = 2
  
  ! Allocate memory and compute matrix entries.
  integer(I32), parameter, public :: CCMASM_ALLOCANDCOMPUTE       = 3
  
  ! Bypass memory allocation for matrices.
  integer(I32), parameter, public :: CMASM_QUICKREFERENCES        = 4
  
  ! Don't clear the old matrix
  integer(i32), parameter, public :: CMASM_NOCLEAR                = 8
  
!</constantblock>

!<constantblock description="Identifiers for the IUPWIND parameter that specifies how to set up the nonlinearity or stabilisation.">

  ! Streamline diffusion; configured by dupsam
  integer, parameter, public :: CCMASM_STAB_STREAMLINEDIFF    = 0

  ! 1st-order upwind; configured by dupsam
  integer, parameter, public :: CCMASM_STAB_UPWIND            = 1
  
  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'
  integer, parameter, public :: CCMASM_STAB_EDGEORIENTED      = 2

  ! Streamline diffusion; configured by dupsam, new implementation
  integer, parameter, public :: CCMASM_STAB_STREAMLINEDIFF2   = 3

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation
  integer, parameter, public :: CCMASM_STAB_EDGEORIENTED2     = 4

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation.
  ! Precomputed matrix.
  integer, parameter, public :: CCMASM_STAB_EDGEORIENTED3     = 5
!</constantblock>

!<constantblock description="Matrix type ID's specifying the general matrix class to set up.">

  ! Standard matrix.
  integer, parameter, public :: CCMASM_MTP_AUTOMATIC         = 0
  
  ! Standard matrix with decoupled velocity blocks
  integer, parameter, public :: CCMASM_MTP_DECOUPLED         = 1
  
  ! Extended 'full-tensor' matrix with submatrices A11, A12, A21, A22, all independent from
  ! each other.
  integer, parameter, public :: CCMASM_MTP_FULLTENSOR        = 2

!</constantblock>

!<constantblock description="Matrix types"

  ! Optimal control matrix (6 equations)
  integer, parameter, public :: MATT_OPTCONTROL = 0
  
  ! Linearised optimal control matrix (6 equations, Newton type)
  integer, parameter, public :: MATT_LINOPTCONTROL = 1

  ! Primal equation matrix (3 equations)
  integer, parameter, public :: MATT_PRIMAL = 2
  
  ! Linearised primal equation matrix (3 equations, Newton type)
  integer, parameter, public :: MATT_LINPRIMAL = 3

  ! Dual equation matrix (3 equations)
  integer, parameter, public :: MATT_DUAL = 4
  
  ! Linearised dual equation matrix (3 equations, Newton type)
  integer, parameter, public :: MATT_LINDUAL = 5


!</constantblock>

!</constants>

  public :: t_nonlinearSpatialMatrix
  public :: t_spatialMatrixNonlinearData
  public :: t_matrixAssemblyFlags
  public :: smva_initNonlinMatrix
  public :: smva_assembleMatrix
  public :: smva_assembleDefect
  public :: smva_projectControlTstepConst
  public :: smva_projectControlTstepVec
  public :: smva_getDiscrData
  public :: smva_initNonlinearData
  public :: smva_addBdEOJvector
  public :: smva_addBdEOJOperator

!<types>

!<typeblock>
  ! Structure that encapsules all nonlinearities that may appear in the nonlinear
  ! spatial matrix.
  type t_spatialMatrixNonlinearData
  
    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'previous' timestep. If there is no previous timestep (e.g.
    ! like in the 0th timestep), the vector can be undefined.
    type(t_vectorBlock), pointer :: p_rvector1 => null()

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'current' timestep. 
    type(t_vectorBlock), pointer :: p_rvector2 => null()

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'next' timestep. If there is no next timestep (e.g.
    ! like in the last timestep), the vector can be undefined.
    type(t_vectorBlock), pointer :: p_rvector3 => null()

  end type
  
!</typeblock>

!<typeblock>

  ! All discretisation and level related data structures that are necessary
  ! to evaluate the matrix or create a defect.
  type t_spatialMatrixDiscrData
  
    ! The physics of the problem
    type(t_settings_physics) :: rphysicsPrimal

    ! Stabilisation parameters for the primal and dual system.
    type(t_settings_stabil) :: rstabilPrimal
    type(t_settings_stabil) :: rstabilDual
    
    ! Structure defining constraints of the problem.
    type(t_optcConstraintsSpace) :: rconstraints

    ! Discretisation of the level, primal space.
    type(t_blockDiscretisation), pointer :: p_rdiscrPrimal => null()

    ! Discretisation of the level, primal/dual space.
    type(t_blockDiscretisation), pointer :: p_rdiscrPrimalDual => null()
    
    ! Assembly template data on that level.
    type(t_staticSpaceAsmTemplates), pointer :: p_rstaticAsmTemplates => null()

    ! Assembly template data on that level for the optimal control problem.
    type(t_staticSpaceAsmTemplatesOptC), pointer :: p_rstaticAsmTemplatesOptC => null()

    ! Pointer to program wide debug flags.
    type(t_optcDebugFlags), pointer :: p_rdebugFlags => null()

  end type

!</typeblock>

  public :: t_spatialMatrixDiscrData

!<typeblock>

!<!--
! The following structure encapsules all coefficients in front of the terms in each 
! submatrix of the global matrix. Here, we don't address each submatrix in the
! global matrix but provide coefficients always for a complete block of all dimensions
! for the primal as well as for the dual term.
!
! Example: The Dalpha, Deta and Dtau coefficients are used in the following way:
!
!    ( Dalpha(1,1)*M   Dalpha(1,1)*M   Deta(1)*B1    Dalpha(1,2)*M   Dalpha(1,2)*M              )
!    ( Dalpha(1,1)*M   Dalpha(1,1)*M   Deta(1)*B2    Dalpha(1,2)*M   Dalpha(1,2)*M              )
!    (  Dtau(1)*B1^t    Dtau(1)*B2^t                                                            )
!    ( Dalpha(2,1)*M   Dalpha(2,1)*M                 Dalpha(2,2)*M   Dalpha(2,2)*M   Deta(2)*B1 )
!    ( Dalpha(2,1)*M   Dalpha(2,1)*M                 Dalpha(2,2)*M   Dalpha(2,2)*M   Deta(2)*B2 )
!    (                                                Dtau(2)*B1^t    Dtau(2)*B2^t              )
!
! The use of the other parameters is similar.
! -->

  ! This structure describes a nonlinear matrix in space.
  ! The nonlinearity usually stems from a dependence of the Navier-Stokes operator
  ! on a solution (Oseen system) and stabilisation operators.
  type t_nonlinearSpatialMatrix
  
    ! IOTA-parameters that switch the identity on/off.
    real(DP), dimension(2,2) :: Diota = 0.0_DP
  
    ! ALPHA-parameters that switch the mass matrix on/off.
    real(DP), dimension(2,2) :: Dalpha = 0.0_DP
  
    ! THETA-parameters that switch the Stokes matrix on/off
    real(DP), dimension(2,2) :: Dtheta = 0.0_DP
    
    ! GAMMA-parameters that switch the nonlinearity (y*grad(.)) on/off
    real(DP), dimension(2,2) :: Dgamma = 0.0_DP
  
    ! NEWTON-parameters that switch the Newton term ((.)*grad(y)) on/off
    real(DP), dimension(2,2) :: Dnewton = 0.0_DP

    ! NEWTON-parameters that switch a 2nd Newton term ((.)*grad(y)) on/off
    real(DP), dimension(2,2) :: Dnewton2 = 0.0_DP

    ! GAMMAT-parameters that switch the transposed nonlinearity (y*grad(.)^T) on/off
    real(DP), dimension(2,2) :: DgammaT = 0.0_DP

    ! GAMMAT-parameters that switch a 2nd transposed nonlinearity (y*grad(.)^T) on/off
    real(DP), dimension(2,2) :: DgammaT2 = 0.0_DP
  
    ! NEWTONT-parameters that switch the transposed Newton term ((.)*grad(y)^T) on/off
    real(DP), dimension(2,2) :: DnewtonT = 0.0_DP
    
    ! ETA-parameters that switch the B-terms on/off.
    real(DP), dimension(2,2) :: Deta = 0.0_DP
    
    ! TAU-parameters that switch the B^T-terms on/off
    real(DP), dimension(2,2) :: Dtau = 0.0_DP
    
    ! KAPPA-parameters that switch the I matrix in the continuity equation
    ! on/off.
    real(DP), dimension(2,2) :: Dkappa = 0.0_DP
    
    ! When evaluating nonlinear terms, the evaluation routine accepts
    ! in timestep i three solution vectors -- one corresponding to
    ! the matrix Aii-1, one corresponding to matrix Aii and one corresponding
    ! to matrix Aii+1. The following id's specify which primal solution is
    ! used for which matrix. A value of 1 correspponds to solution i-1,
    ! a value of 2 to solution i and a value of 3 to solution i+1.
    ! iprimalSol=1 means e.g. that A=A(xi-1). Similarly,
    ! iprimalSol=3 would mean that A=A(xi+1).
    ! =0: undefined.
    integer :: iprimalSol  = 2

    ! The following three variables specify the id of the dual solution to be
    ! evaluated. 1=lambda_{i-1}, 2=lambda_i, 3=lambda_{i+1}.
    ! =0: undefined.
    integer :: idualSol  = 2

    ! The following three variables specify the id of a 2nd dual solution to be
    ! evaluated. 1=lambda_{i-1}, 2=lambda_i, 3=lambda_{i+1}.
    ! =0: undefined.
    ! This is used for the discretisation of the 2nd Newton/GammT-term which
    ! appears on the main diagonal of the system for the dual equation
    ! if Crank-Nicolson is active.
    integer :: idualSol2  = 2

    ! Regularisation parameter ALPHA which controls the transformation
    ! from the dual variable $\lambda$ to the control $u$:
    ! $u=-1/\alpha \lambda$.
    real(DP) :: dalphaC = 0.0_DP
    
    ! Type of this matrix. One of the MATT_xxxx constants.
    integer :: cmatrixType = 0

    ! Discretisation related data.
    type(t_spatialMatrixDiscrData) :: rdiscrData

    ! Pointer to a structure specifying all nonlinearities.
    ! May point to NULL if there is no nonlinearity.
    type(t_spatialMatrixNonlinearData), pointer :: p_rnonlinearity => null()

  end type

!</typeblock>

!<typeblock>
  ! Additional status flags that configure how to assemble a matrix.
  type t_matrixAssemblyFlags
  
    ! MATRIX RESTRICTION: Parameter to activate matrix restriction.
    ! Can be used to generate parts of the matrices on coarse grids where the
    ! aspect ratio of the cells is large. Only applicable for $\tilde Q_1$
    ! discretisations.
    ! Standard = 0 = deactivate matrix restriction
    integer :: iadaptiveMatrices = 0
    
    ! MATRIX RESTRICTION: Threshold parameter for adaptive matrix generation
    ! of coarse grid matrices (here: maximum aspect ratio). 
    ! Only applicable if iadaptiveMatrices <> 0.
    ! Standard = 20.0
    real(DP) :: dadmatthreshold = 20.0_DP
    
    ! ASSEMBLY SPECIALS: If set to TRUE, the D-matrices are assembled as 
    ! virtually transposed B-matrices. With the standard setting FALSE, the
    ! D-matrices are assembled as standard divergence matrices.
    ! WARNING: If matrices are to be assembled as virtually transposed
    ! matrices, they must not be overwritten, as they share their data with
    ! other matrices.
    ! If the matrix is assembled in the standard way without virtual
    ! transposition, separate memory is requested for the matrices and so
    ! they are allowed to be overwritten if necessary. (Neessary for UMFPACK
    ! on a pure Dirichlet problem.)
    logical :: bvirtualTransposedD = .false.

  end type
!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_getDiscrData (rsettings, ilevel, rdiscrData,rphysics)
  
!<description>
  ! Fetches all discretisation data from the main program structure that is
  ! necessary to set up nonlinear matrices.
!</description>

!<input>
  ! The structure of the main solver
  type(t_settings_optflow), intent(in), target :: rsettings
  
  ! Level where the discretisation takes place.
  integer, intent(in) :: ilevel

  ! OPTIONAL: Alternative physics definition to use.
  ! If not present, the standard global physics settings are used.
  type(t_settings_physics), intent(in), optional :: rphysics
!</input>

!<output>
  ! Discretisation related data, can be used to generate nonlinear matrices.
  type(t_spatialMatrixDiscrData), intent(out) :: rdiscrData
!</output>

!</subroutine>

    ! Get level-independent data.
    rdiscrData%rphysicsPrimal = rsettings%rphysicsPrimal
    if (present(rphysics)) then
      rdiscrData%rphysicsPrimal = rphysics
    end if
    rdiscrData%rstabilPrimal = rsettings%rstabilPrimal
    rdiscrData%rstabilDual = rsettings%rstabilDual
    rdiscrData%rconstraints = rdiscrData%rconstraints

    rdiscrData%p_rdiscrPrimal => &
        rsettings%rfeHierPrimal%p_rfeSpaces(ilevel)%p_rdiscretisation

    rdiscrData%p_rdiscrPrimalDual => &
        rsettings%rfeHierPrimalDual%p_rfeSpaces(ilevel)%p_rdiscretisation
    
    rdiscrData%p_rstaticAsmTemplates => &
        rsettings%rspaceAsmHierarchy%p_RasmTemplList(ilevel)

    rdiscrData%p_rstaticAsmTemplatesOptC => &
        rsettings%rspaceAsmHierarchyOptC%p_RasmTemplList(ilevel)

    rdiscrData%p_rdebugFlags => rsettings%rdebugFlags

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_initNonlinMatrix (rnonlinearSpatialMatrix,rdiscrData,rnonlinearity)

!<description>
  ! Initialises the rnonlinearCCMatrix structure.
!</description>

!<input>
  ! Discretisation related data with information about the level, the
  ! matrix acts on.
  type(t_spatialMatrixDiscrData) :: rdiscrData

  ! Structure specifying the nonlinearity evaluation points of the matrix.
  type(t_spatialMatrixNonlinearData), target :: rnonlinearity
!</input>

!<inputoutput>
  ! Nonlinear matrix structure.
  ! Basic parameters in this structure are filled with data.
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>
              
!</subroutine>
    
    ! Initialise the matrix assembly structure  
    ! with basic global information.
    rnonlinearSpatialMatrix%rdiscrData = rdiscrData

    ! Change the sign of dupsam2 for a consistent stabilisation.
    ! Reason: The stablisation is added to the dual operator by the SD/
    ! EOJ stabilisation in the following way:
    !
    !    ... - (u grad lamda + dupsam2*stabilisation) + ... = rhs
    !
    ! We want to *add* the stabilisation, so we have to introduce a "-" sign
    ! in dupsam2 to get
    !
    !    ... - (u grad lamda) - (-dupsam2*stabilisation) + ... = rhs
    ! <=>
    !    ... - (u grad lamda) + dupsam2*stabilisation + ... = rhs
    
    rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam = &
        -rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam
        
    ! Remember the evaluation point of the nonlinearity.
    rnonlinearSpatialMatrix%p_rnonlinearity => rnonlinearity
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_nonlinMatrixSetNonlin (rnonlinearSpatialMatrix,rnonlinearity)

!<description>
  ! Defines the nonlinearity of a nonlinear matrix. Allows to change the
  ! nonlinerity in rnonlinearSpatialMatrix to rnonlinearity.
!</description>

!<input>
  ! Structure specifying the nonlinearity evaluation points of the matrix.
  type(t_spatialMatrixNonlinearData), target :: rnonlinearity
!</input>

!<inputoutput>
  ! Nonlinear matrix structure.
  ! Basic parameters in this structure are filled with data.
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>
              
!</subroutine>
    
    ! Remember the evaluation point of the nonlinearity.
    rnonlinearSpatialMatrix%p_rnonlinearity => rnonlinearity
    
  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_ProjMass (rdiscretisationTrial,rdiscretisationTest,rform, &
      nelements,npointsPerElement,Dpoints, IdofsTrial,IdofsTest,rdomainIntSubset, &
      Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the projective mass matrix.
    !
    ! The coefficients is c=c(lambda_1) or =c(lambda_2) with
    ! c=1/alpha if a < -1/alpha lambda_i < b and c=0 otherwise. This is the derivative
    ! of the projection operator "-P[a,b](-1/alpha lambda_i)" on the left hand
    ! side of the equation.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    type(t_vectorScalar), pointer :: p_rsubvector
    real(dp), dimension(:,:), allocatable :: Dfunc
    integer(I32) :: celement
    real(DP) :: da, db, dalpha, dp1
    integer :: ipt, iel
    
    ! Get the bounds and the multiplier from the collection
    dalpha = rcollection%DquickAccess(3)
    dp1 = rcollection%DquickAccess(4)
    
    ! Scale the bounds by -alpha as we analyse lambda
    ! and not u. Change the role of min/max because of the "-"
    ! sign!
    da = -rcollection%DquickAccess(2)*dalpha
    db = -rcollection%DquickAccess(1)*dalpha
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Do we have to analyse lambda_1 or lambda_2?
    if (rcollection%IquickAccess(1) .eq. 1) then
      p_rsubvector => p_Rvector%RvectorBlock(4)
    else
      p_rsubvector => p_Rvector%RvectorBlock(5)
    end if
  
    ! Allocate memory for the function values in the cubature points:
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the 
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector T we evaluate!
    
    celement = rdiscretisationTrial%RelementDistr(&
        max(1,rdomainIntSubset%ielementDistribution))%celement
    
    call fevl_evaluate_sim (p_rsubvector, &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc)
    
    ! Now check the function values lambda.
    ! If b < -1/alpha lambda < a, return 1/alpha.
    ! Otherwise, return 0.
    do iel = 1,ubound(Dcoefficients,3)
      do ipt = 1,ubound(Dcoefficients,2)
        ! Check if the dual variable is in the bounds for the control.
        if ((Dfunc(ipt,iel) .gt. da) .and. (Dfunc(ipt,iel) .lt. db)) then
          Dcoefficients(1,ipt,iel) = dp1
        else
          Dcoefficients(1,ipt,iel) = 0.0_DP
        end if
      end do
    end do
    
    ! Release memory
    deallocate(Dfunc)

  end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_ProjMassCollect (rdiscretisationTrial,rdiscretisationTest,rform, &
      nelements,npointsPerElement,Dpoints, IdofsTrial,IdofsTest,rdomainIntSubset, &
      Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form
    ! that assembles the projective mass matrix.
    !
    ! The coefficients is c=c(lambda_1) or =c(lambda_2) with
    ! c=1/alpha if a < -1/alpha lambda_i < b and c=0 otherwise. This is the derivative
    ! of the projection operator "-P[a,b](-1/alpha lambda_i)" on the left hand
    ! side of the equation.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; trial space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTrial
    
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.; test space.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisationTest

    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    type(t_vectorScalar), pointer :: p_rsubvector
    real(dp), dimension(:,:), allocatable :: Dfunc
    integer(I32) :: celement
    real(DP) :: da, db, dalpha, dp1
    integer :: ipt, iel
    integer :: nptsInactive
    integer, dimension(:), pointer :: p_IelementList
    
    ! Get the bounds and the multiplier from the collection
    dalpha = rcollection%DquickAccess(3)
    dp1 = rcollection%DquickAccess(4)
    
    ! Scale the bounds by -alpha as we analyse lambda
    ! and not u. Change the role of min/max because of the "-"
    ! sign!
    da = -rcollection%DquickAccess(2)*dalpha
    db = -rcollection%DquickAccess(1)*dalpha
    
    ! Get a pointer to the FE solution from the collection.
    ! The routine below wrote a pointer to the vector T to the
    ! first quick-access vector pointer in the collection.
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Do we have to analyse lambda_1 or lambda_2?
    if (rcollection%IquickAccess(1) .eq. 1) then
      p_rsubvector => p_Rvector%RvectorBlock(4)
    else
      p_rsubvector => p_Rvector%RvectorBlock(5)
    end if
  
    ! Allocate memory for the function values in the cubature points:
    allocate(Dfunc(ubound(Dcoefficients,2),ubound(Dcoefficients,3)))
    
    ! Calculate the function value of the solution vector in all
    ! our cubature points:
    !
    ! Figure out the element type, then call the 
    ! evaluation routine for a prepared element set.
    ! This works only if the trial space of the matrix coincides
    ! with the FE space of the vector T we evaluate!
    
    celement = rdiscretisationTrial%RelementDistr(&
        rdomainIntSubset%ielementDistribution)%celement
    
    call fevl_evaluate_sim (p_rsubvector, &
        rdomainIntSubset%p_revalElementSet, &
        celement, rdomainIntSubset%p_IdofsTrial, DER_FUNC, Dfunc)
    
    call storage_getbase_int (rcollection%IquickAccess(2),p_IelementList)
    
    ! Now check the function values lambda.
    ! If b < -1/alpha lambda < a, return 1/alpha.
    ! Otherwise, return 0.
    !
    ! If we detect an element, which is on the border of the active set,
    ! return 0 everywhere and collect the element to the element list,
    ! so the assembly routine can assemble there on a 2nd pass later.
    do iel = 1,ubound(Dcoefficients,3)
      
      ! Count the inactive points.
      nptsInactive = 0
      do ipt = 1,ubound(Dcoefficients,2)
        ! Check if the dual variable is in the bounds for the control.
        if ((Dfunc(ipt,iel) .gt. da) .and. (Dfunc(ipt,iel) .lt. db)) then
          nptsInactive = nptsInactive + 1
        end if
      end do
      
      ! All points inactive? Ok. 
      ! Partially active? Remember the element.
      ! Completely active? Return 0 everywhere.
      if (nptsInactive .eq.  ubound(Dcoefficients,2)) then
        do ipt = 1,ubound(Dcoefficients,2)
          Dcoefficients(1,ipt,iel) = dp1
        end do
      else
        do ipt = 1,ubound(Dcoefficients,2)
          Dcoefficients(1,ipt,iel) = 0.0_DP
        end do
        if (nptsInactive .gt. 0) then
          rcollection%IquickAccess(3) = rcollection%IquickAccess(3) + 1
          p_IelementList(rcollection%IquickAccess(3)) = &
              rdomainIntSubset%p_Ielements(iel)
        end if
      end if
    end do
    
    ! Release memory
    deallocate(Dfunc)

  end subroutine

  ! -----------------------------------------------------

  subroutine massmatfilter (rmatrix, rvector, dalphaC, dmin, dmax)
      
  ! Filters a mass matrix. The lines in the matrix rmatrix corresponding
  ! to all entries in the (control-)vector violating the constraints
  ! of the problem.
  
  ! Matrix to be filtered
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! Vector containing a dial solution lambda. Whereever -1/alpha*lambda
  ! violates the control constraints given by rnonlinearSpatialMatrix, the corresponding
  ! lines are set to 0.
  type(t_vectorScalar), intent(in) :: rvector
  
  ! ALPHA regularisation parameter from the space-time matrix
  real(dp), intent(in) :: dalphaC
  
  ! minimum bound for the control
  real(dp), intent(in) :: dmin

  ! maximum bound for the control
  real(dp), intent(in) :: dmax
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata
    integer, dimension(:), allocatable :: p_Idofs
    integer :: i,nviolate
    real(dp) :: du
    
    ! Get the vector data
    call lsyssc_getbase_double (rvector,p_Ddata)
    
    ! Figure out the DOF's violating the constraints
    allocate(p_Idofs(rvector%NEQ))
    
    nviolate = 0
    do i=1,rvector%NEQ
      du = -p_Ddata(i)/dalphaC
      if ((du .le. dmin) .or. (du .ge. dmax)) then
        nviolate = nviolate + 1
        p_Idofs(nviolate) = i
      end if
    end do
    
    if (nviolate .gt. 0) then
      ! Filter the matrix
      call mmod_replaceLinesByZero (rmatrix,p_Idofs(1:nviolate))
    end if
    
    deallocate(p_Idofs)

  end subroutine    

  ! -----------------------------------------------------

  subroutine massmatfilterVar (rmatrix, rvector, dalphaC, rvectorMin, rvectorMax)
      
  ! Filters a mass matrix. The lines in the matrix rmatrix corresponding
  ! to all entries in the (control-)vector violating the constraints
  ! of the problem.
  ! Non-constant variant.
  
  ! Matrix to be filtered
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! Vector containing a dial solution lambda. Whereever -1/alpha*lambda
  ! violates the control constraints given by rnonlinearSpatialMatrix, the corresponding
  ! lines are set to 0.
  type(t_vectorScalar), intent(in) :: rvector
  
  ! ALPHA regularisation parameter from the space-time matrix
  real(dp), intent(in) :: dalphaC
  
  ! minimum bound for the control
  type(t_vectorScalar), intent(in) :: rvectorMin

  ! maximum bound for the control
  type(t_vectorScalar), intent(in) :: rvectorMax
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata,p_DdataMin,p_DdataMax
    integer, dimension(:), allocatable :: p_Idofs
    integer :: i,nviolate
    real(dp) :: du
    
    ! Get the vector data
    call lsyssc_getbase_double (rvector,p_Ddata)
    call lsyssc_getbase_double (rvectorMin,p_DdataMin)
    call lsyssc_getbase_double (rvectorMax,p_DdataMax)
    
    ! Figure out the DOF's violating the constraints
    allocate(p_Idofs(rvector%NEQ))
    
    nviolate = 0
    do i=1,rvector%NEQ
      du = -p_Ddata(i)/dalphaC
      if ((du .le. p_DdataMin(i)) .or. (du .ge. p_DdataMax(i))) then
        nviolate = nviolate + 1
        p_Idofs(nviolate) = i
      end if
    end do
    
    if (nviolate .gt. 0) then
      ! Filter the matrix
      call mmod_replaceLinesByZero (rmatrix,p_Idofs(1:nviolate))
    end if
    
    deallocate(p_Idofs)

  end subroutine    

  ! -----------------------------------------------------

  subroutine approxProjectionDerivative (rmatrix, rvector, dalphaC, dmin, dmax, dh)
      
  ! Calculates the approximative Newton matrix of the projection operator
  ! P(-1/alpha lambda) by deriving the operator in a discrete sense.
  
  ! Matrix to be set up
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! Vector containing a dial solution lambda. Whereever -1/alpha*lambda
  ! violates the control constraints given by rnonlinearSpatialMatrix, the corresponding
  ! lines are set to 0.
  type(t_vectorScalar), intent(in) :: rvector
  
  ! ALPHA regularisation parameter from the space-time matrix
  real(dp), intent(in) :: dalphaC
  
  ! minimum bound for the control
  real(dp), intent(in) :: dmin

  ! maximum bound for the control
  real(dp), intent(in) :: dmax
  
  ! Step length for the approximative derivative
  real(dp), intent(in) :: dh
  
    ! The projection operator is given by:
    !
    !          a, if u <= a
    !  P(u) =  u, if a <= u <= b
    !          b, if u >= b
    !
    ! The Frechet derivative of this operator can be calculated by
    !
    !   (P(u+h)-P(u-h))/2 = DP(u)h
    ! 
    ! with h being an arbitrary function <> 0.
    ! Rewriting this in a discrete sense yields
    !
    !   ( P(u + h e_i) - P(u - h e_i) ) / (2h)  =  DP(u) h e_i
    !
    ! with u being the vector of the corresponding FE function.
    ! Let us denote the matrix B:=DP(u), then we obtain by this formula:
    !
    !   B_ij  =  [ ( P(u + h e_j) - P(u - h e_j) ) / (2h) ]_i
    !
    ! If we now treat the P()-operator coponent-wise instead of vector
    ! wise, this means:
    !
    !   B_ij  =  ( P(u_j+h) - P(u_j-h) ) / (2h)   , i=j
    !            0                                , otherwise
      
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata, p_Dmatrix
    integer, dimension(:), pointer :: p_Kdiagonal
    integer :: i
    real(DP) :: du1,du2
    
    ! Get the vector data
    call lsyssc_getbase_double (rvector,p_Ddata)
    call lsyssc_getbase_double (rmatrix,p_Dmatrix)
    call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
    
    call lsyssc_clearMatrix (rmatrix)
    
    ! Loop through the diagonal entries
    do i=1,rmatrix%NEQ

      ! Calculate the diagonal entry, that's it
      
      du1 = -min(dmax,max(dmin,-(p_Ddata(i)/dalphaC) + dh ))
      du2 = -min(dmax,max(dmin,-(p_Ddata(i)/dalphaC) - dh ))
      
      p_Dmatrix(p_Kdiagonal(i)) = (du1-du2)/(2.0_DP*dh)
    
    end do

  end subroutine    

  ! ***************************************************************************

!<subroutine>

  subroutine smva_assembleMatrix (coperation,cmatrixType,rflags,&
      rnonlinearSpatialMatrix,rmatrix,rfineMatrix)

!<description>
  ! This routine assembles a matrix space. The caller must initialise 
  ! the rnonlinearSpatialMatrix according to how the matrix should look like.
  ! The 'coperation' parameter tells the routine what to do.
  ! The destination matrix rmatrix is then set up or updated.
  !
  ! The parameters rvector and rfineMatrix are optional. rvector must be
  ! specified, if the nonlinearity is activated (parameter $\gamma\not=0$ in 
  ! rnonlinearSpatialMatrix). This vector specifies the 'solution' where the nonlinearity 
  ! $u\nabla u$ is evaluated.
  ! rfineMatrix allows to specify a matrix of a 'one level refined mesh'. This
  ! is usually used when setting up preconditioners over multiple levels.
  ! Specifying such a matrix allows the routine (depending on the discretisation)
  ! to include some special stabilisation terms into the matrix rmatrix.
  !
  ! The routine will not include any boundary conditions in the matrix.
!</description>

!<input>

  ! One of the CCMASM_xxxx-constants. This parameter specifies 'what to do'.
  ! Using the CCMASM_ALLOCMEM constant, the routine will allocate memory in
  ! rmatrix for all the matrix blocks but will not compute any entries.
  ! Using the CCMASM_COMPUTE constant, the routine will compute the matrix
  ! entries, while assuming that the memory was allocated previously.
  ! Using the CCMASM_ALLOCANDCOMPUTE constant, the routine will do both,
  ! allocate memory and compute the entries.
  !
  ! If any of the CCMASM_ALLOCxxxx constants is specified here and rmatrix is
  ! already initialised, rmatrix is released and completely rebuild in memory!
  ! Therefore, CCMASM_COMPUTE should be used to update an existing matrix.
  !
  ! The constant CMASM_QUICKREFERENCES may be specified additional to one of
  ! the other constants (e.g. as 'CCMASM_ALLOCANDCOMPUTE+CMASM_QUICKREFERENCES').
  ! If specified, the routine tries to avoid memory allocation. This means e.g.
  ! that references to the original gradient (B-)matrices from rnonlinearSpatialMatrix
  ! are written to rmatrix; matrix entries are not copied!
  ! (This can be used e.g. for setting up a matrix for building a defect
  !  vector without copying matrix data.)
  ! In this case, the caller MUST NOT CHANGE rmatrix in any way, otherwise
  ! the original (template) matrices would be changed!
  integer(I32), intent(IN) :: coperation

  ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
  ! constants.
  ! Usually, CCMASM_MTP_AUTOMATIC is used here. This will automatically determine
  ! the 'smallest possible' matrix structure fitting the needs of the input
  ! parameters. By specifying another matrix type, the caller can explicitly
  ! take influence on the general matrix structure.
  !
  ! If the matrix already exists and only its entries are to be computed,
  ! CCMASM_MTP_AUTOMATIC should be specified here.
  integer, intent(IN) :: cmatrixType
  
  ! Additional matrix assembly flags that configure how to generate the matrix.
  type(t_matrixAssemblyFlags), intent(in) :: rflags

  ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix. 
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, p_rmatrixTemplateXXXX
  ! must be initialised as well as p_rdiscretisation!
  ! The new matrix is created based p_rmatrixTemplateXXXX as well as p_rdiscretisation.
  ! Memory is automatically allocated if it's missing.
  type(t_nonlinearSpatialMatrix), intent(in) :: rnonlinearSpatialMatrix

  ! OPTIONAL: This parameter allows to specify a 'fine grid matrix'. This is 
  ! usually done when assembling matrices on multiple levels. If specified, the
  ! routine will (if possible) try to include a level-dependent stabilisation
  ! term into the matrix (-> e.g. constant matrix restriction for nonparametric
  ! Rannacher-Turek element if cells are too anisotropic).
  type(t_matrixBlock), intent(IN), optional :: rfineMatrix
  
!</input>

!<inputoutput>

  ! The destination matrix which should be set up.
  ! If not initialised, a new matrix is created (as if CCMASM_ALLOCxxxx 
  ! was specified).
  ! If initialised, the existing matrix is updated or recreated, depending on
  ! coperation.
  type(t_matrixBlock), intent(inout) :: rmatrix
  
!</inputoutput>
  
!</subroutine>

    ! local variables
    logical :: ballocate
    type(t_matrixBlock) :: rtempMatrix
    type(t_vectorBlock) :: rtempVector
    type(t_settings_stabil) :: rstabilisation
    integer :: i,j
    type(t_optcoperator) :: roptcoperator
    type(t_vectorBlock), pointer :: p_rprimalSol, p_rdualSol
    type(t_blockDiscretisation) :: rvelDiscr
    real(dp), dimension(:), pointer :: p_Ddata
    real(DP) :: dweightConvection
    
    logical, parameter :: bnewmethod = .false.
    
    ! Debug weight for the convection
    dweightConvection = rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightConvection

    ballocate = .false.
    if ((rmatrix%NEQ .le. 0) .or. &
        iand(coperation,CCMASM_ALLOCMEM) .ne. 0) then
      ballocate = .true.
    end if
    
    ! What should we do? Allocate memory?
    if (ballocate) then
    
      ! Release the matrix if present.
      call lsysbl_releaseMatrix (rmatrix)
    
      ! Allocate the system matrix
      call lsysbl_createMatBlockByDiscr(&
          rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual,rmatrix)
      do j=1,2
        do i=1,2
          call lsysbl_createEmptyMatrix (rtempMatrix,NDIM2D+1,NDIM2D+1)
          call allocSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rtempMatrix,&
              rnonlinearSpatialMatrix%Diota(i,j),rnonlinearSpatialMatrix%Dalpha(i,j),&
              rnonlinearSpatialMatrix%Dtheta(i,j),rnonlinearSpatialMatrix%Dgamma(i,j),&
              rnonlinearSpatialMatrix%Dnewton(i,j),rnonlinearSpatialMatrix%DgammaT(i,j),&
              rnonlinearSpatialMatrix%DnewtonT(i,j),rnonlinearSpatialMatrix%Deta(i,j),&
              rnonlinearSpatialMatrix%Dtau(i,j),rnonlinearSpatialMatrix%Dkappa(i,j),&
              (i .eq. 2) .and. (j .eq. 2))
          call lsysbl_updateMatStrucInfo(rtempMatrix)
          call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,(i-1)*3+1,(j-1)*3+1)
          call lsysbl_releaseMatrix (rtempMatrix)
        end do
      end do
      
      call lsysbl_allocEmptyMatrix (rmatrix,LSYSSC_SETM_UNDEFINED,.true.)

      ! Share the B-matrices of the dual equation with those of the primal
      ! equation; we need them only once since they are the same, even if
      ! boundary conditions are imposed.
      call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,3),rmatrix%RmatrixBlock(4,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(2,3),rmatrix%RmatrixBlock(5,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      rmatrix%RmatrixBlock(4,6)%dscaleFactor = rnonlinearSpatialMatrix%Deta(2,2)
      rmatrix%RmatrixBlock(5,6)%dscaleFactor = rnonlinearSpatialMatrix%Deta(2,2)

    end if   
   
    if (iand(coperation,CCMASM_COMPUTE) .ne. 0) then

      if (.not. bnewmethod) then
        ! The system matrix in the whole looks like:
        !          
        !    ( A11  A12  B1  R1           ) 
        !    ( A21  A22  B2       R1      ) 
        !    ( B1^T B2^T I                ) 
        !    ( R2            A44  A45  B1 ) 
        !    (      R2       A54  A55  B2 ) 
        !    (               B1^T B2^T I  ) 
        !
        ! With some multiplication factors in front of the matrices.
        !
        ! At first switch off all submatrices. Those submatrices that we
        ! assign some data are explicitely switched on again.
        rmatrix%RmatrixBlock(:,:)%dscaleFactor = 0.0_DP
        
        ! Derive a discretisation structure for either primal or
        ! dual variables. Create a temporary matrix based on that.
        call spdiscr_deriveBlockDiscr (&
            rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual, rvelDiscr, 1, 3)
        call lsysbl_createMatBlockByDiscr (rvelDiscr,rtempMatrix)
        
        ! Primal equation
        ! ---------------
        ! In the first step, assemble the linear parts
        ! (Laplace, Mass, B/B^T) on the main diagonal of the primal equation.
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,1,3)
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,&
            rtempMatrix,rnonlinearSpatialMatrix%Diota(1,1),rnonlinearSpatialMatrix%Dalpha(1,1),&
            rnonlinearSpatialMatrix%Dtheta(1,1),&
            rnonlinearSpatialMatrix%Deta(1,1),rnonlinearSpatialMatrix%Dtau(1,1),.false.)


        ! Assemble the nonlinearity u*grad(.) or the Newton nonlinearity
        ! u*grad(.)+grad(u)*(.) to the velocity.
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,&
              rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
              rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
              rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal,&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ1)      
        case (2)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,&
              rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
              rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
              rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal,&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ1)      
        case (3)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,&
              rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
              rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
              rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal,&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ1)      
        end select

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,1,1)

        ! Include the mass matrix blocks
        !
        !    (  .    .    .    R1   .    . ) 
        !    (  .    .    .    .    R1   . ) 
        !    (  .    .    .    .    .    . )
        
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,1,3,4,6)
        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, &
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1, &
              rnonlinearSpatialMatrix%Dalpha(1,2))
        case (2)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, &
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2, &
              rnonlinearSpatialMatrix%Dalpha(1,2))
        case (3)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, &
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,&
              rnonlinearSpatialMatrix%Dalpha(1,2))
        end select

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,1,4)

        ! Dual equation
        ! -------------
        ! In the first step, assemble the linear parts
        ! (Laplace, Mass, B/B^T) on the main diagonal of the primal equation.
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,4,6)
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rtempMatrix,&
            rnonlinearSpatialMatrix%Diota(2,2),rnonlinearSpatialMatrix%Dalpha(2,2),&
            rnonlinearSpatialMatrix%Dtheta(2,2),&
            rnonlinearSpatialMatrix%Deta(2,2),rnonlinearSpatialMatrix%Dtau(2,2),&
            rnonlinearSpatialMatrix%Deta(1,1) .eq. rnonlinearSpatialMatrix%Deta(2,2))

        ! Assemble the nonlinearity u*grad(.) or the Newton nonlinearity
        ! u*grad(.)+grad(u)*(.) to the velocity.
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,&
              rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
              rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
              rnonlinearSpatialMatrix%rdiscrData%rstabilDual,&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ2)      
        case (2)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,&
              rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
              rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
              rnonlinearSpatialMatrix%rdiscrData%rstabilDual,&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ2)      
        case (3)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,&
              rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
              rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
              rnonlinearSpatialMatrix%rdiscrData%rstabilDual,&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ2)      
        end select

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,4,4)

        ! Include the mass matrix blocks
        !
        !    (  R2   .    .    .    .    .  ) 
        !    (  .    R2   .    .    .    .  ) 
        !    (  .    .    .    .    .    .  )
        !
        ! Remember that they are nonlinear if we have a preconditioner with
        ! Newton activated!
        
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,4,6,1,3)
        
        ! Assemble linear parts.
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rtempMatrix,&
            rnonlinearSpatialMatrix%Diota(2,1),rnonlinearSpatialMatrix%Dalpha(2,1),&
            rnonlinearSpatialMatrix%Dtheta(2,1),&
            rnonlinearSpatialMatrix%Deta(2,1),rnonlinearSpatialMatrix%Dtau(2,1),.false.)

        ! Co stabilisation in the convective parts here.            
        ! rstabilisation = t_convecStabilisation(&
        !    rnonlinearSpatialMatrix%iupwind2,rnonlinearSpatialMatrix%dupsam2)
        rstabilisation = t_settings_stabil(0,0.0_DP,1)

        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rtempVector, 4,6,.true.)
        case (2)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rtempVector, 4,6,.true.)
        case (3)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rtempVector, 4,6,.true.)
        end select
            
        call assembleConvection (&
            rnonlinearSpatialMatrix,rtempMatrix,rtempVector,&
            rnonlinearSpatialMatrix%Dgamma(2,1),rnonlinearSpatialMatrix%DgammaT(2,1),&
            rnonlinearSpatialMatrix%Dnewton(2,1),rnonlinearSpatialMatrix%DnewtonT(2,1),&
            rstabilisation)
            
        if (rtempVector%NEQ .ne. 0) &
          call lsysbl_releaseVector (rtempVector)
            
        ! There is probably a 2nd reactive term stemming from the next time step.
        ! Assemble it.
        
        select case (rnonlinearSpatialMatrix%idualSol2)
        case (1)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rtempVector, 4,6,.true.)
        case (2)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rtempVector, 4,6,.true.)
        case (3)
          call lsysbl_deriveSubvector(&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rtempVector, 4,6,.true.)
        end select
            
        call assembleConvection (&
            rnonlinearSpatialMatrix,rtempMatrix,rtempVector,&
            0.0_DP,rnonlinearSpatialMatrix%DgammaT2(2,1),&
            rnonlinearSpatialMatrix%Dnewton2(2,1),0.0_DP,&
            rstabilisation)      

        if (rtempVector%NEQ .ne. 0) &
          call lsysbl_releaseVector (rtempVector)

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,4,1)

        ! Switch the I-matrix in the continuity equation on/off.
        ! The matrix always exists -- but is usually filled with zeroes
        ! or switched off.
        rmatrix%RmatrixBlock(3,3)%dscaleFactor = rnonlinearSpatialMatrix%Dkappa(1,1)
        if (rnonlinearSpatialMatrix%Dkappa(1,1) .ne. 0.0_DP) then
          call lsyssc_initialiseIdentityMatrix (rmatrix%RmatrixBlock(3,3))
        end if
        
        rmatrix%RmatrixBlock(6,6)%dscaleFactor = rnonlinearSpatialMatrix%Dkappa(2,2)
        if (rnonlinearSpatialMatrix%Dkappa(2,2) .ne. 0.0_DP) then
          call lsyssc_initialiseIdentityMatrix (rmatrix%RmatrixBlock(6,6))
        end if

        ! Matrix restriction
        ! ---------------------------------------------------

        ! For the construction of matrices on lower levels, call the matrix
        ! restriction. In case we have a uniform discretisation with Q1~,
        ! iadaptivematrix may be <> 0 and so this will rebuild some matrix entries
        ! by a Galerkin approach using constant prolongation/restriction.
        ! This helps to stabilise the solver if there are elements in the
        ! mesh with high aspect ratio.
        if (present(rfineMatrix)) then
        
          ! Primal system:
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(1,1), &
              rmatrix%RmatrixBlock(1,1), &
              rflags%iadaptiveMatrices, &
              rflags%dadmatthreshold)
              
          if (.not. lsyssc_isMatrixContentShared(&
              rfineMatrix%RmatrixBlock(1,1),rfineMatrix%RmatrixBlock(2,2))) then
            call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(2,2), &
                rmatrix%RmatrixBlock(1,1), &
                rflags%iadaptiveMatrices, &
                rflags%dadmatthreshold)
          end if
            
          ! Dual system
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(4,4), &
              rmatrix%RmatrixBlock(4,4), &
              rflags%iadaptiveMatrices, &
              rflags%dadmatthreshold)
              
          if (.not. lsyssc_isMatrixContentShared(&
            rmatrix%RmatrixBlock(4,4),rmatrix%RmatrixBlock(5,5))) then
            call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(5,5), &
                rmatrix%RmatrixBlock(5,5), &
                rflags%iadaptiveMatrices, &
                rflags%dadmatthreshold)
          end if
          
        end if

        ! Release memory
        call lsysbl_releaseMatrix(rtempMatrix)
        call spdiscr_releaseBlockDiscr(rvelDiscr)

        !call matio_writeBlockMatrixHR (rmatrix, 'matrix',&
        !    .true., 0, 'matrix1.txt', '(E12.5)', 1E-10_DP)
        
      else

        ! New implementation: Use the 'direct' assembly routine.
        
        rmatrix%RmatrixBlock(1:2,1:2)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(4:5,4:5)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(1:2,4:5)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(4:5,1:2)%dscaleFactor = 1.0_DP

        ! Clear the old matrix
        if (iand(coperation,CMASM_NOCLEAR) .eq. 0) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,4))        
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,5))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,4))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,5))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,4))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,4))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,5))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,5))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,2))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,1))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,2))
          
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,3))
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(6,6))
        end if
        
        ! The B1/B2 matrices exist up to now only in rnonlinearSpatialMatrix.
        ! Put a copy of them into the block matrix.
        !
        ! Note that we share the structure of B1/B2 with those B1/B2 of the
        ! block matrix, while we create copies of the entries. The B-blocks
        ! are already prepared and memory for the entries is already allocated;
        ! so we only have to copy the entries.
        !
        ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
        ! memory if necessary.
        if (rnonlinearSpatialMatrix%Deta(1,1) .ne. 0.0_DP) then
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
              rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
              rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        end if
        
        if (rnonlinearSpatialMatrix%Dtau(1,1) .ne. 0.0_DP) then                              
          ! Furthermore, put B1^T and B2^T to the block matrix.
          ! These matrices are always 'shared'.
          if (rflags%bvirtualTransposedD) then
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
                rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
                rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
          else
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
                rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
                rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          end if
        end if
                                        
        if (rnonlinearSpatialMatrix%Dtau(2,2) .ne. 0.0_DP) then
          if (rflags%bvirtualTransposedD) then
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
                rmatrix%RmatrixBlock(6,4),LSYSSC_TR_VIRTUAL)
            call lsyssc_transposeMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
                rmatrix%RmatrixBlock(6,5),LSYSSC_TR_VIRTUAL)
          else
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
                rmatrix%RmatrixBlock(6,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
                rmatrix%RmatrixBlock(6,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          end if
        end if

        rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearSpatialMatrix%Deta(1,1)
        rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearSpatialMatrix%Deta(1,1)

        rmatrix%RmatrixBlock(4,6)%dscaleFactor = rnonlinearSpatialMatrix%Deta(2,2)
        rmatrix%RmatrixBlock(5,6)%dscaleFactor = rnonlinearSpatialMatrix%Deta(2,2)
        
        rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearSpatialMatrix%Dtau(1,1)
        rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearSpatialMatrix%Dtau(1,1)

        rmatrix%RmatrixBlock(6,4)%dscaleFactor = rnonlinearSpatialMatrix%Dtau(2,2)
        rmatrix%RmatrixBlock(6,5)%dscaleFactor = rnonlinearSpatialMatrix%Dtau(2,2)

        ! ---------------------------------------------------
        ! Now a slightly more advanced task for which we use a separate
        ! routine and some submatrices/vectors: The nonlinearity.

        ! Initialise the operator structure for what we need.
        roptcoperator%dupsamPrimal = rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal%dupsam
        roptcoperator%dupsamDual = rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam
        
        ! Timestep-weights
        roptcoperator%dprimalAlpha = rnonlinearSpatialMatrix%Dalpha(1,1)
        roptcoperator%ddualAlpha   = rnonlinearSpatialMatrix%Dalpha(2,2)

        ! Stokes operator
        roptcoperator%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
        roptcoperator%dprimalBeta = rnonlinearSpatialMatrix%Dtheta(1,1)
        roptcoperator%ddualBeta   = rnonlinearSpatialMatrix%Dtheta(2,2)
        
        ! Nonlinearity
        if (rnonlinearSpatialMatrix%Dgamma(1,1) .ne. 0.0_DP) then
          roptcoperator%dprimalDelta = dweightconvection * rnonlinearSpatialMatrix%Dgamma(1,1)
          roptcoperator%ddualDelta   = dweightconvection * rnonlinearSpatialMatrix%Dgamma(2,2)
          roptcoperator%ddualNewtonTrans = dweightconvection * rnonlinearSpatialMatrix%DnewtonT(2,2)
          
          ! Whether or not Newton is active has no influence to the
          ! defect, so the following lines are commented out.
          ! if (rparams%bnewton) then
          roptcoperator%dprimalNewton    = dweightconvection * rnonlinearSpatialMatrix%Dnewton(1,1)
          roptcoperator%ddualRDeltaTrans = dweightconvection * rnonlinearSpatialMatrix%DgammaT(2,1)
          roptcoperator%ddualRNewton     = dweightconvection * rnonlinearSpatialMatrix%Dnewton(2,1)
          ! end if
          
        end if
        
        ! Coupling matrices
        !if (rparams%bdualcoupledtoprimal) then
          roptcoperator%ddualRAlpha = rnonlinearSpatialMatrix%Dalpha(2,1)
        !end if

        !if (rparams%bcontrolactive) then
          roptcoperator%dcontrolWeight = &
              -rnonlinearSpatialMatrix%Dalpha(1,2)*rnonlinearSpatialMatrix%dalphaC
          roptcoperator%dcontrolMultiplier = -1.0_DP/rnonlinearSpatialMatrix%dalphaC
        !end if
        
        if (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints .ne. 0) then

          roptcoperator%ccontrolProjection = &
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints

          roptcoperator%cconstraintsType = &
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType
              
          roptcoperator%dmin1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
          roptcoperator%dmax1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
          roptcoperator%dmin2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
          roptcoperator%dmax2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

          roptcoperator%p_rumin1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin1
          roptcoperator%p_rumax1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax1
          roptcoperator%p_rumin2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin2
          roptcoperator%p_rumax2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax2

          roptcoperator%p_rvectorumin => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin
          roptcoperator%p_rvectorumax => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax
        end if
        
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
        case (2)
          p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
        case (3)
          p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
        end select
        
        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
        case (2)
          p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
        case (3)
          p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
        end select
        
        ! Calculate the velocity-dependent part of the system matrix.
        call conv_strdiffOptC2dgetMatrix (rmatrix,roptcoperator,1.0_DP,&
            p_rprimalSol,p_rdualSol)
            
        !call matio_writeBlockMatrixHR (rmatrix, 'matrix',&
        !    .true., 0, 'matrix2.txt', '(E12.5)', 1E-10_DP)
      
      end if
      
    end if
    
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,1),p_Ddata)
    
  contains
  
    ! -----------------------------------------------------
    
    subroutine allocSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,rsubmatrix,&
        diota,dalpha,dtheta,dgamma,dnewton,dgammaT,dnewtonT,deta,dtau,dkappa,&
        bignoreEta)
        
    ! Allocates memory for a submatrix of the system matrix.
    ! The submatrix can then be inserted into a larger system matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_nonlinearSpatialMatrix), intent(IN), target :: rnonlinearSpatialMatrix

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    integer, intent(IN) :: cmatrixType

    ! Additional matrix assembly flags that configure how to generate the matrix.
    type(t_matrixAssemblyFlags), intent(in) :: rflags

    ! A block matrix that receives the basic system matrix.
    type(t_matrixBlock), intent(INOUT) :: rsubmatrix
    
    ! Coefficients in front of all operators in the system matrix.
    ! Depending on which operators are 'active', memory is allocated.
    real(DP), intent(in) :: diota,dalpha,dtheta,dgamma,dnewton,&
        dgammaT,dnewtonT,deta,dtau
        
    ! Coefficient in front of a pressure matrix at position (3,3).
    ! Switches that pressure-matrix on/off.
    real(DP), intent(in) :: dkappa
       
    ! If set to TRUE, the assembly of the data in the B-matrices will be skipped.
    ! This can be used if the B-matrices of one equation share their information
    ! with B-matrices from another equation and thus do not need to be assembled
    ! twice.
    logical, intent(in) :: bignoreEta

      ! local variables
      logical :: bdecoupled,bfulltensor
       
      ! Determine the shape of the matrix
      bdecoupled = cmatrixType .eq. CCMASM_MTP_DECOUPLED
      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
      
      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
        bfulltensor = dnewton .ne. 0.0_DP
      end if
    
      ! Let's consider the global system in detail. The standard matrix It has 
      ! roughly the following shape:
      !
      !    ( A11  .    B1  ) 
      !    ( .    A22  B2  ) 
      !    ( B1^T B2^T C   ) 
      !
      ! All matices may have multiplication factors in their front.
      !
      ! The structure of the matrices A11 and A22 of the global system matrix 
      ! is governed by the template FEM matrix.
      ! Initialise them with the same structure, i.e. A11, A22 share (!) their
      ! structure (not the entries) with that of the template matrix.
      !
      ! For this purpose, use the "duplicate matrix" routine.
      ! The structure of the matrix is shared with the template FEM matrix.
      ! For the content, a new empty array is allocated which will later receive
      ! the entries.
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rsubmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
          
      if (.not. bfulltensor) then     
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix 
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content.
        call lsyssc_duplicateMatrix (rsubmatrix%RmatrixBlock(1,1),&
                    rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
      else
      
        ! Otherwise, create another copy of the template matrix.
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
            rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                    
      end if
      
      ! A 'full tensor matrix' consists also of blocks A12 and A21.
      if (bfulltensor) then

        ! We have a matrix in the following shape:
        !
        !    ( A11  A12  B1  ) 
        !    ( A21  A22  B2  ) 
        !    ( B1^T B2^T C   ) 
        !
        ! Create A12 and A21.
      
        if (rsubmatrix%RmatrixBlock(1,2)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM, &
              rsubmatrix%RmatrixBlock(1,2), LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! call lsyssc_allocEmptyMatrix (&
          !     rsubmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rsubmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          ! Create a new matrix A21 in memory. create a new matrix
          ! using the template FEM matrix...
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM, &
              rsubmatrix%RmatrixBlock(2,1), LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! call lsyssc_allocEmptyMatrix (&
          !     rsubmatrix%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
            
        end if
        
      end if

      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create empty space for the entries. 
      ! Later, the B-matrices are copied into here and modified for boundary
      ! conditions.
      if ((deta .ne. 0.0_DP) .and. .not. bignoreEta) then
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rsubmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rsubmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      end if
        
      ! In the same manner, insert an identiy matrix for the pressure
      ! to the system matrix. The matrix is zero by default but may
      ! partially or fully be initialised to an identity matrix depending
      ! on the situation. (It's mostly used for direct solvers/UMFPACK)
      if (dkappa .ne. 0.0_DP) then
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEMPressure, &
            rsubmatrix%RmatrixBlock(3,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      end if

      ! Now, prepare B1^T and B2^T. Are these matrices given?
      !
      ! If yes, take the given matrices. If not, create them by
      ! 'virtually transposing' B1 and B2 (i.e. they share the same
      ! data as B1 and B2 but hate the 'transpose'-flag set).

      if (dtau .ne. 0.0_DP) then    
        if (rflags%bvirtualTransposedD) then
          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
              rsubmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)

          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
              rsubmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
        else  
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
              rsubmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
              rsubmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end if
      end if

    end subroutine

    ! -----------------------------------------------------
    
    subroutine assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rflags,&
        rmatrix,diota,dalpha,dtheta,deta,dtau,bignoreEta)
        
    ! Assembles a matrix containing only linear terms. The matrix is a
    ! submatrix of a larger system matrix and can then be inserted 
    ! into a larger system matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
    ! constants.
    integer, intent(IN) :: cmatrixType

    ! Additional matrix assembly flags that configure how to generate the matrix.
    type(t_matrixAssemblyFlags), intent(in) :: rflags

    ! A block matrix that receives the matrix.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Coefficient in front of an identity matrix.
    real(DP), intent(in) :: diota
    
    ! Coefficient in front of a mass matrix.
    real(DP), intent(in) :: dalpha
    
    ! Coefficient in front of the Laplace matrix.
    real(DP), intent(in) :: dtheta
    
    ! Coefficient in front of the B-matrices
    real(DP), intent(in) :: deta
    
    ! Coefficient in front of the B^T-matrices
    real(DP), intent(in) :: dtau
        
    ! If set to TRUE, the assembly of the data in the B-matrices will be skipped.
    ! This can be used if the B-matrices of one equation share their information
    ! with B-matrices from another equation and thus do not need to be assembled
    ! twice.
    logical, intent(in) :: bignoreEta
        
    ! local variables
    logical :: bshared
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2))
                    
      if ((diota .eq. 0.0_DP) .and. (dalpha .eq. 0.0_DP) .and.&
          (dtheta .eq. 0.0_DP)) then
          
        ! Switch off that matrices
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 0.0_DP
          
      else
                    
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
      
        ! Let's consider the global system in detail:
        !
        !    ( A11  A12  B1  ) = ( A11  A12  A13 )
        !    ( A21  A22  B2  )   ( A21  A22  A23 )
        !    ( B1^T B2^T 0   )   ( A31  A32  A33 )
        !
        ! Here, we set up the linear parts in Aij and Bi.
        
        ! ---------------------------------------------------
        ! If diota <> 0, initialise the matrix with the identity,
        ! otherwise with zero.
        if (diota .eq. 0.0_DP) then
          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,1))
          
          if (.not. bshared) then
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,2))
          end if
          
        else
          
          call lsyssc_initialiseIdentityMatrix (&
              rmatrix%RmatrixBlock(1,1))
          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),diota)
          
          if (.not. bshared) then
            call lsyssc_initialiseIdentityMatrix (&
                rmatrix%RmatrixBlock(2,2))
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),diota)
          end if
          
        end if
      
        ! ---------------------------------------------------
        ! Plug in the mass matrix?
        if (dalpha .ne. 0.0_DP) then
         
          call lsyssc_matrixLinearComb (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              dalpha,rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
              
          if (.not. bshared) then

            call lsyssc_matrixLinearComb (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                dalpha,rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if
          
        end if
        
        ! ---------------------------------------------------
        ! Plug in the Stokes matrix?
        if (dtheta .ne. 0.0_DP) then
          call lsyssc_matrixLinearComb (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace,&
              rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu * dtheta,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
              
          if (.not. bshared) then
            call lsyssc_matrixLinearComb (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace,&
                rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu * dtheta,&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if
        end if
        
      end if
      
      ! We exclude the velocity submatrices here, so our system looks like:
      !
      !    (           B1 ) = (           A13 )
      !    (           B2 )   (           A23 )
      !    ( B1^T B2^T    )   ( A31  A32      )

      ! The B1/B2 matrices exist up to now only in rnonlinearSpatialMatrix.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The B-blocks
      ! are already prepared and memory for the entries is already allocated;
      ! so we only have to copy the entries.
      !
      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
      ! memory if necessary.
      if ((deta .ne. 0.0_DP) .and. .not. bignoreEta) then
        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)

        call lsyssc_duplicateMatrix (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
      end if
      
      if (dtau .ne. 0.0_DP) then                 
        ! Furthermore, put B1^T and B2^T to the block matrix.
        if (rflags%bvirtualTransposedD) then
          ! Virtually transposed matrices are always 'shared'.
          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
              rmatrix%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
          call lsyssc_transposeMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
              rmatrix%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
        else
          ! Real matrices are not shared. This allows them to be overwritten
          ! if necessary (e.g. if an UMFPACK solver is used or similar).
          ! In a situation, where the matrix needs to be overwritten,
          ! the caller must therefore always request a not virtually transposed
          ! matix!
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
              rmatrix%RmatrixBlock(3,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPYOVERWRITE)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
              rmatrix%RmatrixBlock(3,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPYOVERWRITE)
        end if
      end if

      rmatrix%RmatrixBlock(1,3)%dscaleFactor = deta
      rmatrix%RmatrixBlock(2,3)%dscaleFactor = deta
      
      rmatrix%RmatrixBlock(3,1)%dscaleFactor = dtau
      rmatrix%RmatrixBlock(3,2)%dscaleFactor = dtau

    end subroutine
    
    ! -----------------------------------------------------

    subroutine assembleConvection (&
        rnonlinearSpatialMatrix,rmatrix,rvector,dgamma,dgammaT,dnewton,dnewtonT,&
        rstabilisation,rmatrixEOJ)
        
    ! Assembles the convection matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dgamma*N(rvector) + dgammaT*N^t(rvector) +
    !            dnewton*N*(rvector) + dnewtonT*N*^t(rvector)
    !
    ! Even if no nonlinearity is present, the routine can be used to
    ! add stabilisation into the matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    type(t_vectorBlock) :: rvector
    
    ! Weight for the nonlinear term u\grad(.)
    real(DP), intent(in) :: dgamma

    ! Weight for the nonlinear term u(\grad(.))^t
    real(DP), intent(in) :: dgammaT

    ! Weight for the nonlinear term (\grad(.))u
    real(DP), intent(in) :: dnewton

    ! Weight for the nonlinear term (\grad(.))^t u
    real(DP), intent(in) :: dnewtonT
    
    ! Stabilisation parameters.
    type(t_settings_stabil), intent(in) :: rstabilisation
    
    ! OPTIONAL: EOJ matrix in case the precomputed EOJ stabilisation is active.
    type(t_matrixScalar), intent(in), optional :: rmatrixEOJ
    
    ! local variables
    logical :: bshared
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_convStreamDiff2) :: rstreamlineDiffusion3
    type(t_jumpStabilisation) :: rjumpStabil
    type(t_convUpwind) :: rupwindStabil
    real(dp), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
    real(DP) :: dweightConvection
    
      ! Debug weight for the convection
      dweightConvection = rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightConvection
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      call lsysbl_getbase_double (rvector,p_Ddata1)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Ddata2)
             
      if ((dgamma .ne. 0.0_DP) .or. (dgammaT .ne. 0.0_DP) .or. &
          (dnewton .ne. 0.0_DP) .or. (dnewtonT .ne. 0.0_DP)) then
                      
        ! Switch on the offdiagonal matrices if necessary
        if ((dgammaT .ne. 0.0_DP) .or. &
            (dnewton .ne. 0.0_DP) .or. (dnewtonT .ne. 0.0_DP)) then
          ! Switch on missing matrices. Initialise with zero when they
          ! are created the first time.
          if (rmatrix%RmatrixBlock(1,2)%dscaleFactor .eq. 0.0_DP) then
            rmatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
            rmatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          end if
        else
          ! If the submatrices A12 and A21 exist, fill them with zero.
          ! If they don't exist, we don't have to do anything.
          if (lsysbl_isSubmatrixPresent (rmatrix,1,2)) then
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,2))
            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,1))
          end if
        end if
                      
        select case (rstabilisation%cupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion%ddelta            = dweightConvection * dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dweightConvection * dgammaT
          rstreamlineDiffusion%dnewton           = dweightConvection * dnewton
          rstreamlineDiffusion%dnewtonTransposed = dweightConvection * dnewtonT
          
          if (rstabilisation%dupsam .ne. 0) then
            call output_line ("assembleConvection: Warning. Please use the")
            call output_line ("alternative SD method for setting up the stabilised operator!!!")
          end if
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dweightConvection * dgamma
          rstreamlineDiffusion2%ddeltaT = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)
          
        case (CCMASM_STAB_UPWIND)
        
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = 0.0_DP
          rstreamlineDiffusion2%ddeltaT = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity for everything except
          ! for the convection. As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)

          ! Prepare the upwind structure for the assembly of the convection.
          ! Note: Stabilisation weight is actually wrong, but it is not possible
          ! to specify in rupwindStabil%dupsam whether the stabilisation
          ! is added or subtracted!
          rupwindStabil%dupsam = abs(rstabilisation%dupsam)
          rupwindStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          rupwindStabil%dtheta = dweightConvection*dgamma
          
          ! Apply the upwind operator.
          call conv_upwind2d (rvector, rvector, 1.0_DP, 0.0_DP,&
              rupwindStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(1,1))

          if (.not. bshared) then
            call conv_upwind2d (rvector, rvector, 1.0_DP, 0.0_DP,&
                rupwindStabil, CONV_MODMATRIX, rmatrix%RmatrixBlock(2,2))
          end if

          ! Prepare the upwind structure for the assembly of the convection.
          !rstreamlineDiffusion3%dupsam = rstabilisation%dupsam
          !rstreamlineDiffusion3%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          !rstreamlineDiffusion3%dtheta = dgamma
          !rstreamlineDiffusion3%ddelta = 1.0_DP
          
          ! Apply the upwind operator.
          !call conv_streamDiff2Blk2dMat (rstreamlineDiffusion3,rmatrix,rvector)
          
          !call output_line ('Upwind not supported.', &
          !                  OU_CLASS_ERROR,OU_MODE_STD,'assembleConvection')
          !call sys_halt()

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta            = dweightConvection * dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dweightConvection * dgammaT
          rstreamlineDiffusion%dnewton           = dweightConvection * dnewton
          rstreamlineDiffusion%dnewtonTransposed = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODMATRIX, &
                              rmatrix)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight. Compensate for any "-" sign in dgamma!
          rjumpStabil%dtheta = dgamma * mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1))
          end if

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))   
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(2,2))
            end if
          end if

        case (CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dweightConvection * dgamma
          rstreamlineDiffusion2%ddeltaT  = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton  = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight. Compensate for any "-" sign in dgamma!
          rjumpStabil%dtheta = dgamma * mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1))
          end if

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(2,2))
            end if
          end if

        case (CCMASM_STAB_EDGEORIENTED3)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dweightConvection * dgamma
          rstreamlineDiffusion2%ddeltaT  = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton  = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1), p_Ddata1)
          call lsyssc_getbase_double (rmatrixEOJ, p_Ddata2)
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)
                              
          ! We use the precomputed EOJ matrix and sum it up to the
          ! existing matrix.
          ! Matrix weight. Compensate for any "-" sign in dgamma!
          call lsyssc_matrixLinearComb(rmatrixEOJ,&
              dgamma*mprim_signum(rstabilisation%dupsam),&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

          if (.not. bshared) then
            ! Also for the Y-velocity.
            call lsyssc_matrixLinearComb(rmatrixEOJ,&
                dgamma*mprim_signum(rstabilisation%dupsam),&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
          end if

        case default
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select
        
      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
      
        if (rstabilisation%dupsam .ne. 0.0_DP) then
          select case (rstabilisation%cupwind)
          case (CCMASM_STAB_EDGEORIENTED,CCMASM_STAB_EDGEORIENTED2)
            
            ! Set up the jump stabilisation structure.
            ! There's not much to do, only initialise the viscosity...
            rjumpStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
            
            ! Set stabilisation parameter
            rjumpStabil%dgamma = abs(rstabilisation%dupsam)
            
            ! Matrix weight. Compensate for any "-" sign in dgamma!
            rjumpStabil%dtheta = dgamma * mprim_signum(rstabilisation%dupsam)

            ! Call the jump stabilisation technique to stabilise that stuff.   
            ! We can assemble the jump part any time as it's independent of any
            ! convective parts...
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(1,1))
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1))
            end if

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                                  rjumpStabil, CONV_MODMATRIX, &
                                  rmatrix%RmatrixBlock(2,2))
              if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
                call smva_addBdEOJOperator (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(2,2))
              end if
            end if
            
          case (CCMASM_STAB_EDGEORIENTED3)
            
            ! We use the precomputed EOJ matrix and sum it up to the
            ! existing matrix.
            call lsyssc_matrixLinearComb(rmatrixEOJ,&
                dgamma*mprim_signum(rstabilisation%dupsam),&
                rmatrix%RmatrixBlock(1,1),1.0_DP,&
                rmatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            if (.not. bshared) then
              ! Also for the Y-velocity.
              call lsyssc_matrixLinearComb(rmatrixEOJ,&
                  dgamma*mprim_signum(rstabilisation%dupsam),&
                  rmatrix%RmatrixBlock(2,2),1.0_DP,&
                  rmatrix%RmatrixBlock(2,2),.false.,.false.,.true.,.true.)
            end if
            
          case default
            ! No stabilisation
          
          end select
        end if
        
      end if
      
    end subroutine  
      
    ! -----------------------------------------------------

    subroutine assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rmatrix, &
        rvector, dweight)
        
    ! Assembles a 2x2 block matrix with
    ! probably nonlinear mass matrices on the diagonal.
    ! The mass matrices are projected according to where the velocity
    ! vector -1/alphaC*rvector%RvectorBlock(4/5) (from rnonlinearSpatialMatrix)
    ! violates the bounds dmin/dmax.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
    
    ! Block matrix where the 2x2-velocity submatrix should be assembled.
    ! The matrix is cleared on call to this routine!
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Vector that specifies where to evaluate nonlinear terms
    type(t_vectorBlock), intent(IN), target :: rvector
    
    ! Weight for the (projected) mass matrices when adding them to the
    ! system.
    real(DP), intent(in) :: dweight
    
      ! local variables
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      type(t_bilinearForm) :: rform
      type(t_collection) :: rcollection
      integer, dimension(:), pointer :: p_IelementList
      integer :: ielemHandle,nelements
      type(t_bilfMatrixAssembly) :: rmatrixAssembly
      integer(I32) :: celement,ccubType
      integer :: ccontrolConstraints

      ! Assemble A14/A25? 
      if (dweight .ne. 0.0_DP) then
          
        ! Switch on these matrices
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
        
        ! Determine how to assemble...
        ccontrolConstraints = rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints
          
        ! Calculate the usual mass matrix if conrol constraints are deactivated
        ! or if Newton is not active.
        if ((ccontrolConstraints .eq. 0) .or. &
            (rnonlinearSpatialMatrix%cmatrixType .eq. 0)) then
      
          ! Copy the entries of the mass matrix. Share the structure.
          ! We must not share the entries as these might be changed by the caller
          ! e.g. due to boundary conditions!
          
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
              
          ! Scale the entries by dweight.
          if (dweight .ne. 1.0_DP) then
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),dweight)
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),dweight)
          end if
          
        else if (rnonlinearSpatialMatrix%cmatrixType .eq. 1) then
          
          select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints)
          
          case (1)
          
            ! Copy the entries of the mass matrix. Share the structure.
            ! We must not share the entries as these might be changed by the caller
            ! e.g. due to boundary conditions!
            
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
          
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
                
            ! Scale the entries by the weight.
            if (dweight .ne. 1.0_DP) then
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),dweight)
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),dweight)
            end if

            ! Filter the matrix. All the rows corresponding to DOF's that violate
            ! the bounds must be set to zero.
            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Constant bounds
              call massmatfilter (rmatrix%RmatrixBlock(1,1),rvector%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%dalphaC,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)
              call massmatfilter (rmatrix%RmatrixBlock(2,2),rvector%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%dalphaC,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)
            case (1)
              ! Variable bounds
              call massmatfilterVar (rmatrix%RmatrixBlock(1,1),rvector%RvectorBlock(4),&
                  rnonlinearSpatialMatrix%dalphaC,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))
              call massmatfilterVar (rmatrix%RmatrixBlock(2,2),rvector%RvectorBlock(5),&
                  rnonlinearSpatialMatrix%dalphaC,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
            case default
              ! Not implemented.
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          case (2)
          
            ! Exact reassembly of the mass matrices.
          
            ! In A11/A22 we have to create a 'projective mass matrix'.
            ! This is the derivative of a projection operator
            ! P[a,b](f)=a if f<a, =b if f>b, =f otherwise.
            ! For a<f<b, this is the mass matrix. Everywhere else, this is =0.
            ! We assemble this matrix just as a standard mass matrix with noconstant
            ! coefficients. Whereever u = -1/alpha * lambda is out of bounds,
            ! we return 0 as coefficient, otherwise 1.
          
            rform%itermCount = 1
            rform%Idescriptors(1,1) = DER_FUNC
            rform%Idescriptors(2,1) = DER_FUNC

            ! In this case, we have nonconstant coefficients.
            rform%ballCoeffConstant = .FALSE.
            rform%BconstantCoeff(:) = .FALSE.

            ! Prepare a collection structure to be passed to the callback
            ! routine. We attach the vector T in the quick-access variables
            ! so the callback routine can access it.
            ! The bounds and the alpha value are passed in the
            ! quickaccess-arrays.
            call collct_init(rcollection)
            rcollection%p_rvectorQuickAccess1 => rvector

            ! Coefficient is dmu1=1/alpha or 0, depending on lambda
            rcollection%DquickAccess(3)  = rnonlinearSpatialMatrix%dalphaC
            rcollection%DquickAccess(4)  = dweight
            
            ! At first, set up A14, depending on lambda_1.
            rcollection%IquickAccess(1) = 1
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
            
            ! Now we can build the matrix entries.
            ! We specify the callback function coeff_Laplace for the coefficients.
            ! As long as we use constant coefficients, this routine is not used.
            ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
            ! the framework will call the callback routine to get analytical
            ! data.
            ! The collection is passed as additional parameter. That's the way
            ! how we get the vector to the callback routine.
            call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(1,1),&
                coeff_ProjMass,rcollection)

            ! Now, set up A25, depending on lambda_2.
            rcollection%IquickAccess(1) = 2
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

            call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(2,2),&
                coeff_ProjMass,rcollection)
            
            ! Now we can forget about the collection again.
            call collct_done (rcollection)
            
          case (3)
          
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
            call lsyssc_duplicateMatrix (&
                rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
                rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Create the matrix
              call approxProjectionDerivative (rmatrix%RmatrixBlock(1,1), &
                  rvector%RvectorBlock(4), rnonlinearSpatialMatrix%dalphaC,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1,0.001_DP)

              call approxProjectionDerivative (rmatrix%RmatrixBlock(2,2), &
                  rvector%RvectorBlock(5), rnonlinearSpatialMatrix%dalphaC,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                  rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2,0.001_DP)
            case default
              ! Not implemented.
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select
          
            ! Scale the entries by the weight if necessary
            if (dweight .ne. 1.0_DP) then
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),dweight)
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),dweight)
            end if
            
          case (4)
          
            ! Exact reassembly of the mass matrices with adaptive integration.
          
            ! Create an array that saves all elements on the border of the active set.  
            nelements = rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%&
                RspatialDiscr(1)%p_rtriangulation%NEL
            call storage_new ('', 'Ielements', nelements, ST_INT, ielemhandle, &
                ST_NEWBLOCK_NOINIT)
          
            ! In A11/A22 we have to create a 'projective mass matrix'.
            ! This is the derivative of a projection operator
            ! P[a,b](f)=a if f<a, =b if f>b, =f otherwise.
            ! For a<f<b, this is the mass matrix. Everywhere else, this is =0.
            ! We assemble this matrix just as a standard mass matrix with noconstant
            ! coefficients. Whereever u = -1/alpha * lambda is out of bounds,
            ! we return 0 as coefficient, otherwise 1.
          
            rform%itermCount = 1
            rform%Idescriptors(1,1) = DER_FUNC
            rform%Idescriptors(2,1) = DER_FUNC

            ! In this case, we have nonconstant coefficients.
            rform%ballCoeffConstant = .FALSE.
            rform%BconstantCoeff(:) = .FALSE.

            ! Prepare a collection structure to be passed to the callback
            ! routine. We attach the vector T in the quick-access variables
            ! so the callback routine can access it.
            ! The bounds and the alpha value are passed in the
            ! quickaccess-arrays.
            call collct_init(rcollection)
            rcollection%p_rvectorQuickAccess1 => rvector
            
            ! Coefficient is dmu1=1/alpha or 0, depending on lambda
            rcollection%DquickAccess(3)  = rnonlinearSpatialMatrix%dalphaC
            rcollection%DquickAccess(4)  = dweight
            
            ! At first, set up A14, depending on lambda_1.
            rcollection%IquickAccess(1) = 1
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
            
            ! The IquickAccess(2) element saves the handle of the element list.
            ! IquickAccess(3) saves how many elements are collected.
            rcollection%IquickAccess(2) = ielemhandle
            rcollection%IquickAccess(3) = 0

            ! Now we can build the matrix entries.
            ! We specify the callback function coeff_Laplace for the coefficients.
            ! As long as we use constant coefficients, this routine is not used.
            ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
            ! the framework will call the callback routine to get analytical
            ! data.
            ! The collection is passed as additional parameter. That's the way
            ! how we get the vector to the callback routine.
            call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(1,1),&
                coeff_ProjMassCollect,rcollection)
                
            ! Assemble a submesh matrix on the elements in the list
            ! with a summed cubature formula.
            ! Note: Up to now, this works only for uniform meshes!
            if (rcollection%IquickAccess(3) .gt. 0) then
              celement = rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(1)%celement
              ccubType = rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(1)%ccubTypeBilForm
              call storage_getbase_int(ielemhandle,p_IelementList)
              call bilf_initAssembly(rmatrixAssembly,rform,celement,celement,&
                  cub_getSummedCubType(ccubType,1))
              call bilf_assembleSubmeshMatrix9(rmatrixAssembly,rmatrix%RmatrixBlock(1,1),&
                  p_IelementList(1:rcollection%IquickAccess(3)),coeff_ProjMass,rcollection)
              call bilf_doneAssembly(rmatrixAssembly)
            end if            

            ! Now, set up A25, depending on lambda_2.
            rcollection%IquickAccess(1) = 2
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

            ! Create a new element list
            rcollection%IquickAccess(3) = 0

            call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(2,2),&
                coeff_ProjMassCollect,rcollection)
            
            select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
            case (0)
              ! Assemble a submesh matrix on the elements in the list
              ! with a summed cubature formula.
              ! Note: Up to now, this works only for uniform meshes!
              if (rcollection%IquickAccess(3) .gt. 0) then
                celement = rmatrix%RmatrixBlock(2,2)%p_rspatialDiscrTest%RelementDistr(1)%celement
                ccubType = rmatrix%RmatrixBlock(2,2)%p_rspatialDiscrTest%RelementDistr(1)%ccubTypeBilForm
                call storage_getbase_int(ielemhandle,p_IelementList)
                call bilf_initAssembly(rmatrixAssembly,rform,celement,celement,&
                    cub_getSummedCubType(ccubType,1))
                call bilf_assembleSubmeshMatrix9(rmatrixAssembly,rmatrix%RmatrixBlock(2,2),&
                    p_IelementList(1:rcollection%IquickAccess(3)),coeff_ProjMass,rcollection)
                call bilf_doneAssembly(rmatrixAssembly)
              end if
            case default
              ! Not implemented.
              call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
              call sys_halt()
            end select

            ! Now we can forget about the collection again.
            call collct_done (rcollection)
            
            ! Release the element set.
            call storage_free (ielemHandle)

          case default
          
            call output_line ('Unsupported ccontrolConstraints flag.', &
                              OU_CLASS_ERROR,OU_MODE_STD,'assembleProjectedMassBlocks')
            call sys_halt()
            
          end select

        else
        
          call output_line ('Unsupported cmatrixType.', &
                            OU_CLASS_ERROR,OU_MODE_STD,'assembleProjectedMassBlocks')
          call sys_halt()
          
        end if
        
      else
      
        ! Switch off these matrices
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 0.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 0.0_DP
        
      end if
            
    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_assembleDefect (rnonlinearSpatialMatrix,rx,rd,cx)

!<description>
  ! This routine assembles the nonlinear defect
  !      rd := rd - cx A(ry) rx
  ! with the system matrix A(.) defined by the configuration in rnonlinearSpatialMatrix.
  ! The caller must initialise the rnonlinearSpatialMatrix according to how the 
  ! matrix should look like.
  !
  ! The parameters rvectorI are optional. If specified, these parameter defines where to
  ! evaluate the nonlinearity (if the system matrix $A$ contains a nonlinearity).
  ! The multiplication factor cx is optional as well. If not specified,
  ! cx=1.0 is assumed.
  !
  ! Note that for evaluating the nonlinearity / Newton matrix in
  ! the dual equation, the primal velocity is used!
  ! This is due to the fact that the 'nonlinearity' in the dual equation
  ! is defined as
  !       $$ y \grad(\lambda)  +   \lambda \grad(y)  $$
  ! so both parts are evaluated at the point y of the primal velocity
  ! before they are multiplied by the dual velocity $\lambda$!
  !
  ! The routine will not include any boundary conditions in the defect.
!</description>

  ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix. 
  !
  ! The caller must provide either p_rmatrixTemplateXXXX in this structure
  ! or set the p_rmatrixTemplateXXXX as well as p_rdiscretisation to
  ! appropriate values. This is necessary for exploiting then structure
  ! of the matrix.
  type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

  ! This vector specifies the 'x' that is multiplied to the matrix.
  type(t_vectorBlock), intent(IN), target :: rx

  ! OPTIONAL: Multiplication factor in front of the term 'A(ry) rx'.
  ! If not specified, cx=1.0 is assumed.
  real(DP), intent(IN), optional :: cx

!</input>

!<inputoutput>
  ! Destination vector. cx*A(ry)*rx is subtracted from this vector.
  type(t_vectorBlock), intent(INOUT) :: rd
!</inputoutput>
  
!</subroutine>

    ! local variables
    real(DP) :: dcx
    type(t_matrixBlock) :: rtempmatrix
    type(t_vectorBlock) :: rtempVectorX,rtempVectorB
    type(t_vectorBlock) :: rvectorPrimal,rvectorDual,rvectorDual2
    type(t_vectorBlock), pointer :: p_rprimalSol, p_rdualSol
    type(t_settings_stabil) :: rstabilisation    
    type(t_optcoperator) :: roptcoperator
    type(t_blockDIscretisation) :: rvelDiscr
    
    logical, parameter :: bnewmethod = .false.
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dd,p_Dx
    
    call lsysbl_getbase_double (rd,p_Dd)
    call lsysbl_getbase_double (rx,p_Dx)
    
    dcx = 1.0_DP
    if (present(cx)) dcx = cx
    if (dcx .eq. 0.0_DP) return
    
    if (.not. bnewmethod) then
    
      ! The system matrix looks like:
      !          
      !    ( A11  A12  B1  M~           ) 
      !    ( A21  A22  B2       M~      ) 
      !    ( B1^T B2^T                  ) 
      !    ( M             A44  A45  B1 ) 
      !    (      M        A54  A55  B2 ) 
      !    (               B1^T B2^T    ) 
      !
      ! We have now to manually perform a matrix-vector multiplication with
      ! all submatrices in that matrix. At first, we start with the
      ! linear part, that's the easiest.
      !
      ! ---------------------------------------------------
      ! 1.) Identity.
      !    ( I                          ) 
      !    (      I                     ) 
      !    (                            ) 
      !    (               I            ) 
      !    (                    I       ) 
      !    (                            ) 

      if (rnonlinearSpatialMatrix%Diota(1,1) .ne. 0.0_DP) then
        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(1), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Diota(1,1)*dcx, 1.0_DP)

        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(2), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Diota(1,1)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%Diota(2,2) .ne. 0.0_DP) then
        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(4), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Diota(2,2)*dcx, 1.0_DP)

        call lsyssc_vectorLinearComb (&
            rx%RvectorBlock(5), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Diota(2,2)*dcx, 1.0_DP)
      end if

      ! ---------------------------------------------------
      ! 2.) Mass matrices
      !    ( M                          ) 
      !    (      M                     ) 
      !    (                            ) 
      !    ( M             M            ) 
      !    (      M             M       ) 
      !    (                            ) 
      if (rnonlinearSpatialMatrix%Dalpha(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(1), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Dalpha(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(2), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Dalpha(1,1)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%Dalpha(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(4), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Dalpha(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(5), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Dalpha(2,2)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%Dalpha(2,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(1), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Dalpha(2,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rx%RvectorBlock(2), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Dalpha(2,1)*dcx, 1.0_DP)
      end if
      
      ! Don't do anything with Dalpha(1,2) -- the mass matrices here
      ! are probably nonlinear! We assemble them later!
      
      ! ---------------------------------------------------
      ! 3.) Stokes matrices
      !    ( L                          ) 
      !    (      L                     ) 
      !    (                            ) 
      !    (               L            ) 
      !    (                    L       ) 
      !    (                            ) 
      if (rnonlinearSpatialMatrix%Dtheta(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
            rx%RvectorBlock(1), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Dtheta(1,1)*dcx*rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu,&
            1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
            rx%RvectorBlock(2), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Dtheta(1,1)*dcx*rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu,&
            1.0_DP)
      end if
            
      if (rnonlinearSpatialMatrix%Dtheta(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
            rx%RvectorBlock(4), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Dtheta(2,2)*dcx*rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu,&
            1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixLaplace, &
            rx%RvectorBlock(5), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Dtheta(2,2)*dcx*rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu,&
            1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! 3.) B-matrices
      !    (           B1               ) 
      !    (           B2               ) 
      !    (                            ) 
      !    (                        B1  ) 
      !    (                        B2  ) 
      !    (                            ) 
      
      if (rnonlinearSpatialMatrix%Deta(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(3), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(3), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Deta(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(6), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Deta(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(6), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Deta(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! 4.) B^T-matrices
      !    (                            ) 
      !    (                            ) 
      !    ( B1^T B2^T                  ) 
      !    (                            ) 
      !    (                            ) 
      !    (              B1^T B2^T     ) 
      
      if (rnonlinearSpatialMatrix%Dtau(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(1), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(2), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Dtau(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(4), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%Dtau(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(5), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%Dtau(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! Now a slightly more advanced task for which we use a separate
      ! routine and some submatrices/vectors: The nonlinearity.
      !
      ! Get a reference to the correct velocity vector we have to use for the
      ! nonlinearity.
      select case (rnonlinearSpatialMatrix%iprimalSol)
      case (1)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorPrimal, 1,2,.true.)
      case (2)
        call lsysbl_deriveSubvector(&
          rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorPrimal, 1,2,.true.)
      case (3)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorPrimal, 1,2,.true.)
      end select

      select case (rnonlinearSpatialMatrix%idualSol)
      case (1)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorDual, 4,5,.true.)
      case (2)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorDual, 4,5,.true.)
      case (3)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorDual, 4,5,.true.)
      end select

      select case (rnonlinearSpatialMatrix%idualSol2)
      case (1)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1,rvectorDual2, 4,5,.true.)
      case (2)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2,rvectorDual2, 4,5,.true.)
      case (3)
        call lsysbl_deriveSubvector(&
            rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3,rvectorDual2, 4,5,.true.)
      end select

      ! Create a block discretisation by deriving it from the 'full' matrix.
      ! This will serve as a local discretisation structure for all
      ! velocity modifications.
      call spdiscr_deriveBlockDiscr (rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual, &
          rvelDiscr, 1, 2)

      ! Create a 2x2 matrix based on the structure of the FE space.
      ! The matrix does not need any entries, we only need the structure.
      call lsysbl_createMatBlockByDiscr (rvelDiscr,rtempMatrix)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEMOffdiag,&
          rtempMatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEMOffdiag,&
          rtempMatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixTemplateFEM,&
          rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      ! 1.) Primal equation, y*grad(.), probably + grad(.)*y
      !
      !    ( N(y) N(y)                  ) 
      !    ( N(y) N(y)                  ) 
      !    (                            ) 
      !    (                            ) 
      !    (                            ) 
      !    (                            ) 
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,1,2,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,1,2,.true.)

      call assembleConvectionDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,rvectorPrimal,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
          rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
          rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal,dcx,&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ1)    
      
      call lsysbl_releaseVector (rtempVectorX)
      call lsysbl_releaseVector (rtempVectorB)
      
      ! 2.) Dual equation, y*grad(.) + grad(.)^T*y
      !
      !    (                            ) 
      !    (                            ) 
      !    (                            ) 
      !    (              N*(y) N*(y)   ) 
      !    (              N*(y) N*(y)   ) 
      !    (                            ) 
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,4,5,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,4,5,.true.)

      call assembleConvectionDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,rvectorPrimal,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
          rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
          rnonlinearSpatialMatrix%rdiscrData%rstabilDual,dcx,&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplatesOptC%rmatrixEOJ2)    
      
      call lsysbl_releaseVector (rtempVectorX)
      call lsysbl_releaseVector (rtempVectorB)
      
      ! 3.) Dual equation, Newton term lambda*grad(.) + grad(.)^T*lambda
      !
      !    (                            ) 
      !    (                            ) 
      !    (                            ) 
      !    (   N*(l) N*(l)              ) 
      !    (   N*(l) N*(l)              ) 
      !    (                            ) 
      
      ! No stabilisation here
      ! rstabilisation = t_convecStabilisation(&
      !     rnonlinearSpatialMatrix%iupwind2,rnonlinearSpatialMatrix%dupsam2)
      rstabilisation = t_settings_stabil(0,0.0_DP,1)
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,1,2,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,4,5,.true.)

      call assembleConvectionDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,rvectorDual,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dgamma(2,1),rnonlinearSpatialMatrix%DgammaT(2,1),&
          rnonlinearSpatialMatrix%Dnewton(2,1),rnonlinearSpatialMatrix%DnewtonT(2,1),&
          rstabilisation,dcx)    
      
      ! There is probably a 2nd reactive term involved stemming from
      ! the next timestep when Crank-Nicolson is used.

      rstabilisation = t_settings_stabil(0,0.0_DP,1)
      
      call assembleConvectionDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,rvectorDual2,rtempVectorX,rtempVectorB,&
          0.0_DP,rnonlinearSpatialMatrix%DgammaT2(2,1),&
          rnonlinearSpatialMatrix%Dnewton2(2,1),0.0_DP,&
          rstabilisation,dcx)    
      
      call lsysbl_releaseVector (rtempVectorX)
      call lsysbl_releaseVector (rtempVectorB)
      
      ! 4.) (Projected) mass matrix
      !
      !    (              PM(l)         ) 
      !    (                    PM(l)   ) 
      !    (                            ) 
      !    (                            ) 
      !    (                            ) 
      !    (                            ) 
      
      ! What's the type of the current matrix? Is this a Newton-matrix?
      if (rnonlinearSpatialMatrix%cmatrixType .eq. 0) then

        ! No, this is a standard matrix. That means, we just have to project
        ! the control u and multiply it with the mass matrix.

        ! Copy our solution vector \lambda. Scale it by -1/alpha.
        call lsysbl_deriveSubvector(rx,rtempVectorX,4,5,.false.)
        call lsysbl_scaleVector(rtempVectorX,-rnonlinearSpatialMatrix%Dalpha(1,2))
        
        ! Project that to the allowed range.
        if (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints .ne. 0) then

          select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
          case (0)
            call smva_projectControlTstepConst (rtempVectorX%RvectorBlock(1),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)
            call smva_projectControlTstepConst (rtempVectorX%RvectorBlock(2),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)
          case (1)
            call smva_projectControlTstepVec (rtempVectorX%RvectorBlock(1),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))
            call smva_projectControlTstepVec (rtempVectorX%RvectorBlock(2),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
          case default
            ! Not implemented.
            call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
            call sys_halt()
          end select
          
        end if

        ! Now carry out MV and include it to the defect.
        ! Note that the multiplication factor is -(-cx) = cx because
        ! it's put on the RHS of the system for creating the defect.
        ! d = b - cx A x = b - ... + \nu Laplace(y) - y\grad(y) - grad(p) + P(-1/alpha lambda)
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rtempVectorX%RvectorBlock(1), &
            rd%RvectorBlock(1), dcx, 1.0_DP)
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity, &
            rtempVectorX%RvectorBlock(2), &
            rd%RvectorBlock(2), dcx, 1.0_DP)

        call lsysbl_releaseVector (rtempVectorX)

      else

        ! Yes, that's a Newton matrix. That means, we have to multiply the
        ! vector with the derivative of the projection operator:
        ! b-(-P[a,b]'(-1/alpha lambda)).
        ! For that purpose, we have to assemble special mass matrices:
        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          call assemblePrimalUConstrMassDefect (&
              rnonlinearSpatialMatrix,rx,&
              rd,dcx*rnonlinearSpatialMatrix%Dalpha(1,2),&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1)
        case (2)
          call assemblePrimalUConstrMassDefect (&
              rnonlinearSpatialMatrix,rx,&
              rd,dcx*rnonlinearSpatialMatrix%Dalpha(1,2),&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2)
        case (3)
          call assemblePrimalUConstrMassDefect (&
              rnonlinearSpatialMatrix,rx,&
              rd,dcx*rnonlinearSpatialMatrix%Dalpha(1,2),&
              rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3)
        end select
        
      end if
      
      ! Release the temp vectors/matrices, that's it.
      call spdiscr_releaseBlockDiscr(rvelDiscr)
      call lsysbl_releaseVector (rvectorDual2)
      call lsysbl_releaseVector (rvectorDual)
      call lsysbl_releaseVector (rvectorPrimal)
      call lsysbl_releaseMatrix (rtempMatrix)

    else
    
      ! ---------------------------------------------------
      ! 3.) B-matrices
      !    (           B1               ) 
      !    (           B2               ) 
      !    (                            ) 
      !    (                        B1  ) 
      !    (                        B2  ) 
      !    (                            ) 
      
      if (rnonlinearSpatialMatrix%Deta(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(3), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(3), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Deta(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB1, &
            rx%RvectorBlock(6), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Deta(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixB2, &
            rx%RvectorBlock(6), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Deta(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! 4.) B^T-matrices
      !    (                            ) 
      !    (                            ) 
      !    ( B1^T B2^T                  ) 
      !    (                            ) 
      !    (                            ) 
      !    (              B1^T B2^T     ) 
      
      if (rnonlinearSpatialMatrix%Dtau(1,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(1), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(2), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Dtau(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD1, &
            rx%RvectorBlock(4), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%Dtau(2,2)*dcx, 1.0_DP)
        call lsyssc_scalarMatVec (&
            rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixD2, &
            rx%RvectorBlock(5), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%Dtau(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! Now a slightly more advanced task for which we use a separate
      ! routine and some submatrices/vectors: The nonlinearity.

      ! Initialise the operator structure for what we need.
      roptcoperator%dupsamPrimal = rnonlinearSpatialMatrix%rdiscrData%rstabilPrimal%dupsam
      roptcoperator%dupsamDual = rnonlinearSpatialMatrix%rdiscrData%rstabilDual%dupsam
      
      ! Timestep-weights
      roptcoperator%dprimalAlpha = rnonlinearSpatialMatrix%Dalpha(1,1)
      roptcoperator%ddualAlpha   = rnonlinearSpatialMatrix%Dalpha(2,2)

      ! Stokes operator
      roptcoperator%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
      roptcoperator%dprimalBeta = rnonlinearSpatialMatrix%Dtheta(1,1)
      roptcoperator%ddualBeta   = rnonlinearSpatialMatrix%Dtheta(2,2)
      
      ! Nonlinearity
      if (rnonlinearSpatialMatrix%Dgamma(1,1) .ne. 0.0_DP) then
        roptcoperator%dprimalDelta = rnonlinearSpatialMatrix%Dgamma(1,1)
        roptcoperator%ddualDelta   = rnonlinearSpatialMatrix%Dgamma(2,2)
        roptcoperator%ddualNewtonTrans = rnonlinearSpatialMatrix%DnewtonT(2,2)
        
        ! Newton implies additional operators.
        if (rnonlinearSpatialMatrix%Dnewton(1,1) .ne. 0.0_DP) then
          roptcoperator%dprimalNewton    = 1.0_DP
          roptcoperator%ddualRDeltaTrans = 1.0_DP
          roptcoperator%ddualRNewton     = -1.0_DP
        end if
        
      end if
      
      ! Coupling matrices
      !if (rparams%bdualcoupledtoprimal) then
        roptcoperator%ddualRAlpha = rnonlinearSpatialMatrix%Dalpha(2,1)
      !end if

      !if (rparams%bcontrolactive) then
        roptcoperator%dcontrolWeight = -rnonlinearSpatialMatrix%Dalpha(1,2)*rnonlinearSpatialMatrix%dalphaC
        roptcoperator%dcontrolMultiplier = -1.0_DP/rnonlinearSpatialMatrix%dalphaC
      !end if
      
      if (rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints .ne. 0) then
        roptcoperator%ccontrolProjection = &
            rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints

        roptcoperator%cconstraintsType = &
            rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType

        roptcoperator%dmin1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
        roptcoperator%dmax1 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
        roptcoperator%dmin2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
        roptcoperator%dmax2 = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

        roptcoperator%p_rumin1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin1
        roptcoperator%p_rumax1 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax1
        roptcoperator%p_rumin2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumin2
        roptcoperator%p_rumax2 => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rumax2

        roptcoperator%p_rvectorumin => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin
        roptcoperator%p_rvectorumax => rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax
      end if
      
      select case (rnonlinearSpatialMatrix%iprimalSol)
      case (1)
        p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
      case (2)
        p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
      case (3)
        p_rprimalSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
      end select
      
      select case (rnonlinearSpatialMatrix%idualSol)
      case (1)
        p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector1
      case (2)
        p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector2
      case (3)
        p_rdualSol => rnonlinearSpatialMatrix%p_rnonlinearity%p_rvector3
      end select
      
      ! Calculate the velocity-dependent part of the system matrix.
      call conv_strdiffOptC2dgetDefect (&
          rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
          roptcoperator,p_rprimalSol,p_rdualSol,dcx,rx,rd)
    
    end if
    
  contains

    ! -----------------------------------------------------

    subroutine assembleConvectionDefect (&
        rnonlinearSpatialMatrix,rmatrix,rvector,rx,rb,dgamma,dgammaT,dnewton,dnewtonT,&
        rstabilisation,dcx,rmatrixEOJ)
        
    ! Assembles the convection matrix in the block matrix rmatrix at position (1,1):
    !
    ! rmatrix := dgamma*N(rvector) + dgammaT*N^t(rvector) +
    !            dnewton*N*(rvector) + dnewtonT*N*^t(rvector)
    !
    ! Even if no nonlinearity is present, the routine can be used to
    ! add stabilisation into the matrix.
    
    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how to set up the matrix. 
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
    
    ! 2X2 block matrix that specifies the structure of the velocity FE space.
    type(t_matrixBlock), intent(INOUT) :: rmatrix
    
    ! Velocity vector for the nonlinearity. Must be specified if
    ! GAMMA <> 0; can be omitted if GAMMA=0.
    type(t_vectorBlock), intent(in) :: rvector
    
    ! The current solution vector for the velocity (x- and y-velocity)
    type(t_vectorBlock), intent(in) :: rx

    ! The RHS vector; a defect will be created in this vector.
    type(t_vectorBlock), intent(inout) :: rb
    
    ! Weight for the nonlinear term u\grad(.)
    real(DP), intent(in) :: dgamma

    ! Weight for the nonlinear term u(\grad(.))^t
    real(DP), intent(in) :: dgammaT

    ! Weight for the nonlinear term (\grad(.))u
    real(DP), intent(in) :: dnewton

    ! Weight for the nonlinear term (\grad(.))^t u
    real(DP), intent(in) :: dnewtonT
    
    ! Weight for the operator when multiplying: d = b + dcx * A x. Standard = -1.0_DP
    real(DP), intent(in) :: dcx
    
    ! Stabilisation parameters
    type(t_settings_stabil), intent(in) :: rstabilisation

    ! OPTIONAL: EOJ matrix in case the precomputed EOJ stabilisation is active.
    type(t_matrixScalar), intent(in), optional :: rmatrixEOJ
    
    ! local variables
    logical :: bshared
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_convStreamDiff2) :: rstreamlineDiffusion3
    type(t_jumpStabilisation) :: rjumpStabil
    type(t_convupwind) :: rupwindStabil
    real(DP) :: dweightConvection
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_Ddata1,p_Ddata2
    call lsysbl_getbase_double (rx,p_Ddata1)
    call lsysbl_getbase_double (rb,p_Ddata2)
    
      ! Debug weight for the convection
      dweightConvection = rnonlinearSpatialMatrix%rdiscrData%p_rdebugFlags%dweightConvection
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      if ((dgamma .ne. 0.0_DP) .or. (dgammaT .ne. 0.0_DP) .or. &
          (dnewton .ne. 0.0_DP) .or. (dnewtonT .ne. 0.0_DP)) then
        select case (rstabilisation%cupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          rstreamlineDiffusion%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion%ddelta            = dweightConvection * dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dweightConvection * dgammaT
          rstreamlineDiffusion%dnewton           = dweightConvection * dnewton
          rstreamlineDiffusion%dnewtonTransposed = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rx,rb)
                              
        case (CCMASM_STAB_STREAMLINEDIFF2)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          rstreamlineDiffusion2%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dweightConvection * dgamma
          rstreamlineDiffusion2%ddeltaT = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvector)

        case (CCMASM_STAB_UPWIND)
        
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          rstreamlineDiffusion2%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = 0.0_DP
          rstreamlineDiffusion2%ddeltaT = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity except for the
          ! convection part. As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvector)

          ! Prepare the upwind structure for the assembly of the convection.
          ! Note: Stabilisation weight is actually wrong, but it is not possible
          ! to specify in rupwindStabil%dupsam whether the stabilisation
          ! is added or subtracted!
          rupwindStabil%dupsam = abs(rstabilisation%dupsam)
          rupwindStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          rupwindStabil%dtheta = dcx*dgamma
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Apply the upwind operator.
          call conv_upwind2d (rvector, rvector, 1.0_DP, 0.0_DP,&
              rupwindStabil, CONV_MODDEFECT, rmatrix%RmatrixBlock(1,1), rx, rb)

          ! Prepare the upwind structure for the assembly of the convection.
          !rstreamlineDiffusion3%dupsam = rstabilisation%dupsam
          !rstreamlineDiffusion3%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          !rstreamlineDiffusion3%dtheta = dgamma*dcx
          !rstreamlineDiffusion3%ddelta = 1.0_DP
          
          ! Apply the upwind operator.
          !call conv_streamDiff2Blk2dDef (rstreamlineDiffusion3,rmatrix,&
          !    rx,rb,rvector)

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          rstreamlineDiffusion%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta            = dweightConvection * dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dweightConvection * dgammaT
          rstreamlineDiffusion%dnewton           = dweightConvection * dnewton
          rstreamlineDiffusion%dnewtonTransposed = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamlineDiffusionBlk2d (&
                              rvector, rvector, &
                              1.0_DP, 0.0_DP,&
                              rstreamlineDiffusion, CONV_MODDEFECT, &
                              rmatrix,rx,rb)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight. Compensate for any "-" sign in dgamma!
          rjumpStabil%dtheta = dcx*dgamma*mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rx,rb)
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJvector (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1),rx,rb)
          end if

!          if (.not. bshared) then
!            call conv_jumpStabilisation2d (&
!                                rjumpStabil, CONV_MODDEFECT, &
!                                rmatrix%RmatrixBlock(2,2),rx,rb)   
!          end if

        case (CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          rstreamlineDiffusion2%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dweightConvection * dgamma
          rstreamlineDiffusion2%ddeltaT  = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton  = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvector)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = abs(rstabilisation%dupsam)
          
          ! Matrix weight
          rjumpStabil%dtheta = dcx*dgamma*mprim_signum(rstabilisation%dupsam)

          ! Call the jump stabilisation technique to stabilise that stuff.   
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          !
          ! The routine processes the X as well as the Y velocity!
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rx,rb)
          if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
            call smva_addBdEOJvector (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1),rx,rb)
          end if

        case (CCMASM_STAB_EDGEORIENTED3)
        
          ! Jump stabilisation, precomputed matrix.
          !
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
          
          rstreamlineDiffusion2%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dweightConvection * dgamma
          rstreamlineDiffusion2%ddeltaT  = dweightConvection * dgammaT
          rstreamlineDiffusion2%dnewton  = dweightConvection * dnewton
          rstreamlineDiffusion2%dnewtonT = dweightConvection * dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvector)
                              
          ! Just do some mat-vec's to compute the defect.
          !
          ! The stabilisation parameter is already incorporated into
          ! the matrix -- only its sign is not incorporated, so we do that here.
          call lsyssc_scalarMatVec(rmatrixEOJ,&
              rx%RvectorBlock(1),rb%RvectorBlock(1),&
              -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

          call lsyssc_scalarMatVec(rmatrixEOJ,&
              rx%RvectorBlock(2),rb%RvectorBlock(2),&
              -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

        case default
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select

      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
        if (rstabilisation%dupsam .ne. 0.0_DP) then
          select case (rstabilisation%cupwind)
          case (CCMASM_STAB_EDGEORIENTED,CCMASM_STAB_EDGEORIENTED2)
            
            ! Set up the jump stabilisation structure.
            ! There's not much to do, only initialise the viscosity...
            rjumpStabil%dnu = rnonlinearSpatialMatrix%rdiscrData%rphysicsPrimal%dnu
            
            ! Set stabilisation parameter
            rjumpStabil%dgamma = abs(rstabilisation%dupsam)
            
            ! Matrix weight. Compensate for any "-" sign in dgamma!
            rjumpStabil%dtheta = dcx*dgamma*mprim_signum(rstabilisation%dupsam)

            ! Call the jump stabilisation technique to stabilise that stuff.   
            ! We can assemble the jump part any time as it's independent of any
            ! convective parts...
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODDEFECT, &
                                rmatrix%RmatrixBlock(1,1),rx,rb)
            if (rstabilisation%ceojStabilOnBoundary .eq. 0) then
              call smva_addBdEOJvector (rjumpStabil,-1.0_DP,rmatrix%RmatrixBlock(1,1),rx,rb)
            end if

!            if (.not. bshared) then
!              call conv_jumpStabilisation2d (&
!                                  rjumpStabil, CONV_MODDEFECT, &
!                                  rmatrix%RmatrixBlock(2,2),rx,rb)   
!            end if
            
          case (CCMASM_STAB_EDGEORIENTED3)
          
            ! Jump stabilisation, precomputed matrix.
            !
            ! Just do some mat-vec's to compute the defect.
            ! The stabilisation parameter is already incorporated into
            ! the matrix -- only its sign is not incorporated, so we do that here.
            call lsyssc_scalarMatVec(rmatrixEOJ,&
                rx%RvectorBlock(1),rb%RvectorBlock(1),&
                -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

            call lsyssc_scalarMatVec(rmatrixEOJ,&
                rx%RvectorBlock(2),rb%RvectorBlock(2),&
                -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

          case default
            ! No stabilisation
          
          end select
        end if
        
      end if
      
    end subroutine  

    subroutine assemblePrimalUConstrMassDefect (&
        rnonlinearSpatialMatrix,rvector,rdefect,dcx,rvelocityVector) !,cprojectionType)
        
    ! Assembles the defect arising from the projective coupling mass
    ! matrices in the primal equation which comes from constraints on u. 
    ! rdefect must have been initialised with the right hand side vector.
    !
    ! Let the mass matrix M~ be the usual mass matrix where u is ok
    ! and the 0-operator where u is out of bounds.
    ! Then, we assemble
    !
    !       rdefect = r(primal)defect - dcx (dmu1 M~ r(dual)vector)

    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
    ! about how the mass matrix part is weighted (dr11, dr12).
    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

    ! Solution vector. 
    type(t_vectorBlock), intent(IN) :: rvector
    
    ! On entry: RHS vector.
    ! Is overwritten by the defect vector in the velocity subsystem.
    type(t_vectorBlock), intent(INOUT) :: rdefect
    
    ! Multiplication factor for the whole operator A*rvector
    real(DP), intent(IN) :: dcx
    
    ! Velocity vector field that should be used for the assembly of the
    ! nonlinearity. Block 1 and 2 in that block vector are used as velocity 
    ! field and shall contain the dual velocity field.
    type(t_vectorBlock), intent(IN), target :: rvelocityVector

    ! Type of projection.
    ! =0: Simple DOF-based projection, rows in the mass matrices are
    !     replaced by zero (faster)
    ! =1: Reassembly of the mass matrices, projection based on cubature points.
    !integer, intent(in) :: cprojectionType

      ! local variables
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
      type(t_collection) :: rcollection
      type(t_bilinearForm) :: rform
      integer, dimension(:), pointer :: p_IelementList
      integer :: ielemHandle,nelements
      type(t_bilfMatrixAssembly) :: rmatrixAssembly
      integer(I32) :: celement,ccubType
      integer :: ccontrolConstraints
      
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if (dcx .ne. 0.0_DP) then

        call lsysbl_createEmptyMatrix (rtempMatrix,2)
        
        ccontrolConstraints = rnonlinearSpatialMatrix%rdiscrData%rconstraints%ccontrolConstraints

        ! The ccontrolConstraints decides on whether we use the 'quick-setup'
        ! method for the mass matrices or rebuild them again.
        ! ccontrolConstraints=0 uses standard mass matrices.
        select case (ccontrolConstraints)
        case (0) 
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
           
          call lsysbl_updateMatStrucInfo (rtempMatrix)
             
        case (1)
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
           
          call lsysbl_updateMatStrucInfo (rtempMatrix)
             
!          ! Scale the entries by dmu1.
!          if (rnonlinearSpatialMatrix%dmu1 .ne. 1.0_DP) then
!            call lsyssc_scaleMatrix (rtempMatrix%RmatrixBlock(1,1),rnonlinearSpatialMatrix%dmu1)
!            call lsyssc_scaleMatrix (rtempMatrix%RmatrixBlock(2,2),rnonlinearSpatialMatrix%dmu1)
!          end if

          ! In the case when we have constraints, filter the matrix. 
          ! All the rows corresponding to DOF's that violate
          ! the bounds must be set to zero.
          select case (rnonlinearSpatialMatrix%rdiscrData%rconstraints%cconstraintsType)
          case (0)
            ! Constant bounds
            call massmatfilter (rtempMatrix%RmatrixBlock(1,1),rvelocityVector%RvectorBlock(4),&
                rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1)
            call massmatfilter (rtempMatrix%RmatrixBlock(2,2),rvelocityVector%RvectorBlock(5),&
                rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2)
          case (1)
            ! Variable bounds
            call massmatfilterVar (rtempMatrix%RmatrixBlock(1,1),rvelocityVector%RvectorBlock(4),&
                rnonlinearSpatialMatrix%dalphaC,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(1),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(1))
            call massmatfilterVar (rtempMatrix%RmatrixBlock(2,2),rvelocityVector%RvectorBlock(5),&
                rnonlinearSpatialMatrix%dalphaC,&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumin%RvectorBlock(2),&
                rnonlinearSpatialMatrix%rdiscrData%rconstraints%p_rvectorumax%RvectorBlock(2))
          case default
            ! Not implemented.
            call output_line("CONSTRAINTS TO BE IMPLEMENTED!!!")
            call sys_halt()
          end select
            
        case (2) 

          ! Create a matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly      
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsysbl_updateMatStrucInfo (rtempMatrix)

          ! Assemble the modified mass matrices.
          
          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC

          ! In this case, we have nonconstant coefficients.
          rform%ballCoeffConstant = .FALSE.
          rform%BconstantCoeff(:) = .FALSE.

          ! Prepare a collection structure to be passed to the callback
          ! routine. We attach the vector T in the quick-access variables
          ! so the callback routine can access it.
          ! The bounds and the alpha value are passed in the
          ! quickaccess-arrays.
          call collct_init(rcollection)
          rcollection%p_rvectorQuickAccess1 => rvelocityVector

          ! Coefficient is dmu1=1/alpha or 0, depending on lambda
          rcollection%DquickAccess(3)  = rnonlinearSpatialMatrix%dalphaC
          rcollection%DquickAccess(4)  = 1.0_DP
          
          ! At first, set up A14, depending on lambda_1.
          rcollection%IquickAccess(1) = 1
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
          
          ! Now we can build the matrix entries.
          ! We specify the callback function coeff_Laplace for the coefficients.
          ! As long as we use constant coefficients, this routine is not used.
          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
          ! the framework will call the callback routine to get analytical
          ! data.
          ! The collection is passed as additional parameter. That's the way
          ! how we get the vector to the callback routine.
          call bilf_buildMatrixScalar (rform,.TRUE.,rtempmatrix%RmatrixBlock(1,1),&
              coeff_ProjMass,rcollection)

          ! Now, set up A22, depending on lambda_2.
          rcollection%IquickAccess(1) = 2
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

          call bilf_buildMatrixScalar (rform,.TRUE.,rtempmatrix%RmatrixBlock(2,2),&
              coeff_ProjMass,rcollection)
          
          ! Now we can forget about the collection again.
          call collct_done (rcollection)
          
        case (3)
        
          ! Create a matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly      
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsysbl_updateMatStrucInfo (rtempMatrix)

          call approxProjectionDerivative (rtempMatrix%RmatrixBlock(1,1), &
              rvelocityVector%RvectorBlock(4),rnonlinearSpatialMatrix%dalphaC,&
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1,&
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1,0.001_DP)

          call approxProjectionDerivative (rtempMatrix%RmatrixBlock(2,2), &
              rvelocityVector%RvectorBlock(5),rnonlinearSpatialMatrix%dalphaC,&
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2,&
              rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2,0.001_DP)
        
          call lsysbl_updateMatStrucInfo (rtempMatrix)
            
        case (4)
        
          ! Exact reassembly of the mass matrices with adaptive integration.

          ! Create a matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly      
          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_duplicateMatrix (&
              rnonlinearSpatialMatrix%rdiscrData%p_rstaticAsmTemplates%rmatrixMassVelocity,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsysbl_updateMatStrucInfo (rtempMatrix)

          ! Create an array that saves all elements on the border of the active set.  
          nelements = rnonlinearSpatialMatrix%rdiscrData%p_rdiscrPrimalDual%&
              RspatialDiscr(1)%p_rtriangulation%NEL
          call storage_new ('', 'Ielements', nelements, ST_INT, ielemhandle, &
              ST_NEWBLOCK_NOINIT)
        
          ! In A11/A22 we have to create a 'projective mass matrix'.
          ! This is the derivative of a projection operator
          ! P[a,b](f)=a if f<a, =b if f>b, =f otherwise.
          ! For a<f<b, this is the mass matrix. Everywhere else, this is =0.
          ! We assemble this matrix just as a standard mass matrix with noconstant
          ! coefficients. Whereever u = -1/alpha * lambda is out of bounds,
          ! we return 0 as coefficient, otherwise 1.
        
          rform%itermCount = 1
          rform%Idescriptors(1,1) = DER_FUNC
          rform%Idescriptors(2,1) = DER_FUNC

          ! In this case, we have nonconstant coefficients.
          rform%ballCoeffConstant = .FALSE.
          rform%BconstantCoeff(:) = .FALSE.

          ! Prepare a collection structure to be passed to the callback
          ! routine. We attach the vector T in the quick-access variables
          ! so the callback routine can access it.
          ! The bounds and the alpha value are passed in the
          ! quickaccess-arrays.
          call collct_init(rcollection)
          rcollection%p_rvectorQuickAccess1 => rvelocityVector

          ! Coefficient is dmu1=1/alpha or 0, depending on lambda
          rcollection%DquickAccess(3)  = rnonlinearSpatialMatrix%dalphaC
          rcollection%DquickAccess(4)  = 1.0_DP
          
          ! At first, set up A14, depending on lambda_1.
          rcollection%IquickAccess(1) = 1
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin1
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax1
          
          ! The IquickAccess(2) element saves the handle of the element list.
          ! IquickAccess(3) saves how many elements are collected.
          rcollection%IquickAccess(2) = ielemhandle
          rcollection%IquickAccess(3) = 0

          ! Now we can build the matrix entries.
          ! We specify the callback function coeff_Laplace for the coefficients.
          ! As long as we use constant coefficients, this routine is not used.
          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
          ! the framework will call the callback routine to get analytical
          ! data.
          ! The collection is passed as additional parameter. That's the way
          ! how we get the vector to the callback routine.
          call bilf_buildMatrixScalar (rform,.TRUE.,rtempmatrix%RmatrixBlock(1,1),&
              coeff_ProjMassCollect,rcollection)
              
          ! Assemble a submesh matrix on the elements in the list
          ! with a summed cubature formula.
          ! Note: Up to now, this works only for uniform meshes!
          if (rcollection%IquickAccess(3) .gt. 0) then
            celement = &
              rtempmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(1)%celement
            ccubType = &
              rtempmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(1)%ccubTypeBilForm
            call storage_getbase_int(ielemhandle,p_IelementList)
            call bilf_initAssembly(rmatrixAssembly,rform,celement,celement,&
                cub_getSummedCubType(ccubType,1))
            call bilf_assembleSubmeshMatrix9(rmatrixAssembly,rtempmatrix%RmatrixBlock(1,1),&
                p_IelementList(1:rcollection%IquickAccess(3)),coeff_ProjMass,rcollection)
            call bilf_doneAssembly(rmatrixAssembly)
          end if            

          ! Now, set up A25, depending on lambda_2.
          rcollection%IquickAccess(1) = 2
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumin2
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%rdiscrData%rconstraints%dumax2

          ! Create a new element list
          rcollection%IquickAccess(3) = 0

          call bilf_buildMatrixScalar (rform,.TRUE.,rtempmatrix%RmatrixBlock(2,2),&
              coeff_ProjMassCollect,rcollection)
          
          ! Assemble a submesh matrix on the elements in the list
          ! with a summed cubature formula.
          ! Note: Up to now, this works only for uniform meshes!
          if (rcollection%IquickAccess(3) .gt. 0) then
            celement = rtempmatrix%RmatrixBlock(2,2)%p_rspatialDiscrTest%RelementDistr(1)%celement
            ccubType = rtempmatrix%RmatrixBlock(2,2)%p_rspatialDiscrTest%RelementDistr(1)%ccubTypeBilForm
            call storage_getbase_int(ielemhandle,p_IelementList)
            call bilf_initAssembly(rmatrixAssembly,rform,celement,celement,&
                cub_getSummedCubType(ccubType,1))
            call bilf_assembleSubmeshMatrix9(rmatrixAssembly,rtempmatrix%RmatrixBlock(2,2),&
                p_IelementList(1:rcollection%IquickAccess(3)),coeff_ProjMass,rcollection)
            call bilf_doneAssembly(rmatrixAssembly)
          end if            

          ! Now we can forget about the collection again.
          call collct_done (rcollection)
          
          ! Release the element set.
          call storage_free (ielemHandle)

        case default
        
          ! Cancel
          call lsysbl_releaseMatrix (rtempMatrix)
          return
          
        end select

        ! Create a temporary block vector that points to the dual velocity.
        ! This has to be evaluated during the assembly.
        call lsysbl_deriveSubvector (rvector,rtempvectorEval,4,5,.true.)
        
        ! Create a temporary block vector for the dual defect.
        ! Matrix*primal velocity is subtracted from this.
        call lsysbl_deriveSubvector (rdefect,rtempvectorDef,1,2,.true.)

        ! Create the defect
        call lsysbl_blockMatVec (rtempmatrix, rtempvectorEval, rtempvectorDef, -dcx, 1.0_DP)
        
        ! Release memory
        call lsysbl_releaseVector (rtempvectorDef)
        call lsysbl_releaseVector (rtempvectorEval)
        call lsysbl_releaseMatrix (rtempMatrix)
      
      end if
      
    end subroutine

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_projectControlTstepConst (rdualSolution,dumin,dumax)

!<description>
  ! Projects a dual solution vector u in such a way, that
  ! dumin <= u <= dumax holds.
!</description>

!<input>
  ! Minimum value for u
  real(DP), intent(in) :: dumin

  ! Maximum value for u
  real(DP), intent(in) :: dumax
!</input>

!<inputoutput>
  ! Vector to be restricted
  type(t_vectorScalar), intent(inout) :: rdualSolution
!</inputoutput>

!</subroutine>
 
    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: i
    
    ! Get the vector array
    call lsyssc_getbase_double (rdualSolution,p_Ddata)
    
    ! Restrict the vector
    do i=1,rdualSolution%NEQ
      p_Ddata(i) = min(max(p_Ddata(i),dumin),dumax)
    end do

  end subroutine   

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_projectControlTstepVec (rdualSolution,rvectorumin,rvectorumax)

!<description>
  ! Projects a dual solution vector u in such a way, that
  ! dumin <= u <= dumax holds.
!</description>

!<input>
  ! Vector specifying the minimum value for u in each DOF
  type(t_vectorScalar), intent(in) :: rvectorumin

  ! Vector specifying the maximum value for u in each DOF
  type(t_vectorScalar), intent(in) :: rvectorumax
!</input>

!<inputoutput>
  ! Vector to be restricted
  type(t_vectorScalar), intent(inout) :: rdualSolution
!</inputoutput>

!</subroutine>
 
    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata, p_DdataMin, p_DdataMax
    integer :: i
    
    ! Get the vector array
    call lsyssc_getbase_double (rdualSolution,p_Ddata)
    call lsyssc_getbase_double (rvectorumin,p_DdataMin)
    call lsyssc_getbase_double (rvectorumax,p_DdataMax)
    
    ! Restrict the vector
    do i=1,rdualSolution%NEQ
      p_Ddata(i) = min(max(p_Ddata(i),p_DdataMin(i)),p_DdataMax(i))
    end do

  end subroutine   
 
  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_initNonlinearData (rnonlinearData,rvector1,rvector2,rvector3)

!<description>
  ! Initialises a nonlinear-data structure that defines the nonlinearity
  ! in a nonlinear matrix.
!</description>

!<input>
    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'previous' timestep. If there is no previous timestep (e.g.
    ! like in the 0th timestep), the vector can be undefined.
    type(t_vectorBlock), intent(in), target :: rvector1

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'current' timestep. 
    type(t_vectorBlock), intent(in), target :: rvector2

    ! Specifies where to evaluate the nonlinearity and must contain the data
    ! for the 'next' timestep. If there is no next timestep (e.g.
    ! like in the last timestep), the vector can be undefined.
    type(t_vectorBlock), intent(in), target :: rvector3
!</input>

!<inputoutput>
  ! Nonlinear-data structure to be created.
  type(t_spatialMatrixNonlinearData), intent(out) :: rnonlinearData
!</inputoutput>

!</subroutine>
 
    ! Set the pointers.
    rnonlinearData%p_rvector1 => rvector1
    rnonlinearData%p_rvector2 => rvector2
    rnonlinearData%p_rvector3 => rvector3

  end subroutine   

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_addBdEOJOperator (rjumpStabil,dweight,rmatrix)

!<description>
  ! Adds the jump stabilisation on all boundary edges to a matrix.
!</description>

!<input>
    ! Jump stabilisation structure that defines the jump stabilisation operator.
    type(t_jumpStabilisation), intent(in) :: rjumpStabil
    
    ! Weight of the operator. Default = 1.0
    real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
    ! Matrix to be modified.
    type(t_matrixScalar), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_jumpStabilisation) :: rjumpStabilLocal
    type(t_boundaryRegion) :: rboundaryRegion
    integer, dimension(:), allocatable :: p_Iedges, p_Itemp
    integer :: iedge,ibc,nedgecount
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    
    p_rtriangulation => rmatrix%p_rspatialDiscrTrial%p_rtriangulation
    p_rboundary => rmatrix%p_rspatialDiscrTrial%p_rboundary
    
    ! prepare a negative operator.
    rjumpStabilLocal = rjumpStabil
    
    rjumpStabilLocal%dtheta = rjumpStabilLocal%dtheta*dweight

    allocate (p_Iedges(p_rtriangulation%NMT))
    allocate (p_Itemp(p_rtriangulation%NMT))

    ! Loop over the edges and boundary components
    do ibc = 1,boundary_igetNBoundComp(p_rboundary)
      do iedge = 1,boundary_igetNsegments(p_rboundary, ibc)

        ! Get a region identifying that boundary edge
        call boundary_createRegion(p_rboundary,ibc,iedge,rboundaryRegion)

        ! Get triangulation edges on that boundary edge
        call bcasm_getEdgesInBdRegion (p_rtriangulation,rboundaryRegion, &
            nedgecount, p_Iedges, p_Itemp)

        ! Subtract the jump
        call conv_JumpStabilisation2d(rjumpStabilLocal,CONV_MODMATRIX,rmatrix, &
            InodeList=p_Iedges(1:nedgecount))

      end do
    end do

    deallocate (p_Iedges,p_Itemp)

  end subroutine   

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_addBdEOJvector (rjumpStabil,dweight,rmatrix,rx,rb)

!<description>
  ! Adds the jump stabilisation on all boundary edges.
!</description>

!<input>
    ! Jump stabilisation structure that defines the jump stabilisation operator.
    type(t_jumpStabilisation), intent(in) :: rjumpStabil
    
    ! Weight of the operator. Default = 1.0.
    real(DP), intent(in) :: dweight

    ! Solution vector.
    type(t_vectorBlock), intent(in) :: rx
    
    ! Matrix that specifies the sparsity pattern of the operator.
    type(t_matrixScalar), intent(inout) :: rmatrix
!</input>

!<inputoutput>
    ! Defect vector to be modified.
    type(t_vectorBlock), intent(inout), target :: rb
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_jumpStabilisation) :: rjumpStabilLocal
    type(t_boundaryRegion) :: rboundaryRegion
    integer, dimension(:), allocatable :: p_Iedges, p_Itemp
    integer :: iedge,ibc,nedgecount
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_boundary), pointer :: p_rboundary
    
    p_rtriangulation => rmatrix%p_rspatialDiscrTrial%p_rtriangulation
    p_rboundary => rmatrix%p_rspatialDiscrTrial%p_rboundary
    
    ! prepare a negative operator.
    rjumpStabilLocal = rjumpStabil
    
    rjumpStabilLocal%dtheta = rjumpStabilLocal%dtheta*dweight

    allocate (p_Iedges(p_rtriangulation%NMT))
    allocate (p_Itemp(p_rtriangulation%NMT))

    ! Loop over the edges and boundary components
    do ibc = 1,boundary_igetNBoundComp(p_rboundary)
      do iedge = 1,boundary_igetNsegments(p_rboundary, ibc)

        ! Get a region identifying that boundary edge
        call boundary_createRegion(p_rboundary,ibc,iedge,rboundaryRegion)

        ! Get triangulation edges on that boundary edge
        call bcasm_getEdgesInBdRegion (p_rtriangulation,rboundaryRegion, &
            nedgecount, p_Iedges, p_Itemp)

        ! Add the jump
        call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix,rx,rb,&
                              InodeList=p_Iedges(1:nedgecount))   

      end do
    end do

    deallocate (p_Iedges,p_Itemp)

  end subroutine   

end module

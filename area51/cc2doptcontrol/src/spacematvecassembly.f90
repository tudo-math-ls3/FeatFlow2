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
!#   $\theta_{ij}$           - weight for the Laplace matrix,
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
!# 2.) cc_assembleDefect
!#     -> Set up a defect vector d:=b-A(x)x
!#
!# 3.) cc_projectControlTimestep
!#     Projects a vector into a range of numbers.
!#
!# 4.) cc_initNonlinMatrix
!#     Initialises a nonlinear-matrix structure with basic parameters.
!# </purpose>
!##############################################################################

module spacematvecassembly

  use fsystem
  use storage
  use genoutput
  use boundary
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
  use scalarpde
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
  use optcontrolconvection
  use convection
  use collection
  use numbersets
  
  use basicstructures
    
  implicit none
  
!<constants>

!<constantblock description="Identifiers for the 'coperation' input parameter of the matrix assembly routine">

  ! Allocate memory if necessary.
  integer(I32), parameter :: CCMASM_ALLOCMEM              = 1
  
  ! Compute all matrix entries.
  integer(I32), parameter :: CCMASM_COMPUTE               = 2
  
  ! Allocate memory and compute matrix entries.
  integer(I32), parameter :: CCMASM_ALLOCANDCOMPUTE       = 3
  
  ! Bypass memory allocation for matrices.
  integer(I32), parameter :: CMASM_QUICKREFERENCES        = 4
  
  ! Don't clear the old matrix
  integer(i32), parameter :: CMASM_NOCLEAR                = 8
  
!</constantblock>

!<constantblock description="Identifiers for the IUPWIND parameter that specifies how to set up the nonlinearity or stabilisation.">

  ! Streamline diffusion; configured by dupsam
  integer, parameter :: CCMASM_STAB_STREAMLINEDIFF    = 0

  ! 1st-order upwind; configured by dupsam
  integer, parameter :: CCMASM_STAB_UPWIND            = 1
  
  ! Edge-oriented stabilisation; configured by dupsam as 'gamma'
  integer, parameter :: CCMASM_STAB_EDGEORIENTED      = 2

  ! Streamline diffusion; configured by dupsam, new implementation
  integer, parameter :: CCMASM_STAB_STREAMLINEDIFF2   = 3

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation
  integer, parameter :: CCMASM_STAB_EDGEORIENTED2     = 4

  ! Edge-oriented stabilisation; configured by dupsam as 'gamma', new implementation.
  ! Precomputed matrix.
  integer, parameter :: CCMASM_STAB_EDGEORIENTED3     = 5
!</constantblock>

!<constantblock description="Matrix type ID's specifying the general matrix class to set up.">

  ! Standard matrix.
  integer, parameter :: CCMASM_MTP_AUTOMATIC         = 0
  
  ! Standard matrix with decoupled velocity blocks
  integer, parameter :: CCMASM_MTP_DECOUPLED         = 1
  
  ! Extended 'full-tensor' matrix with submatrices A11, A12, A21, A22, all independent from
  ! each other.
  integer, parameter :: CCMASM_MTP_FULLTENSOR        = 2

!</constantblock>

!</constants>

!<types>

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
    
  
!    ! IOTA-parameters that switch the identity in the primal/dual equation on/off.
!    real(DP) :: diota1 = 0.0_DP
!    real(DP) :: diota2 = 0.0_DP
!
!    ! KAPPA-parameters that switch the I matrix in the continuity equation
!    ! on/off.
!    real(DP) :: dkappa1 = 0.0_DP
!    real(DP) :: dkappa2 = 0.0_DP
!
!    ! ALPHA-parameters that control the weight of the mass matrix in the
!    ! core equation. =0.0 for stationary simulations.
!    real(DP) :: dalpha1 = 0.0_DP
!    real(DP) :: dalpha2 = 0.0_DP
!
!    ! THETA-parameters that control the weight of the Stokes matrix
!    ! in the core equation. =1.0 for stationary simulations.
!    real(DP) :: dtheta1 = 0.0_DP
!    real(DP) :: dtheta2 = 0.0_DP
!
!    ! GAMMA-parameters that control the weight in front of the
!    ! nonlinearity. =1.0 for Navier-Stokes, =0.0 for Stokes equation.
!    real(DP) :: dgamma1 = 0.0_DP
!    real(DP) :: dgamma2 = 0.0_DP
!
!    ! DNEWTON-parameters that control the weight in font of the
!    ! newton matrix (adjoint of the nonlinearity). =0.0 to dectivate.
!    real(DP) :: dnewton1 = 0.0_DP
!    real(DP) :: dnewton2 = 0.0_DP
!
!    ! ETA-parameters that switch the B-terms on/off.
!    real(DP) :: deta1 = 0.0_DP
!    real(DP) :: deta2 = 0.0_DP
!
!    ! TAU-parameters that switch the B^T-terms on/off
!    real(DP) :: dtau1 = 0.0_DP
!    real(DP) :: dtau2 = 0.0_DP
!
!    ! MU-parameters that weight the coupling mass matrices.
!    real(DP) :: dmu1 = 0.0_DP
!    real(DP) :: dmu2 = 0.0_DP
!
!    ! R-parameter that weight the reactive coupling mass matrix
!    ! in the primal equation.
!    real(DP) :: dr11 = 0.0_DP
!    real(DP) :: dr12 = 0.0_DP
!
!    ! R-parameter that weight the reactive coupling mass matrix
!    ! in the dual equation.
!    real(DP) :: dr21 = 0.0_DP
!    real(DP) :: dr22 = 0.0_DP
    
    ! Regularisation parameter ALPHA which controls the transformation
    ! from the dual variable $\lambda$ to the control $u$:
    ! $u=-1/\alpha \lambda$.
    real(DP) :: dalphaC = 0.0_DP
    
    ! Type of this matrix.
    ! =0: Standard matrix.
    ! =1: Newton matrix.
    integer :: cmatrixType = 0

    ! Type of constraints to apply to the control u.
    ! =0: No constraints.
    ! =1: Constant constraints on u active: dumin1 <= u_1 <= dumax1, dumin2 <= u_2 <= dumax2.
    ! Necessary for applying the correct projection during matrix-vector
    ! multiplication.
    integer :: ccontrolConstraints = 0
    
    ! Minimum and maximum value of u in case projection of u is
    ! active (i.e. dp <> 0). Min/max values for u1.
    real(DP) :: dumin1 = 0.0_DP
    real(DP) :: dumax1 = 0.0_DP
    
    ! Minimum and maximum value of u in case projection of u is
    ! active (i.e. dp <> 0). Min/max values for u2.
    real(DP) :: dumin2 = 0.0_DP
    real(DP) :: dumax2 = 0.0_DP

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

    ! STABILISATION: Parameter that defines how to set up the nonlinearity and
    ! whether to use some kind of stabilisation. One of the CCMASM_STAB_xxxx
    ! constants. Standard is CCMASM_STAB_STREAMLINEDIFF.
    integer :: iupwind1 = CCMASM_STAB_STREAMLINEDIFF2
    integer :: iupwind2 = CCMASM_STAB_STREAMLINEDIFF2
    
    ! STABILISATION: Viscosity parameter. Used for stabilisation schemes when
    ! a nonlinearity is set up.
    real(DP) :: dnu = 0.0_DP
    
    ! STABILISATION: Stabilisation parameter for streamline diffusion, upwind and
    ! edge-oriented stabilisation. If iupwind=CCMASM_STAB_STREAMLINEDIFF, a value of
    ! 0.0 deactivates any stabilisation.
    real(DP) :: dupsam1 = 0.0_DP
    real(DP) :: dupsam2 = 0.0_DP
    
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
    logical :: bvirtualTransposedD = .false.

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...).
    ! Only used during matrix creation.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation => null()

    ! Pointer to static/precalculated information on this level.
    type(t_staticLevelInfo), pointer :: p_rstaticInfo => null()

  end type

!</typeblock>

!<typeblock>

  ! Collects all parameters that specify how to stabilise the nonlinearity
  ! when assembling it.
  type t_convecStabilisation
    ! Stabilisation type
    integer :: iupwind = CCMASM_STAB_STREAMLINEDIFF2
    
    ! Stabilisation parameter
    real(DP) :: dupsam = 0.0_DP
    
    ! Pointer to a precomputed EOJ matrix or NULL()
    type(t_matrixScalar), pointer :: p_rmatrixEOJ
  end type

!</types>

contains
  
  ! ***************************************************************************

  !<subroutine>
  
    subroutine cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
        rdiscretisation,rstaticLevelInfo)
  
  !<description>
    ! Initialises the rnonlinearCCMatrix structure with parameters and pointers
    ! from the main problem and precalculated information.
  !</description>

  !<input>
    ! Global problem structure.
    type(t_problem), intent(inout) :: rproblem
    
    ! Discretisation of the level where the matrix is to be assembled.
    type(t_blockDiscretisation), intent(in), target :: rdiscretisation
    
    ! Core equation structure of one level.
    type(t_staticLevelInfo), intent(in), target :: rstaticLevelInfo
  !</input>
  
  !<inputoutput>
    ! Nonlinear matrix structure.
    ! Basic parameters in this structure are filled with data.
    type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
  !</inputoutput>
               
  !</subroutine>
      
      ! Initialise the matrix assembly structure rnonlinearCCMatrix
      ! with basic global information.
      !
      ! 1.) Model, stabilisation
      rnonlinearSpatialMatrix%dnu = rproblem%rphysicsPrimal%dnu

      ! Get stabilisation parameters
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'iUpwind1',rnonlinearSpatialMatrix%iupwind1,0)
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'iUpwind2',rnonlinearSpatialMatrix%iupwind2,0)
      
      call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'dUpsam1',rnonlinearSpatialMatrix%dupsam1,0.0_DP)
      call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'dUpsam2',rnonlinearSpatialMatrix%dupsam2,0.0_DP)

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
      
      rnonlinearSpatialMatrix%dupsam2 = -rnonlinearSpatialMatrix%dupsam2
      
      ! 2.) Pointers to global precalculated matrices.
      rnonlinearSpatialMatrix%p_rdiscretisation => rdiscretisation
      rnonlinearSpatialMatrix%p_rstaticInfo => rstaticLevelInfo
      
    end subroutine

!  ! ***************************************************************************
!
!  !<subroutine>
!
!    subroutine cc_initMatrixStructure (rnonlinearSpatialMatrix,rmatrix)
!
!  !<description>
!    ! Analyses the matrix configuration in rnonlinearSpatialMatrix and creates
!    ! an empty matrix rmatrix with the structure that allows to assemble
!    ! all operators specified by rnonlinearSpatialMatrix.
!  !</description>
!
!  !<input>
!    ! Nonlinear matrix structure defining the matrix configuration.
!    ! Must have been initialised at least with cc_initNonlinMatrix.
!    type(t_nonlinearSpatialMatrix), intent(i) :: rnonlinearSpatialMatrix
!  !</input>
!
!  !<inputoutput>
!    ! Block matrix receiving the structure for all operators.
!    type(t_matrixBlock), intent(in) :: rmatrix
!  !</inputoutput>
!
!  !</subroutine>
!
!
!
!
!    end subroutine

! ***************************************************************************
  !<subroutine>

  subroutine coeff_ProjMass (rdiscretisationTrial,rdiscretisationTest,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
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
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, &
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
    integer(i32), dimension(:), allocatable :: p_Idofs
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

  subroutine smva_assembleMatrix (coperation,cmatrixType,&
      rmatrix,rnonlinearSpatialMatrix,&
      rvector1,rvector2,rvector3,rfineMatrix)

!<description>
  ! This routine assembles a global matrix. The caller must initialise the
  ! rnonlinearSpatialMatrix according to how the matrix should look like.
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

  ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
  ! about how to set up the matrix.
  !
  ! Note that if coperation=CCMASM_ALLOCxxxx is specified, p_rmatrixTemplateXXXX
  ! must be initialised as well as p_rdiscretisation!
  ! The new matrix is created based p_rmatrixTemplateXXXX as well as p_rdiscretisation.
  ! Memory is automatically allocated if it's missing.
  type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'previous' timestep. If there is no previous timestep (e.g.
  ! like in the 0th timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), target, optional :: rvector1

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'current' timestep.
  type(t_vectorBlock), intent(IN), target, optional :: rvector2

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'next' timestep. If there is no next timestep (e.g.
  ! like in the last timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), target, optional :: rvector3

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
  type(t_matrixBlock), intent(INOUT) :: rmatrix
  
!</inputoutput>
  
!</subroutine>

    ! local variables
    logical :: ballocate
    type(t_matrixBlock) :: rtempMatrix
    type(t_vectorBlock) :: rtempVector
    type(t_convecStabilisation) :: rstabilisation
    integer :: i,j
    type(t_optcoperator) :: roptcoperator
    type(t_vectorBlock), pointer :: p_rprimalSol, p_rdualSol
    type(t_blockDiscretisation) :: rvelDiscr
    
    logical, parameter :: bnewmethod = .false.
    
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
      call lsysbl_createMatBlockByDiscr(rnonlinearSpatialMatrix%p_rdiscretisation,rmatrix)
      do j=1,2
        do i=1,2
          call lsysbl_createEmptyMatrix (rtempMatrix,NDIM2D+1,NDIM2D+1)
          call allocSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rtempMatrix,&
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
        call spdiscr_deriveBlockDiscr (rnonlinearSpatialMatrix%p_rdiscretisation, rvelDiscr, &
                                       1, 3)
        call lsysbl_createMatBlockByDiscr (rvelDiscr,rtempMatrix)
        
        ! Primal equation
        ! ---------------
        ! In the first step, assemble the linear parts
        ! (Laplace, Mass, B/B^T) on the main diagonal of the primal equation.
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,1,3)
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,&
            rtempMatrix,rnonlinearSpatialMatrix%Diota(1,1),rnonlinearSpatialMatrix%Dalpha(1,1),&
            rnonlinearSpatialMatrix%Dtheta(1,1),&
            rnonlinearSpatialMatrix%Deta(1,1),rnonlinearSpatialMatrix%Dtau(1,1),.false.)


        ! Assemble the nonlinearity u*grad(.) or the Newton nonlinearity
        ! u*grad(.)+grad(u)*(.) to the velocity.
        rstabilisation = t_convecStabilisation(&
            rnonlinearSpatialMatrix%iupwind1,rnonlinearSpatialMatrix%dupsam1,&
            rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixEOJ1)
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,rvector1,&
              rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
              rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
              rstabilisation)
        case (2)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,rvector2,&
              rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
              rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
              rstabilisation)
        case (3)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,rvector3,&
              rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
              rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
              rstabilisation)
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
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, rvector1, &
              rnonlinearSpatialMatrix%Dalpha(1,2))
        case (2)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, rvector2, &
              rnonlinearSpatialMatrix%Dalpha(1,2))
        case (3)
          call assembleProjectedMassBlocks (rnonlinearSpatialMatrix,rtempMatrix, rvector3, &
              rnonlinearSpatialMatrix%Dalpha(1,2))
        end select

        ! Reintegrate the computed matrix
        call lsysbl_moveToSubmatrix (rtempMatrix,rmatrix,1,4)

        ! Dual equation
        ! -------------
        ! In the first step, assemble the linear parts
        ! (Laplace, Mass, B/B^T) on the main diagonal of the primal equation.
        call lsysbl_extractSubmatrix (rmatrix,rtempMatrix,4,6)
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rtempMatrix,&
            rnonlinearSpatialMatrix%Diota(2,2),rnonlinearSpatialMatrix%Dalpha(2,2),&
            rnonlinearSpatialMatrix%Dtheta(2,2),&
            rnonlinearSpatialMatrix%Deta(2,2),rnonlinearSpatialMatrix%Dtau(2,2),.true.)

        ! Assemble the nonlinearity u*grad(.) or the Newton nonlinearity
        ! u*grad(.)+grad(u)*(.) to the velocity.
        rstabilisation = t_convecStabilisation(&
            rnonlinearSpatialMatrix%iupwind2,rnonlinearSpatialMatrix%dupsam2,&
            rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixEOJ2)
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,rvector1,&
              rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
              rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
              rstabilisation)
        case (2)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,rvector2,&
              rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
              rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
              rstabilisation)
        case (3)
          call assembleConvection (&
              rnonlinearSpatialMatrix,rtempMatrix,rvector3,&
              rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
              rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
              rstabilisation)
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
        call assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rtempMatrix,&
            rnonlinearSpatialMatrix%Diota(2,1),rnonlinearSpatialMatrix%Dalpha(2,1),&
            rnonlinearSpatialMatrix%Dtheta(2,1),&
            rnonlinearSpatialMatrix%Deta(2,1),rnonlinearSpatialMatrix%Dtau(2,1),.false.)

        ! Co stabilisation in the convective parts here.
        ! rstabilisation = t_convecStabilisation(&
        !    rnonlinearSpatialMatrix%iupwind2,rnonlinearSpatialMatrix%dupsam2)
        rstabilisation = t_convecStabilisation(0,0.0_DP,NULL())

        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          call lsysbl_deriveSubvector(rvector1,rtempVector, 4,6,.true.)
        case (2)
          call lsysbl_deriveSubvector(rvector2,rtempVector, 4,6,.true.)
        case (3)
          call lsysbl_deriveSubvector(rvector3,rtempVector, 4,6,.true.)
        end select
            
        call assembleConvection (&
            rnonlinearSpatialMatrix,rtempMatrix,rtempVector,&
            rnonlinearSpatialMatrix%Dgamma(2,1),rnonlinearSpatialMatrix%DgammaT(2,1),&
            rnonlinearSpatialMatrix%Dnewton(2,1),rnonlinearSpatialMatrix%DnewtonT(2,1),&
            rstabilisation)
            
        ! There is probably a 2nd reactive term stemming from the next time step.
        ! Assemble it.
        
        select case (rnonlinearSpatialMatrix%idualSol2)
        case (1)
          call lsysbl_deriveSubvector(rvector1,rtempVector, 4,6,.true.)
        case (2)
          call lsysbl_deriveSubvector(rvector2,rtempVector, 4,6,.true.)
        case (3)
          call lsysbl_deriveSubvector(rvector3,rtempVector, 4,6,.true.)
        end select
            
        call assembleConvection (&
            rnonlinearSpatialMatrix,rtempMatrix,rtempVector,&
            0.0_DP,rnonlinearSpatialMatrix%DgammaT2(2,1),&
            rnonlinearSpatialMatrix%Dnewton2(2,1),0.0_DP,&
            rstabilisation)

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
              rnonlinearSpatialMatrix%iadaptiveMatrices, &
              rnonlinearSpatialMatrix%dadmatthreshold)
              
          if (.not. lsyssc_isMatrixContentShared(&
              rfineMatrix%RmatrixBlock(1,1),rfineMatrix%RmatrixBlock(2,2))) then
            call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(2,2), &
                rmatrix%RmatrixBlock(1,1), &
                rnonlinearSpatialMatrix%iadaptiveMatrices, &
                rnonlinearSpatialMatrix%dadmatthreshold)
          end if
            
          ! Dual system
          call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(4,4), &
              rmatrix%RmatrixBlock(4,4), &
              rnonlinearSpatialMatrix%iadaptiveMatrices, &
              rnonlinearSpatialMatrix%dadmatthreshold)
              
          if (.not. lsyssc_isMatrixContentShared(&
            rmatrix%RmatrixBlock(4,4),rmatrix%RmatrixBlock(5,5))) then
            call mrest_matrixRestrictionEX3Y (rfineMatrix%RmatrixBlock(5,5), &
                rmatrix%RmatrixBlock(5,5), &
                rnonlinearSpatialMatrix%iadaptiveMatrices, &
                rnonlinearSpatialMatrix%dadmatthreshold)
          end if
          
        end if

        ! Release memory
        call lsysbl_releaseMatrix(rtempMatrix)
        if (rtempVector%NEQ .ne. 0) &
          call lsysbl_releaseVector (rtempVector)
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
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
                                        rmatrix%RmatrixBlock(1,3),&
                                        LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)

          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
                                        rmatrix%RmatrixBlock(2,3),&
                                        LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        end if
        
        if (rnonlinearSpatialMatrix%Dtau(1,1) .ne. 0.0_DP) then
          ! Furthermore, put B1^T and B2^T to the block matrix.
          ! These matrices are always 'shared'.
          if (rnonlinearSpatialMatrix%bvirtualTransposedD) then
            call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
                                          rmatrix%RmatrixBlock(3,1),&
                                          LSYSSC_TR_VIRTUAL)
            call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
                                          rmatrix%RmatrixBlock(3,2),&
                                          LSYSSC_TR_VIRTUAL)
          else
            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
                                          rmatrix%RmatrixBlock(3,1),&
                                          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
                                          rmatrix%RmatrixBlock(3,2),&
                                          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          end if
        end if
                                        
        if (rnonlinearSpatialMatrix%Dtau(2,2) .ne. 0.0_DP) then
          if (rnonlinearSpatialMatrix%bvirtualTransposedD) then
            call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
                                          rmatrix%RmatrixBlock(6,4),&
                                          LSYSSC_TR_VIRTUAL)
            call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
                                          rmatrix%RmatrixBlock(6,5),&
                                          LSYSSC_TR_VIRTUAL)
          else
            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
                                          rmatrix%RmatrixBlock(6,4),&
                                          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
                                          rmatrix%RmatrixBlock(6,5),&
                                          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
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
        roptcoperator%dupsamPrimal = rnonlinearSpatialMatrix%dupsam1
        roptcoperator%dupsamDual = rnonlinearSpatialMatrix%dupsam2
        
        ! Timestep-weights
        roptcoperator%dprimalAlpha = rnonlinearSpatialMatrix%Dalpha(1,1)
        roptcoperator%ddualAlpha   = rnonlinearSpatialMatrix%Dalpha(2,2)

        ! Stokes operator
        roptcoperator%dnu = rnonlinearSpatialMatrix%dnu
        roptcoperator%dprimalBeta = rnonlinearSpatialMatrix%Dtheta(1,1)
        roptcoperator%ddualBeta   = rnonlinearSpatialMatrix%Dtheta(2,2)
        
        ! Nonlinearity
        if (rnonlinearSpatialMatrix%Dgamma(1,1) .ne. 0.0_DP) then
          roptcoperator%dprimalDelta = rnonlinearSpatialMatrix%Dgamma(1,1)
          roptcoperator%ddualDelta   = rnonlinearSpatialMatrix%Dgamma(2,2)
          roptcoperator%ddualNewtonTrans = rnonlinearSpatialMatrix%DnewtonT(2,2)
          
          ! Whether or not Newton is active has no influence to the
          ! defect, so the following lines are commented out.
          ! if (rparams%bnewton) then
          roptcoperator%dprimalNewton    = rnonlinearSpatialMatrix%Dnewton(1,1)
          roptcoperator%ddualRDeltaTrans = rnonlinearSpatialMatrix%DgammaT(2,1)
          roptcoperator%ddualRNewton     = rnonlinearSpatialMatrix%Dnewton(2,1)
          ! end if
          
        end if
        
        ! Coupling matrices
        !if (rparams%bdualcoupledtoprimal) then
          roptcoperator%ddualRAlpha = rnonlinearSpatialMatrix%Dalpha(2,1)
        !end if

        !if (rparams%bcontrolactive) then
          roptcoperator%dcontrolWeight = -rnonlinearSpatialMatrix%Dalpha(1,2)*rnonlinearSpatialMatrix%dalphaC
          roptcoperator%dcontrolMultiplier = -1.0_DP/rnonlinearSpatialMatrix%dalphaC
        !end if
        
        if (rnonlinearSpatialMatrix%ccontrolConstraints .ne. 0) then
          roptcoperator%ccontrolProjection = rnonlinearSpatialMatrix%ccontrolConstraints
          roptcoperator%dmin1 = rnonlinearSpatialMatrix%dumin1
          roptcoperator%dmax1 = rnonlinearSpatialMatrix%dumax1
          roptcoperator%dmin2 = rnonlinearSpatialMatrix%dumin2
          roptcoperator%dmax2 = rnonlinearSpatialMatrix%dumax2
        end if
        
        select case (rnonlinearSpatialMatrix%iprimalSol)
        case (1)
          p_rprimalSol => rvector1
        case (2)
          p_rprimalSol => rvector2
        case (3)
          p_rprimalSol => rvector3
        end select
        
        select case (rnonlinearSpatialMatrix%idualSol)
        case (1)
          p_rdualSol => rvector1
        case (2)
          p_rdualSol => rvector2
        case (3)
          p_rdualSol => rvector3
        end select
        
        ! Calculate the velocity-dependent part of the system matrix.
        call conv_strdiffOptC2dgetMatrix (rmatrix,roptcoperator,1.0_DP,&
            p_rprimalSol,p_rdualSol)
            
        !call matio_writeBlockMatrixHR (rmatrix, 'matrix',&
        !    .true., 0, 'matrix2.txt', '(E12.5)', 1E-10_DP)
      
      end if
      
    end if
    
  contains
  
    ! -----------------------------------------------------
    
    subroutine allocSubmatrix (rnonlinearSpatialMatrix,cmatrixType,rsubmatrix,&
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
      ! All matrices may have multiplication factors in their front.
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
      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixTemplateFEM,&
                  rsubmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
          
      if (.not. bfulltensor) then
           
        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix
        ! A22 is identical to A11! So mirror A11 to A22 sharing the
        ! structure and the content.
        call lsyssc_duplicateMatrix (rsubmatrix%RmatrixBlock(1,1),&
                    rsubmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                    
      else
      
        ! Otherwise, create another copy of the template matrix.
        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixTemplateFEM,&
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
            
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixTemplateFEM, &
            rsubmatrix%RmatrixBlock(1,2), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
          !     rsubmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
            
        end if

        if (rsubmatrix%RmatrixBlock(2,1)%cmatrixFormat &
            .eq. LSYSSC_MATRIXUNDEFINED) then
            
          ! Create a new matrix A21 in memory. create a new matrix
          ! using the template FEM matrix...
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixTemplateFEM, &
            rsubmatrix%RmatrixBlock(2,1), &
            LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
            
          ! Allocate memory for the entries; don't initialise the memory.
          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
          ! zero.
          ! CALL lsyssc_allocEmptyMatrix (&
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
        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
                                      rsubmatrix%RmatrixBlock(1,3),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
                                      rsubmatrix%RmatrixBlock(2,3),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      end if
        
      ! In the same manner, insert an identiy matrix for the pressure
      ! to the system matrix. The matrix is zero by default but may
      ! partially or fully be initialised to an identity matrix depending
      ! on the situation. (It's mostly used for direct solvers/UMFPACK)
      if (dkappa .ne. 0.0_DP) then
        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixTemplateFEMPressure, &
                                      rsubmatrix%RmatrixBlock(3,3),&
                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      end if

      ! Now, prepare B1^T and B2^T. Are these matrices given?
      !
      ! If yes, take the given matrices. If not, create them by
      ! 'virtually transposing' B1 and B2 (i.e. they share the same
      ! data as B1 and B2 but hate the 'transpose'-flag set).

      if (dtau .ne. 0.0_DP) then
        if (rnonlinearSpatialMatrix%bvirtualTransposedD) then
          call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
                                        rsubmatrix%RmatrixBlock(3,1),&
                                        LSYSSC_TR_VIRTUAL)

          call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
                                        rsubmatrix%RmatrixBlock(3,2),&
                                        LSYSSC_TR_VIRTUAL)
        else
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
                                        rsubmatrix%RmatrixBlock(3,1),&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
                                        rsubmatrix%RmatrixBlock(3,2),&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end if
      end if

    end subroutine

    ! -----------------------------------------------------
    
    subroutine assembleLinearSubmatrix (rnonlinearSpatialMatrix,cmatrixType,&
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
              rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass  ,dalpha,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
              
          if (.not. bshared) then

            call lsyssc_matrixLinearComb (&
                rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass     ,dalpha,&
                rmatrix%RmatrixBlock(2,2),1.0_DP,&
                rmatrix%RmatrixBlock(2,2),&
                .false.,.false.,.true.,.true.)
          end if
          
        end if
        
        ! ---------------------------------------------------
        ! Plug in the Stokes matrix?
        if (dtheta .ne. 0.0_DP) then
          call lsyssc_matrixLinearComb (&
              rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes     ,dtheta,&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),&
              .false.,.false.,.true.,.true.)
              
          if (.not. bshared) then
            call lsyssc_matrixLinearComb (&
                rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes   ,dtheta,&
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
        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
                                      rmatrix%RmatrixBlock(1,3),&
                                      LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)

        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
                                      rmatrix%RmatrixBlock(2,3),&
                                      LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
      end if
      
      if (dtau .ne. 0.0_DP) then
        ! Furthermore, put B1^T and B2^T to the block matrix.
        ! These matrices are always 'shared'.
        if (rnonlinearSpatialMatrix%bvirtualTransposedD) then
          call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
                                        rmatrix%RmatrixBlock(3,1),&
                                        LSYSSC_TR_VIRTUAL)
          call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
                                        rmatrix%RmatrixBlock(3,2),&
                                        LSYSSC_TR_VIRTUAL)
        else
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
                                        rmatrix%RmatrixBlock(3,1),&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
                                        rmatrix%RmatrixBlock(3,2),&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
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
        rstabilisation)
        
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
    
    ! Stabilisation parameters
    type(t_convecStabilisation), intent(in) :: rstabilisation
    
    ! local variables
    logical :: bshared
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_jumpStabilisation) :: rjumpStabil
    real(dp), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
    
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
                      
        select case (rstabilisation%iupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion%ddelta            = dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dgammaT
          rstreamlineDiffusion%dnewton           = dnewton
          rstreamlineDiffusion%dnewtonTransposed = dnewtonT
          
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
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%dnu
          
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dgamma
          rstreamlineDiffusion2%ddeltaT = dgammaT
          rstreamlineDiffusion2%dnewton = dnewton
          rstreamlineDiffusion2%dnewtonT = dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)
          
        case (CCMASM_STAB_UPWIND)
        
          call output_line ('Upwind not supported.', &
                            OU_CLASS_ERROR,OU_MODE_STD,'assembleConvection')
          call sys_halt()

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta            = dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dgammaT
          rstreamlineDiffusion%dnewton           = dnewton
          rstreamlineDiffusion%dnewtonTransposed = dnewtonT
          
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
          rjumpStabil%dgamma = rstabilisation%dupsam
          
          ! Matrix weight
          rjumpStabil%dtheta = dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))
          end if

        case (CCMASM_STAB_EDGEORIENTED2)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dgamma
          rstreamlineDiffusion2%ddeltaT  = dgammaT
          rstreamlineDiffusion2%dnewton  = dnewton
          rstreamlineDiffusion2%dnewtonT = dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)
                              
          ! Set up the jump stabilisation structure.
          ! There's not much to do, only initialise the viscosity...
          rjumpStabil%dnu = rstreamlineDiffusion%dnu
          
          ! Set stabilisation parameter
          rjumpStabil%dgamma = rstabilisation%dupsam
          
          ! Matrix weight
          rjumpStabil%dtheta = dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODMATRIX, &
                              rmatrix%RmatrixBlock(1,1))

          if (.not. bshared) then
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(2,2))
          end if

        case (CCMASM_STAB_EDGEORIENTED3)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%dnu
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dgamma
          rstreamlineDiffusion2%ddeltaT  = dgammaT
          rstreamlineDiffusion2%dnewton  = dnewton
          rstreamlineDiffusion2%dnewtonT = dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1), p_Ddata1)
          call lsyssc_getbase_double (rstabilisation%p_rmatrixEOJ, p_Ddata2)
          call conv_streamDiff2Blk2dMat (rstreamlineDiffusion2,rmatrix,rvector)
                              
          ! We use the precomputed EOJ matrix and sum it up to the
          ! existing matrix.
          call lsyssc_matrixLinearComb(rstabilisation%p_rmatrixEOJ,&
              dgamma*mprim_signum(rstabilisation%dupsam),&
              rmatrix%RmatrixBlock(1,1),1.0_DP,&
              rmatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

          if (.not. bshared) then
            ! Also for the Y-velocity.
            call lsyssc_matrixLinearComb(rstabilisation%p_rmatrixEOJ,&
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
          select case (rstabilisation%iupwind)
          case (CCMASM_STAB_EDGEORIENTED,CCMASM_STAB_EDGEORIENTED2)
            
            ! Set up the jump stabilisation structure.
            ! There's not much to do, only initialise the viscosity...
            rjumpStabil%dnu = rnonlinearSpatialMatrix%dnu
            
            ! Set stabilisation parameter
            rjumpStabil%dgamma = rstabilisation%dupsam
            
            ! Matrix weight
            rjumpStabil%dtheta = dgamma

            ! Call the jump stabilisation technique to stabilise that stuff.
            ! We can assemble the jump part any time as it's independent of any
            ! convective parts...
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODMATRIX, &
                                rmatrix%RmatrixBlock(1,1))

            if (.not. bshared) then
              call conv_jumpStabilisation2d (&
                                  rjumpStabil, CONV_MODMATRIX, &
                                  rmatrix%RmatrixBlock(2,2))
            end if
            
          case (CCMASM_STAB_EDGEORIENTED3)
            
            ! We use the precomputed EOJ matrix and sum it up to the
            ! existing matrix.
            call lsyssc_matrixLinearComb(rstabilisation%p_rmatrixEOJ,&
                dgamma*mprim_signum(rstabilisation%dupsam),&
                rmatrix%RmatrixBlock(1,1),1.0_DP,&
                rmatrix%RmatrixBlock(1,1),.false.,.false.,.true.,.true.)

            if (.not. bshared) then
              ! Also for the Y-velocity.
              call lsyssc_matrixLinearComb(rstabilisation%p_rmatrixEOJ,&
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

      ! Assemble A14/A25?
      if (dweight .ne. 0.0_DP) then
          
        ! Switch on these matrices
        rmatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP
        rmatrix%RmatrixBlock(2,2)%dscaleFactor = 1.0_DP
          
        ! Calculate the usual mass matrix if conrol constraints are deactivated
        ! or if Newton is not active.
        if ((rnonlinearSpatialMatrix%ccontrolConstraints .eq. 0) .or. &
            (rnonlinearSpatialMatrix%cmatrixType .eq. 0)) then
      
          ! Copy the entries of the mass matrix. Share the structure.
          ! We must not share the entries as these might be changed by the caller
          ! e.g. due to boundary conditions!
          
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
              
          ! Scale the entries by dweight.
          if (dweight .ne. 1.0_DP) then
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),dweight)
            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),dweight)
          end if
          
        else if (rnonlinearSpatialMatrix%cmatrixType .eq. 1) then
          
          select case (rnonlinearSpatialMatrix%ccontrolConstraints)
          
          case (1)
          
            ! Copy the entries of the mass matrix. Share the structure.
            ! We must not share the entries as these might be changed by the caller
            ! e.g. due to boundary conditions!
            
            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
                rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
          
            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
                rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
                
            ! Scale the entries by the weight.
            if (dweight .ne. 1.0_DP) then
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),dweight)
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),dweight)
            end if

            ! Filter the matrix. All the rows corresponding to DOF's that violate
            ! the bounds must be set to zero.
            call massmatfilter (rmatrix%RmatrixBlock(1,1),rvector%RvectorBlock(4),&
                rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%dumin1,rnonlinearSpatialMatrix%dumax1)
            call massmatfilter (rmatrix%RmatrixBlock(2,2),rvector%RvectorBlock(5),&
                rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%dumin2,rnonlinearSpatialMatrix%dumax2)
                
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
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin1
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax1
            
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
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin2
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax2

            call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(2,2),&
                coeff_ProjMass,rcollection)
            
            ! Now we can forget about the collection again.
            call collct_done (rcollection)
            
          case (3)
          
            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
                rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
                rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

            ! Create the matrix
            call approxProjectionDerivative (rmatrix%RmatrixBlock(1,1), &
                rvector%RvectorBlock(4), rnonlinearSpatialMatrix%dalphaC,&
                rnonlinearSpatialMatrix%dumin1,rnonlinearSpatialMatrix%dumax1,0.001_DP)

            call approxProjectionDerivative (rmatrix%RmatrixBlock(2,2), &
                rvector%RvectorBlock(5), rnonlinearSpatialMatrix%dalphaC,&
                rnonlinearSpatialMatrix%dumin2,rnonlinearSpatialMatrix%dumax2,0.001_DP)
          
            ! Scale the entries by the weight if necessary
            if (dweight .ne. 1.0_DP) then
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,1),dweight)
              call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,2),dweight)
            end if
            
          case (4)
          
            ! Exact reassembly of the mass matrices with adaptive integration.
          
            ! Create an array that saves all elements on the border of the active set.
            nelements = rnonlinearSpatialMatrix%p_rdiscretisation%RspatialDiscr(1)%&
                p_rtriangulation%NEL
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
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin1
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax1
            
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
            rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin2
            rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax2

            ! Create a new element list
            rcollection%IquickAccess(3) = 0

            call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(2,2),&
                coeff_ProjMassCollect,rcollection)
            
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

!    ! -----------------------------------------------------
!
!    subroutine allocMatrix (cmatrixType,rnonlinearSpatialMatrix,rmatrix)
!
!    ! Allocates memory for the system matrix. rnonlinearSpatialMatrix provides information
!    ! about the submatrices that are 'plugged into' rmatrix.
!    ! Therefore, before this routine is called, rnonlinearSpatialMatrix must have been set up.
!
!    ! Type of matrix that should be set up in rmatrix. One of the CCMASM_MTP_xxxx
!    ! constants.
!    integer, intent(IN) :: cmatrixType
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how to set up the matrix.
!    type(t_nonlinearSpatialMatrix), intent(IN), target :: rnonlinearSpatialMatrix
!
!    ! A block matrix that receives the basic system matrix.
!    type(t_matrixBlock), intent(INOUT) :: rmatrix
!
!      ! local variables
!      logical :: bdecoupled, bfulltensor
!
!      ! A pointer to the system matrix and the RHS/solution vectors.
!      type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient
!
!      ! A pointer to the discretisation structure with the data.
!      type(t_blockDiscretisation), pointer :: p_rdiscretisation
!
!      ! Ask the problem structure to give us the discretisation structure
!      p_rdiscretisation => rnonlinearSpatialMatrix%p_rdiscretisation
!
!      ! Get a pointer to the template FEM matrix. If that doesn't exist,
!      ! take the Stokes matrix as template.
!      p_rmatrixTemplateFEM => rnonlinearSpatialMatrix%p_rmatrixTemplateFEM
!      if (.not. associated(p_rmatrixTemplateFEM)) &
!        p_rmatrixTemplateFEM => rnonlinearSpatialMatrix%p_rmatrixStokes
!      if (.not. associated(p_rmatrixTemplateFEM)) then
!        print *,'allocMatrix: Cannot set up A matrices in system matrix!'
!        stop
!      end if
!
!      ! In the global system, there are two gradient matrices B1 and B2.
!      ! Get a pointer to the template structure for these.
!      ! If there is no pointer, try to get use a pointer to one of these
!      ! matrices directly.
!      p_rmatrixTemplateGradient => rnonlinearSpatialMatrix%p_rmatrixTemplateGradient
!      if (.not. associated(p_rmatrixTemplateGradient)) &
!        p_rmatrixTemplateGradient => rnonlinearSpatialMatrix%p_rmatrixB1
!      if (.not. associated(p_rmatrixTemplateGradient)) then
!        print *,'allocMatrix: Cannot set up B matrix in system matrix!'
!        stop
!      end if
!
!      ! Initialise the block matrix with default values based on
!      ! the discretisation.
!      if (associated(p_rdiscretisation)) then
!        call lsysbl_createMatBlockByDiscr (p_rdiscretisation,rmatrix)
!      else
!        ! No discretisation structure; create the matrix directly as 3x3 matrix.
!        call lsysbl_createEmptyMatrix (rmatrix,NDIM2D+1)
!      end if
!
!      ! -------------------------------------------------------------
!      ! Primal equation
!      ! -------------------------------------------------------------
!
!      ! Determine the shape of the matrix
!      bdecoupled = cmatrixType .eq. CCMASM_MTP_DECOUPLED
!      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
!
!      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
!        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
!        bfulltensor = rnonlinearSpatialMatrix%dnewton1 .ne. 0.0_DP
!      end if
!
!      ! Let's consider the global system in detail. The standard matrix It has
!      ! roughly the following shape:
!      !
!      !    ( A11  .    B1  M   .   .   ) = ( A11  .    A13 A14 .   .   )
!      !    ( .    A22  B2  .   M   .   )   ( .    A22  A23 .   A25 .   )
!      !    ( B1^T B2^T .   .   .   .   )   ( A31  A32  .   .   .   .   )
!      !
!      ! All matrices may have multiplication factors in their front.
!      !
!      ! The structure of the matrices A11 and A22 of the global system matrix
!      ! is governed by the template FEM matrix.
!      ! Initialise them with the same structure, i.e. A11, A22 share (!) their
!      ! structure (not the entries) with that of the template matrix.
!      !
!      ! For this purpose, use the "duplicate matrix" routine.
!      ! The structure of the matrix is shared with the template FEM matrix.
!      ! For the content, a new empty array is allocated which will later receive
!      ! the entries.
!      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                  rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      if (.not. bfulltensor) then
!
!        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix
!        ! A22 is identical to A11! So mirror A11 to A22 sharing the
!        ! structure and the content.
!        call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(1,1),&
!                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!      else
!
!        ! Otherwise, create another copy of the template matrix.
!        call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                    rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      end if
!
!      ! Manually change the discretisation structure of the Y-velocity
!      ! matrix to the Y-discretisation structure.
!      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
!      ! so this is not really necessary - we do this for sure...
!      call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(2,2),&
!          p_rdiscretisation%RspatialDiscr(2))
!
!      ! A 'full tensor matrix' consists also of blocks A12 and A21.
!      if (bfulltensor) then
!
!        ! We have a matrix in the following shape:
!        !
!        !    ( A11  A12  B1  R1  .   .   ) = ( A11  A12  A13 A14 .   .   )
!        !    ( A21  A22  B2  .   R1  .   )   ( A21  A22  A23 .   A25 .   )
!        !    ( B1^T B2^T .   .   .   .   )   ( A31  A32  .   .   .   .   )
!        !
!        ! Create A12 and A21.
!
!        if (rmatrix%RmatrixBlock(1,2)%cmatrixFormat &
!            .eq. LSYSSC_MATRIXUNDEFINED) then
!
!          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
!            rmatrix%RmatrixBlock(1,2), &
!            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!          ! Allocate memory for the entries; don't initialise the memory.
!          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
!          ! zero.
!          ! CALL lsyssc_allocEmptyMatrix (&
!          !     rmatrix%RmatrixBlock(1,2),LSYSSC_SETM_UNDEFINED)
!
!        end if
!
!        if (rmatrix%RmatrixBlock(2,1)%cmatrixFormat &
!            .eq. LSYSSC_MATRIXUNDEFINED) then
!
!          ! Create a new matrix A21 in memory. create a new matrix
!          ! using the template FEM matrix...
!          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
!            rmatrix%RmatrixBlock(2,1), &
!            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!          ! Allocate memory for the entries; don't initialise the memory.
!          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
!          ! zero.
!          ! CALL lsyssc_allocEmptyMatrix (&
!          !     p_rmatrixPreconditioner%RmatrixBlock(2,1),LSYSSC_SETM_UNDEFINED)
!
!        end if
!
!      end if
!
!      ! The B1/B2 matrices exist up to now only in our local problem structure.
!      ! Put a copy of them into the block matrix.
!      !
!      ! Note that we share the structure of B1/B2 with those B1/B2 of the
!      ! block matrix, while we create empty space for the entries.
!      ! Later, the B-matrices are copied into here and modified for boundary
!      ! conditions.
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                    rmatrix%RmatrixBlock(1,3),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                    rmatrix%RmatrixBlock(2,3),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      ! In the same manner, insert an identiy matrix for the pressure
!      ! to the system matrix. The matrix is zero by default but may
!      ! partially or fully be initialised to an identity matrix depending
!      ! on the situation. (It's mostly used for direct solvers/UMFPACK)
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixIdentityPressure, &
!                                    rmatrix%RmatrixBlock(3,3),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!      call lsyssc_clearMatrix (rmatrix%RmatrixBlock(3,3))
!
!      ! Now, prepare B1^T and B2^T. Are these matrices given?
!      !
!      ! If yes, take the given matrices. If not, create them by
!      ! 'virtually transposing' B1 and B2 (i.e. they share the same
!      ! data as B1 and B2 but hate the 'transpose'-flag set).
!
!      if (associated(rnonlinearSpatialMatrix%p_rmatrixB1T)) then
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1T, &
!                                      rmatrix%RmatrixBlock(3,1),&
!                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      else
!        call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                      rmatrix%RmatrixBlock(3,1),&
!                                      LSYSSC_TR_VIRTUAL)
!      end if
!
!      if (associated(rnonlinearSpatialMatrix%p_rmatrixB2T)) then
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2T, &
!                                      rmatrix%RmatrixBlock(3,2),&
!                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      else
!        call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                      rmatrix%RmatrixBlock(3,2),&
!                                      LSYSSC_TR_VIRTUAL)
!      end if
!
!      ! Insert free space for the mass matrices to the matrix.
!      !
!      ! Note that we share the structure of M, while we create empty space
!      ! for the entries.
!      ! Later, the M-matrices are copied into here and modified for boundary
!      ! conditions.
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass, &
!                                    rmatrix%RmatrixBlock(1,4),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass, &
!                                    rmatrix%RmatrixBlock(2,5),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      ! -------------------------------------------------------------
!      ! Dual equation
!      ! -------------------------------------------------------------
!
!      ! Determine the shape of the matrix
!      bdecoupled = cmatrixType .eq. CCMASM_MTP_DECOUPLED
!      bfulltensor = cmatrixType .eq. CCMASM_MTP_FULLTENSOR
!
!      if (cmatrixType .eq. CCMASM_MTP_AUTOMATIC) then
!        ! Should we assemble Newton? If yes, we have a full-tensor matrix.
!        bfulltensor = rnonlinearSpatialMatrix%dnewton2 .ne. 0.0_DP
!      end if
!
!      ! Let's consider the global system in detail. The standard matrix It has
!      ! roughly the following shape:
!      !
!      !    ( R2   .    .   A44  .    B1  ) = ( A41  .    .   A44 .   A46 )
!      !    ( .    R2   .   .    A55  B2  )   ( .    A52  .   .   A55 A56 )
!      !    ( .    .    .   B1^T B2^T .   )   ( .    .    .   A64 A65 .   )
!      !
!      ! All matrices may have multiplication factors in their front.
!      !
!      ! The structure of the matrices A44 and A55 of the global system matrix
!      ! is governed by the template FEM matrix.
!      ! Initialise them with the same structure, i.e. A44, A55 share (!) their
!      ! structure (not the entries) with that of the template matrix.
!      !
!      ! For this purpose, use the "duplicate matrix" routine.
!      ! The structure of the matrix is shared with the template FEM matrix.
!      ! For the content, a new empty array is allocated which will later receive
!      ! the entries.
!      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                  rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      if (.not. bdecoupled .and. .not. bfulltensor) then
!
!        ! If X- and Y-velocity is to be treated in a 'coupled' way, the matrix
!        ! A22 is identical to A44! So mirror A44 to A55 sharing the
!        ! structure and the content.
!        call lsyssc_duplicateMatrix (rmatrix%RmatrixBlock(4,4),&
!                    rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!      else
!
!        ! Otherwise, create another copy of the template matrix.
!        call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
!                    rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      end if
!
!      ! Manually change the discretisation structure of the Y-velocity
!      ! matrix to the Y-discretisation structure.
!      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
!      ! so this is not really necessary - we do this for sure...
!      call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(5,5),&
!          p_rdiscretisation%RspatialDiscr(5))
!
!      ! A 'full tensor matrix' consists also of blocks A12 and A21.
!      if (bfulltensor) then
!
!        ! We have a submatrix in the following shape:
!        !
!        !    ( R2   .    .   A44  A45  B1  ) = ( A41  .    .   A44 A45 A46 )
!        !    ( .    R2   .   A54  A55  B2  )   ( .    A52  .   A54 A55 A56 )
!        !    ( .    .    .   B1^T B2^T .   )   ( .    .    .   A64 A65 .   )
!        !
!        ! Create A45 and A54.
!
!        if (rmatrix%RmatrixBlock(4,5)%cmatrixFormat &
!            .eq. LSYSSC_MATRIXUNDEFINED) then
!
!          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
!            rmatrix%RmatrixBlock(4,5), &
!            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!          ! Allocate memory for the entries; don't initialise the memory.
!          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
!          ! zero.
!          ! CALL lsyssc_allocEmptyMatrix (&
!          !     rmatrix%RmatrixBlock(4,5),LSYSSC_SETM_UNDEFINED)
!
!        end if
!
!        if (rmatrix%RmatrixBlock(5,4)%cmatrixFormat &
!            .eq. LSYSSC_MATRIXUNDEFINED) then
!
!          ! Create a new matrix A54 in memory. create a new matrix
!          ! using the template FEM matrix...
!          call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM, &
!            rmatrix%RmatrixBlock(5,4), &
!            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!          ! Allocate memory for the entries; don't initialise the memory.
!          ! Probably possible, but up to now, LSYSSC_DUP_EMPTY above initialises with
!          ! zero.
!          ! CALL lsyssc_allocEmptyMatrix (&
!          !     p_rmatrixPreconditioner%RmatrixBlock(5,4),LSYSSC_SETM_UNDEFINED)
!
!        end if
!
!      end if
!
!      ! The B1/B2 matrices exist up to now only in our local problem structure.
!      ! Put a copy of them into the block matrix.
!      !
!      ! Note that we share the structure of B1/B2 with those B1/B2 of the
!      ! block matrix, while we create empty space for the entries.
!      ! Later, the B-matrices are copied into here and modified for boundary
!      ! conditions.
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                    rmatrix%RmatrixBlock(4,6),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                    rmatrix%RmatrixBlock(5,6),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      ! Now, prepare B1^T and B2^T. Are these matrices given?
!      !
!      ! If yes, take the given matrices. If not, create them by
!      ! 'virtually transposing' B1 and B2 (i.e. they share the same
!      ! data as B1 and B2 but hate the 'transpose'-flag set).
!
!      if (associated(rnonlinearSpatialMatrix%p_rmatrixB1T)) then
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1T, &
!                                      rmatrix%RmatrixBlock(6,4),&
!                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      else
!        call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                      rmatrix%RmatrixBlock(6,4),&
!                                      LSYSSC_TR_VIRTUAL)
!      end if
!
!      if (associated(rnonlinearSpatialMatrix%p_rmatrixB2T)) then
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2T, &
!                                      rmatrix%RmatrixBlock(6,5),&
!                                      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      else
!        call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                      rmatrix%RmatrixBlock(6,5),&
!                                      LSYSSC_TR_VIRTUAL)
!      end if
!
!      ! In the same manner, insert an identiy matrix for the pressure
!      ! to the system matrix; as the enties aren't changed, we can share
!      ! the entries with the original one.
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixIdentityPressure, &
!                                    rmatrix%RmatrixBlock(6,6),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!      call lsyssc_clearMatrix (rmatrix%RmatrixBlock(6,6))
!
!      ! Insert free space for the mass matrices to the matrix.
!      !
!      ! Note that we share the structure of M, while we create empty space
!      ! for the entries.
!      ! Later, the M-matrices are copied into here and modified for boundary
!      ! conditions.
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass, &
!                                    rmatrix%RmatrixBlock(4,1),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass, &
!                                    rmatrix%RmatrixBlock(5,2),&
!                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!      ! That's it, all submatrices are basically set up.
!      !
!      ! Update the structural information of the block matrix, as we manually
!      ! changed the submatrices:
!      call lsysbl_updateMatStrucInfo (rmatrix)
!
!    end subroutine
!
!    ! -----------------------------------------------------
!
!    subroutine assembleVelocityBlocks (bdualEquation,&
!        rnonlinearSpatialMatrix,rmatrix,rvector,dvectorWeight)
!
!    ! Assembles the velocity matrix in the block matrix rmatrix at position (1,1):
!    !
!    ! rmatrix := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
!    !            dnewton*N*(p_rvector)
!
!    ! Whether to set up the primal or the dual equation.
!    ! FALSE=primal, TRUE=dual equation.
!    logical, intent(IN) :: bdualEquation
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how to set up the matrix.
!    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
!
!    ! Block matrix where the 2x2-velocity submatrix should be assembled
!    type(t_matrixBlock), intent(INOUT) :: rmatrix
!
!    ! Velocity vector for the nonlinearity. Must be specified if
!    ! GAMMA <> 0; can be omitted if GAMMA=0.
!    type(t_vectorBlock), optional :: rvector
!
!    ! Weight for the velocity vector; standard = 1.
!    real(DP), intent(IN), optional :: dvectorWeight
!
!    ! local variables
!    logical :: bshared
!    integer :: iupwind
!    real(DP) :: dupsam
!    type(t_convUpwind) :: rupwind
!    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
!    type(t_jumpStabilisation) :: rjumpStabil
!    type(t_matrixBlock) :: rtempMatrix
!    real(DP) :: dvecWeight
!    integer :: imatOffset
!    real(DP) :: diota, dalpha, dtheta, dnewton, dgamma
!
!      if (.not. bdualEquation) then
!        ! Set the weights used here according to the primal equation.
!        ! Set imatOffset=0 so the submatrix at position 1,1 is tackled.
!        imatOffset = 0
!        diota = rnonlinearSpatialMatrix%diota1
!        dalpha = rnonlinearSpatialMatrix%dalpha1
!        dtheta = rnonlinearSpatialMatrix%dtheta1
!        dgamma = rnonlinearSpatialMatrix%dgamma1
!        dnewton = rnonlinearSpatialMatrix%dnewton1
!        iupwind = rnonlinearSpatialMatrix%iupwind1
!        dupsam = rnonlinearSpatialMatrix%dupsam1
!      else
!        ! Set the weights used here according to the primal equation.
!        ! Set imatOffset=3 so the submatrix at position 4,4 is tackled.
!        imatOffset = 3
!        diota = rnonlinearSpatialMatrix%diota2
!        dalpha = rnonlinearSpatialMatrix%dalpha2
!        dtheta = rnonlinearSpatialMatrix%dtheta2
!        dgamma = rnonlinearSpatialMatrix%dgamma2
!        dnewton = rnonlinearSpatialMatrix%dnewton2
!        iupwind = rnonlinearSpatialMatrix%iupwind2
!        dupsam = rnonlinearSpatialMatrix%dupsam2
!      end if
!
!      ! Standard value for dvectorWeight is = -1.
!      dvecWeight = -1.0_DP
!      if (present(dvectorWeight)) dvecWeight = dvectorWeight
!
!      ! Is A11=A22 physically?
!      bshared = lsyssc_isMatrixContentShared(&
!                    rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),&
!                    rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
!
!      ! Allocate memory if necessary. Normally this should not be necessary...
!      ! A11:
!      if (.not. lsyssc_hasMatrixContent (&
!          rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))) then
!        call lsyssc_allocEmptyMatrix (&
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),LSYSSC_SETM_UNDEFINED)
!      end if
!
!      ! A22:
!      if (.not. bshared) then
!        if (.not. lsyssc_hasMatrixContent (&
!            rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))) then
!          call lsyssc_allocEmptyMatrix (&
!              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),LSYSSC_SETM_UNDEFINED)
!        end if
!      end if
!
!      ! A12/ A21:
!      if (lsysbl_isSubmatrixPresent (rmatrix,imatOffset+1,imatOffset+2)) then
!        if (.not. lsyssc_hasMatrixContent (&
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))) then
!          call lsyssc_allocEmptyMatrix (&
!              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2),LSYSSC_SETM_UNDEFINED)
!        end if
!        if (.not. lsyssc_hasMatrixContent (&
!            rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))) then
!          call lsyssc_allocEmptyMatrix (&
!              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1),LSYSSC_SETM_UNDEFINED)
!        end if
!      end if
!
!      ! ---------------------------------------------------
!      ! If diota <> 0, initialise the matrix with the identity,
!      ! otherwise with zero.
!      if (diota .eq. 0.0_DP) then
!        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))
!
!        if (.not. bshared) then
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
!        end if
!
!      else
!
!        call lsyssc_initialiseIdentityMatrix (&
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))
!        call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),diota)
!
!        if (.not. bshared) then
!          call lsyssc_initialiseIdentityMatrix (&
!              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
!          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),diota)
!        end if
!
!      end if
!
!      ! ---------------------------------------------------
!      ! Plug in the mass matrix?
!      if (dalpha .ne. 0.0_DP) then
!
!        call lsyssc_matrixLinearComb (&
!            rnonlinearSpatialMatrix%p_rmatrixMass  ,dalpha,&
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),1.0_DP,&
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),&
!            .false.,.false.,.true.,.true.)
!
!        if (.not. bshared) then
!
!          call lsyssc_matrixLinearComb (&
!              rnonlinearSpatialMatrix%p_rmatrixMass     ,dalpha,&
!              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),1.0_DP,&
!              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),&
!              .false.,.false.,.true.,.true.)
!        end if
!
!      end if
!
!      ! ---------------------------------------------------
!      ! Plug in the Stokes matrix?
!      if (dtheta .ne. 0.0_DP) then
!        call lsyssc_matrixLinearComb (&
!            rnonlinearSpatialMatrix%p_rmatrixStokes     ,dtheta,&
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),1.0_DP,&
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1),&
!            .false.,.false.,.true.,.true.)
!
!        if (.not. bshared) then
!          call lsyssc_matrixLinearComb (&
!              rnonlinearSpatialMatrix%p_rmatrixStokes   ,dtheta,&
!              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),1.0_DP,&
!              rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2),&
!              .false.,.false.,.true.,.true.)
!        end if
!      end if
!
!      ! ---------------------------------------------------
!      ! That was easy -- the adventure begins now... The nonlinearity!
!      if (dgamma .ne. 0.0_DP) then
!
!        if (.not. present(rvector)) then
!          call output_line ('Velocity vector not present!', &
!                             OU_CLASS_ERROR,OU_MODE_STD,'smva_assembleMatrix')
!          call sys_halt()
!        end if
!
!        select case (iupwind)
!        case (CCMASM_STAB_STREAMLINEDIFF)
!          ! Set up the SD structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter
!          rstreamlineDiffusion%dupsam = dupsam
!
!          ! Matrix weight
!          rstreamlineDiffusion%ddelta = dgamma
!
!          ! Weight for the Newton part; =0 deactivates Newton.
!          rstreamlineDiffusion%dnewton = dnewton
!
!          if (dnewton .eq. 0.0_DP) then
!
!            ! If the submatrices A12 and A21 exist, fill them with zero.
!            ! If they don't exist, we don't have to do anything.
!            if (lsysbl_isSubmatrixPresent (rmatrix,imatOffset+1,imatOffset+2)) then
!              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
!              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
!            end if
!
!          else
!
!            ! Clear A12/A21 that may receive parts of the Newton matrix
!            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
!            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
!
!          end if
!
!          ! Create a temporary block matrix only contining the velocity submatrices
!          ! we want to change. Share structure and entries such that changing
!          ! the temporary matrix will also change the original matrix.
!          call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
!                                       LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,&
!                                       imatOffset+1,imatOffset+2)
!
!          ! Call the SD method to calculate the nonlinearity.
!          ! As velocity vector, specify rvector!
!          ! Therefore, the primal velcity is always used for assembling
!          ! that thing!
!          call conv_streamlineDiffusionBlk2d (&
!                              rvector, rvector, &
!                              dvecWeight, 0.0_DP,&
!                              rstreamlineDiffusion, CONV_MODMATRIX, &
!                              rtempMatrix)
!
!          ! Release the temp matrix.
!          call lsysbl_releaseMatrix (rtempMatrix)
!
!        case (CCMASM_STAB_UPWIND)
!          ! Set up the upwind structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rupwind%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter
!          rupwind%dupsam = dupsam
!
!          ! Matrix weight
!          rupwind%dtheta = dgamma
!
!          if (dnewton .ne. 0.0_DP) then
!            call output_line ('Warning: Upwind does not support assembly '&
!                //'of the Newton matrix!',OU_CLASS_TRACE1)
!          end if
!
!          ! Call the upwind method to calculate the nonlinear matrix.
!          ! As velocity vector, specify rvector!
!          ! Therefore, the primal velcity is always used for assembling
!          ! that thing!
!          call conv_upwind2d (rvector, rvector, &
!                              dvecWeight, 0.0_DP,&
!                              rupwind, CONV_MODMATRIX, &
!                              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))
!
!          if (.not. bshared) then
!            ! Modify also the matrix block (2,2)
!            call conv_upwind2d (rvector, rvector, &
!                                dvecWeight, 0.0_DP,&
!                                rupwind, CONV_MODMATRIX, &
!                                rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
!          end if
!
!        case (CCMASM_STAB_EDGEORIENTED)
!          ! Jump stabilisation.
!          ! In the first step, set up the matrix as above with central discretisation,
!          ! i.e. call SD to calculate the matrix without SD stabilisation.
!          ! Set up the SD structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
!          rstreamlineDiffusion%dupsam = 0.0_DP
!
!          ! Matrix weight
!          rstreamlineDiffusion%ddelta = dgamma
!
!          ! Weight for the Newtop part; =0 deactivates Newton.
!          rstreamlineDiffusion%dnewton = dnewton
!
!          if (dnewton .eq. 0.0_DP) then
!
!            ! If the submatrices A12 and A21 exist, fill them with zero.
!            ! If they don't exist, we don't have to do anything.
!            if (lsysbl_isSubmatrixPresent (rmatrix,imatOffset+1,imatOffset+2)) then
!              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
!              call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
!            end if
!
!          else
!
!            ! Clear A12/A21 that receives parts of the Newton matrix
!            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2))
!            call lsyssc_clearMatrix (rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1))
!
!            ! Activate the submatrices A12 and A21 if they aren't.
!            rmatrix%RmatrixBlock(imatOffset+1,imatOffset+2)%dscaleFactor = 1.0_DP
!            rmatrix%RmatrixBlock(imatOffset+2,imatOffset+1)%dscaleFactor = 1.0_DP
!
!          end if
!
!          ! Create a temporary block matrix only contining the velocity submatrices
!          ! we want to change. Share structure and entries such that changing
!          ! the temporary matrix will also change the original matrix.
!          call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
!                                       LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,&
!                                       imatOffset+1,imatOffset+2)
!
!          ! Call the SD method to calculate the nonlinearity.
!          ! As velocity vector, specify rvector!
!          ! Therefore, the primal velcity is always used for assembling
!          ! that thing!
!          call conv_streamlineDiffusionBlk2d (&
!                              rvector, rvector, &
!                              dvecWeight, 0.0_DP,&
!                              rstreamlineDiffusion, CONV_MODMATRIX, &
!                              rtempMatrix)
!
!          ! Release the temp matrix.
!          call lsysbl_releaseMatrix (rtempMatrix)
!
!          ! Set up the jump stabilisation structure.
!          ! There's not much to do, only initialise the viscosity...
!          rjumpStabil%dnu = rstreamlineDiffusion%dnu
!
!          ! Set stabilisation parameter
!          rjumpStabil%dgammastar = dupsam
!          rjumpStabil%dgamma = rjumpStabil%dgammastar
!
!          ! Matrix weight
!          rjumpStabil%dtheta = dgamma
!
!          ! Call the jump stabilisation technique to stabilise that stuff.
!          ! We can assemble the jump part any time as it's independent of any
!          ! convective parts...
!          call conv_jumpStabilisation2d (&
!                              rjumpStabil, CONV_MODMATRIX, &
!                              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))
!
!          if (.not. bshared) then
!            call conv_jumpStabilisation2d (&
!                                rjumpStabil, CONV_MODMATRIX, &
!                                rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
!          end if
!
!        case DEFAULT
!          print *,'Don''t know how to set up nonlinearity!?!'
!          stop
!
!        end select
!
!      else
!
!        ! That's the Stokes-case. Jump stabilisation is possible...
!
!        select case (iupwind)
!        case (CCMASM_STAB_EDGEORIENTED)
!
!          ! Set up the jump stabilisation structure.
!          ! There's not much to do, only initialise the viscosity...
!          rjumpStabil%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter
!          rjumpStabil%dgammastar = dupsam
!          rjumpStabil%dgamma = rjumpStabil%dgammastar
!
!          ! Matrix weight
!          rjumpStabil%dtheta = dgamma
!
!          ! Call the jump stabilisation technique to stabilise that stuff.
!          ! We can assemble the jump part any time as it's independent of any
!          ! convective parts...
!          call conv_jumpStabilisation2d (&
!                              rjumpStabil, CONV_MODMATRIX, &
!                              rmatrix%RmatrixBlock(imatOffset+1,imatOffset+1))
!
!          if (.not. bshared) then
!            call conv_jumpStabilisation2d (&
!                                rjumpStabil, CONV_MODMATRIX, &
!                                rmatrix%RmatrixBlock(imatOffset+2,imatOffset+2))
!          end if
!
!        case DEFAULT
!          ! No stabilisation
!
!        end select
!
!      end if ! gamma <> 0
!
!    end subroutine
!
!    ! -----------------------------------------------------
!
!    subroutine assembleMassBlocks (rnonlinearSpatialMatrix,rmatrix, rvector)
!
!    ! Assembles a 2x2 block matrix with mass matrices and
!    ! probably nonlinear submatrices on the diagonal.
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how to set up the matrix.
!    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
!
!    ! Block matrix where the 2x2-velocity submatrix should be assembled
!    type(t_matrixBlock), intent(INOUT) :: rmatrix
!
!    ! Vector that specifies where to evaluate nonlinear terms
!    type(t_vectorBlock), intent(IN), target :: rvector
!
!      ! local variables
!      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
!      type(t_matrixBlock) :: rtempmatrix
!      type(t_vectorBlock) :: rtempvector
!      type(t_bilinearForm) :: rform
!      type(t_collection) :: rcollection
!
!      ! Assemble A14/A25?
!      if (rnonlinearSpatialMatrix%dmu1 .ne. 0.0_DP) then
!
!        ! Calculate the usual mass matrix if conrol constraints are deactivated
!        ! or if Newton is not active.
!        if ((rnonlinearSpatialMatrix%ccontrolConstraints .eq. 0) .or. &
!            (rnonlinearSpatialMatrix%cmatrixType .eq. 0)) then
!
!          ! Copy the entries of the mass matrix. Share the structure.
!          ! We must not share the entries as these might be changed by the caller
!          ! e.g. due to boundary conditions!
!
!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!              rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!
!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!              rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!
!          ! Scale the entries by dmu2.
!          if (rnonlinearSpatialMatrix%dmu1 .ne. 1.0_DP) then
!            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,4),rnonlinearSpatialMatrix%dmu1)
!            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,5),rnonlinearSpatialMatrix%dmu1)
!          end if
!
!        else if ((rnonlinearSpatialMatrix%ccontrolConstraints .eq. 1) .and. &
!                 (rnonlinearSpatialMatrix%cmatrixType .eq. 1)) then
!
!          ! In A14/A25 we have to create a 'projective mass matrix'.
!          ! This is the derivative of a projection operator
!          ! P[a,b](f)=a if f<a, =b if f>b, =f otherwise.
!          ! For a<f<b, this is the mass matrix. Everywhere else, this is =0.
!          ! We assemble this matrix just as a standard mass matrix with noconstant
!          ! coefficients. Whereever u = -1/alpha * lambda is out of bounds,
!          ! we return 0 as coefficient, otherwise 1.
!
!          rform%itermCount = 1
!          rform%Idescriptors(1,1) = DER_FUNC
!          rform%Idescriptors(2,1) = DER_FUNC
!
!          ! In this case, we have nonconstant coefficients.
!          rform%ballCoeffConstant = .FALSE.
!          rform%BconstantCoeff(:) = .FALSE.
!
!          ! Prepare a collection structure to be passed to the callback
!          ! routine. We attach the vector T in the quick-access variables
!          ! so the callback routine can access it.
!          ! The bounds and the alpha value are passed in the
!          ! quickaccess-arrays.
!          call collct_init(rcollection)
!          rcollection%p_rvectorQuickAccess1 => rvector
!
!          ! Coefficient is dmu1=1/alpha or 0, depending on lambda
!          rcollection%DquickAccess(3)  = rnonlinearSpatialMatrix%dalphaC
!          rcollection%DquickAccess(4)  = rnonlinearSpatialMatrix%dmu1
!
!          ! At first, set up A14, depending on lambda_1.
!          rcollection%IquickAccess(1) = 1
!          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin1
!          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax1
!
!          ! Now we can build the matrix entries.
!          ! We specify the callback function coeff_Laplace for the coefficients.
!          ! As long as we use constant coefficients, this routine is not used.
!          ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!          ! the framework will call the callback routine to get analytical
!          ! data.
!          ! The collection is passed as additional parameter. That's the way
!          ! how we get the vector to the callback routine.
!          call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(1,4),&
!              coeff_ProjMass,rcollection)
!
!          ! Now, set up A25, depending on lambda_2.
!          rcollection%IquickAccess(1)  = 2
!          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin2
!          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax2
!
!          call bilf_buildMatrixScalar (rform,.TRUE.,rmatrix%RmatrixBlock(2,5),&
!              coeff_ProjMass,rcollection)
!
!          ! Now we can forget about the collection again.
!          call collct_done (rcollection)
!
!!          ! Copy the entries of the mass matrix. Share the structure.
!!          ! We must not share the entries as these might be changed by the caller
!!          ! e.g. due to boundary conditions!
!!
!!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!!              rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!!
!!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!!              rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!!
!!          ! Scale the entries by dmu1.
!!          if (rnonlinearSpatialMatrix%dmu1 .ne. 1.0_DP) then
!!            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(1,4),rnonlinearSpatialMatrix%dmu1)
!!            call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(2,5),rnonlinearSpatialMatrix%dmu1)
!!          end if
!!
!!          ! Filter the matrix. All the rows corresponding to DOF's that violate
!!          ! the bounds must be set to zero.
!!          call massmatfilter (rmatrix%RmatrixBlock(1,4),rvector%RvectorBlock(4),&
!!              rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%dumin1,rnonlinearSpatialMatrix%dumax1)
!!          call massmatfilter (rmatrix%RmatrixBlock(2,5),rvector%RvectorBlock(5),&
!!              rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%dumin2,rnonlinearSpatialMatrix%dumax2)
!        end if
!
!        rmatrix%RmatrixBlock(1,4)%dscaleFactor = 1.0_DP
!        rmatrix%RmatrixBlock(2,5)%dscaleFactor = 1.0_DP
!
!        if (rnonlinearSpatialMatrix%dr12 .ne. 0.0_DP) then
!          ! There is some data in A24/A15, so create empty space there
!          ! in case it's missing.
!          if (.not. lsysbl_isSubmatrixPresent(rmatrix,2,4)) then
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(1,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!          end if
!
!          ! Clear the offdiagonal matrices, switch them on
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,5))
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,4))
!
!          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 1.0_DP
!          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 1.0_DP
!
!        else
!
!          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 0.0_DP
!          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 0.0_DP
!
!        end if
!
!      else if ((rnonlinearSpatialMatrix%dr11 .ne. 0.0_DP) .or. &
!               (rnonlinearSpatialMatrix%dr12 .ne. 0.0_DP)) then
!
!        ! There is some data in A14/A25, so create empty space there
!        ! in case it's missing.
!        if (.not. lsysbl_isSubmatrixPresent(rmatrix,1,4)) then
!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!              rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!              rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!        end if
!
!        ! Clear the diagonal matrices, switch them on
!        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,4))
!        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,5))
!
!        rmatrix%RmatrixBlock(1,4)%dscaleFactor = 1.0_DP
!        rmatrix%RmatrixBlock(2,5)%dscaleFactor = 1.0_DP
!
!        if (rnonlinearSpatialMatrix%dr12 .ne. 0.0_DP) then
!          ! There is some data in A42/A51, so create empty space there
!          ! in case it's missing.
!          if (.not. lsysbl_isSubmatrixPresent(rmatrix,2,4)) then
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(2,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(1,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!          end if
!
!          ! Clear the offdiagonal matrices, switch them on
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(1,5))
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(2,4))
!
!          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 1.0_DP
!          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 1.0_DP
!
!        else
!
!          rmatrix%RmatrixBlock(1,5)%dscaleFactor = 0.0_DP
!          rmatrix%RmatrixBlock(2,4)%dscaleFactor = 0.0_DP
!
!        end if
!
!      else
!
!        ! Deactivate this submatrix
!        rmatrix%RmatrixBlock(1,4)%dscaleFactor = 0.0_DP
!        rmatrix%RmatrixBlock(1,5)%dscaleFactor = 0.0_DP
!        rmatrix%RmatrixBlock(2,4)%dscaleFactor = 0.0_DP
!        rmatrix%RmatrixBlock(2,5)%dscaleFactor = 0.0_DP
!
!      end if
!
!      ! If we have a reactive coupling mass matrix, it gets interesting...
!      if ((rnonlinearSpatialMatrix%dr11 .ne. 0.0_DP) .or. &
!          (rnonlinearSpatialMatrix%dr12 .ne. 0.0_DP)) then
!
!        ! The reactive part is: "dr2 * . * grad(lambda)".
!        ! This is exactly the 'Newton' part assembled by the streamline diffusion
!        ! method if we use the dual velocity as velocity field!
!        ! So prepare to call streamline diffusion.
!
!        ! Viscosity; ok, actually not used.
!        rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!        ! Set stabilisation parameter
!        rstreamlineDiffusion%dupsam = rnonlinearSpatialMatrix%dupsam2
!
!        ! Weight dr21 of the convective part.
!        rstreamlineDiffusion%ddelta = rnonlinearSpatialMatrix%dr11
!
!        ! Weight for the Newton part; here, this is the dr22 weight.
!        rstreamlineDiffusion%dnewton = rnonlinearSpatialMatrix%dr12
!
!        ! Create a temporary block matrix only contining the velocity submatrices
!        ! we want to change. Share structure and entries such that changing
!        ! the temporary matrix will also change the original matrix.
!        call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
!                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,1,2,4,5)
!
!        ! Create a temporary block vector that points to the dual velocity.
!        call lsysbl_deriveSubvector (rvector,rtempvector,4,5,.true.)
!
!        ! Call the SD method to calculate the nonlinearity.
!        ! As velocity vector, specify rtempvector, which points to
!        ! the dual velocity.
!        call conv_streamlineDiffusionBlk2d (&
!                            rtempvector, rtempvector, &
!                            1.0_DP, 0.0_DP,&
!                            rstreamlineDiffusion, CONV_MODMATRIX, &
!                            rtempMatrix)
!
!        ! Release the temp matrix and temp vector.
!        call lsysbl_releaseVector (rtempVector)
!        call lsysbl_releaseMatrix (rtempMatrix)
!
!      end if
!
!    end subroutine
!
!    ! -----------------------------------------------------
!
!    subroutine assembleDualMassBlocks (rnonlinearSpatialMatrix,rmatrix,rvector)
!
!    ! Assembles a 2x2 block matrix with mass matrices on the diagonal
!    ! in the dual equation. These matrices consist of a standard mass
!    ! matrix and/or a reactive coupling mass matrix depending on the
!    ! dual velocity.
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how to set up the matrix.
!    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
!
!    ! Block matrix where the 2x2-velocity submatrix should be assembled
!    type(t_matrixBlock), intent(INOUT) :: rmatrix
!
!    ! Vector that specifies where to evaluate nonlinear terms
!    type(t_vectorBlock), intent(IN) :: rvector
!
!      ! local variables
!      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
!      type(t_matrixBlock) :: rtempmatrix
!      type(t_vectorBlock) :: rtempvector
!
!      if (rnonlinearSpatialMatrix%dmu2 .ne. 0.0_DP) then
!
!        ! Copy the entries of the mass matrix. Share the structure.
!        ! We must not share the entries as these might be changed by the caller
!        ! e.g. due to boundary conditions!
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
!
!        ! Scale the entries by dmu2.
!        if (rnonlinearSpatialMatrix%dmu2 .ne. 1.0_DP) then
!          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(4,1),rnonlinearSpatialMatrix%dmu2)
!          call lsyssc_scaleMatrix (rmatrix%RmatrixBlock(5,2),rnonlinearSpatialMatrix%dmu2)
!        end if
!
!        rmatrix%RmatrixBlock(4,1)%dscaleFactor = 1.0_DP
!        rmatrix%RmatrixBlock(5,2)%dscaleFactor = 1.0_DP
!
!        if (rnonlinearSpatialMatrix%dr22 .ne. 0.0_DP) then
!          ! There is some data in A42/A51, so create empty space there
!          ! in case it's missing.
!          if (.not. lsysbl_isSubmatrixPresent(rmatrix,4,2)) then
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(4,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(5,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!          end if
!
!          ! Clear the offdiagonal matrices, switch them on
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,1))
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,2))
!
!          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 1.0_DP
!          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 1.0_DP
!
!        else
!
!          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 0.0_DP
!          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 0.0_DP
!
!        end if
!
!      else if ((rnonlinearSpatialMatrix%dr21 .ne. 0.0_DP) .or. &
!               (rnonlinearSpatialMatrix%dr22 .ne. 0.0_DP)) then
!
!        ! There is some data in A41/A52, so create empty space there
!        ! in case it's missing.
!        if (.not. lsysbl_isSubmatrixPresent(rmatrix,4,1)) then
!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!              rmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!              rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!        end if
!
!        ! Clear the diagonal matrices, switch them on
!        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,1))
!        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,2))
!
!        rmatrix%RmatrixBlock(4,1)%dscaleFactor = 1.0_DP
!        rmatrix%RmatrixBlock(5,2)%dscaleFactor = 1.0_DP
!
!        if (rnonlinearSpatialMatrix%dr22 .ne. 0.0_DP) then
!          ! There is some data in A42/A51, so create empty space there
!          ! in case it's missing.
!          if (.not. lsysbl_isSubmatrixPresent(rmatrix,4,2)) then
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(4,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!            call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!                rmatrix%RmatrixBlock(5,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!          end if
!
!          ! Clear the offdiagonal matrices, switch them on
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(5,1))
!          call lsyssc_clearMatrix (rmatrix%RmatrixBlock(4,2))
!
!          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 1.0_DP
!          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 1.0_DP
!
!        else
!
!          rmatrix%RmatrixBlock(5,1)%dscaleFactor = 0.0_DP
!          rmatrix%RmatrixBlock(4,2)%dscaleFactor = 0.0_DP
!
!        end if
!
!      else
!
!        ! Deactivate this submatrix
!        rmatrix%RmatrixBlock(4,1)%dscaleFactor = 0.0_DP
!        rmatrix%RmatrixBlock(5,1)%dscaleFactor = 0.0_DP
!        rmatrix%RmatrixBlock(4,2)%dscaleFactor = 0.0_DP
!        rmatrix%RmatrixBlock(5,2)%dscaleFactor = 0.0_DP
!
!      end if
!
!      ! If we have a reactive coupling mass matrix, it gets interesting...
!      if ((rnonlinearSpatialMatrix%dr21 .ne. 0.0_DP) .or. &
!          (rnonlinearSpatialMatrix%dr22 .ne. 0.0_DP)) then
!
!        ! The reactive part is: "dr2 * . * grad(lambda)".
!        ! This is exactly the 'Newton' part assembled by the streamline diffusion
!        ! method if we use the dual velocity as velocity field!
!        ! So prepare to call streamline diffusion.
!
!        ! Viscosity; ok, actually not used.
!        rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!        ! Set stabilisation parameter
!        rstreamlineDiffusion%dupsam = rnonlinearSpatialMatrix%dupsam2
!
!        ! Weight dr21 of the convective part.
!        rstreamlineDiffusion%ddelta = rnonlinearSpatialMatrix%dr21
!
!        ! Weight for the Newton part; here, this is the dr22 weight.
!        rstreamlineDiffusion%dnewton = rnonlinearSpatialMatrix%dr22
!
!        ! Create a temporary block matrix only contining the velocity submatrices
!        ! we want to change. Share structure and entries such that changing
!        ! the temporary matrix will also change the original matrix.
!        call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
!                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,4,5,1,2)
!
!        ! Create a temporary block vector that points to the dual velocity.
!        call lsysbl_deriveSubvector (rvector,rtempvector,4,5,.true.)
!
!        ! Call the SD method to calculate the nonlinearity.
!        ! As velocity vector, specify rtempvector, which points to
!        ! the dual velocity.
!        call conv_streamlineDiffusionBlk2d (&
!                            rtempvector, rtempvector, &
!                            1.0_DP, 0.0_DP,&
!                            rstreamlineDiffusion, CONV_MODMATRIX, &
!                            rtempMatrix)
!
!        ! Release the temp matrix and temp vector.
!        call lsysbl_releaseVector (rtempVector)
!        call lsysbl_releaseMatrix (rtempMatrix)
!
!      end if
!
!    end subroutine
!
!    ! -----------------------------------------------------
!
!    subroutine assembleGradientMatrices (bdualEquation,&
!        rnonlinearSpatialMatrix,rmatrix,bsharedMatrix)
!
!    ! Initialises the gradient/divergence matrices with entries from
!    ! the rnonlinearSpatialMatrix structure.
!    !
!    ! The routine copies references from the submatrices tormatrix,
!    ! but it does not initialise any matrix weights / scaling factors.
!    !
!    ! If bsharedMatrix=TRUE, the matrix is created using references to the
!    ! matrix building blocks in rlevelInfo, thus sharing all information
!    ! with those matrices in rnonlinearSpatialMatrix. In this case, the caller must
!    ! not change the matrix entries, because this would change the
!    ! original 'template' matrices!
!    ! (This can be used e.g. for setting up a matrix for building a defect
!    !  vector without copying matrix data.)
!    ! If bsharedMatrix=TRUE on the other hand, the matrix entries of the
!    ! original template (B-) matrices are copied in memory,
!    ! so the new matrix is allowed to be changed!
!
!    ! Whether to set up the proimal or the dual equation.
!    ! FALSE=primal, TRUE=dual equation.
!    logical, intent(IN) :: bdualEquation
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how to set up the matrix.
!    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
!
!    ! Block matrix where the B-matrices should be set up
!    type(t_matrixBlock), intent(INOUT) :: rmatrix
!
!    ! Whether or not the matrix entries of the source gradient-matrices
!    ! should be copied in memory.
!    !
!    ! If set to TRUE, the routine tries to initialise rmatrix
!    ! only with references to the original matrices, thus the caller must not
!    ! change the entries. Nevertheless, if rmatrix is the owner of one of the
!    ! submatrices, the routine will always copy the matrix entries,
!    ! as otherwise memory would have to be deallocated!
!    !
!    ! If set to FALSE, the entries of the source matrices in rnonlinearSpatialMatrix are
!    ! copied, so the caller can change rmatrix afterwards (e.g. to implement
!    ! boundary conditions).
!    logical, intent(IN) :: bsharedMatrix
!
!      ! local variables
!      integer :: idubStructure,idubContent,imatOffset
!
!      if (.not. bdualEquation) then
!        ! Set imatOffset=0 so the submatrix at position 1,1 is tackled.
!        imatOffset = 0
!      else
!        ! Set the weights used here according to the primal equation.
!        imatOffset = 3
!      end if
!
!      ! Initialise a copy flag that tells the duplicateMatrix-routine whether to
!      ! copy the entries or to create references.
!      if (bsharedMatrix) then
!
!        idubContent = LSYSSC_DUP_SHARE
!
!        ! Normally we share entries -- except for if the submatrices belong to
!        ! rmatrix! To avoid memory deallocation in this case, we copy
!        ! the entries.
!        if ((.not. lsyssc_isMatrixContentShared(&
!                Rmatrix%RmatrixBlock(imatOffset+1,imatOffset+3))) .or.&
!            (.not. lsyssc_isMatrixContentShared(&
!                Rmatrix%RmatrixBlock(imatOffset+2,imatOffset+3)))) then
!          idubContent = LSYSSC_DUP_COPY
!        end if
!
!      else
!
!        idubContent = LSYSSC_DUP_COPY
!
!      end if
!
!      idubStructure = LSYSSC_DUP_SHARE
!
!      ! Let's consider the global system in detail:
!      !
!      !    ( A11  A12  B1  ) = ( A11  A12  A13 )
!      !    ( A21  A22  B2  )   ( A21  A22  A23 )
!      !    ( B1^T B2^T 0   )   ( A31  A32  A33 )
!      !
!      ! We exclude the velocity submatrices here, so our system looks like:
!      !
!      !    (           B1 ) = (           A13 )
!      !    (           B2 )   (           A23 )
!      !    ( B1^T B2^T    )   ( A31  A32      )
!
!      ! The B1/B2 matrices exist up to now only in rnonlinearSpatialMatrix.
!      ! Put a copy of them into the block matrix.
!      !
!      ! Note that we share the structure of B1/B2 with those B1/B2 of the
!      ! block matrix, while we create copies of the entries. The B-blocks
!      ! are already prepared and memory for the entries is already allocated;
!      ! so we only have to copy the entries.
!      !
!      ! Note that idubContent = LSYSSC_DUP_COPY will automatically allocate
!      ! memory if necessary.
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                    rmatrix%RmatrixBlock(imatOffset+1,imatOffset+3),&
!                                    idubStructure,idubContent)
!
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                    rmatrix%RmatrixBlock(imatOffset+2,imatOffset+3),&
!                                    idubStructure,idubContent)
!
!      ! Furthermore, put B1^T and B2^T to the block matrix.
!      ! These matrices are always 'shared'.
!      call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                    rmatrix%RmatrixBlock(imatOffset+3,imatOffset+1),&
!                                    LSYSSC_TR_VIRTUAL)
!
!      call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                    rmatrix%RmatrixBlock(imatOffset+3,imatOffset+2),&
!                                    LSYSSC_TR_VIRTUAL)
!
!    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_assembleDefect (rnonlinearSpatialMatrix,rx,rd,cx,rvector1,rvector2,rvector3)

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

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'previous' timestep. If there is no previous timestep (e.g.
  ! like in the 0th timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), target, optional :: rvector1

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'current' timestep.
  type(t_vectorBlock), intent(IN), target, optional :: rvector2

  ! OPTIONAL: If a nonlinearity is to be set up, this vector must be specified.
  ! It specifies where to evaluate the nonlinearity and must contain the data
  ! for the 'next' timestep. If there is no next timestep (e.g.
  ! like in the last timestep), the vector can be undefined.
  type(t_vectorBlock), intent(IN), target, optional :: rvector3

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
    type(t_convecStabilisation) :: rstabilisation
    type(t_optcoperator) :: roptcoperator
    type(t_blockDIscretisation) :: rvelDiscr
    
    logical, parameter :: bnewmethod = .false.
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dd
    
    call lsysbl_getbase_double (rd,p_Dd)
    
    dcx = 1.0_DP
    if (present(cx)) dcx = cx
    
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
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
            rx%RvectorBlock(1), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Dalpha(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
            rx%RvectorBlock(2), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Dalpha(1,1)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%Dalpha(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
            rx%RvectorBlock(4), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Dalpha(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
            rx%RvectorBlock(5), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Dalpha(2,2)*dcx, 1.0_DP)
      end if

      if (rnonlinearSpatialMatrix%Dalpha(2,1) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
            rx%RvectorBlock(1), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Dalpha(2,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
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
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes, &
            rx%RvectorBlock(1), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Dtheta(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes, &
            rx%RvectorBlock(2), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Dtheta(1,1)*dcx, 1.0_DP)
      end if
            
      if (rnonlinearSpatialMatrix%Dtheta(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes, &
            rx%RvectorBlock(4), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Dtheta(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes, &
            rx%RvectorBlock(5), rd%RvectorBlock(5), &
            -rnonlinearSpatialMatrix%Dtheta(2,2)*dcx, 1.0_DP)
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
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
            rx%RvectorBlock(3), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
            rx%RvectorBlock(3), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Deta(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
            rx%RvectorBlock(6), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Deta(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
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
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
            rx%RvectorBlock(1), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
            rx%RvectorBlock(2), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Dtau(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
            rx%RvectorBlock(4), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%Dtau(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
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
        call lsysbl_deriveSubvector(rvector1,rvectorPrimal, 1,2,.true.)
      case (2)
        call lsysbl_deriveSubvector(rvector2,rvectorPrimal, 1,2,.true.)
      case (3)
        call lsysbl_deriveSubvector(rvector3,rvectorPrimal, 1,2,.true.)
      end select

      select case (rnonlinearSpatialMatrix%idualSol)
      case (1)
        call lsysbl_deriveSubvector(rvector1,rvectorDual, 4,5,.true.)
      case (2)
        call lsysbl_deriveSubvector(rvector2,rvectorDual, 4,5,.true.)
      case (3)
        call lsysbl_deriveSubvector(rvector3,rvectorDual, 4,5,.true.)
      end select

      select case (rnonlinearSpatialMatrix%idualSol2)
      case (1)
        call lsysbl_deriveSubvector(rvector1,rvectorDual2, 4,5,.true.)
      case (2)
        call lsysbl_deriveSubvector(rvector2,rvectorDual2, 4,5,.true.)
      case (3)
        call lsysbl_deriveSubvector(rvector3,rvectorDual2, 4,5,.true.)
      end select

      ! Create a block discretisation by deriving it from the 'full' matrix.
      ! This will serve as a local discretisation structure for all
      ! velocity modifications.
      call spdiscr_deriveBlockDiscr (rnonlinearSpatialMatrix%p_rdiscretisation, &
          rvelDiscr, 1, 2)

      ! Create a 2x2 matrix based on the structure of the FE space.
      ! The matrix does not need any entries, we only need the structure.
      call lsysbl_createMatBlockByDiscr (rvelDiscr,rtempMatrix)
      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes,&
          rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes,&
          rtempMatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes,&
          rtempMatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixStokes,&
          rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

      ! 1.) Primal equation, y*grad(.), probably + grad(.)*y
      !
      !    ( N(y) N(y)                  )
      !    ( N(y) N(y)                  )
      !    (                            )
      !    (                            )
      !    (                            )
      !    (                            )
      
      rstabilisation = t_convecStabilisation(&
          rnonlinearSpatialMatrix%iupwind1,rnonlinearSpatialMatrix%dupsam1,&
            rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixEOJ1)
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,1,2,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,1,2,.true.)

      call assembleConvectionDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,rvectorPrimal,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dgamma(1,1),rnonlinearSpatialMatrix%DgammaT(1,1),&
          rnonlinearSpatialMatrix%Dnewton(1,1),rnonlinearSpatialMatrix%DnewtonT(1,1),&
          rstabilisation,dcx)
      
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
      
      rstabilisation = t_convecStabilisation(&
          rnonlinearSpatialMatrix%iupwind2,rnonlinearSpatialMatrix%dupsam2,&
            rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixEOJ2)
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,4,5,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,4,5,.true.)

      call assembleConvectionDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,rvectorPrimal,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dgamma(2,2),rnonlinearSpatialMatrix%DgammaT(2,2),&
          rnonlinearSpatialMatrix%Dnewton(2,2),rnonlinearSpatialMatrix%DnewtonT(2,2),&
          rstabilisation,dcx)
      
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
      rstabilisation = t_convecStabilisation(0,0.0_DP,NULL())
      
      call lsysbl_deriveSubvector(rx,rtempVectorX,1,2,.true.)
      call lsysbl_deriveSubvector(rd,rtempVectorB,4,5,.true.)

      call assembleConvectionDefect (&
          rnonlinearSpatialMatrix,rtempMatrix,rvectorDual,rtempVectorX,rtempVectorB,&
          rnonlinearSpatialMatrix%Dgamma(2,1),rnonlinearSpatialMatrix%DgammaT(2,1),&
          rnonlinearSpatialMatrix%Dnewton(2,1),rnonlinearSpatialMatrix%DnewtonT(2,1),&
          rstabilisation,dcx)
      
      ! There is probably a 2nd reactive term involved stemming from
      ! the next timestep when Crank-Nicolson is used.

      rstabilisation = t_convecStabilisation(0,0.0_DP,NULL())
      
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
        if (rnonlinearSpatialMatrix%ccontrolConstraints .ne. 0) then
          call cc_projectControlTimestep (rtempVectorX%RvectorBlock(1),&
              rnonlinearSpatialMatrix%dumin1,rnonlinearSpatialMatrix%dumax1)
          call cc_projectControlTimestep (rtempVectorX%RvectorBlock(2),&
              rnonlinearSpatialMatrix%dumin2,rnonlinearSpatialMatrix%dumax2)
        end if

        ! Now carry out MV and include it to the defect.
        ! Note that the multiplication factor is -(-cx) = cx because
        ! it's put on the RHS of the system for creating the defect.
        ! d = b - cx A x = b - ... + \nu Laplace(y) - y\grad(y) - grad(p) + P(-1/alpha lambda)
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
            rtempVectorX%RvectorBlock(1), &
            rd%RvectorBlock(1), dcx, 1.0_DP)
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass, &
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
          call assemblePrimalUConstrMassDefect (rnonlinearSpatialMatrix,rx,&
              rd,dcx*rnonlinearSpatialMatrix%Dalpha(1,2),rvector1)
        case (2)
          call assemblePrimalUConstrMassDefect (rnonlinearSpatialMatrix,rx,&
              rd,dcx*rnonlinearSpatialMatrix%Dalpha(1,2),rvector2)
        case (3)
          call assemblePrimalUConstrMassDefect (rnonlinearSpatialMatrix,rx,&
              rd,dcx*rnonlinearSpatialMatrix%Dalpha(1,2),rvector3)
        end select
        
      end if
      
      ! Release the temp vectors/matrices, that's it.
      call spdiscr_releaseBlockDiscr(rvelDiscr)
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
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
            rx%RvectorBlock(3), rd%RvectorBlock(1), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
            rx%RvectorBlock(3), rd%RvectorBlock(2), &
            -rnonlinearSpatialMatrix%Deta(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Deta(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB1, &
            rx%RvectorBlock(6), rd%RvectorBlock(4), &
            -rnonlinearSpatialMatrix%Deta(2,2)*dcx, 1.0_DP)

        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixB2, &
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
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
            rx%RvectorBlock(1), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
            rx%RvectorBlock(2), rd%RvectorBlock(3), &
            -rnonlinearSpatialMatrix%Dtau(1,1)*dcx, 1.0_DP)
      end if
      
      if (rnonlinearSpatialMatrix%Dtau(2,2) .ne. 0.0_DP) then
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD1, &
            rx%RvectorBlock(4), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%Dtau(2,2)*dcx, 1.0_DP)
        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixD2, &
            rx%RvectorBlock(5), rd%RvectorBlock(6), &
            -rnonlinearSpatialMatrix%Dtau(2,2)*dcx, 1.0_DP)
      end if
      
      ! ---------------------------------------------------
      ! Now a slightly more advanced task for which we use a separate
      ! routine and some submatrices/vectors: The nonlinearity.

      ! Initialise the operator structure for what we need.
      roptcoperator%dupsamPrimal = rnonlinearSpatialMatrix%dupsam1
      roptcoperator%dupsamDual = rnonlinearSpatialMatrix%dupsam2
      
      ! Timestep-weights
      roptcoperator%dprimalAlpha = rnonlinearSpatialMatrix%Dalpha(1,1)
      roptcoperator%ddualAlpha   = rnonlinearSpatialMatrix%Dalpha(2,2)

      ! Stokes operator
      roptcoperator%dnu = rnonlinearSpatialMatrix%dnu
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
      
      if (rnonlinearSpatialMatrix%ccontrolConstraints .ne. 0) then
        roptcoperator%ccontrolProjection = rnonlinearSpatialMatrix%ccontrolConstraints
        roptcoperator%dmin1 = rnonlinearSpatialMatrix%dumin1
        roptcoperator%dmax1 = rnonlinearSpatialMatrix%dumax1
        roptcoperator%dmin2 = rnonlinearSpatialMatrix%dumin2
        roptcoperator%dmax2 = rnonlinearSpatialMatrix%dumax2
      end if
      
      select case (rnonlinearSpatialMatrix%iprimalSol)
      case (1)
        p_rprimalSol => rvector1
      case (2)
        p_rprimalSol => rvector2
      case (3)
        p_rprimalSol => rvector3
      end select
      
      select case (rnonlinearSpatialMatrix%idualSol)
      case (1)
        p_rdualSol => rvector1
      case (2)
        p_rdualSol => rvector2
      case (3)
        p_rdualSol => rvector3
      end select
      
      ! Calculate the velocity-dependent part of the system matrix.
      call conv_strdiffOptC2dgetDefect (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
          roptcoperator,p_rprimalSol,p_rdualSol,dcx,rx,rd)
    
    end if
    
!    ! --------------------------------------------------------------------------
!    !
!    ! Create a temporary matrix that covers this structure.
!    call lsysbl_createEmptyMatrix (rmatrix,2*(NDIM2D+1))
!
!    ! Put references to the Stokes- and B-matrices to Aij. assembleVelocityDefect
!    ! needs this template matrix to provide the structure for the stabilisation
!    ! routines! The B-matrices are needed later.
!    ! -----
!    ! Primal equation
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!        rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!        rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    if (rnonlinearSpatialMatrix%dnewton1 .ne. 0.0_DP) then
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!          rmatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!          rmatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!    end if
!
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1,&
!        rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2,&
!        rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                  rmatrix%RmatrixBlock(3,1),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                  rmatrix%RmatrixBlock(3,2),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!        rmatrix%RmatrixBlock(1,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!        rmatrix%RmatrixBlock(2,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!    ! -----
!    ! Dual equation
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!        rmatrix%RmatrixBlock(4,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!        rmatrix%RmatrixBlock(5,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    if (rnonlinearSpatialMatrix%dnewton2 .ne. 0.0_DP) then
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!          rmatrix%RmatrixBlock(4,5),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!      call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixStokes,&
!          rmatrix%RmatrixBlock(5,4),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!    end if
!
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1,&
!        rmatrix%RmatrixBlock(4,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2,&
!        rmatrix%RmatrixBlock(5,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                  rmatrix%RmatrixBlock(6,4),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                  rmatrix%RmatrixBlock(6,5),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!        rmatrix%RmatrixBlock(4,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!        rmatrix%RmatrixBlock(5,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!    ! Update the structural information of the block matrix, as we manually
!    ! changed the submatrices:
!    call lsysbl_updateMatStrucInfo (rmatrix)
!
!    ! Get a reference to the correct velocity vector we have to use for the
!    ! nonlinearity.
!    select case (rnonlinearSpatialMatrix%iprimalSol)
!    case (1)
!      p_rvectorPrimal => rvector1
!    case (2)
!      p_rvectorPrimal => rvector2
!    case (3)
!      p_rvectorPrimal => rvector3
!    end select
!
!    select case (rnonlinearSpatialMatrix%idualSol)
!    case (1)
!      p_rvectorDual => rvector1
!    case (2)
!      p_rvectorDual => rvector2
!    case (3)
!      p_rvectorDual => rvector3
!    end select
!
!    ! OK, LET'S GO...
!    !
!    ! At first, initialise a matrix containing only the B-matrices. That's the easiest
!    ! part.
!    !
!    !    (           B1               )
!    !    (           B2               )
!    !    ( B1^T B2^T                  )
!    !    (                         B1 )
!    !    (                         B2 )
!    !    (               B1^T B2^T    )
!    !
!    call lsysbl_createEmptyMatrix (rmatrix,2*(NDIM2D+1))
!
!    ! Primal
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1,&
!        rmatrix%RmatrixBlock(1,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2,&
!        rmatrix%RmatrixBlock(2,3),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                  rmatrix%RmatrixBlock(3,1),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                  rmatrix%RmatrixBlock(3,2),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    ! Dual
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB1,&
!        rmatrix%RmatrixBlock(4,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!    call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixB2,&
!        rmatrix%RmatrixBlock(5,6),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB1, &
!                                  rmatrix%RmatrixBlock(6,4),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    call lsyssc_transposeMatrix (rnonlinearSpatialMatrix%p_rmatrixB2, &
!                                  rmatrix%RmatrixBlock(6,5),&
!                                  LSYSSC_TR_VIRTUAL)
!
!    ! Initialise the weights for the B/B^T matrices
!    rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearSpatialMatrix%deta1
!    rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearSpatialMatrix%deta1
!
!    rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearSpatialMatrix%dtau1
!    rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearSpatialMatrix%dtau1
!
!    rmatrix%RmatrixBlock(4,6)%dscaleFactor = rnonlinearSpatialMatrix%deta2
!    rmatrix%RmatrixBlock(5,6)%dscaleFactor = rnonlinearSpatialMatrix%deta2
!
!    rmatrix%RmatrixBlock(6,4)%dscaleFactor = rnonlinearSpatialMatrix%dtau2
!    rmatrix%RmatrixBlock(6,5)%dscaleFactor = rnonlinearSpatialMatrix%dtau2
!
!    ! Create the defect
!    call lsysbl_blockMatVec (rmatrix, rx, rd, -dcx, 1.0_DP)
!
!
!
!
!
!
!
!
!    ! In the first step, we assemble the defect that arises in the velocity
!    ! components. This is characterised by the following submatrix:
!    !
!    !    ( A11  A12                   )
!    !    ( A21  A22                   )
!    !    (                            )
!    !    (               A44  A45     )
!    !    (               A54  A55     )
!    !    (                            )
!    !
!    ! assembleVelocityDefect handles exactly these submatrices.
!    ! We call the routine twice -- once for the primal and once for the dual
!    ! equation. In both cases, we specify the primal velocity p_ry
!    ! as velocity field (!).
!
!    select case (rnonlinearSpatialMatrix%iprimalSol)
!    case (1)
!      call assembleVelocityDefect (&
!          .false.,rnonlinearSpatialMatrix,rmatrix,rx,rd,dcx,rvector1,1.0_DP)
!      call assembleVelocityDefect (&
!          .true.,rnonlinearSpatialMatrix,rmatrix,rx,rd,dcx,rvector1,1.0_DP)
!    case (2)
!      call assembleVelocityDefect (&
!          .false.,rnonlinearSpatialMatrix,rmatrix,rx,rd,dcx,rvector2,1.0_DP)
!      call assembleVelocityDefect (&
!          .true.,rnonlinearSpatialMatrix,rmatrix,rx,rd,dcx,rvector2,1.0_DP)
!    case (3)
!      call assembleVelocityDefect (&
!          .false.,rnonlinearSpatialMatrix,rmatrix,rx,rd,dcx,rvector3,1.0_DP)
!      call assembleVelocityDefect (&
!          .true.,rnonlinearSpatialMatrix,rmatrix,rx,rd,dcx,rvector3,1.0_DP)
!    end select
!
!    ! Now, we treat all the remaining blocks. Let's see what is missing:
!    !
!    !    ( .    .    B1  M             )
!    !    ( .    .    B2       M        )
!    !    ( B1^T B2^T .                 )
!    !    ( M             .    .    B1  )
!    !    (      M        .    .    B2  )
!    !    (               B1^T B2^T .   )
!
!    ! To build the appropriate defect, we first remove the velocity blocks:
!
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,1))
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(1,2))
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,1))
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(2,2))
!
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(4,4))
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(4,5))
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(5,4))
!    call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(5,5))
!
!    ! Initialise the weights for the B/B^T matrices
!    rmatrix%RmatrixBlock(1,3)%dscaleFactor = rnonlinearSpatialMatrix%deta1
!    rmatrix%RmatrixBlock(2,3)%dscaleFactor = rnonlinearSpatialMatrix%deta1
!
!    rmatrix%RmatrixBlock(3,1)%dscaleFactor = rnonlinearSpatialMatrix%dtau1
!    rmatrix%RmatrixBlock(3,2)%dscaleFactor = rnonlinearSpatialMatrix%dtau1
!
!    rmatrix%RmatrixBlock(4,6)%dscaleFactor = rnonlinearSpatialMatrix%deta2
!    rmatrix%RmatrixBlock(5,6)%dscaleFactor = rnonlinearSpatialMatrix%deta2
!
!    rmatrix%RmatrixBlock(6,4)%dscaleFactor = rnonlinearSpatialMatrix%dtau2
!    rmatrix%RmatrixBlock(6,5)%dscaleFactor = rnonlinearSpatialMatrix%dtau2
!
!    ! Initialise the weights for the mass matrices
!    !if (rnonlinearSpatialMatrix%ccontrolConstraints .eq. 0) then
!    !  rmatrix%RmatrixBlock(1,4)%dscaleFactor = rnonlinearSpatialMatrix%dmu1
!    !  rmatrix%RmatrixBlock(2,5)%dscaleFactor = rnonlinearSpatialMatrix%dmu1
!    !else
!      rmatrix%RmatrixBlock(1,4)%dscaleFactor = 0.0_DP
!      rmatrix%RmatrixBlock(2,5)%dscaleFactor = 0.0_DP
!    !end if
!
!    rmatrix%RmatrixBlock(4,1)%dscaleFactor = rnonlinearSpatialMatrix%dmu2
!    rmatrix%RmatrixBlock(5,2)%dscaleFactor = rnonlinearSpatialMatrix%dmu2
!
!    ! ------------------------------------------------
!    ! Build the defect by matrix-vector multiplication
!    !
!    ! Note that no time step or whatever is included here; everything
!    ! is initialised with the multiplication factors in the submatrices
!    ! from above!
!    call lsysbl_blockMatVec (rmatrix, rx, rd, -dcx, 1.0_DP)
!
!    ! The coupling of lambda to the primal equation is a little bit tricky.
!    ! Ok, if there are no constraints active, it's easy -- that case was handled
!    ! in the MV above...
!    if ( &!(rnonlinearSpatialMatrix%ccontrolConstraints .eq. 1) .and. &
!        (rnonlinearSpatialMatrix%dmu1 .ne. 0.0_DP)) then
!
!      ! But now it get's interesting. When control constraints are active,
!      ! we have the system
!      !
!      !    ( .    .    B1  PM            )
!      !    ( .    .    B2       PM       )
!      !    ( B1^T B2^T .                 )
!      !    ( M             .    .    B1  )
!      !    (      M        .    .    B2  )
!      !    (               B1^T B2^T .   )
!      !
!      ! Where P is a projection operator onto the allowed space. The dual variable
!      ! \lambda must not be changed, but before adding M\lambda to the RHS,
!      ! we have to apply a projection!
!      !
!      ! That's a little bit ugly because for this we have to carry out the matrix
!      ! vector multiplication 'by hand'.
!
!
!      ! What's the type of the current matrix? Is this a Newton-matrix?
!      if (rnonlinearSpatialMatrix%cmatrixType .eq. 0) then
!
!        ! No, this is a standard matrix. That means, we just have to project
!        ! the control u and multiply it with the mass matrix.
!
!!      Projection of the solution. Test code. Not used as the inequality
!!      defining the coupling must be used in integral form -- for what
!!      the implementation below that is correct.
!
!!      ! At first, create a temporary vector that
!!      ! receives the projected solution.
!!      call lsyssc_duplicateVector (rx%RvectorBlock(4),rtempVector,&
!!        LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
!!
!!      ! Multiply with the mass matrix, correctly scaled.
!!      call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixMass, rx%RvectorBlock(4), &
!!        rtempVector, rnonlinearSpatialMatrix%dmu1, 0.0_DP)
!!
!!      ! Project that to the allowed range.
!!      call cc_projectControlTimestep (rtempVector,&
!!          -rnonlinearSpatialMatrix%dumax1,-rnonlinearSpatialMatrix%dumin1)
!!
!!      ! And then finally, carry our the defect calculation for y_1.
!!      call lsyssc_vectorLinearComb (rtempVector,rd%RvectorBlock(1),-dcx,1.0_DP)
!!
!!      ! The same stuff has to be done for y_1 / lambda_2:
!!      call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixMass, rx%RvectorBlock(5), &
!!        rtempVector, rnonlinearSpatialMatrix%dmu1, 0.0_DP)
!!      call cc_projectControlTimestep (rtempVector,&
!!          -rnonlinearSpatialMatrix%dumax2,-rnonlinearSpatialMatrix%dumin2)
!!      call lsyssc_vectorLinearComb (rtempVector,rd%RvectorBlock(2),-dcx,1.0_DP)
!
!        ! Copy our solution vector \lambda_1. Scale it by -1/alpha.
!        call lsyssc_duplicateVector (rx%RvectorBlock(4),rtempVector,&
!            LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)
!
!        call lsyssc_scaleVector (rtempVector,-rnonlinearSpatialMatrix%dmu1)
!
!        ! Project that to the allowed range to create u_1.
!        if (rnonlinearSpatialMatrix%ccontrolConstraints .eq. 1) then
!          call cc_projectControlTimestep (rtempVector,&
!              rnonlinearSpatialMatrix%dumin1,rnonlinearSpatialMatrix%dumax1)
!        end if
!
!        ! Now carry out MV and include it to the defect.
!        ! Note that the multiplication factor is -(-cx) = cx because
!        ! it's put on the RHS of the system for creating the defect.
!        ! d = b - cx A x = b - ... + \nu Laplace(y) - y\grad(y) - grad(p) + P(-1/alpha lambda)
!        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixMass, rtempVector, &
!            rd%RvectorBlock(1), dcx, 1.0_DP)
!
!        ! The same stuff has to be done for y_2 / lambda_2:
!        call lsyssc_duplicateVector (rx%RvectorBlock(5),rtempVector,&
!            LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)
!        call lsyssc_scaleVector (rtempVector,-rnonlinearSpatialMatrix%dmu1)
!
!        if (rnonlinearSpatialMatrix%ccontrolConstraints .eq. 1) then
!          call cc_projectControlTimestep (rtempVector,&
!              rnonlinearSpatialMatrix%dumin2,rnonlinearSpatialMatrix%dumax2)
!        end if
!
!        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixMass, rtempVector, &
!          rd%RvectorBlock(2), dcx, 1.0_DP)
!
!        ! Release the temp vector, that's it.
!        call lsyssc_releaseVector (rtempVector)
!
!      else
!
!        ! Yes, that's a Newton matrix. That means, we have to multiply the
!        ! vector with the derivative of the projection operator:
!        ! b-P[a,b]'(lambda).
!        ! For that purpose, we have to assemble special mass matrices:
!
!        select case (rnonlinearSpatialMatrix%idualSol)
!        case (1)
!          call assemblePrimalUConstrMassDefect (rnonlinearSpatialMatrix,rx,&
!              rd,dcx,rvector1)
!        case (2)
!          call assemblePrimalUConstrMassDefect (rnonlinearSpatialMatrix,rx,&
!              rd,dcx,rvector2)
!        case (3)
!          call assemblePrimalUConstrMassDefect (rnonlinearSpatialMatrix,rx,&
!              rd,dcx,rvector3)
!        end select
!
!      end if
!
!    end if
!
!    ! If we have a reactive coupling mass matrix, it gets interesting...
!    !
!    ! Primal equation
!    if ((rnonlinearSpatialMatrix%dr11 .ne. 0.0_DP) .or. &
!        (rnonlinearSpatialMatrix%dr12 .ne. 0.0_DP)) then
!
!      ! Assemble the defect of the reactive coupling mass matrices
!      select case (rnonlinearSpatialMatrix%idualSol)
!      case (1)
!        call assemblePrimalReactMassDefect (rnonlinearSpatialMatrix,rx,&
!            rd,dcx,rvector1,1.0_DP)
!      case (2)
!        call assemblePrimalReactMassDefect (rnonlinearSpatialMatrix,rx,&
!            rd,dcx,rvector2,1.0_DP)
!      case (3)
!        call assemblePrimalReactMassDefect (rnonlinearSpatialMatrix,rx,&
!            rd,dcx,rvector3,1.0_DP)
!      end select
!
!    end if
!
!    ! Dual equation
!    if ((rnonlinearSpatialMatrix%dr21 .ne. 0.0_DP) .or. &
!        (rnonlinearSpatialMatrix%dr22 .ne. 0.0_DP)) then
!
!      ! Assemble the defect of the reactive coupling mass matrices
!      select case (rnonlinearSpatialMatrix%idualSol)
!      case (1)
!        call assembleDualReactMassDefect (rnonlinearSpatialMatrix,rx,&
!            rd,dcx,rvector1,1.0_DP)
!      case (2)
!        call assembleDualReactMassDefect (rnonlinearSpatialMatrix,rx,&
!            rd,dcx,rvector2,1.0_DP)
!      case (3)
!        call assembleDualReactMassDefect (rnonlinearSpatialMatrix,rx,&
!            rd,dcx,rvector3,1.0_DP)
!      end select
!
!    end if
!
!    ! Release the temporary matrix, we don't need it anymore.
!    call lsysbl_releaseMatrix (rmatrix)

  contains

    ! -----------------------------------------------------

    subroutine assembleConvectionDefect (&
        rnonlinearSpatialMatrix,rmatrix,rvector,rx,rb,dgamma,dgammaT,dnewton,dnewtonT,&
        rstabilisation,dcx)
        
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
    type(t_convecStabilisation), intent(in) :: rstabilisation
    
    ! local variables
    logical :: bshared
    type(t_convUpwind) :: rupwind
    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
    type(t_convStreamDiff2) :: rstreamlineDiffusion2
    type(t_jumpStabilisation) :: rjumpStabil
    
    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_Ddata1,p_Ddata2
    call lsysbl_getbase_double (rx,p_Ddata1)
    call lsysbl_getbase_double (rb,p_Ddata2)
    
    
      ! Is A11=A22 physically?
      bshared = lsyssc_isMatrixContentShared(&
                    rmatrix%RmatrixBlock(1,1),&
                    rmatrix%RmatrixBlock(2,2))
                    
      if ((dgamma .ne. 0.0_DP) .or. (dgammaT .ne. 0.0_DP) .or. &
          (dnewton .ne. 0.0_DP) .or. (dnewtonT .ne. 0.0_DP)) then
        select case (rstabilisation%iupwind)
        case (CCMASM_STAB_STREAMLINEDIFF)
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
          
          rstreamlineDiffusion%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion%ddelta            = dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dgammaT
          rstreamlineDiffusion%dnewton           = dnewton
          rstreamlineDiffusion%dnewtonTransposed = dnewtonT
          
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
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%dnu
          
          rstreamlineDiffusion2%dtheta = dcx
                    
          ! Set stabilisation parameter
          rstreamlineDiffusion2%dupsam = rstabilisation%dupsam
          
          ! Matrix weights
          rstreamlineDiffusion2%ddelta  = dgamma
          rstreamlineDiffusion2%ddeltaT = dgammaT
          rstreamlineDiffusion2%dnewton = dnewton
          rstreamlineDiffusion2%dnewtonT = dnewtonT
          
          ! Call the SD method to calculate the nonlinearity.
          ! As velocity vector, specify rvector!
          ! Therefore, the primal velcity is always used for assembling
          ! that thing!
          call conv_streamDiff2Blk2dDef (rstreamlineDiffusion2,rmatrix,&
              rx,rb,rvector)

        case (CCMASM_STAB_UPWIND)
        
          call output_line ('Upwind not supported.', &
                            OU_CLASS_ERROR,OU_MODE_STD,'assembleConvection')
          call sys_halt()

        case (CCMASM_STAB_EDGEORIENTED)
          ! Jump stabilisation.
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
          
          rstreamlineDiffusion%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion%ddelta            = dgamma
          rstreamlineDiffusion%ddeltaTransposed  = dgammaT
          rstreamlineDiffusion%dnewton           = dnewton
          rstreamlineDiffusion%dnewtonTransposed = dnewtonT
          
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
          rjumpStabil%dgamma = rstabilisation%dupsam
          
          ! Matrix weight
          rjumpStabil%dtheta = dcx*dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rx,rb)

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
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%dnu
          
          rstreamlineDiffusion2%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dgamma
          rstreamlineDiffusion2%ddeltaT  = dgammaT
          rstreamlineDiffusion2%dnewton  = dnewton
          rstreamlineDiffusion2%dnewtonT = dnewtonT
          
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
          rjumpStabil%dgamma = rstabilisation%dupsam
          
          ! Matrix weight
          rjumpStabil%dtheta = dcx*dgamma

          ! Call the jump stabilisation technique to stabilise that stuff.
          ! We can assemble the jump part any time as it's independent of any
          ! convective parts...
          !
          ! The routine processes the X as well as the Y velocity!
          call conv_jumpStabilisation2d (&
                              rjumpStabil, CONV_MODDEFECT, &
                              rmatrix%RmatrixBlock(1,1),rx,rb)

        case (CCMASM_STAB_EDGEORIENTED3)
        
          ! Jump stabilisation, precomputed matrix.
          !
          ! In the first step, set up the matrix as above with central discretisation,
          ! i.e. call SD to calculate the matrix without SD stabilisation.
          ! Set up the SD structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rstreamlineDiffusion2%dnu = rnonlinearSpatialMatrix%dnu
          
          rstreamlineDiffusion2%dtheta = dcx
          
          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
          rstreamlineDiffusion2%dupsam = 0.0_DP
          
          ! Matrix weight
          rstreamlineDiffusion2%ddelta   = dgamma
          rstreamlineDiffusion2%ddeltaT  = dgammaT
          rstreamlineDiffusion2%dnewton  = dnewton
          rstreamlineDiffusion2%dnewtonT = dnewtonT
          
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
          call lsyssc_scalarMatVec(rstabilisation%p_rmatrixEOJ,&
              rx%RvectorBlock(1),rb%RvectorBlock(1),&
              -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

          call lsyssc_scalarMatVec(rstabilisation%p_rmatrixEOJ,&
              rx%RvectorBlock(2),rb%RvectorBlock(2),&
              -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

        case default
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select

      else
      
        ! That's the Stokes-case. Jump stabilisation is possible...
        if (rstabilisation%dupsam .ne. 0.0_DP) then
          select case (rstabilisation%iupwind)
          case (CCMASM_STAB_EDGEORIENTED,CCMASM_STAB_EDGEORIENTED2)
            
            ! Set up the jump stabilisation structure.
            ! There's not much to do, only initialise the viscosity...
            rjumpStabil%dnu = rnonlinearSpatialMatrix%dnu
            
            ! Set stabilisation parameter
            rjumpStabil%dgamma = rstabilisation%dupsam
            
            ! Matrix weight
            rjumpStabil%dtheta = dcx*dgamma

            ! Call the jump stabilisation technique to stabilise that stuff.
            ! We can assemble the jump part any time as it's independent of any
            ! convective parts...
            call conv_jumpStabilisation2d (&
                                rjumpStabil, CONV_MODDEFECT, &
                                rmatrix%RmatrixBlock(1,1),rx,rb)

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
            call lsyssc_scalarMatVec(rstabilisation%p_rmatrixEOJ,&
                rx%RvectorBlock(1),rb%RvectorBlock(1),&
                -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

            call lsyssc_scalarMatVec(rstabilisation%p_rmatrixEOJ,&
                rx%RvectorBlock(2),rb%RvectorBlock(2),&
                -dcx*dgamma*mprim_signum(rstabilisation%dupsam),1.0_DP)

          case default
            ! No stabilisation
          
          end select
        end if
        
      end if
      
    end subroutine

    subroutine assemblePrimalUConstrMassDefect (rnonlinearSpatialMatrix,&
        rvector,rdefect,dcx,rvelocityVector) !,cprojectionType)
        
    ! Assembles the defect arising from the projective coupling mass
    ! matrices in the primal equation which comes from constraints on u.
    ! rdefect must have been initialised with the right hand side vector.
    !
    ! Let the mass matrix M~ be the usual mass matrix where u is ok
    ! and the 0-operator where u is out of bounds.
    ! Then, we assemble
    !
    !       rdefect = r(primal)defect - dcx (dmu1 M~ r(dual)vector)
    !
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
      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
      type(t_matrixBlock) :: rtempmatrix
      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
      type(t_collection) :: rcollection
      type(t_bilinearForm) :: rform
      integer, dimension(:), pointer :: p_IelementList
      integer :: ielemHandle,nelements
      type(t_bilfMatrixAssembly) :: rmatrixAssembly
      integer(I32) :: celement,ccubType
      
      ! If we have a reactive coupling mass matrix, it gets interesting...
      if (dcx .ne. 0.0_DP) then

        call lsysbl_createEmptyMatrix (rtempMatrix,2)

        ! The ccontrolConstraints decides on whether we use the 'quick-setup'
        ! method for the mass matrices or rebuild them again.
        ! ccontrolConstraints=0 uses standard mass matrices.
        select case (rnonlinearSpatialMatrix%ccontrolConstraints)
        case (0)
        
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
           
          call lsysbl_updateMatStrucInfo (rtempMatrix)
             
        case (1)
        
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
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
          call massmatfilter (rtempMatrix%RmatrixBlock(1,1),rvelocityVector%RvectorBlock(4),&
              rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%dumin1,rnonlinearSpatialMatrix%dumax1)
          call massmatfilter (rtempMatrix%RmatrixBlock(2,2),rvelocityVector%RvectorBlock(5),&
              rnonlinearSpatialMatrix%dalphaC,rnonlinearSpatialMatrix%dumin2,rnonlinearSpatialMatrix%dumax2)
            
        case (2)

          ! Create a matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly
          CALL lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          CALL lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          CALL lsysbl_updateMatStrucInfo (rtempMatrix)

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
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin1
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax1
          
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
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin2
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax2

          call bilf_buildMatrixScalar (rform,.TRUE.,rtempmatrix%RmatrixBlock(2,2),&
              coeff_ProjMass,rcollection)
          
          ! Now we can forget about the collection again.
          call collct_done (rcollection)
          
        case (3)
        
          ! Create a matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly
          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          call lsysbl_updateMatStrucInfo (rtempMatrix)

          call approxProjectionDerivative (rtempMatrix%RmatrixBlock(1,1), &
              rvelocityVector%RvectorBlock(4),rnonlinearSpatialMatrix%dalphaC,&
              rnonlinearSpatialMatrix%dumin1,rnonlinearSpatialMatrix%dumax1,0.001_DP)

          call approxProjectionDerivative (rtempMatrix%RmatrixBlock(2,2), &
              rvelocityVector%RvectorBlock(5),rnonlinearSpatialMatrix%dalphaC,&
              rnonlinearSpatialMatrix%dumin2,rnonlinearSpatialMatrix%dumax2,0.001_DP)
        
          call lsysbl_updateMatStrucInfo (rtempMatrix)
            
        case (4)
        
          ! Exact reassembly of the mass matrices with adaptive integration.

          ! Create a matrix with the structure we need. Share the structure
          ! of the mass matrix. Entries are not necessary for the assembly
          CALL lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          CALL lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rstaticInfo%rmatrixMass,&
              rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

          CALL lsysbl_updateMatStrucInfo (rtempMatrix)

          ! Create an array that saves all elements on the border of the active set.
          nelements = rnonlinearSpatialMatrix%p_rdiscretisation%RspatialDiscr(1)%&
              p_rtriangulation%NEL
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
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin1
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax1
          
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
            celement = rtempmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(1)%celement
            ccubType = rtempmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%RelementDistr(1)%ccubTypeBilForm
            call storage_getbase_int(ielemhandle,p_IelementList)
            call bilf_initAssembly(rmatrixAssembly,rform,celement,celement,&
                cub_getSummedCubType(ccubType,1))
            call bilf_assembleSubmeshMatrix9(rmatrixAssembly,rtempmatrix%RmatrixBlock(1,1),&
                p_IelementList(1:rcollection%IquickAccess(3)),coeff_ProjMass,rcollection)
            call bilf_doneAssembly(rmatrixAssembly)
          end if

          ! Now, set up A25, depending on lambda_2.
          rcollection%IquickAccess(1) = 2
          rcollection%DquickAccess(1) = rnonlinearSpatialMatrix%dumin2
          rcollection%DquickAccess(2) = rnonlinearSpatialMatrix%dumax2

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








!    subroutine assembleVelocityDefect (bdualEquation,rnonlinearSpatialMatrix,&
!        rmatrix,rvector,rdefect,dcx,rvelocityVector,dvectorWeight)
!
!    ! Assembles the velocity defect in the block matrix rmatrix at position
!    ! itop..itop+1 in the velocity vector. rdefect must have been initialised
!    ! with the right hand side vector.
!    !
!    ! With a matrix 'A' of the theoretical form
!    !
!    !       A := dalpha*M + dtheta*Laplace + dgamma*N(p_rvector) +
!    !            dnewton*N*(p_rvector)
!    !
!    ! the routine will construct
!    !
!    !       rdefect = rdefect - dcx * (dtheta A rvector)
!
!    ! Whether to set up the primal or the dual equation.
!    ! FALSE=primal, TRUE=dual equation.
!    logical, intent(IN) :: bdualEquation
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how to set up the matrix.
!    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
!
!    ! Reference to the system matrix. Only the structure of the matrix
!    ! is used to reconstruct the structure of the discretisation.
!    ! The content of the matrix is not changed or used.
!    type(t_matrixBlock), intent(INOUT) :: rmatrix
!
!    ! Solution vector.
!    type(t_vectorBlock), intent(IN) :: rvector
!
!    ! On entry: RHS vector.
!    ! Is overwritten by the defect vector in the velocity subsystem.
!    type(t_vectorBlock), intent(INOUT) :: rdefect
!
!    ! Multiplication factor for the whole operator A*rvector
!    real(DP), intent(IN) :: dcx
!
!    ! Weight for the velocity vector rvelocityVector; usually = 1.0
!    real(DP), intent(IN) :: dvectorWeight
!
!    ! Velocity vector field that should be used for the assembly of the
!    ! nonlinearity. The first two blocks in that block vector are
!    ! used as velocity field.
!    type(t_vectorBlock), intent(IN) :: rvelocityVector
!
!    ! local variables
!    logical :: bshared
!    integer :: iupwind,dupsam
!    type(t_convUpwind) :: rupwind
!    type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
!    type(T_jumpStabilisation) :: rjumpStabil
!    type(t_vectorBlock) :: rtempVector,rtempDefect
!    type(t_matrixBlock) :: rtempMatrix
!    integer :: imatOffset
!    real(DP) :: dalpha, dtheta, dnewton, diota, dgamma
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Dd
!
!    call lsysbl_getbase_double (rdefect,p_Dd)
!
!      if (.not. bdualEquation) then
!        ! Set the weights used here according to the primal equation.
!        ! Set imatOffset=0 so the submatrix at position 1,1 is tackled.
!        imatOffset = 0
!        dalpha = rnonlinearSpatialMatrix%dalpha1
!        dtheta = rnonlinearSpatialMatrix%dtheta1
!        diota  = rnonlinearSpatialMatrix%diota1
!        dgamma = rnonlinearSpatialMatrix%dgamma1
!        dnewton = rnonlinearSpatialMatrix%dnewton1
!        iupwind = rnonlinearSpatialMatrix%iupwind1
!        dupsam = rnonlinearSpatialMatrix%dupsam1
!      else
!        ! Set the weights used here according to the primal equation.
!        ! Set imatOffset=3 so the submatrix at position 4,4 is tackled.
!        imatOffset = 3
!        dalpha = rnonlinearSpatialMatrix%dalpha2
!        dtheta = rnonlinearSpatialMatrix%dtheta2
!        diota  = rnonlinearSpatialMatrix%diota2
!        dgamma = rnonlinearSpatialMatrix%dgamma2
!        dnewton = rnonlinearSpatialMatrix%dnewton2
!        iupwind = rnonlinearSpatialMatrix%iupwind2
!        dupsam = rnonlinearSpatialMatrix%dupsam2
!      end if
!
!      ! Derive a temporary vector that contains only those velocity
!      ! subvectors that might affect the matrix.
!      call lsysbl_deriveSubvector(rvector,rtempVector, &
!          imatOffset+1,imatOffset+2,.true.)
!      call lsysbl_deriveSubvector(rdefect,rtempDefect, &
!          imatOffset+1,imatOffset+2,.true.)
!
!      ! Create a temporary block matrix only contining the velocity submatrices
!      ! we want to change. Share structure and entries such that changing
!      ! the temporary matrix will also change the original matrix.
!      call lsysbl_deriveSubmatrix (rmatrix,rtempMatrix,&
!                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,&
!                                    imatOffset+1,imatOffset+2)
!
!      ! Is A11=A22 physically?
!      bshared = lsyssc_isMatrixContentShared(&
!                    rtempmatrix%RmatrixBlock(1,1),&
!                    rtempmatrix%RmatrixBlock(2,2))
!
!      ! ---------------------------------------------------
!      ! Subtract Ix ?
!      if (diota*dcx .ne. 0.0_DP) then
!        call lsyssc_vectorLinearComb (&
!            rvector%RvectorBlock(imatOffset+1), &
!            rdefect%RvectorBlock(imatOffset+1), &
!            -diota*dcx, 1.0_DP)
!
!        call lsyssc_vectorLinearComb (&
!            rvector%RvectorBlock(imatOffset+2), &
!            rdefect%RvectorBlock(imatOffset+2), &
!            -diota*dcx, 1.0_DP)
!      end if
!
!      ! ---------------------------------------------------
!      ! Subtract the mass matrix stuff?
!      if (dalpha*dcx .ne. 0.0_DP) then
!        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixMass, &
!            rvector%RvectorBlock(imatOffset+1), &
!            rdefect%RvectorBlock(imatOffset+1), &
!            -dalpha*dcx, 1.0_DP)
!
!        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixMass, &
!            rvector%RvectorBlock(imatOffset+2), &
!            rdefect%RvectorBlock(imatOffset+2), &
!            -dalpha*dcx, 1.0_DP)
!      end if
!
!      ! ---------------------------------------------------
!      ! Subtract the Stokes matrix stuff?
!      if (dtheta*dcx .ne. 0.0_DP) then
!        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixStokes, &
!            rvector%RvectorBlock(imatOffset+1), &
!            rdefect%RvectorBlock(imatOffset+1), &
!            -dtheta*dcx, 1.0_DP)
!
!        call lsyssc_scalarMatVec (rnonlinearSpatialMatrix%p_rmatrixStokes, &
!            rvector%RvectorBlock(imatOffset+2), &
!            rdefect%RvectorBlock(imatOffset+2), &
!            -dtheta*dcx, 1.0_DP)
!      end if
!
!      ! ---------------------------------------------------
!      ! That was easy -- the adventure begins now... The nonlinearity!
!      if (dgamma*dcx .ne. 0.0_DP) then
!
!        select case (iupwind)
!        case (0)
!          ! Set up the SD structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter
!          rstreamlineDiffusion%dupsam = dupsam
!
!          ! Matrix weight
!          rstreamlineDiffusion%ddelta = dgamma*dcx
!
!          ! Weight for the Newtop part; =0 deactivates Newton.
!          rstreamlineDiffusion%dnewton = dnewton*dcx
!
!          ! Call the SD method to calculate the defect of the nonlinearity.
!          ! As rrhsTemp shares its entries with rdefect, the result is
!          ! directly written to rdefect!
!          ! As velocity field, we specify rvelocityVector here -- the primal
!          ! velocity!. The first two subvectors are used as velocity field.
!
!          call conv_streamlineDiffusionBlk2d (&
!                              rvelocityVector, rvelocityVector, &
!                              dvectorWeight, 0.0_DP,&
!                              rstreamlineDiffusion, CONV_MODDEFECT, &
!                              rtempMatrix,rsolution=rtempVector,rdefect=rtempDefect)
!
!        case (1)
!          ! Set up the upwind structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rupwind%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter
!          rupwind%dupsam = dupsam
!
!          ! Matrix weight
!          rupwind%dtheta = dgamma*dcx
!
!          ! Call the upwind method to calculate the nonlinear defect.
!          call conv_upwind2d (rtempVector, rtempVector, &
!                              dvectorWeight, 0.0_DP,&
!                              rupwind, CONV_MODDEFECT, &
!                              rtempMatrix%RmatrixBlock(1,1),&
!                              rtempVector,rtempDefect)
!
!          if (.not. bshared) then
!            print *,'Upwind does not support independent A11/A22!'
!            stop
!          end if
!
!        case (2)
!          ! Jump stabilisation.
!          ! In the first step, set up the matrix as above with central discretisation,
!          ! i.e. call SD to calculate the matrix without SD stabilisation.
!          ! Set up the SD structure for the creation of the defect.
!          ! There's not much to do, only initialise the viscosity...
!          rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter to 0 to deactivate the stabilisation.
!          rstreamlineDiffusion%dupsam = 0.0_DP
!
!          ! Matrix weight
!          rstreamlineDiffusion%ddelta = dgamma*dcx
!
!          ! Weight for the Newtop part; =0 deactivates Newton.
!          rstreamlineDiffusion%dnewton = dnewton*dcx
!
!          if (dnewton .eq. 0.0_DP) then
!
!            ! Deactivate the matrices A12 and A21 by setting the multiplicators
!            ! to 0.0. Whatever the content is (if there's content at all),
!            ! these matrices are ignored then by the kernel.
!
!            rtempMatrix%RmatrixBlock(1,2)%dscaleFactor = 0.0_DP
!            rtempMatrix%RmatrixBlock(2,1)%dscaleFactor = 0.0_DP
!
!          else
!
!            ! Clear A12/A21 that receives parts of the Newton matrix
!            call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(1,2))
!            call lsyssc_clearMatrix (rtempMatrix%RmatrixBlock(2,1))
!
!            ! Activate the submatrices A12 and A21 if they aren't.
!            rtempMatrix%RmatrixBlock(1,2)%dscaleFactor = 1.0_DP
!            rtempMatrix%RmatrixBlock(2,1)%dscaleFactor = 1.0_DP
!
!          end if
!
!          ! Call the SD method to calculate the nonlinearity.
!          call conv_streamlineDiffusionBlk2d (&
!                              rtempVector, rtempVector, &
!                              dvectorWeight, 0.0_DP,&
!                              rstreamlineDiffusion, CONV_MODDEFECT, &
!                              rtempMatrix,rsolution=rtempVector,rdefect=rtempDefect)
!
!          ! Set up the jump stabilisation structure.
!          ! There's not much to do, only initialise the viscosity...
!          rjumpStabil%dnu = rstreamlineDiffusion%dnu
!
!          ! Set stabilisation parameter
!          rjumpStabil%dgammastar = dupsam
!          rjumpStabil%dgamma = rjumpStabil%dgammastar
!
!          ! Matrix weight
!          rjumpStabil%dtheta = dgamma * dcx
!
!          ! Call the jump stabilisation technique to stabilise that stuff.
!          ! We can assemble the jump part any time as it's independent of any
!          ! convective parts...
!          call conv_jumpStabilisation2d (&
!                              rjumpStabil, CONV_MODDEFECT, &
!                              rtempMatrix%RmatrixBlock(1,1),&
!                              rsolution=rtempVector,rdefect=rtempDefect)
!
!          if (.not. bshared) then
!            print *,'Edge oriented stabilisation does not support independent A11/A22!'
!            stop
!          end if
!
!        case DEFAULT
!          print *,'Don''t know how to set up nonlinearity!?!'
!          stop
!
!        end select
!
!      else
!
!        ! That's the Stokes-case. Jump stabilisation is possible...
!
!        select case (iupwind)
!        case (2)
!          ! Jump stabilisation.
!
!          ! Set up the jump stabilisation structure.
!          ! There's not much to do, only initialise the viscosity...
!          rjumpStabil%dnu = rnonlinearSpatialMatrix%dnu
!
!          ! Set stabilisation parameter
!          rjumpStabil%dgammastar = dupsam
!          rjumpStabil%dgamma = rjumpStabil%dgammastar
!
!          ! Matrix weight
!          rjumpStabil%dtheta = dgamma * dcx
!
!          ! Call the jump stabilisation technique to stabilise that stuff.
!          ! We can assemble the jump part any time as it's independent of any
!          ! convective parts...
!          call conv_jumpStabilisation2d (&
!                              rjumpStabil, CONV_MODDEFECT, &
!                              rtempMatrix%RmatrixBlock(1,1),&
!                              rsolution=rtempVector,rdefect=rtempDefect)
!
!          if (.not. bshared) then
!            call conv_jumpStabilisation2d (&
!                                rjumpStabil, CONV_MODDEFECT, &
!                                rtempMatrix%RmatrixBlock(2,2),&
!                                rsolution=rtempVector,rdefect=rtempDefect)
!          end if
!
!        case DEFAULT
!          ! No stabilisation
!
!        end select
!
!      end if ! gamma <> 0
!
!      ! Release the temp matrix
!      call lsysbl_releaseMatrix (rtempMatrix)
!
!      ! Release the temp vector if allocated.
!      ! Derive a temporary vector that contains only those velocity
!      ! subvectors that might affect the matrix.
!      call lsysbl_releaseVector (rtempVector)
!      call lsysbl_releaseVector (rtempDefect)
!
!    end subroutine

!    subroutine assemblePrimalReactMassDefect (rnonlinearSpatialMatrix,&
!        rvector,rdefect,dcx,rvelocityVector,dvectorWeight)
!
!    ! Assembles the defect arising from the reactive coupling mass
!    ! matrices in the primal equation. rdefect must have been initialised
!    ! with the right hand side vector.
!    !
!    ! This special matrix is added in case that Newton is active. It has the
!    ! form
!    !        $$ M~ = dr_{11} N(\lambda) + dr_{12}   N*(\lambda) $$
!    !
!    ! the routine will construct
!    !
!    !       rdefect = r(primal)defect - dcx * (dtheta M~ r(dual)vector)
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how the mass matrix part is weighted (dr11, dr12).
!    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
!
!    ! Solution vector. Provides the dual equation for the assembly.
!    type(t_vectorBlock), intent(IN) :: rvector
!
!    ! On entry: RHS vector.
!    ! Is overwritten by the defect vector in the velocity subsystem.
!    type(t_vectorBlock), intent(INOUT) :: rdefect
!
!    ! Multiplication factor for the whole operator A*rvector
!    real(DP), intent(IN) :: dcx
!
!    ! Weight for the velocity vector rvelocityVector; usually = 1.0
!    real(DP), intent(IN) :: dvectorWeight
!
!    ! Velocity vector field that should be used for the assembly of the
!    ! nonlinearity. Block 4 and 5 in that block vector are used as velocity
!    ! field.
!    type(t_vectorBlock), intent(IN) :: rvelocityVector
!
!      ! local variables
!      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
!      type(t_matrixBlock) :: rtempmatrix
!      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
!
!      ! If we have a reactive coupling mass matrix, it gets interesting...
!      if ((rnonlinearSpatialMatrix%dr11 .ne. 0.0_DP) .or. &
!          (rnonlinearSpatialMatrix%dr12 .ne. 0.0_DP)) then
!
!        call lsysbl_createEmptyMatrix (rtempMatrix,2)
!
!        ! Create a matrix with the structure we need. Share the structure
!        ! of the mass matrix. Entries are not necessary for the assembly
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        ! Switch the matrices on.
!        rtempMatrix%RmatrixBlock(1:2,1:2)%dscaleFactor = 1.0_DP
!
!        call lsysbl_updateMatStrucInfo (rtempMatrix)
!
!        ! The reactive part is: "dr1 * . * grad(lambda)".
!        ! This is exactly the 'Newton' part assembled by the streamline diffusion
!        ! method if we use the dual velocity as velocity field!
!        ! So prepare to call streamline diffusion.
!
!        ! Viscosity; ok, actually not used.
!        rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!        ! Set stabilisation parameter
!        rstreamlineDiffusion%dupsam = rnonlinearSpatialMatrix%dupsam2
!
!        ! Weight dr21 of the convective part.
!        rstreamlineDiffusion%ddelta = rnonlinearSpatialMatrix%dr11*dcx
!
!        ! Weight for the Newton part; here, this is the dr22 weight.
!        rstreamlineDiffusion%dnewton = rnonlinearSpatialMatrix%dr12*dcx
!
!        ! Create a temporary block vector that points to the dual velocity.
!        ! This has to be evaluated during the assembly.
!        call lsysbl_deriveSubvector (rvelocityVector,rtempvectorEval,4,5,.true.)
!
!        ! Create a temporary block vector for the dual defect.
!        ! Matrix*primal velocity is subtracted from this.
!        call lsysbl_deriveSubvector (rdefect,rtempvectorDef,1,2,.true.)
!
!        ! Call the SD method to calculate the nonlinearity.
!        ! As velocity vector, specify rtempvector, which points to
!        ! the dual velocity.
!        call conv_streamlineDiffusionBlk2d (&
!                            rtempvectorEval, rtempvectorEval, &
!                            dvectorWeight, 0.0_DP,&
!                            rstreamlineDiffusion, CONV_MODDEFECT, &
!                            rtempMatrix,rsolution=rvector,rdefect=rtempvectorDef)
!
!        ! Release the temp matrix and temp vector.
!        call lsysbl_releaseVector (rtempVectorEval)
!        call lsysbl_releaseVector (rtempVectorDef)
!        call lsysbl_releaseMatrix (rtempMatrix)
!
!      end if
!
!    end subroutine
!
!    subroutine assembleDualReactMassDefect (rnonlinearSpatialMatrix,&
!        rvector,rdefect,dcx,rvelocityVector,dvectorWeight)
!
!    ! Assembles the defect arising from the reactive coupling mass
!    ! matrices. rdefect must have been initialised with the right hand side
!    ! vector.
!    !
!    ! This special matrix is added in case that Newton is active. It has the
!    ! form
!    !        $$ M~ = dr_{21} N(\lambda) + dr_{22}   N*(\lambda) $$
!    !
!    ! the routine will construct
!    !
!    !       rdefect = r(dual)defect - dcx * (dtheta M~ r(primal)vector)
!
!    ! A t_nonlinearSpatialMatrix structure providing all necessary 'source' information
!    ! about how the mass matrix part is weighted (dr21, dr22).
!    type(t_nonlinearSpatialMatrix), intent(IN) :: rnonlinearSpatialMatrix
!
!    ! Solution vector. Provides the dual equation for the assembly.
!    type(t_vectorBlock), intent(IN) :: rvector
!
!    ! On entry: RHS vector.
!    ! Is overwritten by the defect vector in the velocity subsystem.
!    type(t_vectorBlock), intent(INOUT) :: rdefect
!
!    ! Multiplication factor for the whole operator A*rvector
!    real(DP), intent(IN) :: dcx
!
!    ! Weight for the velocity vector rvelocityVector; usually = 1.0
!    real(DP), intent(IN) :: dvectorWeight
!
!    ! Velocity vector field that should be used for the assembly of the
!    ! nonlinearity. The first two blocks in that block vector are
!    ! used as velocity field.
!    type(t_vectorBlock), intent(IN) :: rvelocityVector
!
!      ! local variables
!      type(t_convStreamlineDiffusion) :: rstreamlineDiffusion
!      type(t_matrixBlock) :: rtempmatrix
!      type(t_vectorBlock) :: rtempvectorEval,rtempVectorDef
!
!      ! If we have a reactive coupling mass matrix, it gets interesting...
!      if ((rnonlinearSpatialMatrix%dr21 .ne. 0.0_DP) .or. &
!          (rnonlinearSpatialMatrix%dr22 .ne. 0.0_DP)) then
!
!        call lsysbl_createEmptyMatrix (rtempMatrix,2)
!
!        ! Create a matrix with the structure we need. Share the structure
!        ! of the mass matrix. Entries are not necessary for the assembly
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        call lsyssc_duplicateMatrix (rnonlinearSpatialMatrix%p_rmatrixMass,&
!            rtempMatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!
!        ! Switch the matrices on.
!        rtempMatrix%RmatrixBlock(1:2,1:2)%dscaleFactor = 1.0_DP
!
!        call lsysbl_updateMatStrucInfo (rtempMatrix)
!
!        ! The reactive part is: "dr2 * . * grad(lambda)".
!        ! This is exactly the 'Newton' part assembled by the streamline diffusion
!        ! method if we use the dual velocity as velocity field!
!        ! So prepare to call streamline diffusion.
!
!        ! Viscosity; ok, actually not used.
!        rstreamlineDiffusion%dnu = rnonlinearSpatialMatrix%dnu
!
!        ! Set stabilisation parameter
!        rstreamlineDiffusion%dupsam = rnonlinearSpatialMatrix%dupsam2
!
!        ! Weight dr21 of the convective part.
!        rstreamlineDiffusion%ddelta = rnonlinearSpatialMatrix%dr21*dcx
!
!        ! Weight for the Newton part; here, this is the dr22 weight.
!        rstreamlineDiffusion%dnewton = rnonlinearSpatialMatrix%dr22*dcx
!
!        ! Create a temporary block vector that points to the dual velocity.
!        ! This has to be evaluated during the assembly.
!        call lsysbl_deriveSubvector (rvelocityVector,rtempvectorEval,4,5,.true.)
!
!        ! Create a temporary block vector for the dual defect.
!        ! Matrix*primal velocity is subtracted from this.
!        call lsysbl_deriveSubvector (rdefect,rtempvectorDef,4,5,.true.)
!
!        ! Call the SD method to calculate the nonlinearity.
!        ! As velocity vector, specify rtempvector, which points to
!        ! the dual velocity.
!        call conv_streamlineDiffusionBlk2d (&
!                            rtempvectorEval, rtempvectorEval, &
!                            dvectorWeight, 0.0_DP,&
!                            rstreamlineDiffusion, CONV_MODDEFECT, &
!                            rtempMatrix,rsolution=rvector,rdefect=rtempvectorDef)
!
!        ! Release the temp matrix and temp vector.
!        call lsysbl_releaseVector (rtempVectorEval)
!        call lsysbl_releaseVector (rtempVectorDef)
!        call lsysbl_releaseMatrix (rtempMatrix)
!
!      end if
!
!    end subroutine

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_projectControlTimestep (rdualSolution,dumin,dumax)

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

end module

!##############################################################################
!# ****************************************************************************
!# <name> afcstabilisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for performing
!# discrete stabilisation by means of algebraic flux correction
!#
!# The following routines are available:
!#
!# 1.) afcstab_initFromParameterlist
!#     -> Creates a stabilisation structure and initialize
!#        it from the values of a given parameter list
!#
!# 2.) afcstab_releaseStabilisation
!#     -> Releases a stabilisation structure
!#
!# 3.) afcstab_resizeStabilisation = afcstab_resizeStabDirect /
!#                                   afcstab_resizeStabIndScalar
!#                                   afcstab_resizeStabIndBlock
!#     -> Resizes a stabilisation structure
!#
!# 4.) afcstab_copyStabilisation
!#     -> Copies a stabilisation structure to another structure
!#
!# 5.) afcstab_duplicateStabilisation
!#     -> Duplicate (parts of) a stabilisation structure
!#
!# 6.) afcstab_initMatrixCoeffs
!#     -> Initialises the auxiliary matrix coefficients
!#
!# 7.) afcstab_copyMatrixCoeffs
!#     -> Copies auxiliary matrix coefficients into stabilisation structure
!#
!# 8.) afcstab_isMatrixCompatible = afcstab_isMatrixCompatibleSc /
!#                                  afcstab_isMatrixCompatibleBl
!#     -> Checks wether a matrix and a stabilisation structure are compatible
!#
!# 9.) afcstab_isVectorCompatible = afcstab_isVectorCompatibleSc /
!#                                  afcstab_isVectorCompatibleBl
!#     -> Checks wether a vector and a stabilisation structure are compatible
!#
!# 10.) afcstab_getbase_array = afcstab_getbase_arrayScalar /
!#                              afcstab_getbase_arrayBlock
!#      -> Returns the array of pointers to a given block matrix
!#
!# 11.) afcstab_getbase_IvertAtEdgeIdx
!#      -> Returns pointer to the index pointer for the
!#         vertices at edge structure
!#
!# 12.) afcstab_getbase_IverticesAtEdge
!#      -> Returns pointer to the vertices at edge structure
!#
!# 13.) afcstab_getbase_IsupdiagEdgeIdx
!#      -> Returns pointer to the index pointer for the
!#         superdiagonal edge numbers
!#
!# 14.) afcstab_getbase_IsubdiagEdgeIdx
!#      -> Returns pointer to the index pointer for the
!#         subdiagonal edge numbers
!#
!# 15.) afcstab_getbase_IsubdiagEdge
!#      -> Returns pointer to the subdiagonal edge numbers
!#
!# 16.) afcstab_getbase_DcoeffsAtEdge = afcstab_getbase_DcoeffsAtEdge /
!#                                      afcstab_getbase_FcoeffsAtEdge
!#      -> Returns pointer to edge data
!#
!# 17.) afcstab_getbase_DmatCoeffAtNode = afcstab_getbase_DmatCoeffAtNode /
!#                                        afcstab_getbase_FmatCoeffAtNode
!#      -> Returns pointer to the diagonal entries
!#         of the auxiiliary constant matrix coefficients
!#
!# 18.) afcstab_getbase_DmatCoeffAtEdge = afcstab_getbase_DmatCoeffAtEdge /
!#                                        afcstab_getbase_FmatCoeffAtEdge
!#      -> Returns pointer to the off-diagonal entries
!#         of the auxiliary constant matrix coefficients
!#
!# 19.) afcstab_getbase_DboundsAtEdge = afcstab_getbase_DboundsAtEdge /
!#                                      afcstab_getbase_FboundsAtEdge
!#      -> Returns pointer to the bounds at edges
!#
!# 20.) afcstab_generateVerticesAtEdge
!#      -> Generates the standard edge data structure
!#
!# 21.) afcstab_generateOffdiagEdges
!#      -> Generates the subdiagonal edge data structure
!#
!# 22.) afcstab_generateExtSparsity
!#      -> Generates the extended sparsity pattern
!#
!# 23.) afcstab_failsafeLimiting = afcstab_failsafeLimitingBlock /
!#                                 afcstab_failsafeLimitingArray
!#      -> Perform failsafe flux correction
!#
!# 24.) afcstab_copyH2D_IverticesAtEdge
!#      -> Copies the vertices at edge structure from the host memory
!#         to the device memory.
!#
!# 25.) afcstab_copyD2H_IverticesAtEdge
!#      -> Copies the vertices at edge structure from the device memory
!#         to the host memory.
!#
!# 26.) afcstab_copyH2D_DmatCoeffAtEdge
!# -> Copies the off-diagonal entries of the auxiliary constant matrix
!#         coefficients from the host memory to the device memory.
!#
!# 27.) afcstab_copyD2H_DmatCoeffAtEdge
!#      -> Copies the off-diagonal entries of the auxiliary constant matrix
!#         coefficients from the device memory to the host memory.
!#
!# 28.) afcstab_allocEdgeStructure
!#      -> Allocates the edge data structure
!#
!# 29.) afcstab_allocCoeffsAtEdge
!#      -> Allocates the coefficients at edge data structure
!#
!# 30.) afcstab_allocVectorsPQR
!#      -> Allocates the nodal vectors P, Q, and R each for '+' and '-'
!#
!# 31.) afcstab_allocFlux
!#      -> Allocates the edge-wise flux vector flux
!#
!# 32.) afcstab_allocFlux0
!#      -> Allocates the edge-wise flux vector flux0
!#
!# 33.) afcstab_allocFluxPrel
!#      -> Allocates the edge-wise flux vector fluxPrel
!#
!# 34.) afcstab_allocAlpha
!#      -> Allocates the edge-wise correction factors alpha
!#
!# 35.) afcstab_allocBoundsAtEdge
!#      -> Allocates the bounds at edge data structure
!#
!# 36.) afcstab_limit = afcstab_limit_unbounded/
!#                      afcstab_limit_bounded
!#      -> Compute the nodal correction factors, i.e., the ratio of
!#         admissible solution increments and raw antidiffusion
!#
!# 37.) afcstab_buildBoundsLPT = afcstab_buildBoundsLPT1D /
!#                               afcstab_buildBoundsLPT2D /
!#                               afcstab_buildBoundsLPT3D
!#          
!#
!# </purpose>
!##############################################################################
module afcstabilisation

  use basicgeometry
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none
  
  private
  public :: t_afcstab, t_array

  public :: afcstab_initFromParameterlist
  public :: afcstab_releaseStabilisation
  public :: afcstab_resizeStabilisation
  public :: afcstab_copyStabilisation
  public :: afcstab_duplicateStabilisation
  public :: afcstab_initMatrixCoeffs
  public :: afcstab_copyMatrixCoeffs
  public :: afcstab_isMatrixCompatible
  public :: afcstab_isVectorCompatible

  public :: afcstab_getbase_array
  public :: afcstab_getbase_IvertAtEdgeIdx
  public :: afcstab_getbase_IverticesAtEdge
  public :: afcstab_getbase_IsupdiagEdgeIdx
  public :: afcstab_getbase_IsubdiagEdgeIdx
  public :: afcstab_getbase_IsubdiagEdge
  public :: afcstab_getbase_DcoeffsAtEdge
  public :: afcstab_getbase_DmatCoeffAtNode
  public :: afcstab_getbase_DmatCoeffAtEdge
  public :: afcstab_getbase_DboundsAtEdge

  public :: afcstab_generateVerticesAtEdge
  public :: afcstab_generateOffdiagEdges
  public :: afcstab_generateExtSparsity

  public :: afcstab_copyD2H_IverticesAtEdge
  public :: afcstab_copyD2H_DmatCoeffAtEdge
  public :: afcstab_copyH2D_IverticesAtEdge
  public :: afcstab_copyH2D_DmatCoeffAtEdge

  public :: afcstab_failsafeLimiting
  public :: afcstab_buildBoundsLPT

  public :: afcstab_allocEdgeStructure
  public :: afcstab_allocCoeffsAtEdge
  public :: afcstab_allocBoundsAtEdge
  public :: afcstab_allocVectorsPQR
  public :: afcstab_allocFluxPrel
  public :: afcstab_allocFlux0
  public :: afcstab_allocFlux
  public :: afcstab_allocAlpha
  
  public :: afcstab_limit
 
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
!<constantblock description="Global format flags for AFC stabilisation">

  ! ***** Stabilisation group 0x: no high-resolution stabilisation *****

  ! No stabilisation: use standard high-order Galerkin discretisation
  integer, parameter, public :: AFCSTAB_GALERKIN              = 0
  
  ! Stabilisation of discrete upwind type for non-symmetric operators
  integer, parameter, public :: AFCSTAB_UPWIND                = 1

  ! Stabilisation of discrete maximum principle preserving type for
  ! symmetric operators
  integer, parameter, public :: AFCSTAB_DMP                   = 2


  ! ***** Stabilisation group 1x: nonlinear FEM-FCT stabilisation *****

  ! Stabilisation of semi-explicit FEM-FCT type
  integer, parameter, public :: AFCSTAB_NLINFCT_EXPLICIT      = 10

  ! Stabilisation of semi-implicit FEM-FCT type
  integer, parameter, public :: AFCSTAB_NLINFCT_IMPLICIT      = 11

  ! Stabilisation of iterative FEM-FCT type
  integer, parameter, public :: AFCSTAB_NLINFCT_ITERATIVE     = 12


  ! ***** Stabilisation group 2x: linearised FEM-FCT stabilisation *****

  ! Stabilisation of linearised FEM-FCT type
  integer, parameter, public :: AFCSTAB_LINFCT                = 20

  ! Stabilisation of linearised FEM-FCT type for mass antidiffusion
  integer, parameter, public :: AFCSTAB_LINFCT_MASS           = 21


  ! ***** Stabilisation group 3x: upwind-biased FEM-TVD stabilisation *****
  
  ! Stabilisation of FEM-TVD type
  integer, parameter, public :: AFCSTAB_TVD                   = 30
  
  ! Stabilisation of general purpose type
  integer, parameter, public :: AFCSTAB_GP                    = 31

  ! Stabilisation of symmetric type for diffusion operators
  integer, parameter, public :: AFCSTAB_SYMMETRIC             = 32

  
  ! ***** Stabilisation group 4x: nonlinear FEM-LP stabilisation *****

  ! Stabilisation of linearity-preserving flux correction type
  ! for symmetric mass antidiffusion
  integer, parameter, public :: AFCSTAB_NLINLPT_MASS          = 40

  ! Stabilisation of linearity-preserving flux correction type
  ! for upwind-biased antidiffusion (e.g., convection)
  integer, parameter, public :: AFCSTAB_NLINLPT_UPWINDBIASED  = 41

  ! Stabilisation of linearity-preserving flux correction type
  ! for symmetric antidiffusion (e.g., diffusion)
  integer, parameter, public :: AFCSTAB_NLINLPT_SYMMETRIC     = 42
  
  
  ! ***** Stabilisation group 5x: linearised FEM-LP stabilisation *****

  ! Stabilisation of linearised linearity-preserving flux correction
  ! type for symmetric mass antidiffusion
  integer, parameter, public :: AFCSTAB_LINLPT_MASS           = 50

  ! Stabilisation of linearised linearity-preserving flux correction
  ! type for upwind-biased antidiffusion (e.g., convection)
  integer, parameter, public :: AFCSTAB_LINLPT_UPWINDBIASED   = 51

  ! Stabilisation of linearised linearity-preserving flux correction
  ! type for symmetric antidiffusion (e.g., diffusion)
  integer, parameter, public :: AFCSTAB_LINLPT_SYMMETRIC      = 52
  
!</constantblock>


!<constantblock description="Global format flags for prelimiting">

  ! No prelimiting
  integer, parameter, public :: AFCSTAB_PRELIMITING_NONE      = 0

  ! Standard prelimiting
  integer, parameter, public :: AFCSTAB_PRELIMITING_STD       = 1

  ! Minmod prelimiting
  integer, parameter, public :: AFCSTAB_PRELIMITING_MINMOD    = 2

!</constantblock>


!<constantblock description="Global format flags for limiting">

  ! No limiting at all
  integer, parameter, public :: AFCSTAB_LIMITING_NONE             = 0

  ! Symmetric limiting
  integer, parameter, public :: AFCSTAB_LIMITING_SYMMETRIC        = 1

  ! Upwind-biased limiting (implies edge-orientation)
  integer, parameter, public :: AFCSTAB_LIMITING_UPWINDBIASED     = 2

  ! Characteristic limiting
  integer, parameter, public :: AFCSTAB_LIMITING_CHARACTERISTIC   = 3

!</constantblock>


!<constantblock description="Bitfield identifiers for properties of stabilisation">

  ! Stabilisation is undefined
  integer(I32), parameter, public :: AFCSTAB_UNDEFINED            = 2_I32**0

  ! Stabilisation has been initialised
  integer(I32), parameter, public :: AFCSTAB_INITIALISED          = 2_I32**1

  ! Edge-based structure has been generated: IverticesAtEdge
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGESTRUCTURE    = 2_I32**2

  ! Edge-based structure has been oriented: IverticesAtEdge
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEORIENTATION  = 2_I32**3

  ! Edge-based values have been computed from matrix: DcoefficientsAtEdge
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEVALUES       = 2_I32**4

  ! Subdiagonal edge-based structure has been generated
  integer(I32), parameter, public :: AFCSTAB_HAS_OFFDIAGONALEDGES = 2_I32**5
  
  ! Antidiffusive fluxes have been precomputed
  integer(I32), parameter, public :: AFCSTAB_HAS_ADFLUXES         = 2_I32**6
  
  ! Nodal sums of antidiffusive increments have been computed: PP, PM
  integer(I32), parameter, public :: AFCSTAB_HAS_ADINCREMENTS     = 2_I32**7

  ! Nodal upper/lower bounds have been computed: QP, QM
  integer(I32), parameter, public :: AFCSTAB_HAS_NODEBOUNDS       = 2_I32**8
  
  ! Nodal correction factors have been computed: RP, RM
  integer(I32), parameter, public :: AFCSTAB_HAS_NODELIMITER      = 2_I32**9
  
  ! Edge-wise correction factors have been computed: ALPHA
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGELIMITER      = 2_I32**10

  ! Low-order predictor has been computed
  integer(I32), parameter, public :: AFCSTAB_HAS_PREDICTOR        = 2_I32**11

  ! Auxiliary matrix coefficients have been attached
  integer(I32), parameter, public :: AFCSTAB_HAS_MATRIXCOEFFS     = 2_I32**12

  ! Transformed nodal solution values have been computed
  integer(I32), parameter, public :: AFCSTAB_HAS_NODEVALUES       = 2_I32**13

  ! Slope-based bounds have been computed
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEBOUNDS       = 2_I32**14

!</constantblock>


!<constantblock description="Duplication flags. Specifies which information is duplicated">
  
  ! Duplicate atomic stabilisation structure
  integer(I32), parameter, public :: AFCSTAB_DUP_STRUCTURE        = 2_I32**1

  ! Duplicate edge-based structure: IverticesAtEdge
  integer(I32), parameter, public :: AFCSTAB_DUP_EDGESTRUCTURE    = AFCSTAB_HAS_EDGESTRUCTURE

  ! Duplicate edge-based values: DcoefficientsAtEdge
  integer(I32), parameter, public :: AFCSTAB_DUP_EDGEVALUES       = AFCSTAB_HAS_EDGEVALUES 

  ! Duplicate subdiagonal edge-based
  integer(I32), parameter, public :: AFCSTAB_DUP_OFFDIAGONALEDGES = AFCSTAB_HAS_OFFDIAGONALEDGES
  
  ! Duplicate antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_DUP_ADFLUXES         = AFCSTAB_HAS_ADFLUXES 
  
  ! Duplicate nodal sums of antidiffusive increments: PP, PM
  integer(I32), parameter, public :: AFCSTAB_DUP_ADINCREMENTS     = AFCSTAB_HAS_ADINCREMENTS

  ! Duplicate nodal upper/lower bounds: QP, QM
  integer(I32), parameter, public :: AFCSTAB_DUP_NODEBOUNDS       = AFCSTAB_HAS_NODEBOUNDS
  
  ! Duplicate nodal correction factors: RP, RM
  integer(I32), parameter, public :: AFCSTAB_DUP_NODELIMITER      = AFCSTAB_HAS_NODELIMITER
  
  ! Duplicate edge-wise correction factors: ALPHA
  integer(I32), parameter, public :: AFCSTAB_DUP_EDGELIMITER      = AFCSTAB_HAS_EDGELIMITER

  ! Duplicate low-order predictor
  integer(I32), parameter, public :: AFCSTAB_DUP_PREDICTOR        = AFCSTAB_HAS_PREDICTOR
    
  ! Duplicate auxiliary matrix coefficients
  integer(I32), parameter, public :: AFCSTAB_DUP_MATRIXCOEFFS     = AFCSTAB_HAS_MATRIXCOEFFS

  ! Duplicate transformed nodal solution values
  integer(I32), parameter, public :: AFCSTAB_DUP_NODEVALUES       = AFCSTAB_HAS_NODEVALUES

  ! Duplicate slope-based bounds
  integer(I32), parameter, public :: AFCSTAB_DUP_EDGEBOUNDS       = AFCSTAB_HAS_EDGEBOUNDS

!</constantblock>


!<constantblock description="Duplication flags. Specifies which information is shared \
!                            between stabilisation structures">
  
  ! Duplicate atomic stabilisation structure
  integer(I32), parameter, public :: AFCSTAB_SHARE_STRUCTURE        = AFCSTAB_DUP_STRUCTURE

  ! Share edge-based structure: IverticesAtEdge
  integer(I32), parameter, public :: AFCSTAB_SHARE_EDGESTRUCTURE    = AFCSTAB_DUP_EDGESTRUCTURE

  ! Share edge-based values: DcoefficientsAtEdge
  integer(I32), parameter, public :: AFCSTAB_SHARE_EDGEVALUES       = AFCSTAB_DUP_EDGEVALUES

  ! Share subdiagonal edge-based structure
  integer(I32), parameter, public :: AFCSTAB_SHARE_OFFDIAGONALEDGES = AFCSTAB_DUP_OFFDIAGONALEDGES
  
  ! Share antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_SHARE_ADFLUXES         = AFCSTAB_DUP_ADFLUXES
  
  ! Share nodal sums of antidiffusive increments: PP, PM
  integer(I32), parameter, public :: AFCSTAB_SHARE_ADINCREMENTS     = AFCSTAB_DUP_ADINCREMENTS

  ! Share nodal upper/lower bounds: QP, QM
  integer(I32), parameter, public :: AFCSTAB_SHARE_NODEBOUNDS       = AFCSTAB_DUP_NODEBOUNDS
  
  ! Share nodal correction factors: RP, RM
  integer(I32), parameter, public :: AFCSTAB_SHARE_NODELIMITER      = AFCSTAB_DUP_NODELIMITER
  
  ! Share edge-wise correction factors: ALPHA
  integer(I32), parameter, public :: AFCSTAB_SHARE_EDGELIMITER      = AFCSTAB_DUP_EDGELIMITER

  ! Share low-order predictor
  integer(I32), parameter, public :: AFCSTAB_SHARE_PREDICTOR        = AFCSTAB_DUP_PREDICTOR

  ! Share auxiliary matrix coefficients
  integer(I32), parameter, public :: AFCSTAB_SHARE_MATRIXCOEFFS     = AFCSTAB_DUP_MATRIXCOEFFS

  ! Share transformed nodal solution values
  integer(I32), parameter, public :: AFCSTAB_SHARE_NODEVALUES       = AFCSTAB_DUP_NODEVALUES

  ! Share slope-based bounds
  integer(I32), parameter, public :: AFCSTAB_SHARE_EDGEBOUNDS       = AFCSTAB_DUP_EDGEBOUNDS

!</constantblock>


!<constantblock description="Bitfield identifiers for TVD-algorithm">

  ! Compute the raw-antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_ADFLUXES     = 2_I32**0

  ! Compute the sums of antidiffusive increments
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_ADINCREMENTS = 2_I32**1

  ! Compute the local solution bounds
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_BOUNDS       = 2_I32**2
  
  ! Compute the nodal correction factors
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_LIMIT        = 2_I32**3

  ! Correct raw antidiffusive fluxes and apply them
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_CORRECT      = 2_I32**4

  ! FEM-TVD algorithm without application of the corrected fluxes
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_PREPARE   = AFCSTAB_TVDALGO_ADFLUXES +&
                                                                 AFCSTAB_TVDALGO_ADINCREMENTS +&
                                                                 AFCSTAB_TVDALGO_BOUNDS +&
                                                                 AFCSTAB_TVDALGO_LIMIT

  ! Standard FEM-TVD algorithm
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_STANDARD = AFCSTAB_TVDALGO_PREPARE +&
                                                                AFCSTAB_TVDALGO_CORRECT

!</constantblock>


!<constantblock description="Bitfield identifiers for FCT-fluxes">

   ! Compute explicit part of raw-antidiffusive fluxes
   integer(I32), parameter, public :: AFCSTAB_FCTFLUX_EXPLICIT    = 2_I32**0

   ! Compute implicit part of raw-antidiffusive fluxes
   integer(I32), parameter, public :: AFCSTAB_FCTFLUX_IMPLICIT    = 2_I32**1

   ! Apply rejected antidiffusives fluxes to the implicit part
   integer(I32), parameter, public :: AFCSTAB_FCTFLUX_REJECTED    = 2_I32**2

!</constantblock>


!<constantblock description="Bitfield identifiers for FCT-algorithm">

  ! Initialize the edgewise correction factors by unity
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_INITALPHA    = 2_I32**0

  ! Prelimit the raw antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_PRELIMIT     = 2_I32**1
  
  ! Compute the sums of antidiffusive increments
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_ADINCREMENTS = 2_I32**2

  ! Compute the distances to a local extremum
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_BOUNDS       = 2_I32**3

  ! Compute the nodal correction factors
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_LIMITNODAL   = 2_I32**4

  ! Compute edgewise correction factors
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_LIMITEDGE    = 2_I32**5
  
  ! Correct raw antidiffusive fluxes and apply them
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_CORRECT      = 2_I32**6

  ! Scale corrected antidiffusive fluxes by the inverse of the
  ! lumped mass matrix prior to applying it to the residual/solution
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_SCALEBYMASS  = 2_I32**7

  ! Constrain the raw antidiffusive fluxes by the size of limited
  ! limited rather than by the correction factors directly
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_CONSTRAIN    = 2_I32**8
  
  ! FEM-FCT algorithm without application of the corrected fluxes
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_PREPARE  = AFCSTAB_FCTALGO_INITALPHA +&
                                                                AFCSTAB_FCTALGO_PRELIMIT +&
                                                                AFCSTAB_FCTALGO_ADINCREMENTS +&
                                                                AFCSTAB_FCTALGO_BOUNDS +&
                                                                AFCSTAB_FCTALGO_LIMITNODAL +&
                                                                AFCSTAB_FCTALGO_LIMITEDGE

  ! Standard FEM-FCT algorithm
  ! (also used in the very first iteration of all other variants)
  integer(I32), parameter, public :: AFCSTAB_FCTALGO_STANDARD = AFCSTAB_FCTALGO_PREPARE +&
                                                                AFCSTAB_FCTALGO_CORRECT
  
!</constantblock>


!<constantblock description="Bitfield identifiers for LPT-fluxes">

  ! Compute raw-antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_LPTFLUX          = 2_I32**0

  ! Compute bounds
  integer(I32), parameter, public :: AFCSTAB_LPTFLUX_BOUNDS   = 2_I32**1

!</constantblock>


!<constantblock description="Bitfield identifiers for LPT-algorithm">

  ! Initialize the edgewise correction factors by unity
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_INITALPHA    = 2_I32**0

  ! Compute the sums of antidiffusive increments
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_ADINCREMENTS = 2_I32**1

  ! Compute the distances to a local extremum
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_BOUNDS       = 2_I32**2

  ! Compute the nodal correction factors
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_LIMITNODAL   = 2_I32**3

  ! Compute edgewise correction factors
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_LIMITEDGE    = 2_I32**4

  ! Correct raw antidiffusive fluxes and apply them
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_CORRECT      = 2_I32**5

  ! Scale corrected antidiffusive fluxes by the inverse of the
  ! lumped mass matrix prior to applying it to the residual/solution
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_SCALEBYMASS  = 2_I32**6

  ! FEM-LPT algorithm without application of the corrected fluxes
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_PREPARE  = AFCSTAB_LPTALGO_INITALPHA +&
                                                                AFCSTAB_LPTALGO_ADINCREMENTS +&
                                                                AFCSTAB_LPTALGO_BOUNDS +&
                                                                AFCSTAB_LPTALGO_LIMITNODAL +&
                                                                AFCSTAB_LPTALGO_LIMITEDGE

  ! Standard FEM-LPT algorithm
  integer(I32), parameter, public :: AFCSTAB_LPTALGO_STANDARD = AFCSTAB_LPTALGO_PREPARE +&
                                                                AFCSTAB_LPTALGO_CORRECT

!</constantblock>


!<constantblock description="Default tolerances for stabilisation">
  
  ! Absolute tolerance for prelimiting of antidiffusive fluxes
#ifndef AFCSTAB_PRELIMABS
  real(DP), parameter, public :: AFCSTAB_PRELIMABS = 1e-16
#endif

  ! Relative tolerance for prelimiting of antidiffusive fluxes
#ifndef AFCSTAB_PRELIMREL
  real(DP), parameter, public :: AFCSTAB_PRELIMREL = 1e-6
#endif

  ! Absolute tolerance for stabilisation
#ifndef AFCSTAB_EPSABS
  real(DP), parameter, public :: AFCSTAB_EPSABS = 1e-16
#endif 

  ! Relative tolerance for stabilisation
#ifndef AFCSTAB_EPSREL
  real(DP), parameter, public :: AFCSTAB_EPSREL = 1e-3
#endif

!</constantblock>

!<constantblock>
  ! Minimum number of edges for OpenMP parallelisation: If the number of
  ! edges is below this value, then no parallelisation is performed.
#ifndef AFCSTAB_NEDGEMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: AFCSTAB_NEDGEMIN_OMP = 1000
#else
  integer, public            :: AFCSTAB_NEDGEMIN_OMP = 1000
#endif
#endif
  
!</constantblock>

!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<types>
!<typeblock>

  ! This structures holds all required information for stabilisation
  ! of algebraic flux correction type. It is applicable to scalar
  ! problems and systems of equations alike. Depending on the
  ! stabilisation strategy some internal structures will be generated.
  
  type t_afcstab
    
    ! Format Tag: Identifies the type of stabilisation
    integer :: ctypeAFCstabilisation = AFCSTAB_GALERKIN

    ! Format Tag: Identifies the type of prelimiting
    ! Can take one of the AFCSTAB_PRELIMITING_xxxx flags.
    integer :: ctypePrelimiting = AFCSTAB_PRELIMITING_NONE

    ! Format Tag: Identifies the type of limiting
    integer :: ctypeLimiting = AFCSTAB_LIMITING_NONE

    ! Format Tag: Identifies the format of the scalar template matrix
    ! underlying the edge-based sparsity structure (if any)
    integer :: cmatrixFormat = LSYSSC_MATRIXUNDEFINED

    ! Format Tag: Identifies the data type
    integer :: cdataType = ST_DOUBLE
    
    ! Duplication Flag: This is a bitfield coming from an OR
    ! combination of different AFCSTAB_SHARE_xxxx constants and
    ! specifies which parts of the stabilisation are shared with
    ! another stabilisation structure.
    integer(I32) :: iduplicationFlag = AFCSTAB_UNDEFINED

    ! Specification Flag: Specifies the stabilisation. This is a bitfield
    ! coming from an OR combination of different AFCSTAB_HAS_xxxx
    ! constants and specifies various properties of the stabilisation
    ! structure.
    integer(I32) :: istabilisationSpec = AFCSTAB_UNDEFINED

    ! Number of equations of the sparsity pattern
    integer :: NEQ = 0

    ! Number of local variables; in general scalar solution vectors of
    ! size NEQ posses NEQ entries. However, scalar vectors can be
    ! interleaved, that is, each of the NEQ entries stores NVAR local
    ! variables. In this case, NEQ remains unmodified but NVAR>1 such
    ! that the physical length of the vector is NEQ*NVAR.
    integer :: NVAR = 1

    ! Number of local variables after transformation; this variable is
    ! similar to NVAR except for the fact that it corresponds to the
    ! transformed variable. Flux correction for systems of equations
    ! may require some sort of synchronisation which may be performed
    ! in terms of, e.g., primitive variables whereas the system
    ! matrices and the residual vector/right-hand side are assembled
    ! in terms of conservative variables. 
    integer :: NVARtransformed = 1

    ! Number of edges of the sparsity pattern
    integer :: NEDGE = 0

    ! Maximum number of edges adjacent to one vertex. This
    ! corresponds to the maximum number of nonzero row entries.
    integer :: NNVEDGE = 0

    ! Number of different matrix coefficients stored in
    ! DmatrixCoeffsAtNode and DmatrixCoeffsAtEdge, respectively.
    integer :: nmatrixCoeffs = 0

    ! Number of different coefficients stored in DCoefficientsAtEdge
    integer :: ncoeffs = 0

    ! Handle to index pointer for edge structure
    ! The numbers IverticesAtEdgeIdx(k):IverticesAtEdgeIdx(k+1)-1
    ! denote the edge numbers of the k-th group of edges.
    integer :: h_IverticesAtEdgeIdx = ST_NOHANDLE

    ! Handle to vertices at edge structure
    ! IverticesAtEdge(1:2,1:NEDGE) : the two end-points of the edge
    ! IverticesAtEdge(3:4,1:NEDGE) : the two matrix position that
    !                                correspond to the edge
    integer :: h_IverticesAtEdge = ST_NOHANDLE

    ! Handle to index pointer for superdiagonal edge numbers
    ! The numbers IsuperdiagEdgesIdx(i):IsuperdiagEdgesIdx(i+1)-1
    ! denote the edge numbers of the ith vertex which are located in
    ! the upper right triangular matrix.
    integer :: h_IsuperdiagEdgesIdx = ST_NOHANDLE

    ! Handle to index pointer for subdiagonal edge numbers
    integer :: h_IsubdiagEdgesIdx = ST_NOHANDLE

    ! Handle to the subdiagonal edge numbers
    integer :: h_IsubdiagEdges = ST_NOHANDLE

    ! Handle to coefficient at edge structure
    integer :: h_DcoefficientsAtEdge = ST_NOHANDLE

    ! Handle to auxiliary matrix data at nodes (i.e. diagonal entries)
    integer :: h_DmatrixCoeffsAtNode = ST_NOHANDLE

    ! Handle to auxiliary matrix data at edges (i.e. off-diagonal entries)
    integer :: h_DmatrixCoeffsAtEdge = ST_NOHANDLE

    ! Handle to bounds at edges (i.e. off-diagonal entries)
    integer :: h_DboundsAtEdge = ST_NOHANDLE

    ! Pointer to the vector of correction factors
    type(t_vectorScalar), pointer :: p_rvectorAlpha => null()

    ! Pointer to the vector of explicit antidiffusive fluxes
    type(t_vectorScalar), pointer :: p_rvectorFlux0 => null()

    ! Pointer to the vector of raw antidiffusive fluxes
    type(t_vectorScalar), pointer :: p_rvectorFlux => null()

    ! Pointer to the vector of prelimiting antidiffusive fluxes
    type(t_vectorScalar), pointer :: p_rvectorFluxPrel => null()

    ! Pointer to the vector of nodal solution values
    type(t_vectorScalar), pointer :: p_rvectorDx => null()

    ! Pointers to the vectors of antidiffusive contributions
    type(t_vectorScalar), pointer :: p_rvectorPp => null()
    type(t_vectorScalar), pointer :: p_rvectorPm => null()

    ! Pointers to the vectors of local solution bounds
    type(t_vectorScalar), pointer :: p_rvectorQ  => null()
    type(t_vectorScalar), pointer :: p_rvectorQp => null()
    type(t_vectorScalar), pointer :: p_rvectorQm => null()

    ! Pointers to the vectors of nodal correction factors
    type(t_vectorScalar), pointer :: p_rvectorRp => null()
    type(t_vectorScalar), pointer :: p_rvectorRm => null()

    ! Pointer to the low-order predictor
    type(t_vectorBlock), pointer :: p_rvectorPredictor => null()

  end type t_afcstab
!</typeblock>

!<typeblock>

  ! This structure can be used to realise arrays-of-pointers which is
  ! necessary to address the content of multiple scalar submatrices
  ! simultaneously. 

  type t_array

    ! Type of data associated to the handle ST_DOUBLE, ST_SINGLE)
    integer :: idataType = ST_NOHANDLE

    ! Pointer to the double-valued matrix data
    real(DP), dimension(:), pointer :: p_Ddata => null()

    ! Pointer to the single-valued matrix data
    real(SP), dimension(:), pointer :: p_Fdata => null()

  end type t_array
!</typeblock>

!</types>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface afcstab_resizeStabilisation
    module procedure afcstab_resizeStabDirect
    module procedure afcstab_resizeStabIndScalar
    module procedure afcstab_resizeStabIndBlock
  end interface

  interface afcstab_isMatrixCompatible
    module procedure afcstab_isMatrixCompatibleSc
    module procedure afcstab_isMatrixCompatibleBl
  end interface
  
  interface afcstab_isVectorCompatible
    module procedure afcstab_isVectorCompatibleSc
    module procedure afcstab_isVectorCompatibleBl
  end interface

  interface afcstab_failsafeLimiting
    module procedure afcstab_failsafeLimitingBlock
    module procedure afcstab_failsafeLimitingArray
  end interface

  interface afcstab_limit
    module procedure afcstab_limit_unbounded
    module procedure afcstab_limit_bounded
  end interface

  interface afcstab_getbase_array
    module procedure afcstab_getbase_arrayScalar
    module procedure afcstab_getbase_arrayBlock
  end interface

  interface afcstab_buildBoundsLPT
    module procedure afcstab_buildBoundsLPT1D
    module procedure afcstab_buildBoundsLPT2D
    module procedure afcstab_buildBoundsLPT3D
  end interface

  interface afcstab_getbase_DcoeffsAtEdge
    module procedure afcstab_getbase_DcoeffsAtEdge
    module procedure afcstab_getbase_FcoeffsAtEdge
  end interface

  interface afcstab_getbase_DmatCoeffAtNode
    module procedure afcstab_getbase_DmatCoeffAtNode
    module procedure afcstab_getbase_FmatCoeffAtNode
  end interface

  interface afcstab_getbase_DmatCoeffAtEdge
    module procedure afcstab_getbase_DmatCoeffAtEdge
    module procedure afcstab_getbase_FmatCoeffAtEdge
  end interface

  interface afcstab_getbase_DboundsAtEdge
    module procedure afcstab_getbase_DboundsAtEdge
    module procedure afcstab_getbase_FboundsAtEdge
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine afcstab_initFromParameterlist(rparlist, ssectionName, rafcstab)

!<description>
    ! This subroutine creates a stabilisation structure and initializes
    ! its values from a given parameter list
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in)    :: rparlist

    ! Section name of the parameter list
    character(LEN=*), intent(in)   :: ssectionName
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: istabilisation,iprelimiting

    ! First, we retrieve all information of the stabilisation
    ! structure specified in the parameter file. 

    ! Get type of stabilisation from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "istabilisation", rafcstab%ctypeAFCstabilisation, AFCSTAB_GALERKIN)

    ! Get type of prelimiting from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "iprelimiting", rafcstab%ctypePrelimiting, AFCSTAB_PRELIMITING_NONE)


    ! In a second step, special settings particular to individual
    ! stabilisation techniques are made `by hand'.
    
    ! Check if stabilisation should be applied
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)
      ! high-order Galerkin scheme, low-order scheme:
      ! no prelimiting, no limiting
      rafcstab%ctypePrelimiting = AFCSTAB_PRELIMITING_NONE
      rafcstab%ctypeLimiting    = AFCSTAB_LIMITING_NONE

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_LINFCT)
      ! Symmetric flux limiting with prelimiting as given by the user
      rafcstab%ctypeLimiting = AFCSTAB_LIMITING_SYMMETRIC

    case (AFCSTAB_LINFCT_MASS,&
          AFCSTAB_SYMMETRIC,&
          AFCSTAB_NLINLPT_MASS,&
          AFCSTAB_NLINLPT_SYMMETRIC,&
          AFCSTAB_LINLPT_MASS,&
          AFCSTAB_LINLPT_SYMMETRIC)
      ! Symmetric flux limiting without prelimiting
      rafcstab%ctypePrelimiting = AFCSTAB_PRELIMITING_NONE
      rafcstab%ctypeLimiting    = AFCSTAB_LIMITING_SYMMETRIC

    case (AFCSTAB_LINLPT_UPWINDBIASED)
      ! For this special case, symmetric flux limiting must be
      ! performed as stated in the conclusion of the paper:
      ! D. Kuzmin, Linearity-preserving flux correction and
      ! convergence acceleration for constrained Galerin schemes,
      ! Technical Report 421, TU Dortmund, 2011.
      rafcstab%ctypePrelimiting = AFCSTAB_PRELIMITING_NONE
      rafcstab%ctypeLimiting    = AFCSTAB_LIMITING_SYMMETRIC

    case (AFCSTAB_TVD)
      ! Upwind-biased flux limiting with standard prelimiting
      rafcstab%ctypePrelimiting = AFCSTAB_PRELIMITING_STD
      rafcstab%ctypeLimiting    = AFCSTAB_LIMITING_UPWINDBIASED

    case (AFCSTAB_GP,&
          AFCSTAB_NLINLPT_UPWINDBIASED)
      ! Upwind-biased flux limiting with minmod prelimiting
      rafcstab%ctypePrelimiting = AFCSTAB_PRELIMITING_MINMOD
      rafcstab%ctypeLimiting    = AFCSTAB_LIMITING_UPWINDBIASED

    case default
      call output_line('Invalid AFC type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_initFromParameterlist')
      call sys_halt()
    end select

  end subroutine afcstab_initFromParameterlist

  !*****************************************************************************

!<subroutine>
  
  subroutine afcstab_releaseStabilisation(rafcstab)

!<description>
    ! This subroutine releases a stabilisation structure
!</description>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Release edge structure
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGESTRUCTURE)) then
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdge)
      if (rafcstab%h_IverticesAtEdgeIdx .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IverticesAtEdgeIdx)
    end if
    rafcstab%h_IverticesAtEdge = ST_NOHANDLE
    rafcstab%h_IverticesAtEdgeIdx = ST_NOHANDLE

    ! Release edge values
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEVALUES) .and.&
        (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE))&
        call storage_free(rafcstab%h_DcoefficientsAtEdge)
    rafcstab%h_DcoefficientsAtEdge = ST_NOHANDLE

    ! Release edge bounds
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEBOUNDS) .and.&
        (rafcstab%h_DboundsAtEdge .ne. ST_NOHANDLE))&
        call storage_free(rafcstab%h_DboundsAtEdge)
    rafcstab%h_DboundsAtEdge = ST_NOHANDLE

    ! Release off-diagonal edges
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES)) then
      if (rafcstab%h_IsuperdiagEdgesIdx .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IsuperdiagEdgesIdx)
      if (rafcstab%h_IsubdiagEdgesIdx .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IsubdiagEdgesIdx)
      if (rafcstab%h_IsubdiagEdges .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IsubdiagEdges)
    end if
    rafcstab%h_IsuperdiagEdgesIdx = ST_NOHANDLE
    rafcstab%h_IsubdiagEdgesIdx = ST_NOHANDLE
    rafcstab%h_IsubdiagEdges = ST_NOHANDLE
    
    ! Release matrix data
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_MATRIXCOEFFS)) then
      if (rafcstab%h_DmatrixCoeffsAtNode .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_DmatrixCoeffsAtNode)
      if (rafcstab%h_DmatrixCoeffsAtEdge .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_DmatrixCoeffsAtEdge)
    end if
    rafcstab%h_DmatrixCoeffsAtNode = ST_NOHANDLE
    rafcstab%h_DmatrixCoeffsAtEdge = ST_NOHANDLE

    ! Release antidiffusive fluxes
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES)) then
      if (associated(rafcstab%p_rvectorFlux0)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorFlux0)
        deallocate(rafcstab%p_rvectorFlux0)
      end if
      if (associated(rafcstab%p_rvectorFlux)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorFlux)
        deallocate(rafcstab%p_rvectorFlux)
      end if
      if (associated(rafcstab%p_rvectorFluxPrel)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorFluxPrel)
        deallocate(rafcstab%p_rvectorFluxPrel)
      end if
    end if
    rafcstab%p_rvectorFlux      => null()
    rafcstab%p_rvectorFlux0     => null()
    rafcstab%p_rvectorFluxPrel  => null()

    ! Release transformed nodal solution values
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEVALUES)) then
      if (associated(rafcstab%p_rvectorDx)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorDx)
        deallocate(rafcstab%p_rvectorDx)
      end if
    end if
    rafcstab%p_rvectorDx => null()

    ! Release antidiffusive increments
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADINCREMENTS)) then
      if (associated(rafcstab%p_rvectorPp)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorPp)
        deallocate(rafcstab%p_rvectorPp)
      end if
      if (associated(rafcstab%p_rvectorPm)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorPm)
        deallocate(rafcstab%p_rvectorPm)
      end if
    end if
    rafcstab%p_rvectorPp => null()
    rafcstab%p_rvectorPm => null()

    ! Release upper/lower bounds
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS)) then
      if (associated(rafcstab%p_rvectorQ)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorQ)
        deallocate(rafcstab%p_rvectorQ)
      end if
      if (associated(rafcstab%p_rvectorQp)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorQp)
        deallocate(rafcstab%p_rvectorQp)
      end if
      if (associated(rafcstab%p_rvectorQm)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorQm)
        deallocate(rafcstab%p_rvectorQm)
      end if
    end if
    rafcstab%p_rvectorQ  => null()
    rafcstab%p_rvectorQp => null()
    rafcstab%p_rvectorQm => null()
    
    ! Release nodal limiting factors
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODELIMITER)) then
      if (associated(rafcstab%p_rvectorRp)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorRp)
        deallocate(rafcstab%p_rvectorRp)
      end if
      if (associated(rafcstab%p_rvectorRm)) then
        call lsyssc_releaseVector(rafcstab%p_rvectorRm)
        deallocate(rafcstab%p_rvectorRm)
      end if
    end if
    rafcstab%p_rvectorRp => null()
    rafcstab%p_rvectorRm => null()

    ! Release edge-wise limiting coefficients
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGELIMITER) .and.&
        (associated(rafcstab%p_rvectorAlpha))) then
      call lsyssc_releaseVector(rafcstab%p_rvectorAlpha)
      deallocate(rafcstab%p_rvectorAlpha)
    end if
    rafcstab%p_rvectorAlpha => null()

    ! Release low-order predictor
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_PREDICTOR) .and.&
        (associated(rafcstab%p_rvectorPredictor))) then
      call lsysbl_releaseVector(rafcstab%p_rvectorPredictor)
      deallocate(rafcstab%p_rvectorPredictor)
    end if
    rafcstab%p_rvectorPredictor => null()

    ! Reset atomic data
    rafcstab%ctypeAFCstabilisation = AFCSTAB_GALERKIN
    rafcstab%ctypePrelimiting      = AFCSTAB_PRELIMITING_NONE
    rafcstab%ctypeLimiting         = AFCSTAB_LIMITING_NONE
    rafcstab%cmatrixFormat         = LSYSSC_MATRIXUNDEFINED
    rafcstab%cdataType             = ST_DOUBLE
    rafcstab%istabilisationSpec    = AFCSTAB_UNDEFINED
    rafcstab%iduplicationFlag      = 0
    rafcstab%NEQ                   = 0
    rafcstab%NVAR                  = 1
    rafcstab%NVARtransformed       = 1
    rafcstab%NEDGE                 = 0
    rafcstab%NNVEDGE               = 0
    rafcstab%nmatrixCoeffs         = 0
    rafcstab%ncoeffs               = 0

  contains

    !**************************************************************
    ! Checks if bitfield ibitfield in idupFlag is not set.

    pure function check(idupFlag, ibitfield)

      integer(I32), intent(in) :: idupFlag,ibitfield
      
      logical :: check
      
      check = (iand(idupFlag,ibitfield) .ne. ibitfield)

    end function check

  end subroutine afcstab_releaseStabilisation

   !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabDirect(rafcstab, neq, nedge)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure to the new values NEQ and NEDGE
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! number of equations
    integer, intent(in) :: neq

    ! number of edges
    integer, intent(in) :: nedge   
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>
    
    
    ! REMARK: The handle IverticesAtEdgeIdx is not modified here due
    ! to the following reasons: If the edges are not reordered into
    ! independent groups then IverticesAtEdgeIdx(1:2)=(/1,NEDGE+1/).
    ! Otherwise, the edge-data structure needs to be regenerated anyway.
    
    ! Resize nodal quantities
    if (rafcstab%NEQ .ne. neq) then

      ! Set new number of nodes
      rafcstab%NEQ = neq
      
      ! Resize edge index vector and clear specification 
      if (rafcstab%h_IsuperdiagEdgesIdx .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEQ+1, rafcstab%h_IsuperdiagEdgesIdx,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_IsuperdiagEdgesIdx '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize subdiagonal edge index vector and clear specification
      if (rafcstab%h_IsubdiagEdgesIdx .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEQ+1, rafcstab%h_IsubdiagEdgesIdx,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_IsubdiagEdgesIdx '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

      ! Resize matrix data at nodes and clear specification
      if (rafcstab%h_DmatrixCoeffsAtNode .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_MATRIXCOEFFS)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEQ, rafcstab%h_DmatrixCoeffsAtNode,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_DmatrixCoeffsAtNode '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

      ! Resize transformed nodal vector
      if (associated(rafcstab%p_rvectorDx)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEVALUES)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorDx,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorDx '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

      ! Resize nodal vectors P+
      if (associated(rafcstab%p_rvectorPp)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADINCREMENTS)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorPp,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorPp '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize nodal vectors P-
      if (associated(rafcstab%p_rvectorPm)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADINCREMENTS)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorPm,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorPm '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize nodal vectors Q
      if (associated(rafcstab%p_rvectorQ)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorQ,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorQ '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

      ! Resize nodal vectors Q+
      if (associated(rafcstab%p_rvectorQp)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorQp,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorQp '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize nodal vectors Q-
      if (associated(rafcstab%p_rvectorQm)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorQm,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorQm '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize nodal vectors R+
      if (associated(rafcstab%p_rvectorRp)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODELIMITER)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorRp,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorRp '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize nodal vectors R-
      if (associated(rafcstab%p_rvectorRm)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODELIMITER)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorRm,&
          rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorRm '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize nodal vector for the low-order predictor
      if (associated(rafcstab%p_rvectorPredictor)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_PREDICTOR)) then
          call lsysbl_resizeVectorBlock(rafcstab%p_rvectorPredictor,&
              rafcstab%NEQ, .false., .false.)
        else
          call output_line('Vector p_rvectorPredictor '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

    end if
      
    
    ! Resize edge quantities
    if (rafcstab%NEDGE .ne. nedge) then

      ! Set new number of edges
      rafcstab%NEDGE = nedge

      ! Resize array of edges
      if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGESTRUCTURE)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_IverticesAtEdge,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_IverticesAtEdge '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize array of subdiagonal edges
      if (rafcstab%h_IsubdiagEdges .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_IsubdiagEdges,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_IsubdiagEdges '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

      ! Resize array of edge data
      if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEVALUES)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_DcoefficientsAtEdge,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_DcoefficientsAtEdge '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize array of edge bounds
      if (rafcstab%h_DboundsAtEdge .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEBOUNDS)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_DboundsAtEdge,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_DboundsAtEdge '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      ! Resize matrix data at edges
      if (rafcstab%h_DmatrixCoeffsAtEdge .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_MATRIXCOEFFS)) then
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_DmatrixCoeffsAtEdge,&
              ST_NEWBLOCK_NOINIT, .false.)
        else
          call output_line('Handle h_DmatrixCoeffsAtEdge '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

      ! Resize edge vectors for the correction factor and the fluxes
      if(associated(rafcstab%p_rvectorAlpha)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGELIMITER)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorAlpha,&
              rafcstab%NEDGE, .false., .false.)
        else
          call output_line('Vector p_rvectorAlpha '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      if(associated(rafcstab%p_rvectorFlux0)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorFlux0,&
              rafcstab%NEDGE, .false., .false.)
        else
          call output_line('Vector p_rvectorFlux0 '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

      if(associated(rafcstab%p_rvectorFlux)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorFlux,&
              rafcstab%NEDGE, .false., .false.)
        else
          call output_line('Vector p_rvectorFlux '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if
      
      if(associated(rafcstab%p_rvectorFluxPrel)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES)) then
          call lsyssc_resizeVector(rafcstab%p_rvectorFluxPrel,&
              rafcstab%NEDGE, .false., .false.)
        else
          call output_line('Vector p_rvectorFluxPrel '//&
              'is shared and cannot be resized!',&
              OU_CLASS_WARNING,OU_MODE_STD,'afcstab_resizeStabDirect')
        end if
      end if

    end if
    
    ! Set state of stabilisation to either undefined or initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      rafcstab%istabilisationSpec = AFCSTAB_UNDEFINED
    else
      rafcstab%istabilisationSpec = AFCSTAB_INITIALISED
    end if

  contains

    !**************************************************************
    ! Checks if bitfield ibitfield in idupFlag is not set.

    pure function check(idupFlag, ibitfield)

      integer(I32), intent(in) :: idupFlag,ibitfield
      
      logical :: check
      
      check = (iand(idupFlag,ibitfield) .ne. ibitfield)

    end function check

  end subroutine afcstab_resizeStabDirect

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabIndScalar(rafcstab, rmatrixTemplate)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure so that they are compatible to the template matrix.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(in) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: neq, nedge

    
    ! Determine number of equations and edges
    neq   = rmatrixTemplate%NEQ
    nedge = (rmatrixTemplate%NA-rmatrixTemplate%NEQ)/2

    ! Call resize routine directly
    call afcstab_resizeStabDirect(rafcstab, neq, nedge)

  end subroutine afcstab_resizeStabIndScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabIndBlock(rafcstab, rmatrixBlockTemplate)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure so that they are compatible to the template matrix.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! template matrix
    type(t_matrixBlock), intent(in) :: rmatrixBlockTemplate
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: neq, nedge

    
    ! Check if block matrix has only one block
    if ((rmatrixBlockTemplate%nblocksPerCol .eq. 1) .and.&
        (rmatrixBlockTemplate%nblocksPerRow .eq. 1)) then
      call afcstab_resizeStabIndScalar(rafcstab,&
          rmatrixBlockTemplate%RmatrixBlock(1,1))

      ! That is it
      return
    end if
    

    ! Determine number of equations and edges
    neq   = rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ
    nedge = (rmatrixBlockTemplate%RmatrixBlock(1,1)%NA-&
             rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ)/2

    ! Call resize routine directly
    call afcstab_resizeStabDirect(rafcstab, neq, nedge)

  end subroutine afcstab_resizeStabIndBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyStabilisation(rafcstabSrc, rafcstabDest, idupFlag)

!<description>
    ! This subroutine selectively copies data from the source
    ! stabilisation structure rafcstabScr to the destination
    ! stabilisation structure rsfcatsbDest
!</description>

!<input>
    ! Source stabilisation structure
    type(t_afcstab), intent(in) :: rafcstabSrc
    
    ! Duplication flag that decides on how to set up the structure
    integer(I32), intent(in) :: idupFlag
!</input>

!<inputoutput>
    ! Destination stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstabDest
!</inputoutput>
!</subroutine>
   
    ! Copy structural data
    if (check(idupFlag, AFCSTAB_DUP_STRUCTURE) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_INITIALISED)) then
      rafcstabDest%ctypeAFCstabilisation = rafcstabSrc%ctypeAFCstabilisation
      rafcstabDest%ctypePrelimiting      = rafcstabSrc%ctypePrelimiting
      rafcstabDest%ctypeLimiting         = rafcstabSrc%ctypeLimiting
      rafcstabDest%cmatrixFormat         = rafcstabSrc%cmatrixFormat
      rafcstabDest%cdataType             = rafcstabSrc%cdataType
      rafcstabDest%NEQ                   = rafcstabSrc%NEQ
      rafcstabDest%NVAR                  = rafcstabSrc%NVAR
      rafcstabDest%NVARtransformed       = rafcstabSrc%NVARtransformed
      rafcstabDest%NEDGE                 = rafcstabSrc%NEDGE
      rafcstabDest%NNVEDGE               = rafcstabSrc%NNVEDGE
      rafcstabDest%istabilisationSpec    = ior(rafcstabDest%istabilisationSpec,&
                                               AFCSTAB_INITIALISED)
    end if

    ! Copy edge structre
    if (check(idupFlag, AFCSTAB_DUP_EDGESTRUCTURE) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)) then
      ! Copy content from source to destination structure
      call storage_copy(rafcstabSrc%h_IverticesAtEdgeIdx,&
          rafcstabDest%h_IverticesAtEdgeIdx)
      call storage_copy(rafcstabSrc%h_IverticesAtEdge,&
          rafcstabDest%h_IverticesAtEdge)
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE))
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION))
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_EDGESTRUCTURE))
    end if

    ! Copy edge values
    if (check(idupFlag, AFCSTAB_DUP_EDGEVALUES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)) then
      ! Copy content from source to destination structure
      call storage_copy(rafcstabSrc%h_DcoefficientsAtEdge,&
          rafcstabDest%h_DcoefficientsAtEdge)
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES))
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_EDGEVALUES))
    end if
    
    ! Copy off-diagonal edges
    if (check(idupFlag, AFCSTAB_DUP_OFFDIAGONALEDGES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES)) then
      ! Copy content from source to destination structure
      call storage_copy(rafcstabSrc%h_IsuperdiagEdgesIdx,&
          rafcstabDest%h_IsuperdiagEdgesIdx)
      call storage_copy(rafcstabSrc%h_IsubdiagEdgesIdx,&
          rafcstabDest%h_IsubdiagEdgesIdx)
      call storage_copy(rafcstabSrc%h_IsubdiagEdges,&
          rafcstabDest%h_IsubdiagEdges)
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES))
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_OFFDIAGONALEDGES))
    end if

    ! Copy matrix data at nodes and edges
    if (check(idupFlag, AFCSTAB_DUP_MATRIXCOEFFS) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)) then
      ! Copy content from source to destination structure
      call storage_copy(rafcstabSrc%h_DmatrixCoeffsAtNode,&
          rafcstabDest%h_DmatrixCoeffsAtNode)
      call storage_copy(rafcstabSrc%h_DmatrixCoeffsAtEdge,&
          rafcstabDest%h_DmatrixCoeffsAtEdge)
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS))
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_MATRIXCOEFFS))
    end if

    ! Copy transformed nodal solution values
    if (check(idupFlag, AFCSTAB_DUP_NODEVALUES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEVALUES)) then
      ! Unlink pointers to vectors which are shared by the structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_NODEVALUES))&
          nullify(rafcstabDest%p_rvectorDx)
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_NODEVALUES))
      ! Copy transformed nodal solution values (if any)
      if (associated(rafcstabSrc%p_rvectorDx)) then
        if (.not.(associated(rafcstabDest%p_rvectorDx)))&
            allocate(rafcstabDest%p_rvectorDx)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorDx,&
            rafcstabDest%p_rvectorDx)
      end if
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEVALUES))
    end if

    ! Copy antidiffusive fluxes
    if (check(idupFlag, AFCSTAB_DUP_ADFLUXES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)) then
      ! Unlink pointers to vectors which are shared by the structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES))&
          nullify(rafcstabDest%p_rvectorFlux0, rafcstabDest%p_rvectorFlux)
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_ADFLUXES))
      ! Copy explicit raw antidiffusive flux (if any)
      if (associated(rafcstabSrc%p_rvectorFlux0)) then
        if (.not.(associated(rafcstabDest%p_rvectorFlux0)))&
            allocate(rafcstabDest%p_rvectorFlux0)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorFlux0,&
            rafcstabDest%p_rvectorFlux0)
      end if
      ! Copy raw antidiffusive flux (if any)
      if (associated(rafcstabSrc%p_rvectorFlux)) then
        if (.not.(associated(rafcstabDest%p_rvectorFlux)))&
            allocate(rafcstabDest%p_rvectorFlux)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorFlux,&
            rafcstabDest%p_rvectorFlux)
      end if
      ! Copy prelimited antidiffusive flux (if any)
      if (associated(rafcstabSrc%p_rvectorFluxPrel)) then
        if (.not.(associated(rafcstabDest%p_rvectorFluxPrel)))&
            allocate(rafcstabDest%p_rvectorFluxPrel)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorFluxPrel,&
            rafcstabDest%p_rvectorFluxPrel)
      end if

      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADFLUXES))
    end if

    ! Copy antidiffusive increments
    if (check(idupFlag, AFCSTAB_DUP_ADINCREMENTS) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)) then
      ! Unlink pointers to vectors which are shared by the structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_ADINCREMENTS))&
          nullify(rafcstabDest%p_rvectorPp, rafcstabDest%p_rvectorPm)
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_ADINCREMENTS))
      ! Copy antidiffusive increment Pp
      if (associated(rafcstabSrc%p_rvectorPp)) then
        if (.not.(associated(rafcstabDest%p_rvectorPp)))&
            allocate(rafcstabDest%p_rvectorPp)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorPp,&
            rafcstabDest%p_rvectorPp)
      end if
      ! Copy antidiffusive increment Pm
      if (associated(rafcstabSrc%p_rvectorPm)) then
        if (.not.(associated(rafcstabDest%p_rvectorPm)))&
            allocate(rafcstabDest%p_rvectorPm)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorPm,&
            rafcstabDest%p_rvectorPm)
      end if
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS))
    end if
    
    ! Copy upper/lower bounds
    if (check(idupFlag, AFCSTAB_DUP_NODEBOUNDS) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)) then
      ! Unlink pointers to vectors which are shared by the structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS))&
          nullify(rafcstabDest%p_rvectorQp,&
                  rafcstabDest%p_rvectorQm,&
                  rafcstabDest%p_rvectorQ)
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_NODEBOUNDS))
      ! Copy upper bounds Q
      if (associated(rafcstabSrc%p_rvectorQ)) then
        if (.not.(associated(rafcstabDest%p_rvectorQ)))&
            allocate(rafcstabDest%p_rvectorQ)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorQ,&
            rafcstabDest%p_rvectorQ)
      end if
      ! Copy upper bounds Qp
      if (associated(rafcstabSrc%p_rvectorQp)) then
        if (.not.(associated(rafcstabDest%p_rvectorQp)))&
            allocate(rafcstabDest%p_rvectorQp)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorQp,&
            rafcstabDest%p_rvectorQp)
      end if
      ! Copy lower bounds Qm
      if (associated(rafcstabSrc%p_rvectorQm)) then
        if (.not.(associated(rafcstabDest%p_rvectorQm)))&
            allocate(rafcstabDest%p_rvectorQm)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorQm,&
            rafcstabDest%p_rvectorQm)
      end if
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS))
    end if

    ! Copy nodal limiting coefficients
    if (check(idupFlag, AFCSTAB_DUP_NODELIMITER) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)) then
      ! Unlink pointers to vectors which are shared by the structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_NODELIMITER))&
          nullify(rafcstabDest%p_rvectorRp, rafcstabDest%p_rvectorRm)
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_NODELIMITER))
      ! Copy nodal limiting coefficients Rp
      if (associated(rafcstabSrc%p_rvectorRp)) then
        if (.not.(associated(rafcstabDest%p_rvectorRp)))&
            allocate(rafcstabDest%p_rvectorRp)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorRp,&
            rafcstabDest%p_rvectorRp)
      end if
      ! Copy nodal limiting coefficients Rm
      if (associated(rafcstabSrc%p_rvectorRm)) then
        if (.not.(associated(rafcstabDest%p_rvectorRm)))&
            allocate(rafcstabDest%p_rvectorRm)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorRm,&
            rafcstabDest%p_rvectorRm)
      end if
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODELIMITER))
    end if
    
    ! Copy edge-wise limiting coefficients
    if (check(idupFlag, AFCSTAB_DUP_EDGELIMITER) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)) then
      ! Unlink pointers to vectors which are shared by the structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_EDGELIMITER))&
          nullify(rafcstabDest%p_rvectorAlpha)
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_EDGELIMITER))
      ! Copy edge-wise limiting coefficients
      if (associated(rafcstabSrc%p_rvectorAlpha)) then
        if (.not.(associated(rafcstabDest%p_rvectorAlpha)))&
            allocate(rafcstabDest%p_rvectorAlpha)
        call lsyssc_copyVector(rafcstabSrc%p_rvectorAlpha,&
            rafcstabDest%p_rvectorAlpha)
      end if
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER))
    end if
    
    ! Copy low-order predictor
    if (check(idupFlag, AFCSTAB_DUP_PREDICTOR) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_PREDICTOR)) then
      ! Unlink pointers to vectors which are shared by the structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_PREDICTOR))&
          nullify(rafcstabDest%p_rvectorPredictor)
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_PREDICTOR))
      ! Copy edge-wise limiting coefficients
      if (associated(rafcstabSrc%p_rvectorPredictor)) then
        if (.not.(associated(rafcstabDest%p_rvectorPredictor)))&
            allocate(rafcstabDest%p_rvectorPredictor)
        call lsysbl_copyVector(rafcstabSrc%p_rvectorPredictor,&
            rafcstabDest%p_rvectorPredictor)
      end if
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_PREDICTOR))
    end if

  contains

    !**************************************************************
    ! Checks if idupFlag has all bits ibitfield set.

    pure function check(idupFlag, ibitfield)
      
      integer(I32), intent(in) :: idupFlag,ibitfield
      
      logical :: check
      
      check = (iand(idupFlag,ibitfield) .eq. ibitfield)

    end function check

  end subroutine afcstab_copyStabilisation

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_duplicateStabilisation(rafcstabSrc, rafcstabDest, idupFlag)

!<description>
    ! This subroutine duplicated parts of the stabilisation structure
    ! rafcstabScr in the stabilisation structure rsfcatsbDest. Note
    ! that rafcstabScr is still the owner of the duplicated content.
!</description>

!<input>
    ! Source stabilisation structure
    type(t_afcstab), intent(in) :: rafcstabSrc
    
    ! Duplication flag that decides on how to set up the structure
    integer(I32), intent(in) :: idupFlag
!</input>

!<inputoutput>
    ! Destination stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstabDest
!</inputoutput>
!</subroutine>

    ! Duplicate structural data
    if (check(idupFlag, AFCSTAB_DUP_STRUCTURE) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_INITIALISED)) then
      rafcstabDest%ctypeAFCstabilisation = rafcstabSrc%ctypeAFCstabilisation
      rafcstabDest%ctypePrelimiting      = rafcstabSrc%ctypePrelimiting
      rafcstabDest%ctypeLimiting         = rafcstabSrc%ctypeLimiting
      rafcstabDest%cmatrixFormat         = rafcstabSrc%cmatrixFormat
      rafcstabDest%cdataType             = rafcstabSrc%cdataType
      rafcstabDest%NEQ                   = rafcstabSrc%NEQ
      rafcstabDest%NVAR                  = rafcstabSrc%NVAR
      rafcstabDest%NVARtransformed       = rafcstabSrc%NVARtransformed
      rafcstabDest%NEDGE                 = rafcstabSrc%NEDGE
      rafcstabDest%NNVEDGE               = rafcstabSrc%NNVEDGE
      rafcstabDest%istabilisationSpec    = ior(rafcstabDest%istabilisationSpec,&
                                               AFCSTAB_INITIALISED)
    end if

    ! Duplicate edge structre
    if (check(idupFlag, AFCSTAB_DUP_EDGESTRUCTURE) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_EDGESTRUCTURE))) then
        call storage_free(rafcstabDest%h_IverticesAtEdgeIdx)
        call storage_free(rafcstabDest%h_IverticesAtEdge)
      end if
      ! Copy handle from source to destination structure
      rafcstabDest%h_IverticesAtEdgeIdx = rafcstabSrc%h_IverticesAtEdgeIdx
      rafcstabDest%h_IverticesAtEdge = rafcstabSrc%h_IverticesAtEdge
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE))
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_EDGESTRUCTURE)
    end if

    ! Duplicate edge values
    if (check(idupFlag, AFCSTAB_DUP_EDGEVALUES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_EDGEVALUES)))&
          call storage_free(rafcstabDest%h_DcoefficientsAtEdge)
      ! Copy handle from source to destination structure
      rafcstabDest%h_DcoefficientsAtEdge = rafcstabSrc%h_DcoefficientsAtEdge
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_EDGEVALUES)
    end if

    ! Duplicate edge bounds
    if (check(idupFlag, AFCSTAB_DUP_EDGEBOUNDS) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_EDGEBOUNDS)))&
          call storage_free(rafcstabDest%h_DboundsAtEdge)
      ! Copy handle from source to destination structure
      rafcstabDest%h_DboundsAtEdge = rafcstabSrc%h_DboundsAtEdge
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_EDGEBOUNDS)
    end if

    ! Duplicate off-diagonal edges
    if (check(idupFlag, AFCSTAB_DUP_OFFDIAGONALEDGES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES))) then
        call storage_free(rafcstabDest%h_IsuperdiagEdgesIdx)
        call storage_free(rafcstabDest%h_IsubdiagEdgesIdx)
        call storage_free(rafcstabDest%h_IsubdiagEdges)
      end if
      ! Copy handles from source to destination structure
      rafcstabDest%h_IsuperdiagEdgesIdx = rafcstabSrc%h_IsuperdiagEdgesIdx
      rafcstabDest%h_IsubdiagEdgesIdx   = rafcstabSrc%h_IsubdiagEdgesIdx
      rafcstabDest%h_IsubdiagEdges      = rafcstabSrc%h_IsubdiagEdges
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_OFFDIAGONALEDGES)
    end if

    ! Duplicate matrix data
    if (check(idupFlag, AFCSTAB_DUP_MATRIXCOEFFS) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_MATRIXCOEFFS))) then
        call storage_free(rafcstabDest%h_DmatrixCoeffsAtNode)
        call storage_free(rafcstabDest%h_DmatrixCoeffsAtEdge)
      end if
      ! Copy handles from source to destination structure
      rafcstabDest%h_DmatrixCoeffsAtNode = rafcstabSrc%h_DmatrixCoeffsAtNode
      rafcstabDest%h_DmatrixCoeffsAtEdge = rafcstabSrc%h_DmatrixCoeffsAtEdge
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_MATRIXCOEFFS)
    end if

    ! Duplicate transformed nodal solution values
    if (check(idupFlag, AFCSTAB_DUP_NODEVALUES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEVALUES)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_NODEVALUES))) then
        call lsyssc_releaseVector(rafcstabDest%p_rvectorDx)
        deallocate(rafcstabDest%p_rvectorDx)
      end if
      ! Copy pointers from source to destination structure
      rafcstabDest%p_rvectorDx => rafcstabSrc%p_rvectorDx
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEVALUES))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_NODEVALUES)
    end if

    ! Duplicate antidiffusive fluxes
    if (check(idupFlag, AFCSTAB_DUP_ADFLUXES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES))) then
        call lsyssc_releaseVector(rafcstabDest%p_rvectorFlux0)
        call lsyssc_releaseVector(rafcstabDest%p_rvectorFlux)
        deallocate(rafcstabDest%p_rvectorFlux0, rafcstabDest%p_rvectorFlux)
        if (associated(rafcstabDest%p_rvectorFluxPrel)) then
          call lsyssc_releaseVector(rafcstabDest%p_rvectorFluxPrel)
          deallocate(rafcstabDest%p_rvectorFluxPrel)
        end if
      end if
      ! Copy pointers from source to destination structure
      rafcstabDest%p_rvectorFlux0 => rafcstabSrc%p_rvectorFlux0
      rafcstabDest%p_rvectorFlux  => rafcstabSrc%p_rvectorFlux
      if (associated(rafcstabSrc%p_rvectorFluxPrel))&
          rafcstabDest%p_rvectorFluxPrel => rafcstabSrc%p_rvectorFluxPrel
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADFLUXES))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_ADFLUXES)
    end if

    ! Duplicate antidiffusive incremens
    if (check(idupFlag, AFCSTAB_DUP_ADINCREMENTS) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_ADINCREMENTS))) then
        call lsyssc_releaseVector(rafcstabDest%p_rvectorPp)
        call lsyssc_releaseVector(rafcstabDest%p_rvectorPm)
      end if
      ! Copy pointers from source to destination structure
      rafcstabDest%p_rvectorPp => rafcstabSrc%p_rvectorPp
      rafcstabDest%p_rvectorPm => rafcstabSrc%p_rvectorPm
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_ADINCREMENTS)
    end if

    ! Duplicate upper/lower bounds
    if (check(idupFlag, AFCSTAB_DUP_NODEBOUNDS) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS))) then
        call lsyssc_releaseVector(rafcstabDest%p_rvectorQp)
        call lsyssc_releaseVector(rafcstabDest%p_rvectorQm)
        if (associated(rafcstabDest%p_rvectorQ))&
            call lsyssc_releaseVector(rafcstabDest%p_rvectorQ)
      end if
      ! Copy pointers from source to destination structure
      rafcstabDest%p_rvectorQp => rafcstabSrc%p_rvectorQp
      rafcstabDest%p_rvectorQm => rafcstabSrc%p_rvectorQm
      if (associated(rafcstabSrc%p_rvectorQ))&
          rafcstabDest%p_rvectorQ => rafcstabSrc%p_rvectorQ
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_NODEBOUNDS)
    end if

    ! Duplicate nodal limiting coefficients
    if (check(idupFlag, AFCSTAB_DUP_NODELIMITER) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)) then
      ! Remove existing data owned by the destination structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_NODELIMITER)) then
        call lsyssc_releaseVector(rafcstabDest%p_rvectorRp)
        call lsyssc_releaseVector(rafcstabDest%p_rvectorRm)
      end if
      ! Copy pointers from source to destination structure
      rafcstabDest%p_rvectorRp => rafcstabSrc%p_rvectorRp
      rafcstabDest%p_rvectorRm => rafcstabSrc%p_rvectorRm
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_NODELIMITER))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_NODELIMITER)
    end if

    ! Duplicate edge-wise limiting coefficients
    if (check(idupFlag, AFCSTAB_DUP_EDGELIMITER) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)) then
      ! Remove existing data owned by the destination structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_EDGELIMITER)) then
        call lsyssc_releaseVector(rafcstabDest%p_rvectorAlpha)
      end if
      ! Copy pointers from source to destination structure
      rafcstabDest%p_rvectorAlpha => rafcstabSrc%p_rvectorAlpha
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_EDGELIMITER)
    end if

    ! Duplicate low-order predictor
    if (check(idupFlag, AFCSTAB_DUP_PREDICTOR) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_PREDICTOR)) then
      ! Remove existing data owned by the destination structure
      if (check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_PREDICTOR)) then
        call lsysbl_releaseVector(rafcstabDest%p_rvectorPredictor)
      end if
      ! Copy pointers from source to destination structure
      rafcstabDest%p_rvectorPredictor => rafcstabSrc%p_rvectorPredictor
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_PREDICTOR))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_PREDICTOR)
    end if
    
  contains

    !**************************************************************
    ! Checks if idupFlag has all bits ibitfield set.

    pure function check(idupFlag, ibitfield)
      
      integer(I32), intent(in) :: idupFlag,ibitfield
      
      logical :: check
      
      check = (iand(idupFlag,ibitfield) .eq. ibitfield)

    end function check

  end subroutine afcstab_duplicateStabilisation

  ! ***************************************************************************

!<subroutine>

  subroutine afcstab_initMatrixCoeffs(rafcstab, n)

!<description>
    ! This subroutine initialises the arrays DmatrixCoeffsAtNode and
    ! DmatrixCoeffsAtEdge to store auxiliary data for n matrices.
!</description>

!<input>
    ! Number of auxliary matrices to store at most
    integer, intent(in) :: n
!</input>

!<inputoutpu>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    integer, dimension(2) :: Isize2D
    integer, dimension(3) :: Isize3D

    if (n .le. 0) then
      call output_line('Number of matrices to be stored must be positive!',&
          OU_CLASS_WARNING,OU_MODE_STD,'afcstab_initMatrixCoeffs')
      return
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_initMatrixCoeffs')
      call sys_halt()
    end if

    ! Remove auxiliary matrix data if present and not shared with others
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS) .ne. 0) then
      if (iand(rafcstab%iduplicationFlag, AFCSTAB_SHARE_MATRIXCOEFFS) .eq. 0) then
        call storage_free(rafcstab%h_DmatrixCoeffsAtNode)
        call storage_free(rafcstab%h_DmatrixCoeffsAtEdge)
      end if
    end if

    ! Initialise auxiliary coefficients at nodes for n matrices
    if (rafcstab%NEQ .gt. 0) then
      Isize2D = (/n, rafcstab%NEQ/)   
      call storage_new('afcstab_initMatrixCoeffs', 'DmatrixCoeffsAtNode', Isize2D,&
          ST_DOUBLE, rafcstab%h_DmatrixCoeffsAtNode, ST_NEWBLOCK_ZERO)
    end if

    ! Initialise auxiliary coefficients at edges for n matrices
    if (rafcstab%NEDGE .gt. 0) then
      Isize3D = (/n, 2, rafcstab%NEDGE/)  
      call storage_new('afcstab_initMatrixCoeffs', 'DmatrixCoeffsAtEdge', Isize3D,&
          ST_DOUBLE, rafcstab%h_DmatrixCoeffsAtEdge, ST_NEWBLOCK_ZERO)
    end if

    ! Store number of auxiliary matrices
    rafcstab%nmatrixCoeffs = n

  end subroutine afcstab_initMatrixCoeffs

  ! ***************************************************************************

!<subroutine>

  subroutine afcstab_CopyMatrixCoeffs(rafcstab, Rmatrices, Iposition)

!<description>
    ! This subroutine copies the data of the given matrices into the
    ! arrays DmatrixCoeffsAtNode and DmatrixCoeffsAtEdge of the
    ! stabilisation structure. If these arrays are not allocated then
    ! they are allocated by this routine in correct size.
!</description>

!<input>
    ! Array of scalar coefficient matrices
    type(t_matrixScalar), dimension(:), intent(in) :: Rmatrices

    ! OPTIOANL: Array of integers which indicate the positions of the
    ! given matrices. If this parameter is not given, then the
    ! matrices are stored starting at position one.
    integer, dimension(:), intent(in), optional, target :: Iposition
!</input>

!<inputoutpu>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Iposition, p_Kdiagonal
    integer, dimension(2) :: Isize2D
    integer, dimension(3) :: Isize3D
    integer :: i,ij,ji,iedge,ipos,imatrix,nmatrices,nmaxpos


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_CopyMatrixCoeffs')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) then
      call afcstab_generateVerticesAtEdge(Rmatrices(1), rafcstab)
    end if

    ! Set number of matrices
    nmatrices = size(Rmatrices)

    ! Set pointer to matrix positions
    if (present(Iposition)) then
      p_Iposition => Iposition
    else
      allocate(p_Iposition(size(Rmatrices)))
      do imatrix = 1, nmatrices
        p_Iposition(imatrix) = imatrix
      end do
    end if

    ! Set maximum position
    nmaxpos = maxval(p_Iposition)

    ! Check if array Rmatrices and Iposition have the same size
    if (nmatrices .ne. size(p_Iposition)) then
      call output_line('size(Rmatrices) /= size(Iposition)!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_CopyMatrixCoeffs')
      call sys_halt()
    end if

    ! Note that only double precision matrices are supported
    if (Rmatrices(1)%cdataType .ne. ST_DOUBLE) then
      call output_line('Only double precision matrices are supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_CopyMatrixCoeffs')
      call sys_halt()
    end if
    
    ! Check if first matrix is compatible with the stabilisation structure
    if ((Rmatrices(1)%NEQ .ne. rafcstab%NEQ) .or.&
        ((Rmatrices(1)%NA-Rmatrices(1)%NEQ)/2 .ne. rafcstab%NEDGE) .or.&
        (Rmatrices(1)%cmatrixFormat .ne. rafcstab%cmatrixFormat)) then
      call output_line('Matrix is not compatible with stabilisation structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_CopyMatrixCoeffs')
      call sys_halt()
    end if

    ! Check if all matrices are compatible to the first one
    do imatrix = 2, nmatrices
      call lsyssc_isMatrixCompatible(Rmatrices(1), Rmatrices(imatrix))
    end do
    
    ! Check if auxiliary data arrays for matrix coefficients are initialised
    if ((rafcstab%h_DmatrixCoeffsAtNode .eq. ST_NOHANDLE) .or.&
        (rafcstab%h_DmatrixCoeffsAtEdge .eq. ST_NOHANDLE)) then
      call afcstab_initMatrixCoeffs(rafcstab, nmaxpos)
    else
      call storage_getsize(rafcstab%h_DmatrixCoeffsAtNode, Isize2D)
      call storage_getsize(rafcstab%h_DmatrixCoeffsAtEdge, Isize3D)
      if ((nmaxpos .gt. Isize2D(2)) .or. (nmaxpos .gt. Isize3D(3))) then
        call output_line('Number of matrices exceeds pre-initialised arrays!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_CopyMatrixCoeffs')
        call sys_halt()
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call storage_getbase_double2D(rafcstab%h_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtNode)
    call storage_getbase_double3D(rafcstab%h_DmatrixCoeffsAtEdge, p_DmatrixCoeffsAtEdge)

    ! Associate data of each matrix separateky
    do imatrix = 1, nmatrices
      
      ! Get matrix position
      ipos = p_Iposition(imatrix)

      ! What kind of matrix are we?
      select case(Rmatrices(imatrix)%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-------------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-------------------------------------------------------------------------
        
        ! Set pointer to matrix data
        call lsyssc_getbase_double(Rmatrices(imatrix), p_Ddata)
        
        ! Set diagonal pointer
        if (Rmatrices(imatrix)%cmatrixFormat .eq. LSYSSC_MATRIX7) then
          call lsyssc_getbase_Kld(Rmatrices(imatrix), p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(Rmatrices(imatrix), p_Kdiagonal)
        end if

        ! Loop over all rows and copy diagonal entries
        do i = 1, rafcstab%NEQ
          p_DmatrixCoeffsAtNode(ipos,i) = p_Ddata(p_Kdiagonal(i))
        end do
        
        ! Loop over all edges and copy off-diagonal entries
        do iedge = 1, rafcstab%NEDGE
          ij = p_IverticesAtEdge(3,iedge)
          ji = p_IverticesAtEdge(4,iedge)
          p_DmatrixCoeffsAtEdge(ipos,1,iedge) = p_Ddata(ij)
          p_DmatrixCoeffsAtEdge(ipos,2,iedge) = p_Ddata(ji)
        end do

      case(LSYSSC_MATRIX1)
        !-------------------------------------------------------------------------
        ! Matrix format 1
        !-------------------------------------------------------------------------
        
        ! Set pointer to matrix data
        call lsyssc_getbase_double(Rmatrices(imatrix), p_Ddata)

        ! Loop over all rows and copy diagonal entries
        do i = 1, rafcstab%NEQ
          p_DmatrixCoeffsAtNode(ipos,i) = p_Ddata(rafcstab%NEQ*(i-1)+i)
        end do
        
        ! Loop over all edges and copy off-diagonal entries
        do iedge = 1, rafcstab%NEDGE
          ij = p_IverticesAtEdge(3,iedge)
          ji = p_IverticesAtEdge(4,iedge)
          p_DmatrixCoeffsAtEdge(ipos,1,iedge) = p_Ddata(ij)
          p_DmatrixCoeffsAtEdge(ipos,2,iedge) = p_Ddata(ji)
        end do

      end select
    end do

    ! Set specifiert for auxiliary matrix coefficients
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)

  end subroutine afcstab_CopyMatrixCoeffs

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_isMatrixCompatibleSc(rafcstab, rmatrix, bcompatible)

!<description>
    ! This subroutine checks if a scalar matrix and a discrete 
    ! stabilisation structure are compatible to each other, 
    ! i.e. if they share the same structure, size and so on.
!</description>

!<input>
    ! The scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! The stabilisation structure
    type(t_afcstab), intent(in)      :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE
    ! depending on whether matrix and stabilisation are compatible or
    ! not.  If not given, an error will inform the user if the
    ! matrix/operator are not compatible and the program will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>

    ! Matrix/operator must have the same size
    if (rafcstab%NEQ   .ne. rmatrix%NEQ  .or.&
        rafcstab%NVAR  .ne. rmatrix%NVAR .or.&
        rafcstab%NEDGE .ne. (rmatrix%NA-rmatrix%NEQ)/2) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_isMatrixCompatibleSc')
        call sys_halt()
      end if
    end if
  end subroutine afcstab_isMatrixCompatibleSc
  
  ! *****************************************************************************

!<subroutine>

  subroutine afcstab_isMatrixCompatibleBl(rafcstab, rmatrixBlock, bcompatible)

!<description>
    ! This subroutine checks whether a block matrix and a discrete
    ! stabilisation structure are compatible to each other,
    ! i.e. if they share the same structure, size and so on.
    !
    ! If the matrix has only one block, then the scalar counterpart of this
    ! subroutine is called with the corresponding scalar submatrix.
    ! Otherwise, the matrix is required to possess group structure.
!</description>

!<input>
    ! block matrix
    type(t_matrixBlock), intent(in) :: rmatrixBlock

    ! stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>


    ! Check if matrix has only one block
    if ((rmatrixBlock%nblocksPerCol .eq. 1) .and.&
        (rmatrixBlock%nblocksPerRow .eq. 1)) then
      call afcstab_isMatrixCompatible(rafcstab,&
          rmatrixBlock%RmatrixBlock(1,1), bcompatible)
      return
    end if

    ! Check if number of columns equans number of rows
    if (rmatrixBlock%nblocksPerCol .ne.&
        rmatrixBlock%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_isMatrixCompatibleBl')
      call sys_halt()
    end if

    ! Check if matrix exhibits group structure
    if (rmatrixBlock%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Block matrix must have group structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_isMatrixCompatibleBl')
        call sys_halt()
      end if
    end if

    ! Matrix/operator must have the same size
    if (rafcstab%NVAR  .ne. rmatrixBlock%nblocksPerCol .or.&
        rafcstab%NEQ   .ne. rmatrixBlock%RmatrixBlock(1,1)%NEQ  .or.&
        rafcstab%NEDGE .ne. (rmatrixBlock%RmatrixBlock(1,1)%NA-&
                             rmatrixBlock%RmatrixBlock(1,1)%NEQ)/2) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_isMatrixCompatibleBl')
        call sys_halt()
      end if
    end if
  end subroutine afcstab_isMatrixCompatibleBl

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_isVectorCompatibleSc(rafcstab, rvector, bcompatible)

!<description>
    ! This subroutine checks if a vector and a stabilisation
    ! structure are compatible to each other, i.e., share the
    ! same structure, size and so on.
!</description>

!<input>
    ! The scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! Teh stabilisation structure
    type(t_afcstab), intent(in)      :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE
    ! depending on whether matrix and stabilisation are compatible or
    ! not. If not given, an error will inform the user if the
    ! matrix/operator are not compatible and the program will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>

    ! Matrix/operator must have the same size
    if (rafcstab%NEQ   .ne. rvector%NEQ .or.&
        rafcstab%NVAR  .ne. rvector%NVAR) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Vector/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_isVectorCompatibleSc')
        call sys_halt()
      end if
    end if
  end subroutine afcstab_isVectorCompatibleSc

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_isVectorCompatibleBl(rafcstab, rvectorBlock, bcompatible)

!<description>
    ! This subroutine checks whether a block vector and a stabilisation
    ! structure are compatible to each other, i.e., share the same
    ! structure, size and so on.
    !
    ! If the vectors has only one block, then the scalar counterpart of
    ! this subroutine is called with the corresponding scalar subvector.
!</description>

!<input>
    ! block vector
    type(t_vectorBlock), intent(in) :: rvectorBlock

    ! stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
    ! whether matrix and stabilisation are compatible or not.
    ! If not given, an error will inform the user if the matrix/operator are
    ! not compatible and the program will halt.
    logical, intent(out), optional :: bcompatible
!</output>
!</subroutine>

    ! Check if block vectors has just one block
    if (rvectorBlock%nblocks .eq. 1) then
      call afcstab_isVectorCompatible(rafcstab,&
          rvectorBlock%RvectorBlock(1), bcompatible)
      return
    end if

    ! Vector/operator must have the same size
    if (rafcstab%NEQ*rafcstab%NVAR .ne. rvectorBlock%NEQ) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Vector/Operator not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_isVectorCompatibleBl')
        call sys_halt()
      end if
    end if
  end subroutine afcstab_isVectorCompatibleBl

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_arrayBlock(rmatrix, rarray, bisFullMatrix)

!<description>
    ! This subroutine assigns the pointers of the array to the scalar
    ! submatrices of the given block matrix rmatrix.
    ! If the optional parameter bisFullMatrix is given, then this routine
    ! returns bisFullMatrix = .TRUE. if all blocks of rmatrix are associated.
    ! Otherwise, bisFullMatrix = .FALSE. is returned if only the diagonal
    ! blocks of the block matrix are associated.
!</description>

!<input>
    ! The block matrix
    type(t_matrixBlock), intent(in) :: rmatrix
!</input>

!<output>
    ! The array
    type(t_array), dimension(:,:), intent(out) :: rarray

    ! OPTIONAL: indicator for full block matrix
    logical, intent(out), optional :: bisFullMatrix
!</output>
!</subroutine>

    ! local variables
    integer :: iblock,jblock
    logical :: bisFull

    ! Check if array is compatible
    if (rmatrix%nblocksPerCol .ne. size(rarray,1) .or.&
        rmatrix%nblocksPerRow .ne. size(rarray,2)) then
      call output_line('Block matrix and array are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_getbase_arrayBlock')
      call sys_halt()
    end if

    ! Assign pointers
    bisFull = .true.
    do iblock = 1, rmatrix%nblocksPerCol
      do jblock = 1, rmatrix%nblocksPerRow
        if (lsyssc_isExplicitMatrix1D(rmatrix%RmatrixBlock(iblock,jblock))) then
          select case(rmatrix%RmatrixBlock(iblock,jblock)%cdataType)
          case (ST_DOUBLE)
            rarray(iblock,jblock)%idataType = ST_DOUBLE
            call lsyssc_getbase_double(&
                rmatrix%RmatrixBlock(iblock,jblock), rarray(iblock,jblock)%p_Ddata)
          case (ST_SINGLE)
            rarray(iblock,jblock)%idataType = ST_SINGLE
            call lsyssc_getbase_single(&
                rmatrix%RmatrixBlock(iblock,jblock), rarray(iblock,jblock)%p_Fdata)
          case default
            call output_line('Unsupported data type!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_getbase_arrayBlock')
            call sys_halt()
          end select
        else
          nullify(rarray(iblock,jblock)%p_Ddata)
          nullify(rarray(iblock,jblock)%p_Fdata)
          bisFull = .false.
        end if
      end do
    end do

    if (present(bisFullMatrix)) bisFullMatrix = bisFull

  end subroutine afcstab_getbase_arrayBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_arrayScalar(Rmatrices, rarray)

!<description>
    ! This subroutine assigns the pointers of the array to the array
    ! of scalar matrices.
!</description>

!<input>
    ! The array of scalar matrices
    type(t_matrixScalar), dimension(:), intent(in) :: Rmatrices
!</input>

!<output>
    ! The array
    type(t_array), dimension(:), intent(out) :: rarray
!</output>
!</subroutine>

    ! local variables
    integer :: i

    ! Check if array is compatible
    if (size(Rmatrices) .ne. size(rarray)) then
      call output_line('Array of matrices and array are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_getbase_arrayScalar')
      call sys_halt()
    end if

    ! Assing pointers
    do i = 1, size(Rmatrices)
      if (lsyssc_isExplicitMatrix1D(Rmatrices(i))) then
        select case(Rmatrices(i)%cdataType)
        case (ST_DOUBLE)
          rarray(i)%idataType = ST_DOUBLE
          call lsyssc_getbase_double(&
              Rmatrices(i), rarray(i)%p_Ddata)
        case (ST_SINGLE)
          rarray(i)%idataType = ST_SINGLE
          call lsyssc_getbase_single(&
              Rmatrices(i), rarray(i)%p_Fdata)
        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcstab_getbase_arrayScalar')
          call sys_halt()
        end select
      else
        nullify(rarray(i)%p_Ddata)
        nullify(rarray(i)%p_Fdata)
      end if
    end do

  end subroutine afcstab_getbase_arrayScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IvertAtEdgeIdx(rafcstab, p_IverticesAtEdgeIdx)

!<description>
    ! Returns a pointer to the index pointer for the vertices at edge structure
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the vertices at edge index structure
    ! NULL() if the stabilisation structure does not provide it.
    integer, dimension(:), pointer :: p_IverticesAtEdgeIdx
!</output>
!</subroutine>

    ! Do we vertices at edge structure at all?
    if (rafcstab%h_IverticesAtEdgeIdx .eq. ST_NOHANDLE) then
      nullify(p_IverticesAtEdgeIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IverticesAtEdgeIdx,&
        p_IverticesAtEdgeIdx)

  end subroutine afcstab_getbase_IvertAtEdgeIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)

!<description>
    ! Returns a pointer to the vertices at edge structure
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the vertices at edge structure
    ! NULL() if the stabilisation structure does not provide it.
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
!</output>
!</subroutine>

    ! Do we vertices at edge structure at all?
    if ((rafcstab%h_IverticesAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE             .eq. 0)) then
      nullify(p_IverticesAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_int2D(rafcstab%h_IverticesAtEdge,&
        p_IverticesAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_IverticesAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for vertices at edge structure
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for superdiagonal edge numbers.
    ! NULL() if the stabilisation structure does not provide it.
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsuperdiagEdgesIdx .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ .eq. 0)) then
      nullify(p_IsuperdiagEdgesIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsuperdiagEdgesIdx,&
        p_IsuperdiagEdgesIdx,rafcstab%NEQ+1)

  end subroutine afcstab_getbase_IsupdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for 
    ! subdiagonal edge numbers
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the index pointer for subdiagonal edge numbers
    ! NULL() if the stabilisation structure does not provide it.
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
!</output>
!</subroutine>

    ! Do we have index pointer for subdiagonal edge numbers at all?
    if ((rafcstab%h_IsubdiagEdgesIdx .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ                    .eq. 0)) then
      nullify(p_IsubdiagEdgesIdx)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsubdiagEdgesIdx,&
        p_IsubdiagEdgesIdx,rafcstab%NEQ+1)

  end subroutine afcstab_getbase_IsubdiagEdgeIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagEdges)

!<description>
    ! Returns a pointer to the subdiagonal edge number
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the subdiagonal edge numbers
    ! NULL() if the stabilisation structure does not provide it.
    integer, dimension(:), pointer :: p_IsubdiagEdges
!</output>
!</subroutine>

    ! Do we have an edge separator at all?
    if ((rafcstab%h_IsubdiagEdges .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE               .eq. 0)) then
      nullify(p_IsubdiagEdges)
      return
    end if
    
    ! Get the array
    call storage_getbase_int(rafcstab%h_IsubdiagEdges,&
        p_IsubdiagEdges,rafcstab%NEDGE)

  end subroutine afcstab_getbase_IsubdiagEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DcoeffsAtEdge(rafcstab,p_DcoefficientsAtEdge)

!<description>
    ! Returns a pointer to the double-valued edge data
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the double-valued edge data
    ! NULL() if the stabilisation structure does not provide it.
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
!</output>
!</subroutine>

    ! Do we have double-valued edge data at all?
    if ((rafcstab%h_DcoefficientsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                 .eq. 0)) then
      nullify(p_DcoefficientsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_double2D(rafcstab%h_DcoefficientsAtEdge,&
        p_DcoefficientsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_DcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_FcoeffsAtEdge(rafcstab,p_FcoefficientsAtEdge)

!<description>
    ! Returns a pointer to the single-valued edge data
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the double-valued edge data
    ! NULL() if the stabilisation structure does not provide it.
    real(SP), dimension(:,:), pointer :: p_FcoefficientsAtEdge
!</output>
!</subroutine>

    ! Do we have double-valued edge data at all?
    if ((rafcstab%h_DcoefficientsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                 .eq. 0)) then
      nullify(p_FcoefficientsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_single2D(rafcstab%h_DcoefficientsAtEdge,&
        p_FcoefficientsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_FcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DmatCoeffAtNode(rafcstab,p_DmatrixCoeffsAtNode)

!<description>
    ! Returns a pointer to the matrix coefficients at nodes
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the matrix coefficients at nodes
    ! NULL() if the stabilisation structure does not provide it.
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
!</output>
!</subroutine>

    ! Do we have matrix coefficients at nodes data at all?
    if ((rafcstab%h_DmatrixCoeffsAtNode .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ                   .eq. 0)) then
      nullify(p_DmatrixCoeffsAtNode)
      return
    end if
    
    ! Get the array
    call storage_getbase_double2D(rafcstab%h_DmatrixCoeffsAtNode,&
        p_DmatrixCoeffsAtNode,rafcstab%NEQ)

  end subroutine afcstab_getbase_DmatCoeffAtNode

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_FmatCoeffAtNode(rafcstab,p_FmatrixCoeffsAtNode)

!<description>
    ! Returns a pointer to the matrix coefficients at nodes
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the matrix coefficients at nodes
    ! NULL() if the stabilisation structure does not provide it.
    real(SP), dimension(:,:), pointer :: p_FmatrixCoeffsAtNode
!</output>
!</subroutine>

    ! Do we have matrix coefficients at nodes data at all?
    if ((rafcstab%h_DmatrixCoeffsAtNode .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEQ                   .eq. 0)) then
      nullify(p_FmatrixCoeffsAtNode)
      return
    end if
    
    ! Get the array
    call storage_getbase_single2D(rafcstab%h_DmatrixCoeffsAtNode,&
        p_FmatrixCoeffsAtNode,rafcstab%NEQ)

  end subroutine afcstab_getbase_FmatCoeffAtNode

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DmatCoeffAtEdge(rafcstab,p_DmatrixCoeffsAtEdge)

!<description>
    ! Returns a pointer to the matrix coefficients at edges
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the matrix coefficients at edges
    ! NULL() if the stabilisation structure does not provide it.
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
!</output>
!</subroutine>

    ! Do we have matrix coefficients at edges data at all?
    if ((rafcstab%h_DmatrixCoeffsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                 .eq. 0)) then
      nullify(p_DmatrixCoeffsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_double3D(rafcstab%h_DmatrixCoeffsAtEdge,&
        p_DmatrixCoeffsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_DmatCoeffAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_FmatCoeffAtEdge(rafcstab,p_FmatrixCoeffsAtEdge)

!<description>
    ! Returns a pointer to the matrix coefficients at edges
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the matrix coefficients at edges
    ! NULL() if the stabilisation structure does not provide it.
    real(SP), dimension(:,:,:), pointer :: p_FmatrixCoeffsAtEdge
!</output>
!</subroutine>

    ! Do we have matrix coefficients at edges data at all?
    if ((rafcstab%h_DmatrixCoeffsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                 .eq. 0)) then
      nullify(p_FmatrixCoeffsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_single3D(rafcstab%h_DmatrixCoeffsAtEdge,&
        p_FmatrixCoeffsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_FmatCoeffAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DboundsAtEdge(rafcstab,p_DboundsAtEdge)

!<description>
    ! Returns a pointer to the bounds at edges
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the bounds at edges
    ! NULL() if the stabilisation structure does not provide it.
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
!</output>
!</subroutine>

    ! Do we have bounds at edges data at all?
    if ((rafcstab%h_DboundsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE           .eq. 0)) then
      nullify(p_DboundsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_double2D(rafcstab%h_DboundsAtEdge,&
        p_DboundsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_DboundsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_FboundsAtEdge(rafcstab,p_FboundsAtEdge)

!<description>
    ! Returns a pointer to the bounds at edges
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the bounds at edges
    ! NULL() if the stabilisation structure does not provide it.
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
!</output>
!</subroutine>

    ! Do we have bounds at edges data at all?
    if ((rafcstab%h_DboundsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE           .eq. 0)) then
      nullify(p_FboundsAtEdge)
      return
    end if
    
    ! Get the array
    call storage_getbase_single2D(rafcstab%h_DboundsAtEdge,&
        p_FboundsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_FboundsAtEdge
 
  !*****************************************************************************

!<subroutine>

  subroutine afcstab_generateVerticesAtEdge(rmatrixTemplate, rafcstab)

!<description>
    ! This subroutine generates the list of edges which are
    ! characterized by their two endpoints (i,j) and the absolute
    ! position of matrix entries ij and ji. 
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(in) :: rmatrixTemplate
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge,p_IverticesAtEdgeAux
    integer, dimension(:), pointer :: p_IverticesAtEdgeIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal,p_Ksep
    integer :: h_Ksep,h_IverticesAtEdgeAux
    !$ integer :: nmaxcolors

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_generateVerticesAtEdge')
      call sys_halt()
    end if

    ! Check if edge structure is owned by the stabilisation structure
    if (iand(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGESTRUCTURE) .eq.&
        AFCSTAB_SHARE_EDGESTRUCTURE) then
      call output_line('Edge structure is not owned by stabilisation and '//&
          'therefore cannot be generated',&
          OU_CLASS_WARNING,OU_MODE_STD,'afcstab_generateVerticesAtEdge')
      return
    end if
    
    ! Set pointer to edge structure
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IvertAtEdgeIdx(rafcstab, p_IverticesAtEdgeIdx)

    ! If no OpenMP is used, then all edges belong to the same
    ! group. Otherwise, the edges will be reordered below.
    p_IverticesAtEdgeIdx     = rafcstab%NEDGE+1
    p_IverticesAtEdgeIdx(1)  = 1
    
    ! What kind of matrix are we?
    select case(rmatrixTemplate%cmatrixFormat)
    case(LSYSSC_MATRIX1)
      ! Generate edge structure for dense matrix
      call genEdgesMat1(rmatrixTemplate%NEQ, p_IverticesAtEdge)

      ! OpenMP-Extension: Perform edge-coloring to find groups of
      ! edges which can be processed in parallel, that is, the
      ! vertices of the edges in the group are all distinct

      !$ call storage_free(rafcstab%h_IverticesAtEdgeIdx)
      !$ call storage_new('afcstab_generateVerticesAtEdge','IverticesAtEdgeIdx',&
      !$     2*(rmatrixTemplate%NEQ-1)+1, ST_INT, rafcstab%h_IverticesAtEdgeIdx,&
      !$     ST_NEWBLOCK_ZERO)
      !$ call afcstab_getbase_IvertAtEdgeIdx(rafcstab, p_IverticesAtEdgeIdx)
      !$ h_IverticesAtEdgeAux = ST_NOHANDLE
      !$ call storage_copy(rafcstab%h_IverticesAtEdge, h_IverticesAtEdgeAux)
      !$ call storage_getbase_int2d(h_IverticesAtEdgeAux, p_IverticesAtEdgeAux)
      !$ call genEdgeIndex(rmatrixTemplate%NEQ, rmatrixTemplate%NEQ+2,&
      !$     p_IverticesAtEdgeIdx, p_IverticesAtEdgeAux, p_IverticesAtEdge)
      !$ call storage_free(h_IverticesAtEdgeAux)

      ! Set matrix format
      rafcstab%cmatrixFormat = LSYSSC_MATRIX1

    case(LSYSSC_MATRIX7,&
         LSYSSC_MATRIX7INTL)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrixTemplate, p_Kld)
      call lsyssc_getbase_Kcol(rmatrixTemplate, p_Kcol)

      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrixTemplate%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrixTemplate%NEQ+1)

      ! Generate edge structure
      call genEdgesMat7(p_Kld, p_Kcol, p_Ksep, p_IverticesAtEdge)

      ! Release diagonal separator
      call storage_free(h_Ksep)

      ! OpenMP-Extension: Perform edge-coloring to find groups of
      ! edges which can be processed in parallel, that is, the
      ! vertices of the edges in the group are all distinct
      
      !$ call storage_free(rafcstab%h_IverticesAtEdgeIdx)
      !$ nmaxcolors = computeMaxEdgeColors(rmatrixTemplate%NEQ, p_IverticesAtEdge)
      !$ call storage_new('afcstab_generateVerticesAtEdge','IverticesAtEdgeIdx',&
      !$     nmaxcolors+1, ST_INT, rafcstab%h_IverticesAtEdgeIdx, ST_NEWBLOCK_ZERO)
      !$ call afcstab_getbase_IvertAtEdgeIdx(rafcstab, p_IverticesAtEdgeIdx)
      !$ h_IverticesAtEdgeAux = ST_NOHANDLE
      !$ call storage_copy(rafcstab%h_IverticesAtEdge, h_IverticesAtEdgeAux)
      !$ call storage_getbase_int2d(h_IverticesAtEdgeAux, p_IverticesAtEdgeAux)
      !$ call genEdgeIndex(rmatrixTemplate%NEQ, nmaxcolors,&
      !$     p_IverticesAtEdgeIdx, p_IverticesAtEdgeAux, p_IverticesAtEdge)
      !$ call storage_free(h_IverticesAtEdgeAux)

      ! Set matrix format
      rafcstab%cmatrixFormat = LSYSSC_MATRIX7

      
    case(LSYSSC_MATRIX9,&
         LSYSSC_MATRIX9INTL)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrixTemplate, p_Kld)
      call lsyssc_getbase_Kcol(rmatrixTemplate, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrixTemplate, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrixTemplate%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrixTemplate%NEQ+1)

      ! Generate edge structure
      call genEdgesMat9(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_IverticesAtEdge)

      ! Release diagonal separator
      call storage_free(h_Ksep)

      ! OpenMP-Extension: Perform edge-coloring to find groups of
      ! edges which can be processed in parallel, that is, the
      ! vertices of the edges in the group are all distinct
      
      !$ call storage_free(rafcstab%h_IverticesAtEdgeIdx)
      !$ nmaxcolors = computeMaxEdgeColors(rmatrixTemplate%NEQ, p_IverticesAtEdge)
      !$ call storage_new('afcstab_generateVerticesAtEdge','IverticesAtEdgeIdx',&
      !$     nmaxcolors+1, ST_INT, rafcstab%h_IverticesAtEdgeIdx, ST_NEWBLOCK_ZERO)
      !$ call afcstab_getbase_IvertAtEdgeIdx(rafcstab, p_IverticesAtEdgeIdx)
      !$ h_IverticesAtEdgeAux = ST_NOHANDLE
      !$ call storage_copy(rafcstab%h_IverticesAtEdge, h_IverticesAtEdgeAux)
      !$ call storage_getbase_int2d(h_IverticesAtEdgeAux, p_IverticesAtEdgeAux)
      !$ call genEdgeIndex(rmatrixTemplate%NEQ, nmaxcolors,&
      !$     p_IverticesAtEdgeIdx, p_IverticesAtEdgeAux, p_IverticesAtEdge)
      !$ call storage_free(h_IverticesAtEdgeAux)

      ! Set matrix format
      rafcstab%cmatrixFormat = LSYSSC_MATRIX9


    case(LSYSSC_MATRIXD)
      call output_line('Edge-bases data structure cannot be generated&
         &from diagonal matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_generateVerticesAtEdge')
      call sys_halt()


    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_generateVerticesAtEdge')
      call sys_halt()
    end select

    ! Set state of stabiliation
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)
         
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Generate edge data structure for dense matrices in format 1

    subroutine genEdgesMat1(neq, IverticesAtEdge)

      integer, intent(in) :: neq
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      ! local variables
      integer :: i,j,iedge
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, neq

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do j = i+1, neq

          ! Increatse edge counter
          iedge = iedge+1
          
          ! Set node numbers i and j
          IverticesAtEdge(1, iedge) = i
          IverticesAtEdge(2, iedge) = j
          
          ! Set matrix positions ij and ji
          IverticesAtEdge(3, iedge) = neq*(i-1)+j
          IverticesAtEdge(4, iedge) = neq*(j-1)+i
          
        end do
      end do

    end subroutine genEdgesMat1

    !**************************************************************
    ! Generate edge data structure for matrices in format 7

    subroutine genEdgesMat7(Kld, Kcol, Ksep, IverticesAtEdge)

      integer, dimension(:), intent(in) :: Kld, Kcol
      integer, dimension(:), intent(inout) :: Ksep
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      ! local variables
      integer :: i,j,ij,ji,iedge
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, size(Kld)-1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1
       
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)
          
          ! Increatse edge counter
          iedge = iedge+1

          ! Set node numbers i and j
          IverticesAtEdge(1, iedge) = i
          IverticesAtEdge(2, iedge) = j
          
          ! Set matrix positions ij and ji
          IverticesAtEdge(3, iedge) = ij
          IverticesAtEdge(4, iedge) = ji
          
        end do
      end do

    end subroutine genEdgesMat7

    !**************************************************************
    ! Generate edge data structure for matrices in format 9

    subroutine genEdgesMat9(Kld, Kcol, Kdiagonal, Ksep, IverticesAtEdge)

      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      integer, dimension(:), intent(inout) :: Ksep
      integer, dimension(:,:), intent(inout) :: IverticesAtEdge

      ! local variables
      integer :: i,j,ij,ji,iedge
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, size(Kld)-1
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
       
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Increatse edge counter
          iedge = iedge+1

          ! Set node numbers i and j
          IverticesAtEdge(1, iedge) = i
          IverticesAtEdge(2, iedge) = j
          
          ! Set matrix positions ij and ji
          IverticesAtEdge(3, iedge) = ij
          IverticesAtEdge(4, iedge) = ji
          
        end do
      end do

    end subroutine genEdgesMat9

    !**************************************************************
    ! Compute the maximum number of colors required to color the
    ! edges of the matrix such that no two edges with the same color
    ! share a common node (This is the classical edge-coloring)

    function computeMaxEdgeColors(neq, IverticesAtEdge) result(ncolors)
      
      integer, intent(in) :: neq
      integer, dimension(:,:), intent(in) :: IverticesAtEdge

      integer :: ncolors

      integer, dimension(:), pointer :: p_IedgesAtNode
      integer :: h_IedgesAtNode
      integer :: i,j,iedge

      ! Allocate temporal memory
      h_IedgesAtNode = ST_NOHANDLE
      call storage_new('computeMaxEdgeColors', 'IedgesAtNode', neq,&
          ST_INT, h_IedgesAtNode, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(h_IedgesAtNode, p_IedgesAtNode)

      ! First shuffle pass: 
      ! Compute the number of edges at each node
      do iedge = 1, size(IverticesAtEdge,2)
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)
        p_IedgesAtNode(i) = p_IedgesAtNode(i)+1
        p_IedgesAtNode(j) = p_IedgesAtNode(j)+1
      end do

      ! Second shuffle pass:
      ! Compute the maximum number of edges at nodes
      ncolors = 0
      do i = 1, neq
        ncolors = max(ncolors, p_IedgesAtNode(i))
      end do

      ! Deallocate temporal memory
      call storage_free(h_IedgesAtNode)

      ! At the moment the greedy edge coloring algorithm is 
      ! used which may require 2x chromatic index of colors
      ncolors = 2*ncolors

    end function computeMaxEdgeColors

    !**************************************************************
    ! Generate index for edge data structure 

    subroutine genEdgeIndex(neq, nmaxcolors, IverticesAtEdgeIdx,&
        IverticesAtEdge, IverticesAtEdgeGrouped)

      integer, intent(in) :: neq,nmaxcolors
      integer, dimension(:), intent(inout) :: IverticesAtEdgeIdx
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:,:), intent(inout) :: IverticesAtEdgeGrouped

      integer, dimension(size(IverticesAtEdgeIdx)) :: Iaux
      
      integer, dimension(:), pointer :: p_InodeColor
      integer :: h_InodeColor
      integer :: i,j,iedge,jedge,icolor

      ! Allocate temporal memory
      h_InodeColor = ST_NOHANDLE
      call storage_new('genEdgeIndex', 'InodeColor', neq,&
          ST_INT, h_InodeColor, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(h_InodeColor, p_InodeColor)

      ! First shuffle pass: Loop over all edges and determine the
      ! next free color to be used for that edge
      IverticesAtEdgeIdx = 0
      call lalg_clearVector(p_InodeColor)

      edge1: do iedge = 1, size(IverticesAtEdge,2)
        ! Get node numbers
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Find next unused color for nodes I and J
        color1: do icolor = 0, nmaxcolors
          
          if (.not.btest(p_InodeColor(i), icolor) .and.&
              .not.btest(p_InodeColor(j), icolor)) then

            ! Set color marker for nodes I and J
            p_InodeColor(i) = ibset(p_InodeColor(i), icolor)
            p_InodeColor(j) = ibset(p_InodeColor(j), icolor)

            ! Increase number of edges in this color group by one.
            IverticesAtEdgeIdx(icolor+2) = IverticesAtEdgeIdx(icolor+2)+1

            ! That is it
            exit color1
          end if
        end do color1
      end do edge1

      ! Second shuffle pass: Compute starting position of each group
      IverticesAtEdgeIdx(1) = 1
      do icolor = 2, size(IverticesAtEdgeIdx)
        IverticesAtEdgeIdx(icolor) = IverticesAtEdgeIdx(icolor)&
                                   + IverticesAtEdgeIdx(icolor-1)
      end do

      ! Third shuffle pass: Reorder the edges physically
      call lalg_clearVector(p_InodeColor)

       edge2: do iedge = 1, size(IverticesAtEdge,2)
        ! Get node numbers
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Find next unused color for nodes I and J
        color2: do icolor = 0, nmaxcolors
          
          if (.not.btest(p_InodeColor(i), icolor) .and.&
              .not.btest(p_InodeColor(j), icolor)) then

            ! Set color marker for nodes I and J
            p_InodeColor(i) = ibset(p_InodeColor(i), icolor)
            p_InodeColor(j) = ibset(p_InodeColor(j), icolor)

            ! Get next free position
            jedge = IverticesAtEdgeIdx(icolor+1)
            
            ! Increase number of edges in this color group by one
            IverticesAtEdgeIdx(icolor+1) = IverticesAtEdgeIdx(icolor+1)+1
           
            ! Copy edge data
            IverticesAtEdgeGrouped(:,jedge) = IverticesAtEdge(:,iedge)

            ! That is it
            exit color2
          end if
        end do color2
      end do edge2
      
      ! Fourth shuffle pass:
      ! Adjust the starting index which was tainted in the third pass
      do icolor = size(IverticesAtEdgeIdx), 2, -1
        IverticesAtEdgeIdx(icolor) = IverticesAtEdgeIdx(icolor-1)
      end do
      IverticesAtEdgeIdx(1) = 1

      ! Deallocate temporal memory
      call storage_free(h_InodeColor)

    end subroutine genEdgeIndex
    
  end subroutine afcstab_generateVerticesAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_generateOffdiagEdges(rafcstab)

!<description>
    ! This subroutine generates the edge data structure
    ! (superdiagonal separator and subdiagonal edges)
    ! based on a given edge data structure.
!</description>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer   :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer   :: p_IsubdiagEdges
    integer, dimension(:), pointer   :: p_IsubdiagEdgesIdx
    integer :: iedge,nedge,istor
    integer :: ieq,jeq,neq
    integer :: isize

    ! Check if edge-based data structure is prepared
    if (iand(rafcstab%istabilisationSpec,AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) then
      call output_line('Stabilisation structure does not provide required ' //&
          'edge-based data structure',OU_CLASS_ERROR,OU_MODE_STD,&
          'afcstab_generateOffdiagEdge')
      call sys_halt()
    end if

    ! Store dimensions of stabilisation structure
    neq   = rafcstab%NEQ
    nedge = rafcstab%NEDGE

    ! Allocate memory (if required)
    if (rafcstab%h_IsuperdiagEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateOffdiagEdges','IsuperdiagEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsuperdiagEdgesIdx,ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rafcstab%h_IsuperdiagEdgesIdx,isize)
      if (isize < neq+1) then
        call storage_realloc('afcstab_generateOffdiagEdges',neq+1,&
            rafcstab%h_IsuperdiagEdgesIdx,ST_NEWBLOCK_ZERO,.false.)
      else
        call storage_clear(rafcstab%h_IsuperdiagEdgesIdx)
      end if
    end if

    ! Allocate memory (if required)
    if (rafcstab%h_IsubdiagEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateOffdiagEdges','IsubdiagEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsubdiagEdgesIdx,ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rafcstab%h_IsubdiagEdgesIdx,isize)
      if (isize < neq+1) then
        call storage_realloc('afcstab_generateOffdiagEdges',neq+1,&
            rafcstab%h_IsubdiagEdgesIdx,ST_NEWBLOCK_ZERO,.false.)
      else
        call storage_clear(rafcstab%h_IsubdiagEdgesIdx)
      end if
    end if

    ! Allocate memory (if required)
    if (rafcstab%h_IsubdiagEdges .eq. ST_NOHANDLE) then
      call storage_new('afcstab_generateOffdiagEdges','IsubdiagEdges',&
          nedge,ST_INT,rafcstab%h_IsubdiagEdges,ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rafcstab%h_IsubdiagEdges,isize)
      if (isize < nedge) then
        call storage_realloc('afcstab_generateOffdiagEdges',nedge,&
            rafcstab%h_IsubdiagEdges,ST_NEWBLOCK_NOINIT,.false.)
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab,p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IverticesAtEdge(rafcstab,p_IverticesAtEdge)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab,p_IsubdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagEdges)
    
    ! Count the number of edges in each row
    do iedge = 1, nedge
      
       ! Determine the start-point of the current edge
      ieq = min(p_IverticesAtEdge(1,iedge),&
                p_IverticesAtEdge(2,iedge))

      ! Determine the end-point of the current edge
      jeq = max(p_IverticesAtEdge(1,iedge),&
                p_IverticesAtEdge(2,iedge))

      ! Increase number of edges connected to these points by one
      p_IsuperdiagEdgesIdx(ieq+1) = p_IsuperdiagEdgesIdx(ieq+1)+1
      p_IsubdiagEdgesIdx(jeq+1)   = p_IsubdiagEdgesIdx(jeq+1)+1
    end do
    
    ! Reshuffle pass 1
    rafcstab%NNVEDGE = 0
    p_IsuperdiagEdgesIdx(1) = 1
    p_IsubdiagEdgesIdx(1) = 1
    
    do ieq = 2, neq+1
      ! Compute the maximum number of edges adjacent to one vertex.
      rafcstab%NNVEDGE = max(rafcstab%NNVEDGE,&
          p_IsuperdiagEdgesIdx(ieq)+p_IsubdiagEdgesIdx(ieq))
      
      ! Compute absolute starting positions of edge numbers connected
      ! to the current node IEQ
      p_IsuperdiagEdgesIdx(ieq) =&
          p_IsuperdiagEdgesIdx(ieq) + p_IsuperdiagEdgesIdx(ieq-1)
      p_IsubdiagEdgesIdx(ieq) =&
          p_IsubdiagEdgesIdx(ieq) + p_IsubdiagEdgesIdx(ieq-1)
    end do
    
    ! Store the subdiagonal edge numbers
    do ieq = 1, neq
      ! Loop over the edges located in the upper right triangular matrix
      do iedge = p_IsuperdiagEdgesIdx(ieq),&
                 p_IsuperdiagEdgesIdx(ieq+1)-1
        
        ! Determine the end-point of the edge, i.e. 
        ! the node which is not equal to IEQ
        jeq = p_IverticesAtEdge(1,iedge) + p_IverticesAtEdge(2,iedge) - ieq
        
        ! Get and update next free position
        istor = p_IsubdiagEdgesIdx(jeq)
        p_IsubdiagEdgesIdx(jeq) = istor+1
        
        ! Store global edge number
        p_IsubdiagEdges(istor)  = iedge
      end do
    end do
    
    ! Reshuffle pass 2:
    ! Adjust the index vector which was tainted in the above loop
    do ieq = neq+1, 2, -1
      p_IsubdiagEdgesIdx(ieq) = p_IsubdiagEdgesIdx(ieq-1)
    end do
    p_IsubdiagEdgesIdx(1) = 1

    ! Set specifier for extended edge structure
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES)

  end subroutine afcstab_generateOffdiagEdges
  
  !*****************************************************************************

!<subroutine>

  subroutine afcstab_generateExtSparsity(rmatrixSrc, rmatrixExtended)

!<description>
    ! This subroutine generates the extended sparsity pattern
    ! required to assemble the Jacobian matrix from a sparse finite
    ! element matrix stored in CRS format. The basic idea is as
    ! follows. Let A ={a_ij} denote the adjacency matrix which
    ! represents the undirected connectivity graph of the finite
    ! element matrix (rmatrix). The matrix coefficients are
    ! henceforth given by
    !   a_ij=1 if there exists some edge ij, a_ij=0 otherwise
    ! Now, let <tex>$Z=A^2$</tex>. Then z_ij>0 if and only if there
    ! exists a path of length two connecting nodes i and j. This is due
    !   z_ij=sum_k(a_ik*a_kj)>0 <=> ex. k : a_ik=1 and a_kj=1.
!</description>

!<input>
    ! scalar source matrix with standard sparsity pattern
    type(t_matrixScalar), intent(in)    :: rmatrixSrc
!</input>

!<inputoutput>
    ! scalar destination matrix with extended sparsity pattern
    type(t_matrixScalar), intent(inout) :: rmatrixExtended
!</inputoutput>
!</subroutine>
    
    ! Clear output matrix
    if (lsyssc_hasMatrixStructure(rmatrixExtended) .or.&
        lsyssc_hasMatrixContent(rmatrixExtended)) then
      call lsyssc_releaseMatrix(rmatrixExtended)
    end if

    ! Compute Z=A*A and let the connectivity graph of Z be the
    ! extended sparsity pattern of the Jacobian matrix
    call lsyssc_multMatMat(rmatrixSrc, rmatrixSrc, rmatrixExtended,&
                           .true., .true., .false.)

  end subroutine afcstab_generateExtSparsity

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_failsafeLimitingBlock(rafcstab, rlumpedMassMatrix,&
      Cvariables, dscale, nsteps, fcb_extractVariableBlock, rvector, rvectorTmp)

!<description>
    ! This subroutine performs failsafe flux limiting as described in
    ! the paper by Kuzmin, Moeller, Shadid, and Shashkov: "Failsafe
    ! flux limiting and constrained data projection for equations of
    ! gas dynamics" Journal of Computational Physics, vol. 229,
    ! Nov. 2010, p. 8766-8779.
!</description>

!<input>
    ! stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

    ! control variable names
    character(len=*), dimension(:), intent(in) :: Cvariables
    
    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! number of failsafe steps to be performed
    integer, intent(in) :: nsteps

    ! callback function to extract variables
    include 'intf_extractVariableBlock.inc'
!</input>

!<inputoutput>
    ! vector to be corrected
    type(t_vectorBlock), intent(inout) :: rvector

    ! temporal vector
    type(t_vectorBlock), intent(inout), target, optional :: rvectorTmp
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvectorLbound, rvectorUbound, rvectorControl
    type(t_vectorBlock), pointer :: p_rvectorTmp
    type(t_vectorScalar) :: rbeta
    real(DP), dimension(:), pointer :: p_Dlbound, p_Dubound, p_Dcontrol, p_Dflux
    real(DP), dimension(:), pointer :: p_DlumpedMassMatrix, p_Dalpha, p_Dbeta
    real(DP), dimension(:), pointer :: p_Ddata, p_DdataTmp
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_IverticesAtEdgeIdx
    integer :: istep,ivariable, nvariable

    ! Get number of control variables
    nvariable = size(Cvariables)
    
    ! Create temporal block vectors
    call lsysbl_createVectorBlock(rvectorControl, rafcstab%NEQ, nvariable, .false.)
    call lsysbl_createVectorBlock(rvectorLbound, rafcstab%NEQ, nvariable, .false.)
    call lsysbl_createVectorBlock(rvectorUbound, rafcstab%NEQ, nvariable, .false.)
    call lsyssc_createVector(rbeta, rafcstab%NEDGE, .false.)

    ! Set pointer to temporal vector or create new one
    if (present(rvectorTmp)) then
      p_rvectorTmp => rvectorTmp
    else
      allocate(p_rvectorTmp)
    end if

    ! Make a copy of the initial vector
    call lsysbl_copyVector(rvector, p_rvectorTmp)
    
    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)
    call afcstab_getbase_IvertAtEdgeIdx(rafcstab, p_IverticesAtEdgeIdx)
    call lsyssc_getbase_double(rlumpedMassMatrix, p_DlumpedMassMatrix)
    call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsysbl_getbase_double(rvectorControl, p_Dcontrol)
    call lsysbl_getbase_double(rvectorLbound, p_Dlbound)
    call lsysbl_getbase_double(rvectorUbound, p_Dubound)
    call lsysbl_getbase_double(p_rvectorTmp, p_DdataTmp)
    call lsysbl_getbase_double(rvector, p_Ddata)
    call lsyssc_getbase_double(rbeta, p_Dbeta)
    
    ! Initialise the nodal vectors by the given solution
    do ivariable = 1, nvariable
      call fcb_extractVariableBlock(rvector, trim(Cvariables(ivariable)),&
          rvectorControl%RvectorBlock(ivariable))
    end do
    
    ! Compute upper and lower nodal bounds
    call afcstab_computeBounds(p_IverticesAtEdgeIdx, p_IverticesAtEdge,&
        rafcstab%NEQ, nvariable, p_Dcontrol, p_Dlbound, p_Dubound)

    ! Initialize correction factors to "one time dscale"
    call lalg_setVector(p_Dbeta, dscale)
    
    ! Perform prescribed failsafe steps
    do istep = 1, nsteps+1

      ! Restore the initial vector
      call lalg_copyVector(p_DdataTmp, p_Ddata)

      ! Apply correction factors to solution vector
      if (rvector%nblocks .eq. 1) then
        call afcstab_applyCorrectionDim2(p_IverticesAtEdgeIdx,&
            p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            p_DlumpedMassMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Ddata)
      else
        call afcstab_applyCorrectionDim1(p_IverticesAtEdgeIdx,&
            p_IverticesAtEdge, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            p_DlumpedMassMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Ddata)
      end if

      ! For the last step, no failsafe checks are required
      if (istep .gt. nsteps) exit

      ! Convert solution to control variables
      do ivariable = 1, nvariable
        call fcb_extractVariableBlock(rvector, trim(Cvariables(ivariable)),&
            rvectorControl%RvectorBlock(ivariable))
      end do

      ! Compute failsafe correction factors and  exit if no further
      ! failsafe correction is required
      if (afcstab_computeFailsafeFactors(p_IverticesAtEdge, rafcstab%NEQ,&
          nvariable, dscale*real(nsteps-istep, DP)/real(nsteps, DP),&
          p_Dcontrol, p_Dlbound, p_Dubound, 1e-8_DP, p_Dbeta)) exit      
      
    end do

    ! Release temporal block vectors
    call lsysbl_releaseVector(rvectorControl)
    call lsysbl_releaseVector(rvectorLbound)
    call lsysbl_releaseVector(rvectorUbound)
    call lsyssc_releaseVector(rbeta)

    ! Release temporal vector
    if (.not.present(rvectorTmp)) then
      call lsysbl_releaseVector(p_rvectorTmp)
      deallocate(p_rvectorTmp)
    end if

  end subroutine afcstab_failsafeLimitingBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_failsafeLimitingArray(Rafcstab, rlumpedMassMatrix,&
      Cvariables, dscale, nsteps, fcb_extractVariableArray, Rvector, RvectorTmp)

!<description>
    ! This subroutine performs failsafe flux limiting as described in
    ! the paper by Kuzmin, Moeller, Shadid, and Shashkov: "Failsafe
    ! flux limiting and constrained data projection for equations of
    ! gas dynamics" Journal of Computational Physics, vol. 229,
    ! Nov. 2010, p. 8766-8779.
    !
    ! This subroutine works for arrays of block vectors, that is, is
    ! handles coupled problems consisting of individual subproblems
!</description>

!<input>
    ! array of stabilisation structures
    type(t_afcstab), dimension(:), intent(in) :: Rafcstab

    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rlumpedMassMatrix

    ! control variable names
    character(len=*), dimension(:), intent(in) :: Cvariables
    
    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! number of failsafe steps to be performed
    integer, intent(in) :: nsteps

    ! callback function to extract variables
    include 'intf_extractVariableArray.inc'
!</input>

!<inputoutput>
    ! array of vectors to be corrected
    type(t_vectorBlock), dimension(:), intent(inout) :: Rvector

    ! array of temporal vectors
    type(t_vectorBlock), dimension(:), intent(inout), target, optional :: rvectorTmp
!</inputoutput>
!</subroutine>

  end subroutine afcstab_failsafeLimitingArray
  
  !*****************************************************************************

!<function>
  
  elemental function afcstab_limit_unbounded(p, q, default) result(r)

!<description>
    ! This function computes the ratio Q/P. If the denominator is
    ! too small, then the default value is applied.
!</description>

!<input>
    ! (de)nominator
    real(DP), intent(in) :: p,q

    ! default value
    real(DP), intent(in) :: default
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>

    if (abs(p) > 1e-12_DP) then
      r = q/p
    else
      r = default
    end if
  end function afcstab_limit_unbounded

  !*****************************************************************************

!<function>
  
  elemental function afcstab_limit_bounded(p, q, default, dbound) result(r)

!<description>
    ! This function computes the limited ratio Q/P and bounds the
    ! result by the size of dbound. If the denominator is too small
    ! then the default value is applied.
!</description>

!<input>
    ! (de)nominator
    real(DP), intent(in) :: p,q
    
    ! default value
    real(DP), intent(in) :: default

    ! upper bound
    real(DP), intent(in) :: dbound
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>
    
    if (abs(p) > 1e-12_DP) then
      r = min(q/p, dbound)
    else
      r = default
    end if
  end function afcstab_limit_bounded

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_computeBounds(IverticesAtEdgeIdx, IverticesAtEdge,&
      neq, nvar, Dx, Dlbound, Dubound)
      
!<description>
    ! This subroutine computes the local upper and lower bounds based on
    ! the solution vector Dx evaluated at the neighbouring nodes
!</description>

!<input>
    ! Index pointer to the edge structure
    integer, dimension(:), intent(in) :: IverticesAtEdgeIdx

    ! Nodal numbers of edge-neighbours
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! Solution vector
    real(DP), dimension(neq,nvar), intent(in) :: Dx

    ! Number of equations
    integer, intent(in) :: neq

    ! Number of variables
    integer, intent(in) :: nvar
!</input>

!<inputoutput>
    ! Vector of upper and lower bounds
    real(DP), dimension(neq,nvar), intent(inout) :: Dlbound, Dubound
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: i,iedge,igroup,j
    
    call lalg_copyVector(Dx, Dlbound)
    call lalg_copyVector(Dx, Dubound)

    !$omp parallel default(shared) private(i,j)&
    !$omp if (IverticesAtEdgeIdx(size(IverticesAtEdgeIdx)) > AFCSTAB_NEDGEMIN_OMP)
    
    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronize memory access
    do igroup = 1, size(IverticesAtEdgeIdx)-1
      
      ! Do nothing for empty groups
      if (IverticesAtEdgeIdx(igroup+1)-IverticesAtEdgeIdx(igroup) .le. 0) cycle
      
      ! Loop over all edges
      !$omp do
      do iedge = IverticesAtEdgeIdx(igroup), IverticesAtEdgeIdx(igroup+1)-1
          
        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Compute minimum/maximum value of neighboring nodes
        Dlbound(i,:) = min(Dlbound(i,:), Dx(j,:))
        Dlbound(j,:) = min(Dlbound(j,:), Dx(i,:))
        Dubound(i,:) = max(Dubound(i,:), Dx(j,:))
        Dubound(j,:) = max(Dubound(j,:), Dx(i,:))
      end do
      !$omp end do

    end do ! igroup
    !$omp end parallel
    
  end subroutine afcstab_computeBounds

  !*****************************************************************************

!<function>

  function afcstab_computeFailsafeFactors(IverticesAtEdge, neq, nvar,&
      dscale, Dx, Dlbound, Dubound, dtolerance, Dbeta) result(baccept)

!<description>
    ! This function computes the failsafe correction factors
!</description>

!<input>
    ! Nodal numbers of edge-neighbours
    integer, dimension(:,:), intent(in) :: IverticesAtEdge
    
    ! Solution vector
    real(DP), dimension(neq,nvar), intent(in) :: Dx

    ! Vectors of upper and lower bounds
    real(DP), dimension(neq,nvar), intent(in) :: Dlbound, Dubound
    
    ! Scaling parameter
    real(DP), intent(in) :: dscale

    ! Tolerance parameter
    real(DP), intent(in) :: dtolerance

    ! Number of equations
    integer, intent(in) :: neq

    ! Number of variables
    integer, intent(in) :: nvar
!</input>

!<inputoutput>
    ! Failsafe correction factors
    real(DP), dimension(:), intent(inout) :: Dbeta
!</inputoutput>

!<result>
    ! TRUE if no additional failsafe correction step is required
    logical :: baccept
!</result>
!</function>
      
    ! local variables
    integer :: iedge,i,j,ivar
    
    ! Initialisation
    baccept = .true.
    
    ! Loop over all edges
    !$omp parallel do default(shared) private(i,j,ivar)&
    !$omp if(size(IverticesAtEdge,2) > AFCSTAB_NEDGEMIN_OMP)&
    !$omp reduction(.and.:baccept)
    do iedge = 1, size(IverticesAtEdge,2)
      
      ! Get node numbers
      i  = IverticesAtEdge(1, iedge)
      j  = IverticesAtEdge(2, iedge)
      
      ! Loop over all variables
      do ivar = 1, nvar
        
        if ((Dx(i,ivar) .lt. Dlbound(i,ivar)-dtolerance) .or.&
            (Dx(j,ivar) .lt. Dlbound(j,ivar)-dtolerance) .or.&
            (Dx(i,ivar) .gt. Dubound(i,ivar)+dtolerance) .or.&
            (Dx(j,ivar) .gt. Dubound(j,ivar)+dtolerance)) then
          Dbeta(iedge) = dscale
          baccept = .false.
        end if
      end do
    end do
    !$omp end parallel do
    
  end function afcstab_computeFailsafeFactors
  
  !*****************************************************************************

!<subroutine>

  subroutine afcstab_applyCorrectionDim1(IverticesAtEdgeIdx, IverticesAtEdge,&
      NEDGE, NEQ, NVAR, ML, Dalpha, Dbeta, Dflux, Dx)

!<description>
    ! This subroutine applies the failsafe correction factors Dbeta and 
    ! the regular correction factors Dalpha to the raw antidiffusive
    ! fluxes Dflux and adds the result to the solution vectors Dx
    ! scaled by the lumped mass matrix ML.
    ! Solution vector is stored with leading dimension.
!</description>
   
!<input>
    ! Index pointer to the edge structure
    integer, dimension(:), intent(in) :: IverticesAtEdgeIdx

    ! Nodal numbers of edge-neighbours
    integer, dimension(:,:), intent(in) :: IverticesAtEdge
    
    ! Lumped mass matrix
    real(DP), dimension(:), intent(in) :: ML

    ! Regular correction factors
    real(DP), dimension(:), intent(in) :: Dalpha
    
    ! Failsafe correction factors
    real(DP), dimension(:), intent(in) :: Dbeta

    ! Raw antidiffusive fluxes
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux

    ! Number of edges, euqations, and variables
    integer, intent(in) :: NEDGE, NEQ, NVAR
!</input>

!<inputoutput>
    ! Solution vector to be corrected
    real(DP), dimension(NEQ,NVAR), intent(inout) :: Dx
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(NVAR) :: F_ij
    integer :: i,iedge,igroup,j

    !$omp parallel default(shared) private(i,j,F_ij)&
    !$omp if (NEDGE > AFCSTAB_NEDGEMIN_OMP)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronize memory access
    do igroup = 1, size(IverticesAtEdgeIdx)-1
      
      ! Do nothing for empty groups
      if (IverticesAtEdgeIdx(igroup+1)-IverticesAtEdgeIdx(igroup) .le. 0) cycle
      
      ! Loop over all edges
      !$omp do
      do iedge = IverticesAtEdgeIdx(igroup), IverticesAtEdgeIdx(igroup+1)-1
        
        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Compute portion of corrected antidiffusive flux
        F_ij = Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)
        
        ! Remove flux from solution
        Dx(i,:) = Dx(i,:) + F_ij/ML(i)
        Dx(j,:) = Dx(j,:) - F_ij/ML(j)
      end do
      !$omp end do

    end do ! igroup
    !$omp end parallel
    
  end subroutine afcstab_applyCorrectionDim1

   !*****************************************************************************

!<subroutine>

  subroutine afcstab_applyCorrectionDim2(IverticesAtEdgeIdx, IverticesAtEdge,&
      NEDGE, NEQ, NVAR, ML, Dalpha, Dbeta, Dflux, Dx)

!<description>
    ! This subroutine applies the failsafe correction factors Dbeta and 
    ! the regular correction factors Dalpha to the raw antidiffusive
    ! fluxes Dflux and adds the result to the solution vectors Dx
    ! scaled by the lumped mass matrix ML.
    ! Solution vector is stored with trailing dimension.
!</description>
   
!<input>
    ! Index pointer to the edge structure
    integer, dimension(:), intent(in) :: IverticesAtEdgeIdx

    ! Nodal numbers of edge-neighbours
    integer, dimension(:,:), intent(in) :: IverticesAtEdge
    
    ! Lumped mass matrix
    real(DP), dimension(:), intent(in) :: ML

    ! Regular correction factors
    real(DP), dimension(:), intent(in) :: Dalpha
    
    ! Failsafe correction factors
    real(DP), dimension(:), intent(in) :: Dbeta

    ! Raw antidiffusive fluxes
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux

    ! Number of edges, euqations, and variables
    integer, intent(in) :: NEDGE, NEQ, NVAR
!</input>

!<inputoutput>
    ! Solution vector to be corrected
    real(DP), dimension(NVAR,NEQ), intent(inout) :: Dx
!</inputoutput>
!</subroutine> 

    ! local variables
    real(DP), dimension(NVAR) :: F_ij
    integer :: i,iedge,igroup,j

    !$omp parallel default(shared) private(i,j,F_ij)&
    !$omp if (NEDGE > AFCSTAB_NEDGEMIN_OMP)

    ! Loop over the edge groups and process all edges of one group
    ! in parallel without the need to synchronize memory access
    do igroup = 1, size(IverticesAtEdgeIdx)-1
      
      ! Do nothing for empty groups
      if (IverticesAtEdgeIdx(igroup+1)-IverticesAtEdgeIdx(igroup) .le. 0) cycle
      
      ! Loop over all edges
      !$omp do
      do iedge = IverticesAtEdgeIdx(igroup), IverticesAtEdgeIdx(igroup+1)-1
        
        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Compute portion of corrected antidiffusive flux
        F_ij = Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)
        
        ! Remove flux from solution
        Dx(:,i) = Dx(:,i) + F_ij/ML(i)
        Dx(:,j) = Dx(:,j) - F_ij/ML(j)
      end do
      !$omp end do

    end do ! igroup
    !$omp end parallel

  end subroutine afcstab_applyCorrectionDim2

  

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyH2D_IverticesAtEdge(rafcstab, btranspose)

!<description>
    ! This subroutine copies the vertices at edge structure from the
    ! host memory to the memory of the coprocessor device. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>


    if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemory(rafcstab%h_IverticesAtEdge,&
        ST_SYNCBLOCK_COPY_H2D, btranspose)

  end subroutine afcstab_copyH2D_IverticesAtEdge

!*****************************************************************************

!<subroutine>

  subroutine afcstab_copyD2H_IverticesAtEdge(rafcstab, btranspose)

!<description>
    ! This subroutine copies the vertices at edge structure from the
    ! memory of the coprocessor device to the host memory. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>

    if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemory(rafcstab%h_IverticesAtEdge,&
        ST_SYNCBLOCK_COPY_D2H, btranspose)

  end subroutine afcstab_copyD2H_IverticesAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyH2D_DmatCoeffAtEdge(rafcstab, btranspose)

!<description>
    ! This subroutine copies the off-diagonal entries of the auxiliary
    ! constant matrix coefficients from the host memory to the memory
    ! of the coprocessor device. If no device is available, then an
    ! error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>


    if (rafcstab%h_DmatrixCoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemory(rafcstab%h_DmatrixCoeffsAtEdge,&
        ST_SYNCBLOCK_COPY_H2D, btranspose)

  end subroutine afcstab_copyH2D_DmatCoeffAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyD2H_DmatCoeffAtEdge(rafcstab, btranspose)

!<description>
    ! This subroutine copies the off-diagonal entries of the auxiliary
    ! constant matrix coefficients from the memory of the coprocessor
    ! device to the host memory. If no device is available, then an
    ! error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! If true then the memory is transposed.
    logical, intent(in) :: btranspose
!</input>
!</subroutine>

    if (rafcstab%h_DmatrixCoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemory(rafcstab%h_DmatrixCoeffsAtEdge,&
        ST_SYNCBLOCK_COPY_D2H, btranspose)

  end subroutine afcstab_copyD2H_DmatCoeffAtEdge

  !*****************************************************************************

  !<subroutine>

  subroutine afcstab_allocEdgeStructure(rafcstab, ndata)

!<description>
    ! This subroutine allocates the edge data structure
!</description>

!<input>
    ! Number of data items per edge
    integer, intent(in) :: ndata
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize

    ! Handle for IverticesAtEdgeIdx
    if (rafcstab%h_IverticesAtEdgeIdx .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IverticesAtEdgeIdx)
    call storage_new('afcstab_allocEdgeStructure', 'IverticesAtEdgeIdx',&
        2, ST_INT, rafcstab%h_IverticesAtEdgeIdx, ST_NEWBLOCK_NOINIT)
    
    ! Handle for IverticesAtEdge: (/i,j,ij,ji/)
    Isize = (/ndata, rafcstab%NEDGE/)
    if (rafcstab%h_IverticesAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IverticesAtEdge)
    call storage_new('afcstab_allocEdgeStructure', 'IverticesAtEdge',&
        Isize, ST_INT, rafcstab%h_IverticesAtEdge, ST_NEWBLOCK_NOINIT)
    
    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_EDGESTRUCTURE))

  end subroutine afcstab_allocEdgeStructure

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_allocCoeffsAtEdge(rafcstab, ncoeffs, cdataType)

!<description>
    ! This subroutine allocates the coefficients at edge structuree
!</description>

!<input>
    ! Number of coefficients at edge
    integer, intent(in) :: ncoeffs
    
    ! OPTIONAL: data type which overwrites internal setting
    integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer, dimension(2) :: Isize
    integer :: ctype

    ctype = rafcstab%cdataType
    if (present(cdataType)) ctype = cdataType

    rafcstab%ncoeffs = ncoeffs
    Isize = (/rafcstab%ncoeffs, rafcstab%NEDGE/)
    if (rafcstab%h_DcoefficientsAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_DcoefficientsAtEdge)
    call storage_new('afcstab_allocCoeffsAtEdge', 'DcoefficientsAtEdge',&
        Isize, cType, rafcstab%h_DcoefficientsAtEdge, ST_NEWBLOCK_NOINIT)
    
    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_EDGEVALUES))

  end subroutine afcstab_allocCoeffsAtEdge

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_allocBoundsAtEdge(rafcstab, cdataType)

!<description>
    ! This subroutine allocates the bounds at edge structuree
!</description>

!<input>
    ! OPTIONAL: data type which overwrites internal setting
    integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer, dimension(2) :: Isize
    integer :: ctype

    ctype = rafcstab%cdataType
    if (present(cdataType)) ctype = cdataType
    
    Isize = (/2, rafcstab%NEDGE/)
    if (rafcstab%h_DboundsAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_DboundsAtEdge)
    call storage_new('afcstab_allocBoundsAtEdge', 'DboundsAtEdge',&
        Isize, cType, rafcstab%h_DboundsAtEdge, ST_NEWBLOCK_NOINIT)
    
    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_EDGEBOUNDS))

  end subroutine afcstab_allocBoundsAtEdge

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_allocVectorsPQR(rafcstab, cdataType, ballocCommonQ)

!<description>
    ! This subroutine allocates the nodal vectors P, Q, and R for '+'
    ! and '-'. If ballocCommonQ = .TRUE., then an additional nodal
    ! vector Q is allocated commonly used by Q '+' and '-'
!</description>

!<input>
    ! OPTIONAL: data type which overwrites internal setting
    integer, intent(in), optional :: cdataType

    ! OPTIONAL: switch for common nodal vector Q
    logical, intent(in) , optional :: ballocCommonQ
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ctype
    logical :: ballocQ

    ctype = rafcstab%cdataType
    if (present(cdataType)) ctype = cdataType

    ballocQ = .false.
    if (present(ballocCommonQ)) ballocQ = ballocCommonQ

    if (.not.associated(rafcstab%p_rvectorPp)) allocate(rafcstab%p_rvectorPp)
    if (.not.associated(rafcstab%p_rvectorPm)) allocate(rafcstab%p_rvectorPm)
    if (.not.associated(rafcstab%p_rvectorQp)) allocate(rafcstab%p_rvectorQp)
    if (.not.associated(rafcstab%p_rvectorQm)) allocate(rafcstab%p_rvectorQm)
    if (.not.associated(rafcstab%p_rvectorRp)) allocate(rafcstab%p_rvectorRp)
    if (.not.associated(rafcstab%p_rvectorRm)) allocate(rafcstab%p_rvectorRm)
    
    if (rafcstab%p_rvectorPp%NEQ .ne. 0)&
        call lsyssc_releaseVector(rafcstab%p_rvectorPp)
    if (rafcstab%p_rvectorPm%NEQ .ne. 0)&
        call lsyssc_releaseVector(rafcstab%p_rvectorPm)
    if (rafcstab%p_rvectorQp%NEQ .ne. 0)&
        call lsyssc_releaseVector(rafcstab%p_rvectorQp)
    if (rafcstab%p_rvectorQm%NEQ .ne. 0)&
        call lsyssc_releaseVector(rafcstab%p_rvectorQm)
    if (rafcstab%p_rvectorRp%NEQ .ne. 0)&
        call lsyssc_releaseVector(rafcstab%p_rvectorRp)
    if (rafcstab%p_rvectorRm%NEQ .ne. 0)&
        call lsyssc_releaseVector(rafcstab%p_rvectorRm)

    if (ballocQ) then
      if (.not.associated(rafcstab%p_rvectorQ)) allocate(rafcstab%p_rvectorQ)
      if (rafcstab%p_rvectorQ%NEQ .ne. 0)&
          call lsyssc_releaseVector(rafcstab%p_rvectorQ)
    end if

    if (rafcstab%NVARtransformed .ne. 1) then
     
      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ctype)
    
      if (ballocQ)&
          call lsyssc_createVector(rafcstab%p_rvectorQ, rafcstab%NEQ,&
          rafcstab%NVARtransformed, .false., ctype)
  
    else

      call lsyssc_createVector(rafcstab%p_rvectorPp, rafcstab%NEQ,&
          .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorPm, rafcstab%NEQ,&
          .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorQp, rafcstab%NEQ,&
          .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorQm, rafcstab%NEQ,&
          .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorRp, rafcstab%NEQ,&
          .false., ctype)
      call lsyssc_createVector(rafcstab%p_rvectorRm, rafcstab%NEQ,&
          .false., ctype)

      if (ballocQ)&
          call lsyssc_createVector(rafcstab%p_rvectorQ, rafcstab%NEQ,&
          .false., ctype)

    end if

    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_ADINCREMENTS))
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_NODEBOUNDS))
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_NODELIMITER))

  end subroutine afcstab_allocVectorsPQR

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_allocFlux0(rafcstab, cdataType)

!<description>
    ! This subroutine allocates the edgewise flux vector flux0
!</description>

!<input>
    ! OPTIONAL: data type which overwrites internal setting
    integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ctype

    ctype = rafcstab%cdataType
    if (present(cdataType)) ctype = cdataType

    if (.not.associated(rafcstab%p_rvectorFlux0))&
        allocate(rafcstab%p_rvectorFlux0)
    
    if (rafcstab%NVAR .ne. 1) then

      call lsyssc_createVector(rafcstab%p_rvectorFlux0, rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ctype)
      
    else

      call lsyssc_createVector(rafcstab%p_rvectorFlux0, rafcstab%NEDGE,&
          .false., ctype)

    end if
    
    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_ADFLUXES))

  end subroutine afcstab_allocFlux0
  
  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_allocFlux(rafcstab, cdataType)

!<description>
    ! This subroutine allocates the edgewise flux vector flux
!</description>

!<input>
    ! OPTIONAL: data type which overwrites internal setting
    integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ctype

    ctype = rafcstab%cdataType
    if (present(cdataType)) ctype = cdataType

    if (.not.associated(rafcstab%p_rvectorFlux))&
        allocate(rafcstab%p_rvectorFlux)

    if (rafcstab%NVAR .ne. 1) then

      call lsyssc_createVector(rafcstab%p_rvectorFlux,  rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ctype)
      
    else

      call lsyssc_createVector(rafcstab%p_rvectorFlux,  rafcstab%NEDGE,&
          .false., ctype)

    end if

    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_ADFLUXES))
    
  end subroutine afcstab_allocFlux

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_allocFluxPrel(rafcstab, cdataType)

!<description>
    ! This subroutine allocates the edgewise flux vector fluxPrel
!</description>

!<input>
    ! OPTIONAL: data type which overwrites internal setting
    integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ctype

    ctype = rafcstab%cdataType
    if (present(cdataType)) ctype = cdataType

    if (.not.associated(rafcstab%p_rvectorFluxPrel))&
        allocate(rafcstab%p_rvectorFluxPrel)

    if (rafcstab%NVAR .ne. 1) then

      call lsyssc_createVector(rafcstab%p_rvectorFluxPrel,  rafcstab%NEDGE,&
          rafcstab%NVAR, .false., ctype)
      
    else

      call lsyssc_createVector(rafcstab%p_rvectorFluxPrel,  rafcstab%NEDGE,&
          .false., ctype)

    end if

    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_ADFLUXES))
    
  end subroutine afcstab_allocFluxPrel

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_allocAlpha(rafcstab, cdataType)

!<description>
    ! This subroutine allocates the edgewise correction factors alpha
!</description>

!<input>
    ! OPTIONAL: data type which overwrites internal setting
    integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ctype

    ctype = rafcstab%cdataType
    if (present(cdataType)) ctype = cdataType

    if (.not.associated(rafcstab%p_rvectorAlpha))&
        allocate(rafcstab%p_rvectorAlpha)

    call lsyssc_createVector(rafcstab%p_rvectorAlpha,  rafcstab%NEDGE,&
        .false., ctype)
    
    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_EDGELIMITER))

  end subroutine afcstab_allocAlpha

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_buildBoundsLPT1D(rafcstab, rmatrixCx, rmatrixM, dscale)

!<description>
    ! This subroutine computes the solution bounds for
    ! linearity-preserving flux correction in 1D
!</description>

!<input>
    ! Coefficient matrix
    type(t_matrixScalar), intent(in) :: rmatrixCx

    ! Diagonal scaling matrix
    type(t_matrixScalar), intent(in) :: rmatrixM

    ! Global scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
    real(DP), dimension(:), pointer :: p_DdataM, p_DdataCx
    real(SP), dimension(:), pointer :: p_FdataM, p_FdataCx
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal, p_Kld, p_Kcol

    ! Check spatial discretisation structure
    if (.not. rmatrixCx%bidenticalTrialAndTest) then
      call output_line('Test and trial functions must be identical!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Set pointer to triangulation; note that for more general
    ! discretisations the coordinates of the degrees of freedom must
    ! be computed and cannot be adopted from the vertex coordinates
    p_rtriangulation => rmatrixCx%p_rspatialDiscrTrial%p_rtriangulation

    ! Check if this is a 1D triangulation
    if (p_rtriangulation%ndim .ne. NDIM1D) then
      call output_line('Only 1D triangulations are supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if
    
    ! Get vertex coordinates
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(rmatrixCx, rafcstab)
    end if

    ! Check if stabilisation structure provides array to store bounds;
    ! otherwise allocate new memory using the predefined data type
    if (rafcstab%h_DboundsAtEdge .eq. ST_NOHANDLE)&
        call afcstab_allocBoundsAtEdge(rafcstab, rafcstab%cdataType)

    ! Check if source matrix is compatible with stabilisation structure
    if ((rmatrixCx%NA-rmatrixCx%NEQ .ne. rafcstab%NEDGE * 2) .or.&
        (rmatrixCx%cdataType .ne. rafcstab%cdataType)) then
      call output_line('Matrix/stabilisation structure is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if
    
    ! Check if scaling matrix rmatrixM is diagonal and 
    ! has the same data type as the source matrix
    if ((rmatrixM%cmatrixFormat .ne. LSYSSC_MATRIXD)    .or.&
        (rmatrixM%cdataType .ne. rmatrixCx%cdataType) .or.&
        (rmatrixM%NEQ .ne. rmatrixCx%NEQ)) then
      call output_line('Coefficient matrices are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)   

    ! What type of matrix are we?
    select case(rmatrixCx%cmatrixFormat)
    case (LSYSSC_MATRIX7)
      call lsyssc_getbase_Kld(rmatrixCx, p_Kld)
      
      ! What type of data format are we?
      select case(rmatrixCx%cdatatype)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixM, p_DdataM)
        call lsyssc_getbase_double (rmatrixCx, p_DdataCx)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
        
        call doBoundsMat7Dble(p_IverticesAtEdge, p_Kld,&
            p_DdataCx, p_DdataM, p_DvertexCoords, dscale, p_DboundsAtEdge)
        
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat7Sngl(p_IverticesAtEdge, p_Kld,&
            p_FdataCx, p_FdataM, p_DvertexCoords, dscale, p_FboundsAtEdge)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
        call sys_halt()
      end select

    case (LSYSSC_MATRIX9)
      call lsyssc_getbase_Kdiagonal(rmatrixCx, p_Kdiagonal)
      call lsyssc_getbase_Kld(rmatrixCx, p_Kld)
      
      ! What type of data format are we?
      select case(rmatrixCx%cdatatype)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixM, p_DdataM)
        call lsyssc_getbase_double (rmatrixCx, p_DdataCx)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
        
        call doBoundsMat9Dble(p_IverticesAtEdge, p_Kdiagonal, p_Kld,&
            p_DdataCx, p_DdataM, p_DvertexCoords, dscale, p_DboundsAtEdge)
        
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat9Sngl(p_IverticesAtEdge, p_Kdiagonal, p_Kld,&
            p_FdataCx, p_FdataM, p_DvertexCoords, dscale, p_FboundsAtEdge)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end select

  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7Dble(IverticesAtEdge, Kld,&
        DdataCx, DdataM, DvertexCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DvertexCoords(1,i) - DvertexCoords(1,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*diff)
        end do
        
        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*diff)
        end do
        
        ! Store result into second entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat7Dble

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9Dble(IverticesAtEdge, Kdiagonal, Kld,&
        DdataCx, DdataM, DvertexCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DvertexCoords(1,i) - DvertexCoords(1,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(DdataCx(ik)*diff)
        end do
        
        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*diff)
        end do
        
        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(DdataCx(jk)*diff)
        end do
        
        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*diff)
        end do
        
        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat9Dble

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7Sngl(IverticesAtEdge, Kld,&
        FdataCx, FdataM, DvertexCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DvertexCoords(1,i) - DvertexCoords(1,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*diff)
        end do
        
        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*diff)
        end do
        
        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat7Sngl

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9Sngl(IverticesAtEdge, Kdiagonal, Kld,&
        FdataCx, FdataM, DvertexCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DvertexCoords(1,i) - DvertexCoords(1,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(FdataCx(ik)*diff)
        end do
        
        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*diff)
        end do
        
        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(FdataCx(jk)*diff)
        end do
        
        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*diff)
        end do
        
        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat9Sngl

  end subroutine afcstab_buildBoundsLPT1D

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_buildBoundsLPT2D(rafcstab, rmatrixCx, rmatrixCy,&
      rmatrixM, dscale)

!<description>
    ! This subroutine computes the solution bounds for
    ! linearity-preserving flux correction in 2D
!</description>

!<input>
    ! Coefficient matrices
    type(t_matrixScalar), intent(in) :: rmatrixCx, rmatrixCy

    ! Diagonal scaling matrix
    type(t_matrixScalar), intent(in) :: rmatrixM

    ! Global scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>
  
    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
    real(DP), dimension(:), pointer :: p_DdataM, p_DdataCx, p_DdataCy
    real(SP), dimension(:), pointer :: p_FdataM, p_FdataCx, p_FdataCy
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal, p_Kld, p_Kcol

    ! Check spatial discretisation structure
    if (.not. rmatrixCx%bidenticalTrialAndTest) then
      call output_line('Test and trial functions must be identical!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
      call sys_halt()
    end if

    ! Set pointer to triangulation; note that for more general
    ! discretisations the coordinates of the degrees of freedom must
    ! be computed and cannot be adopted from the vertex coordinates
    p_rtriangulation => rmatrixCx%p_rspatialDiscrTrial%p_rtriangulation

    ! Check if this is a 2D triangulation
    if (p_rtriangulation%ndim .ne. NDIM2D) then
      call output_line('Only 3D triangulations are supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
      call sys_halt()
    end if
    
    ! Get vertex coordinates
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(rmatrixCx, rafcstab)
    end if

    ! Check if stabilisation structure provides array to store bounds;
    ! otherwise allocate new memory using the predefined data type
    if (rafcstab%h_DboundsAtEdge .eq. ST_NOHANDLE)&
        call afcstab_allocBoundsAtEdge(rafcstab, rafcstab%cdataType)

    ! Check if source matrix is compatible with stabilisation structure
    if ((rmatrixCx%NA-rmatrixCx%NEQ .ne. rafcstab%NEDGE * 2) .or.&
        (rmatrixCx%cdataType .ne. rafcstab%cdataType)) then
      call output_line('Matrix/stabilisation structure is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
      call sys_halt()
    end if
    
    ! Check if coefficient matrices are compatible
    call lsyssc_isMatrixCompatible(rmatrixCx, rmatrixCy)

    ! Check if scaling matrix rmatrixM is diagonal and 
    ! has the same data type as the source matrix
    if ((rmatrixM%cmatrixFormat .ne. LSYSSC_MATRIXD)    .or.&
        (rmatrixM%cdataType .ne. rmatrixCx%cdataType) .or.&
        (rmatrixM%NEQ .ne. rmatrixCx%NEQ)) then
      call output_line('Coefficient matrices are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
      call sys_halt()
    end if

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)   

    ! What type of matrix are we?
    select case(rmatrixCx%cmatrixFormat)
    case (LSYSSC_MATRIX7)
      call lsyssc_getbase_Kld(rmatrixCx, p_Kld)
      
      ! What type of data format are we?
      select case(rmatrixCx%cdatatype)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixM, p_DdataM)
        call lsyssc_getbase_double (rmatrixCx, p_DdataCx)
        call lsyssc_getbase_double (rmatrixCy, p_DdataCy)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
        
        call doBoundsMat7Dble(p_IverticesAtEdge, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataM, p_DvertexCoords,&
            dscale, p_DboundsAtEdge)
        
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat7Sngl(p_IverticesAtEdge, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataM, p_DvertexCoords,&
            dscale, p_FboundsAtEdge)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
        call sys_halt()
      end select

    case (LSYSSC_MATRIX9)
      call lsyssc_getbase_Kdiagonal(rmatrixCx, p_Kdiagonal)
      call lsyssc_getbase_Kld(rmatrixCx, p_Kld)
      
      ! What type of data format are we?
      select case(rmatrixCx%cdatatype)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixM, p_DdataM)
        call lsyssc_getbase_double (rmatrixCx, p_DdataCx)
        call lsyssc_getbase_double (rmatrixCy, p_DdataCy)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
        
        call doBoundsMat9Dble(p_IverticesAtEdge, p_Kdiagonal, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataM, p_DvertexCoords,&
            dscale, p_DboundsAtEdge)
        
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat9Sngl(p_IverticesAtEdge, p_Kdiagonal, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataM, p_DvertexCoords,&
            dscale, p_FboundsAtEdge)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
      call sys_halt()
    end select

  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7Dble(IverticesAtEdge, Kld,&
        DdataCx, DdataCy, DdataM, DvertexCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM2D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:2,i) - DvertexCoords(1:2,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*Diff(1)+&
                            DdataCy(ik)*Diff(2))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*Diff(1)+&
                            DdataCy(jk)*Diff(2))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat7Dble

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9Dble(IverticesAtEdge, Kdiagonal, Kld,&
        DdataCx, DdataCy, DdataM, DvertexCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM2D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:2,i) - DvertexCoords(1:2,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(DdataCx(ik)*Diff(1)+&
                            DdataCy(ik)*Diff(2))
        end do
        
        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*Diff(1)+&
                            DdataCy(ik)*Diff(2))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(DdataCx(jk)*Diff(1)+&
                            DdataCy(jk)*Diff(2))
        end do
        
        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*Diff(1)+&
                            DdataCy(jk)*Diff(2))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat9Dble

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7Sngl(IverticesAtEdge, Kld,&
        FdataCx, FdataCy, FdataM, DvertexCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM2D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:2,i) - DvertexCoords(1:2,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*Diff(1)+&
                            FdataCy(ik)*Diff(2))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*Diff(1)+&
                            FdataCy(jk)*Diff(2))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat7Sngl

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9Sngl(IverticesAtEdge, Kdiagonal, Kld,&
        FdataCx, FdataCy, FdataM, DvertexCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM2D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:2,i) - DvertexCoords(1:2,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(FdataCx(ik)*Diff(1)+&
                            FdataCy(ik)*Diff(2))
        end do
        
        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*Diff(1)+&
                            FdataCy(ik)*Diff(2))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(FdataCx(jk)*Diff(1)+&
                            FdataCy(jk)*Diff(2))
        end do
        
        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*Diff(1)+&
                            FdataCy(jk)*Diff(2))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat9Sngl

  end subroutine afcstab_buildBoundsLPT2D

  !*****************************************************************************
  
!<subroutine>

  subroutine afcstab_buildBoundsLPT3D(rafcstab, rmatrixCx, rmatrixCy,&
      rmatrixCz, rmatrixM, dscale)

!<description>
    ! This subroutine computes the solution bounds for
    ! linearity-preserving flux correction in 3D
!</description>

!<input>
    ! Coefficient matrices
    type(t_matrixScalar), intent(in) :: rmatrixCx, rmatrixCy, rmatrixCz

    ! Diagonal scaling matrix
    type(t_matrixScalar), intent(in) :: rmatrixM

    ! Global scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>
  
    ! pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
    real(DP), dimension(:), pointer :: p_DdataM, p_DdataCx, p_DdataCy, p_DdataCz
    real(SP), dimension(:), pointer :: p_FdataM, p_FdataCx, p_FdataCy, p_FdataCz
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:), pointer :: p_Kdiagonal, p_Kld, p_Kcol

    ! Check spatial discretisation structure
    if (.not. rmatrixCx%bidenticalTrialAndTest) then
      call output_line('Test and trial functions must be identical!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
      call sys_halt()
    end if

    ! Set pointer to triangulation; note that for more general
    ! discretisations the coordinates of the degrees of freedom must
    ! be computed and cannot be adopted from the vertex coordinates
    p_rtriangulation => rmatrixCx%p_rspatialDiscrTrial%p_rtriangulation

    ! Check if this is a 3D triangulation
    if (p_rtriangulation%ndim .ne. NDIM3D) then
      call output_line('Only 3D triangulations are supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
      call sys_halt()
    end if
    
    ! Get vertex coordinates
    call storage_getbase_double2D(&
        p_rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(rmatrixCx, rafcstab)
    end if

    ! Check if stabilisation structure provides array to store bounds;
    ! otherwise allocate new memory using the predefined data type
    if (rafcstab%h_DboundsAtEdge .eq. ST_NOHANDLE)&
        call afcstab_allocBoundsAtEdge(rafcstab, rafcstab%cdataType)

    ! Check if source matrix is compatible with stabilisation structure
    if ((rmatrixCx%NA-rmatrixCx%NEQ .ne. rafcstab%NEDGE * 2) .or.&
        (rmatrixCx%cdataType .ne. rafcstab%cdataType)) then
      call output_line('Matrix/stabilisation structure is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
      call sys_halt()
    end if
    
    ! Check if coefficient matrices are compatible
    call lsyssc_isMatrixCompatible(rmatrixCx, rmatrixCy)
    call lsyssc_isMatrixCompatible(rmatrixCx, rmatrixCz)

    ! Check if scaling matrix rmatrixM is diagonal and 
    ! has the same data type as the source matrix
    if ((rmatrixM%cmatrixFormat .ne. LSYSSC_MATRIXD)    .or.&
        (rmatrixM%cdataType .ne. rmatrixCx%cdataType) .or.&
        (rmatrixM%NEQ .ne. rmatrixCx%NEQ)) then
      call output_line('Coefficient matrices are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
      call sys_halt()
    end if

    ! Set pointers
    call afcstab_getbase_IverticesAtEdge(rafcstab, p_IverticesAtEdge)   

    ! What type of matrix are we?
    select case(rmatrixCx%cmatrixFormat)
    case (LSYSSC_MATRIX7)
      call lsyssc_getbase_Kld(rmatrixCx, p_Kld)
      
      ! What type of data format are we?
      select case(rmatrixCx%cdatatype)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixM, p_DdataM)
        call lsyssc_getbase_double (rmatrixCx, p_DdataCx)
        call lsyssc_getbase_double (rmatrixCy, p_DdataCy)
        call lsyssc_getbase_double (rmatrixCz, p_DdataCz)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
        
        call doBoundsMat7Dble(p_IverticesAtEdge, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataCz, p_DdataM,&
            p_DvertexCoords, dscale, p_DboundsAtEdge)
        
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call lsyssc_getbase_single (rmatrixCz, p_FdataCz)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat7Sngl(p_IverticesAtEdge, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataCz, p_FdataM,&
            p_DvertexCoords, dscale, p_FboundsAtEdge)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
        call sys_halt()
      end select

    case (LSYSSC_MATRIX9)
      call lsyssc_getbase_Kdiagonal(rmatrixCx, p_Kdiagonal)
      call lsyssc_getbase_Kld(rmatrixCx, p_Kld)
      
      ! What type of data format are we?
      select case(rmatrixCx%cdatatype)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixM, p_DdataM)
        call lsyssc_getbase_double (rmatrixCx, p_DdataCx)
        call lsyssc_getbase_double (rmatrixCy, p_DdataCy)
        call lsyssc_getbase_double (rmatrixCz, p_DdataCz)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
        
        call doBoundsMat9Dble(p_IverticesAtEdge, p_Kdiagonal, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataCz, p_DdataM,&
            p_DvertexCoords, dscale, p_DboundsAtEdge)
        
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call lsyssc_getbase_single (rmatrixCz, p_FdataCz)
        call afcstab_getbase_DboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat9Sngl(p_IverticesAtEdge, p_Kdiagonal, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataCz, p_FdataM,&
            p_DvertexCoords, dscale, p_FboundsAtEdge)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
      call sys_halt()
    end select

  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7Dble(IverticesAtEdge, Kld,&
        DdataCx, DdataCy, DdataCz, DdataM, DvertexCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataCz, DdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM3D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:3,i) - DvertexCoords(1:3,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*Diff(1)+&
                            DdataCy(ik)*Diff(2)+&
                            DdataCz(ik)*Diff(3))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*Diff(1)+&
                            DdataCy(jk)*Diff(2)+&
                            DdataCz(jk)*Diff(3))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat7Dble

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9Dble(IverticesAtEdge, Kdiagonal, Kld,&
        DdataCx, DdataCy, DdataCz, DdataM, DvertexCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataCz, DdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM3D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:3,i) - DvertexCoords(1:3,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(DdataCx(ik)*Diff(1)+&
                            DdataCy(ik)*Diff(2)+&
                            DdataCz(ik)*Diff(3))
        end do
        
        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*Diff(1)+&
                            DdataCy(ik)*Diff(2)+&
                            DdataCz(ik)*Diff(3))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(DdataCx(jk)*Diff(1)+&
                            DdataCy(jk)*Diff(2)+&
                            DdataCz(jk)*Diff(3))
        end do
        
        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*Diff(1)+&
                            DdataCy(jk)*Diff(2)+&
                            DdataCz(jk)*Diff(3))
        end do
        
        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat9Dble

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7Sngl(IverticesAtEdge, Kld,&
        FdataCx, FdataCy, FdataCz, FdataM, DvertexCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataCz, FdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM3D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:3,i) - DvertexCoords(1:3,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*Diff(1)+&
                            FdataCy(ik)*Diff(2)+&
                            FdataCz(ik)*Diff(3))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*Diff(1)+&
                            FdataCy(jk)*Diff(2)+&
                            FdataCz(jk)*Diff(3))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat7Sngl

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9Sngl(IverticesAtEdge, Kdiagonal, Kld,&
        FdataCx, FdataCy, FdataCz, FdataM, DvertexCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataCz, FdataM
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP), dimension(NDIM3D) :: Diff
      real(DP) :: daux
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IverticesAtEdge,2)

        ! Get the numbers of the two endpoints
        i = IverticesAtEdge(1,iedge)
        j = IverticesAtEdge(2,iedge)

        ! Compute difference (x_i-x_j)
        Diff = DvertexCoords(1:3,i) - DvertexCoords(1:3,j)
        
        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(FdataCx(ik)*Diff(1)+&
                            FdataCy(ik)*Diff(2)+&
                            FdataCz(ik)*Diff(3))
        end do
        
        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*Diff(1)+&
                            FdataCy(ik)*Diff(2)+&
                            FdataCz(ik)*Diff(3))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(FdataCx(jk)*Diff(1)+&
                            FdataCy(jk)*Diff(2)+&
                            FdataCz(jk)*Diff(3))
        end do
        
        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*Diff(1)+&
                            FdataCy(jk)*Diff(2)+&
                            FdataCz(jk)*Diff(3))
        end do
        
        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat9Sngl

  end subroutine afcstab_buildBoundsLPT3D

end module afcstabilisation

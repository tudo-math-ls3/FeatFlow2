!##############################################################################
!# ****************************************************************************
!# <name> afcstabbase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for performing
!# discrete stabilisation by means of algebraic flux correction
!#
!# The following routines are available:
!#
!# 1.) afcstab_initFromParameterlist
!#     -> Creates a stabilisation structure and initialise
!#        it from the values of a given parameter list
!#
!# 2.) afcstab_releaseStabilisation
!#     -> Releases a stabilisation structure
!#
!# 3.) afcstab_resizeStabilisation = afcstab_resizeStabDirect /
!#                                   afcstab_resizeStabIndScalar
!#                                   afcstab_resizeStabIndBlock /
!#                                   afcstab_resizeStabIndGFEM
!#     -> Resizes a stabilisation structure
!#
!# 4.) afcstab_copyStabilisation
!#     -> Copies a stabilisation structure to another structure
!#
!# 5.) afcstab_duplicateStabilisation
!#     -> Duplicate (parts of) a stabilisation structure
!#
!# 6.) afcstab_isMatrixCompatible = afcstab_isMatrixCompatibleSc /
!#                                  afcstab_isMatrixCompatibleBl
!#     -> Checks whether a matrix and a stabilisation structure are compatible
!#
!# 7.) afcstab_isVectorCompatible = afcstab_isVectorCompatibleSc /
!#                                  afcstab_isVectorCompatibleBl
!#     -> Checks whether a vector and a stabilisation structure are compatible
!#
!# 8.) afcstab_getbase_IedgeListIdx
!#     -> Returns pointer to the index pointer for the edge structure
!#
!# 9.) afcstab_getbase_IedgeList
!#     -> Returns pointer to the edge structure
!#
!# 10.) afcstab_getbase_IsupdiagEdgeIdx
!#      -> Returns pointer to the index pointer for the superdiagonal edge numbers
!#
!# 11.) afcstab_getbase_IsubdiagEdgeIdx
!#      -> Returns pointer to the index pointer for the
!#         subdiagonal edge numbers
!#
!# 12.) afcstab_getbase_IsubdiagEdge
!#      -> Returns pointer to the subdiagonal edge numbers
!#
!# 13.) afcstab_getbase_QcoeffsAtEdge
!#      afcstab_getbase_DcoeffsAtEdge / afcstab_getbase_FcoeffsAtEdge
!#      -> Returns pointer to edge data
!#
!# 14.) afcstab_getbase_QboundsAtEdge
!#      afcstab_getbase_DboundsAtEdge / afcstab_getbase_FboundsAtEdge
!#      -> Returns pointer to the bounds at edges
!#
!# 15.) afcstab_genEdgeList
!#      -> Generates the standard edge data structure
!#
!# 16.) afcstab_genOffdiagEdges
!#      -> Generates the subdiagonal edge data structure
!#
!# 17.) afcstab_genExtSparsity
!#      -> Generates the extended sparsity pattern
!#
!# 18.) afcstab_allocInternalData
!#      -> Allocates internal data of the stabilisation structure
!#
!# 19.) afcstab_allocEdgeStructure
!#      -> Allocates the edge data structure
!#
!# 20.) afcstab_allocCoeffsAtEdge
!#      -> Allocates the coefficients at edge data structure
!#
!# 21.) afcstab_allocVectorsPQR
!#      -> Allocates the nodal vectors P, Q, and R each for '+' and '-'
!#
!# 22.) afcstab_allocFlux
!#      -> Allocates the edge-wise flux vector flux
!#
!# 23.) afcstab_allocFlux0
!#      -> Allocates the edge-wise flux vector flux0
!#
!# 24.) afcstab_allocFluxPrel
!#      -> Allocates the edge-wise flux vector fluxPrel
!#
!# 25.) afcstab_allocAlpha
!#      -> Allocates the edge-wise correction factors alpha
!#
!# 26.) afcstab_allocBoundsAtEdge
!#      -> Allocates the bounds at edge data structure
!#
!# 27.) afcstab_buildBoundsLPT = afcstab_buildBoundsLPT1D /
!#                               afcstab_buildBoundsLPT2D /
!#                               afcstab_buildBoundsLPT3D
!#
!# 28.) afcstab_limit = afcstab_limitUnboundedQP /
!#                      afcstab_limitUnboundedDP /
!#                      afcstab_limitUnboundedSP /
!#                      afcstab_limitBoundedQP /
!#                      afcstab_limitBoundedDP /
!#                      afcstab_limitBoundedSP
!#      -> Compute the nodal correction factors, i.e., the ratio of
!#         admissible solution increments and raw antidiffusion
!#
!# 29.) afcstab_combineFluxes = afcstab_combFluxesQP /
!#                              afcstab_combFluxesDP /
!#                              afcstab_combFluxesSP
!#      -> Linear combination of the vectors of fluxes
!#
!# 30.) afcstab_upwindOrientation = afcstab_upwindOrientationDP1 /
!#                                  afcstab_upwindOrientationDP2 /
!#                                  afcstab_upwindOrientationSP1 /
!#                                  afcstab_upwindOrientationSP2
!#      -> Swap edge orientation so that the starting edge is located upwind
!#
!# 31.) afcstab_infoStabilisation
!#      -> Outputs information about the stabilisation structure
!#
!# 32.) afcstab_copyH2D_IedgeListIdx
!#      -> Copies the edge index structure from the host memory
!#         to the device memory.
!#
!# 33.) afcstab_copyD2H_IedgeListIdx
!#      -> Copies the edge index structure from the device memory
!#         to the host memory.
!#
!# 34.) afcstab_copyH2D_IedgeList
!#      -> Copies the edge structure from the host memory
!#         to the device memory.
!#
!# 35.) afcstab_copyD2H_IedgeList
!#      -> Copies the edge structure from the device memory
!#         to the host memory.
!#
!# 36.) afcstab_copyH2D_CoeffsAtEdge
!#      -> Copies the edge data from the host memory
!#         to the device memory.
!#
!# 37.) afcstab_copyD2H_CoeffsAtEdge
!#      -> Copies the edge structure from the device memory
!#         to the host memory.
!#
!# 38.) afcstab_setDataType
!#      -> Converts the data type of the stabilisation structure
!#
!# </purpose>
!##############################################################################
module afcstabbase

!$ use omp_lib
  use basicgeometry
  use fsystem
  use genoutput
  use groupfembase
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use perfconfig
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none

  private
  public :: t_afcstab

  public :: afcstab_initPerfConfig
  public :: afcstab_initFromParameterlist
  public :: afcstab_releaseStabilisation
  public :: afcstab_resizeStabilisation
  public :: afcstab_copyStabilisation
  public :: afcstab_duplicateStabilisation
  public :: afcstab_isMatrixCompatible
  public :: afcstab_isVectorCompatible
  public :: afcstab_infoStabilisation

  public :: afcstab_getbase_IedgeListIdx
  public :: afcstab_getbase_IedgeList
  public :: afcstab_getbase_IsupdiagEdgeIdx
  public :: afcstab_getbase_IsubdiagEdgeIdx
  public :: afcstab_getbase_IsubdiagEdge
  public :: afcstab_getbase_QcoeffsAtEdge
  public :: afcstab_getbase_DcoeffsAtEdge
  public :: afcstab_getbase_FcoeffsAtEdge
  public :: afcstab_getbase_QboundsAtEdge
  public :: afcstab_getbase_DboundsAtEdge
  public :: afcstab_getbase_FboundsAtEdge

  public :: afcstab_genEdgeList
  public :: afcstab_genOffdiagEdges
  public :: afcstab_genExtSparsity
  public :: afcstab_buildBoundsLPT

  public :: afcstab_allocInternalData
  public :: afcstab_allocEdgeStructure
  public :: afcstab_allocCoeffsAtEdge
  public :: afcstab_allocBoundsAtEdge
  public :: afcstab_allocVectorsPQR
  public :: afcstab_allocFluxPrel
  public :: afcstab_allocFlux0
  public :: afcstab_allocFlux
  public :: afcstab_allocAlpha

  public :: afcstab_limit
  public :: afcstab_combineFluxes
  public :: afcstab_combineFluxesQP
  public :: afcstab_combineFluxesDP
  public :: afcstab_combineFluxesSP
  public :: afcstab_upwindOrientation

  public :: afcstab_copyD2H_IedgeListIdx
  public :: afcstab_copyH2D_IedgeListIdx
  public :: afcstab_copyD2H_IedgeList
  public :: afcstab_copyH2D_IedgeList
  public :: afcstab_copyD2H_CoeffsAtEdge
  public :: afcstab_copyH2D_CoeffsAtEdge

  public :: afcstab_setDataType

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

  ! Upwind-biased limiting (implies edge-orientation)
  integer, parameter, public :: AFCSTAB_LIMITING_UPWINDBIASED     = 1

  ! Symmetric limiting (implies no edge-orientation)
  integer, parameter, public :: AFCSTAB_LIMITING_SYMMETRIC        = 2

  ! Characteristic limiting
  integer, parameter, public :: AFCSTAB_LIMITING_CHARACTERISTIC   = 3

!</constantblock>


!<constantblock description="Bitfield identifiers for properties of stabilisation">

  ! Stabilisation is undefined
  integer(I32), parameter, public :: AFCSTAB_NONE                 = 0_I32

  ! Stabilisation has been initialised
  integer(I32), parameter, public :: AFCSTAB_INITIALISED          = 2_I32**0

  ! Edge-based structure has been generated: IedgeList
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGELIST         = 2_I32**1

  ! Edge-based structure has been oriented: IedgeList
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEORIENTATION  = 2_I32**2

  ! Edge-based values have been computed from matrix: DcoefficientsAtEdge
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEVALUES       = 2_I32**3

  ! Subdiagonal edge-based structure has been generated
  integer(I32), parameter, public :: AFCSTAB_HAS_OFFDIAGONALEDGES = 2_I32**4

  ! Antidiffusive fluxes have been precomputed
  integer(I32), parameter, public :: AFCSTAB_HAS_ADFLUXES         = 2_I32**5

  ! Nodal sums of antidiffusive increments have been computed: PP, PM
  integer(I32), parameter, public :: AFCSTAB_HAS_ADINCREMENTS     = 2_I32**6

  ! Nodal upper/lower bounds have been computed: QP, QM
  integer(I32), parameter, public :: AFCSTAB_HAS_NODEBOUNDS       = 2_I32**7

  ! Nodal correction factors have been computed: RP, RM
  integer(I32), parameter, public :: AFCSTAB_HAS_NODELIMITER      = 2_I32**8

  ! Edge-wise correction factors have been computed: ALPHA
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGELIMITER      = 2_I32**9

  ! Low-order predictor has been computed
  integer(I32), parameter, public :: AFCSTAB_HAS_PREDICTOR        = 2_I32**10

  ! Transformed nodal solution values have been computed
  integer(I32), parameter, public :: AFCSTAB_HAS_NODEVALUES       = 2_I32**11

  ! Slope-based bounds have been computed
  integer(I32), parameter, public :: AFCSTAB_HAS_EDGEBOUNDS       = 2_I32**12

!</constantblock>


!<constantblock description="Duplication flags. Specifies which information is duplicated">

  ! Duplicate atomic stabilisation structure
  integer(I32), parameter, public :: AFCSTAB_DUP_STRUCTURE        = 2_I32**1

  ! Duplicate edge-based structure: IedgeList
  integer(I32), parameter, public :: AFCSTAB_DUP_EDGELIST         = AFCSTAB_HAS_EDGELIST

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

  ! Duplicate transformed nodal solution values
  integer(I32), parameter, public :: AFCSTAB_DUP_NODEVALUES       = AFCSTAB_HAS_NODEVALUES

  ! Duplicate slope-based bounds
  integer(I32), parameter, public :: AFCSTAB_DUP_EDGEBOUNDS       = AFCSTAB_HAS_EDGEBOUNDS

!</constantblock>


!<constantblock description="Duplication flags. Specifies which information is shared \
!                            between stabilisation structures">

  ! Duplicate atomic stabilisation structure
  integer(I32), parameter, public :: AFCSTAB_SHARE_STRUCTURE        = AFCSTAB_DUP_STRUCTURE

  ! Share edge-based structure: IedgeList
  integer(I32), parameter, public :: AFCSTAB_SHARE_EDGELIST         = AFCSTAB_DUP_EDGELIST

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

  ! Share transformed nodal solution values
  integer(I32), parameter, public :: AFCSTAB_SHARE_NODEVALUES       = AFCSTAB_DUP_NODEVALUES

  ! Share slope-based bounds
  integer(I32), parameter, public :: AFCSTAB_SHARE_EDGEBOUNDS       = AFCSTAB_DUP_EDGEBOUNDS

!</constantblock>


!<constantblock description="Bitfield identifiers for TVD-algorithm">

  ! Initialise the edgewise correction factors by unity
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_INITALPHA    = 2_I32**0

  ! Compute the raw-antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_ADFLUXES     = 2_I32**1

  ! Compute the sums of antidiffusive increments
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_ADINCREMENTS = 2_I32**2

  ! Compute the local solution bounds
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_BOUNDS       = 2_I32**3

  ! Compute the nodal correction factors
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_LIMITNODAL   = 2_I32**4

  ! Compute edgewise correction factors
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_LIMITEDGE    = 2_I32**5

  ! Correct raw antidiffusive fluxes and apply them
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_CORRECT      = 2_I32**6

  ! FEM-TVD algorithm without application of the corrected fluxes
  integer(I32), parameter, public :: AFCSTAB_TVDALGO_PREPARE   = AFCSTAB_TVDALGO_INITALPHA +&
                                                                 AFCSTAB_TVDALGO_ADFLUXES +&
                                                                 AFCSTAB_TVDALGO_ADINCREMENTS +&
                                                                 AFCSTAB_TVDALGO_BOUNDS +&
                                                                 AFCSTAB_TVDALGO_LIMITNODAL +&
                                                                 AFCSTAB_TVDALGO_LIMITEDGE

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

  ! Initialise the edgewise correction factors by unity
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

  ! Initialise the edgewise correction factors by unity
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


!<constantblock description="Bitfield identifiers for failsafe algorithm">

  ! Initialise the edgewise correction factors by unity
  integer(I32), parameter, public :: AFCSTAB_FAILSAFEALGO_INITBETA = 2_I32**0

  ! Compute the distances to a local extremum
  integer(I32), parameter, public :: AFCSTAB_FAILSAFEALGO_BOUNDS   = 2_I32**1

  ! Perform failsafe correction of antidiffusive fluxes
  integer(I32), parameter, public :: AFCSTAB_FAILSAFEALGO_LIMIT    = 2_I32**2

  ! Apply the failsafe correction
  integer(I32), parameter, public :: AFCSTAB_FAILSAFEALGO_CORRECT  = 2_I32**3

  ! Standard failsafe FCT algorithm
  integer(I32), parameter, public :: AFCSTAB_FAILSAFEALGO_STANDARD = AFCSTAB_FAILSAFEALGO_INITBETA +&
                                                                     AFCSTAB_FAILSAFEALGO_BOUNDS +&
                                                                     AFCSTAB_FAILSAFEALGO_LIMIT +&
                                                                     AFCSTAB_FAILSAFEALGO_CORRECT

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
  real(DP), parameter, public :: AFCSTAB_EPSABS = 1e-10
#endif

  ! Relative tolerance for stabilisation
#ifndef AFCSTAB_EPSREL
  real(DP), parameter, public :: AFCSTAB_EPSREL = 1e-3
#endif

!</constantblock>

!</constants>

  ! *****************************************************************************

!<types>
!<typeblock>

  ! This structures holds all required information for stabilisation
  ! of algebraic flux correction type. It is applicable to scalar
  ! problems and systems of equations alike. Depending on the
  ! stabilisation strategy some internal structures will be generated.

  type t_afcstab

    ! Flag which indicates if the underlying operator is symmetric
    logical :: bisSymmetricOperator = .false.

    ! Format Tag: Identifies the type of stabilisation
    integer :: cafcstabType = AFCSTAB_GALERKIN

    ! Format Tag: Identifies the type of prelimiting
    ! Can take one of the AFCSTAB_PRELIMITING_xxxx flags.
    integer :: cprelimitingType = AFCSTAB_PRELIMITING_NONE

    ! Format Tag: Identifies the type of limiting
    integer :: climitingType = AFCSTAB_LIMITING_NONE

    ! Format Tag: Identifies the data type
    integer :: cdataType = ST_DOUBLE

    ! Duplication Flag: This is a bitfield coming from an OR
    ! combination of different AFCSTAB_SHARE_xxxx constants and
    ! specifies which parts of the stabilisation are shared with
    ! another stabilisation structure.
    integer(I32) :: iduplicationFlag = AFCSTAB_NONE

    ! Specification Flag: Specifies the stabilisation. This is a bitfield
    ! coming from an OR combination of different AFCSTAB_HAS_xxxx
    ! constants and specifies various properties of the stabilisation
    ! structure.
    integer(I32) :: istabilisationSpec = AFCSTAB_NONE

    ! Number of equations of the sparsity pattern
    integer :: NEQ = 0

    ! Number of edges of the sparsity pattern
    integer :: NEDGE = 0

    ! Maximum number of edges adjacent to one vertex. This
    ! corresponds to the maximum number of nonzero row entries.
    integer :: NNVEDGE = 0

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

    ! Number of different coefficients stored in CoefficientsAtEdge
    integer :: ncoeffsAtEdge = 0

    ! Handle to index pointer for edge structure
    ! The numbers IedgeListIdx(k):IedgeListIdx(k+1)-1
    ! denote the edge numbers of the k-th group of edges.
    integer :: h_IedgeListIdx = ST_NOHANDLE

    ! Handle to edge structure
    ! IedgeList(1:2,1:NEDGE) : the two end-points i and j of the edge (ij)
    ! IedgeList(3:4,1:NEDGE) : the two matrix position ij and ji that
    !                          correspond to the edge (ij)
    ! IedgeList(5:6,1:NEDGE) : the two matrix position ii and jj that
    !                          correspond to the diagonal entries
    integer :: h_IedgeList = ST_NOHANDLE

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
    integer :: h_CoeffsAtEdge = ST_NOHANDLE

    ! Handle to bounds at edges (i.e. off-diagonal entries)
    integer :: h_BoundsAtEdge = ST_NOHANDLE

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
!</types>

  !*****************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: afcstab_perfconfig

  ! *****************************************************************************

  interface afcstab_resizeStabilisation
    module procedure afcstab_resizeStabDirect
    module procedure afcstab_resizeStabIndScalar
    module procedure afcstab_resizeStabIndBlock
    module procedure afcstab_resizeStabIndGFEM
  end interface

  interface afcstab_isMatrixCompatible
    module procedure afcstab_isMatrixCompatibleSc
    module procedure afcstab_isMatrixCompatibleBl
  end interface

  interface afcstab_isVectorCompatible
    module procedure afcstab_isVectorCompatibleSc
    module procedure afcstab_isVectorCompatibleBl
  end interface

  interface afcstab_buildBoundsLPT
    module procedure afcstab_buildBoundsLPT1D
    module procedure afcstab_buildBoundsLPT2D
    module procedure afcstab_buildBoundsLPT3D
  end interface

  interface afcstab_limit
#ifdef ENABLE_QUADPREC
    module procedure afcstab_limitUnboundedQP
    module procedure afcstab_limitBoundedQP
#endif
    module procedure afcstab_limitUnboundedDP
    module procedure afcstab_limitBoundedDP
    module procedure afcstab_limitBoundedSP
    module procedure afcstab_limitUnboundedSP
  end interface

  interface afcstab_combineFluxes
#ifdef ENABLE_QUADPREC
    module procedure afcstab_combFluxesQP
#endif
    module procedure afcstab_combFluxesDP
    module procedure afcstab_combFluxesSP
  end interface

  interface afcstab_upwindOrientation
    module procedure afcstab_upwindOrientationDP1
    module procedure afcstab_upwindOrientationDP2
    module procedure afcstab_upwindOrientationSP1
    module procedure afcstab_upwindOrientationSP2
  end interface

contains

  !****************************************************************************

!<subroutine>

  subroutine afcstab_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      afcstab_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(afcstab_perfconfig)
    end if

  end subroutine afcstab_initPerfConfig

  ! ***************************************************************************

!<subroutine>

  subroutine afcstab_initFromParameterlist(rparlist, ssectionName, rafcstab)

!<description>
    ! This subroutine creates a stabilisation structure and initialises
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


    ! First, we retrieve all information of the stabilisation
    ! structure specified in the parameter file.

    ! Get type of stabilisation from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "istabilisation", rafcstab%cafcstabType, AFCSTAB_GALERKIN)

    ! Get type of prelimiting from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
        "iprelimiting", rafcstab%cprelimitingType, AFCSTAB_PRELIMITING_NONE)


    ! In a second step, special settings particular to individual
    ! stabilisation techniques are made `by hand'.

    ! Check if stabilisation should be applied
    select case(rafcstab%cafcstabType)

    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)
      ! high-order Galerkin scheme, low-order scheme:
      ! no prelimiting, no limiting
      rafcstab%cprelimitingType = AFCSTAB_PRELIMITING_NONE
      rafcstab%climitingType    = AFCSTAB_LIMITING_NONE

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_LINFCT)
      ! Symmetric flux limiting with prelimiting as given by the user
      rafcstab%climitingType = AFCSTAB_LIMITING_SYMMETRIC

    case (AFCSTAB_LINFCT_MASS,&
          AFCSTAB_SYMMETRIC,&
          AFCSTAB_NLINLPT_MASS,&
          AFCSTAB_NLINLPT_SYMMETRIC,&
          AFCSTAB_LINLPT_MASS,&
          AFCSTAB_LINLPT_SYMMETRIC)
      ! Symmetric flux limiting without prelimiting
      rafcstab%cprelimitingType = AFCSTAB_PRELIMITING_NONE
      rafcstab%climitingType    = AFCSTAB_LIMITING_SYMMETRIC

    case (AFCSTAB_LINLPT_UPWINDBIASED)
      ! For this special case, symmetric flux limiting must be
      ! performed as stated in the conclusion of the paper:
      ! D. Kuzmin, Linearity-preserving flux correction and
      ! convergence acceleration for constrained Galerin schemes,
      ! Technical Report 421, TU Dortmund, 2011.
      rafcstab%cprelimitingType = AFCSTAB_PRELIMITING_NONE
      rafcstab%climitingType    = AFCSTAB_LIMITING_SYMMETRIC

    case (AFCSTAB_TVD)
      ! Upwind-biased flux limiting with standard prelimiting
      rafcstab%cprelimitingType = AFCSTAB_PRELIMITING_STD
      rafcstab%climitingType    = AFCSTAB_LIMITING_UPWINDBIASED

    case (AFCSTAB_GP,&
          AFCSTAB_NLINLPT_UPWINDBIASED)
      ! Upwind-biased flux limiting with minmod prelimiting
      rafcstab%cprelimitingType = AFCSTAB_PRELIMITING_MINMOD
      rafcstab%climitingType    = AFCSTAB_LIMITING_UPWINDBIASED

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
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGELIST)) then
      if (rafcstab%h_IedgeList .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IedgeList)
      if (rafcstab%h_IedgeListIdx .ne. ST_NOHANDLE)&
          call storage_free(rafcstab%h_IedgeListIdx)
    end if
    rafcstab%h_IedgeList = ST_NOHANDLE
    rafcstab%h_IedgeListIdx = ST_NOHANDLE

    ! Release edge values
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEVALUES) .and.&
        (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE))&
        call storage_free(rafcstab%h_CoeffsAtEdge)
    rafcstab%h_CoeffsAtEdge = ST_NOHANDLE

    ! Release edge bounds
    if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEBOUNDS) .and.&
        (rafcstab%h_BoundsAtEdge .ne. ST_NOHANDLE))&
        call storage_free(rafcstab%h_BoundsAtEdge)
    rafcstab%h_BoundsAtEdge = ST_NOHANDLE

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
    rafcstab%bisSymmetricOperator = .false.
    rafcstab%cafcstabType         = AFCSTAB_GALERKIN
    rafcstab%cprelimitingType     = AFCSTAB_PRELIMITING_NONE
    rafcstab%climitingType        = AFCSTAB_LIMITING_NONE
    rafcstab%cdataType            = ST_DOUBLE
    rafcstab%istabilisationSpec   = AFCSTAB_NONE
    rafcstab%iduplicationFlag     = 0
    rafcstab%NEQ                  = 0
    rafcstab%NVAR                 = 1
    rafcstab%NVARtransformed      = 1
    rafcstab%NEDGE                = 0
    rafcstab%NNVEDGE              = 0
    rafcstab%ncoeffsAtEdge        = 0

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

  subroutine afcstab_resizeStabDirect(rafcstab, NEQ, NEDGE)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure to the new values NEQ and NEDGE
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! number of equations
    integer, intent(in) :: NEQ

    ! number of edges
    integer, intent(in) :: NEDGE
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(2) :: Isize2D
    integer :: isize

    ! REMARK: The handle IedgeListIdx is not modified here due
    ! to the following reasons: If the edges are not reordered into
    ! independent groups then IedgeListIdx(1:2)=(/1,NEDGE+1/).
    ! Otherwise, the edge-data structure needs to be regenerated anyway.

    !---------------------------------------------------------------------------
    ! Resize nodal quantities
    !---------------------------------------------------------------------------
    if (rafcstab%NEQ .ne. NEQ) then

      ! Set new number of nodes
      rafcstab%NEQ = NEQ

      ! Resize edge index vector and clear specification
      !-------------------------------------------------------------------------
      if (rafcstab%h_IsuperdiagEdgesIdx .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES)) then
          call storage_getsize(rafcstab%h_IsuperdiagEdgesIdx, isize)
          if (rafcstab%NEQ+1 .ne. isize) then
            call output_line('Handle h_IsuperdiagEdgesIdx '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEQ+1, rafcstab%h_IsuperdiagEdgesIdx,&
              ST_NEWBLOCK_NOINIT, .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_OFFDIAGONALEDGES))
        end if
      end if

      ! Resize subdiagonal edge index vector and clear specification
      !-------------------------------------------------------------------------
      if (rafcstab%h_IsubdiagEdgesIdx .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES)) then
          call storage_getsize(rafcstab%h_IsubdiagEdgesIdx, isize)
          if (rafcstab%NEQ+1 .ne. isize) then
            call output_line('Handle h_IsubdiagEdgesIdx '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEQ+1, rafcstab%h_IsubdiagEdgesIdx,&
              ST_NEWBLOCK_NOINIT, .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_OFFDIAGONALEDGES))
        end if
      end if

      ! Resize transformed nodal vector and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorDx)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEVALUES)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorDx%NEQ) then
            call output_line('Vector p_rvectorDx '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorDx,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_NODEVALUES))
        end if
      end if

      ! Resize nodal vectors P+ and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorPp)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADINCREMENTS)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorPp%NEQ) then
            call output_line('Vector p_rvectorPp '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorPp,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_ADINCREMENTS))
        end if
      end if

      ! Resize nodal vectors P- and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorPm)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADINCREMENTS)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorPm%NEQ) then
            call output_line('Vector p_rvectorPm '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorPm,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_ADINCREMENTS))
        end if
      end if

      ! Resize nodal vectors Q and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorQ)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorQ%NEQ) then
            call output_line('Vector p_rvectorQ '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorQ,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_NODEBOUNDS))
        end if
      end if

      ! Resize nodal vectors Q+ and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorQp)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorQp%NEQ) then
            call output_line('Vector p_rvectorQp '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorQp,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_NODEBOUNDS))
        end if
      end if

      ! Resize nodal vectors Q- and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorQm)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODEBOUNDS)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorQm%NEQ) then
            call output_line('Vector p_rvectorQm '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorQm,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_NODEBOUNDS))
        end if
      end if

      ! Resize nodal vectors R+ and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorRp)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODELIMITER)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorRp%NEQ) then
            call output_line('Vector p_rvectorRp '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorRp,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_NODELIMITER))
        end if
      end if

      ! Resize nodal vectors R- and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorRm)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_NODELIMITER)) then
          if (rafcstab%NEQ .ne. rafcstab%p_rvectorRm%NEQ) then
            call output_line('Vector p_rvectorRm '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorRm,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_NODELIMITER))
        end if
      end if

      ! Resize nodal vector for the low-order predictor and clear specification
      !-------------------------------------------------------------------------
      if (associated(rafcstab%p_rvectorPredictor)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_PREDICTOR)) then
          if (rafcstab%NEQ*rafcstab%NVAR .ne. rafcstab%p_rvectorPredictor%NEQ) then
            call output_line('Vector p_rvectorPredictor '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsysbl_resizeVector(rafcstab%p_rvectorPredictor,&
              rafcstab%NEQ, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_PREDICTOR))
        end if
      end if

    end if

    !---------------------------------------------------------------------------
    ! Resize edge quantities
    !---------------------------------------------------------------------------
    if (rafcstab%NEDGE .ne. NEDGE) then

      ! Set new number of edges
      rafcstab%NEDGE = NEDGE

      ! Resize array of edges and clear specification
      !-------------------------------------------------------------------------
      if (rafcstab%h_IedgeList .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGELIST)) then
          call storage_getsize(rafcstab%h_IedgeList, Isize2D)
          if (rafcstab%NEDGE .ne. Isize2D(2)) then
            call output_line('Handle h_IedgeList '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_IedgeList,&
              ST_NEWBLOCK_NOINIT, .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_EDGELIST))
        end if
      end if

      ! Resize array of subdiagonal edges and clear specification
      !-------------------------------------------------------------------------
      if (rafcstab%h_IsubdiagEdges .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_OFFDIAGONALEDGES)) then
          call storage_getsize(rafcstab%h_IsubdiagEdges, isize)
          if (rafcstab%NEDGE .ne. isize) then
            call output_line('Handle h_IsubdiagEdges '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_IsubdiagEdges,&
              ST_NEWBLOCK_NOINIT, .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_OFFDIAGONALEDGES))
        end if
      end if

      ! Resize array of edge data and clear specification
      !-------------------------------------------------------------------------
      if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEVALUES)) then
          call storage_getsize(rafcstab%h_CoeffsAtEdge, Isize2D)
          if (rafcstab%NEDGE .ne. Isize2D(2)) then
            call output_line('Handle h_CoeffsAtEdge '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_CoeffsAtEdge,&
              ST_NEWBLOCK_NOINIT, .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_EDGEVALUES))
        end if
      end if

      ! Resize array of edge bounds and clear specification
      !-------------------------------------------------------------------------
      if (rafcstab%h_BoundsAtEdge .ne. ST_NOHANDLE) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGEBOUNDS)) then
          call storage_getsize(rafcstab%h_BoundsAtEdge, Isize2D)
          if (rafcstab%NEDGE .ne. Isize2D(2)) then
            call output_line('Handle h_BoundsAtEdge '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call storage_realloc('afcstab_resizeStabDirect',&
              rafcstab%NEDGE, rafcstab%h_BoundsAtEdge,&
              ST_NEWBLOCK_NOINIT, .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_EDGEBOUNDS))
        end if
      end if

      ! Resize edge vectors for the correction factor and clear specification
      !-------------------------------------------------------------------------
      if(associated(rafcstab%p_rvectorAlpha)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGELIMITER)) then
          if (rafcstab%NEDGE .ne. rafcstab%p_rvectorAlpha%NEQ) then
            call output_line('Vector p_rvectorAlpha '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorAlpha,&
              rafcstab%NEDGE, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_EDGELIMITER))
        end if
      end if

      ! Resize edge vectors for the explicit flux and clear specification
      !-------------------------------------------------------------------------
      if(associated(rafcstab%p_rvectorFlux0)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES)) then
          if (rafcstab%NEDGE .ne. rafcstab%p_rvectorFlux0%NEQ) then
            call output_line('Vector p_rvectorFlux0 '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorFlux0,&
              rafcstab%NEDGE, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_ADFLUXES))
        end if
      end if

      ! Resize edge vectors for the implicit flux and clear specification
      !-------------------------------------------------------------------------
      if(associated(rafcstab%p_rvectorFlux)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES)) then
          if (rafcstab%NEDGE .ne. rafcstab%p_rvectorFlux%NEQ) then
            call output_line('Vector p_rvectorFlux '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorFlux,&
              rafcstab%NEDGE, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_ADFLUXES))
        end if
      end if

      ! Resize edge vectors for the prelimiting flux and clear specification
      !-------------------------------------------------------------------------
      if(associated(rafcstab%p_rvectorFluxPrel)) then
        if (check(rafcstab%iduplicationFlag, AFCSTAB_SHARE_ADFLUXES)) then
          if (rafcstab%NEDGE .ne. rafcstab%p_rvectorFluxPrel%NEQ) then
            call output_line('Vector p_rvectorFluxPrel '//&
                'is shared and cannot be resized!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabDirect')
            call sys_halt()
          end if
        else
          call lsyssc_resizeVector(rafcstab%p_rvectorFluxPrel,&
              rafcstab%NEDGE, .false., .false.)

          ! Reset specifier
          rafcstab%istabilisationSpec = iand(rafcstab%istabilisationSpec,&
                                             not(AFCSTAB_HAS_ADFLUXES))
        end if
      end if

    end if

  contains

    !**************************************************************
    ! Checks if bitfield ibitfield in iflag is set.

    pure function check(iflag, ibitfield)

      integer(I32), intent(in) :: iflag,ibitfield

      logical :: check

      check = (iand(iflag,ibitfield) .eq. ibitfield)

    end function check

  end subroutine afcstab_resizeStabDirect

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabIndScalar(rafcstab, rmatrix)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure so that they are compatible to the template matrix.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: neq, nedge


    ! Determine number of equations and edges
    neq   = rmatrix%NEQ
    nedge = (rmatrix%NA-rmatrix%NEQ)/2

    ! Call resize routine directly
    call afcstab_resizeStabDirect(rafcstab, neq, nedge)

    ! Regenerate edge list
    call afcstab_genEdgeList(rmatrix, rafcstab)

  end subroutine afcstab_resizeStabIndScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabIndBlock(rafcstab, rmatrixBlock)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure so that they are compatible to the template matrix.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! template matrix
    type(t_matrixBlock), intent(in) :: rmatrixBlock
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: neq, nedge


    ! Check if block matrix has only one block
    if ((rmatrixBlock%nblocksPerCol .eq. 1) .and.&
        (rmatrixBlock%nblocksPerRow .eq. 1)) then
      call afcstab_resizeStabIndScalar(rafcstab,&
          rmatrixBlock%RmatrixBlock(1,1))

      ! That is it
      return
    end if


    ! Determine number of equations and edges
    neq   = rmatrixBlock%RmatrixBlock(1,1)%NEQ
    nedge = (rmatrixBlock%RmatrixBlock(1,1)%NA-&
             rmatrixBlock%RmatrixBlock(1,1)%NEQ)/2

    ! Call resize routine directly
    call afcstab_resizeStabDirect(rafcstab, neq, nedge)

    ! Regenerate edge list
    call afcstab_genEdgeList(rmatrixBlock%RmatrixBlock(1,1), rafcstab)

  end subroutine afcstab_resizeStabIndBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_resizeStabIndGFEM(rafcstab, rgroupFEMSet)

!<description>
    ! This subroutine resizes all vectors of the stabilisation
    ! structure so that they are compatible to the template matrix.
    !
    ! NOTE: Only those vectors are resized which are actually present.
!</description>

!<input>
    ! group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout)   :: rafcstab
!</inputoutput>
!</subroutine>

    ! Call resize routine directly
    call afcstab_resizeStabDirect(rafcstab, rgroupFEMSet%NEQ, rgroupFEMSet%NEDGE)

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    if (iand(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGELIST) .ne. 0) then
      ! The edge list is shared with the group finite element
      ! set. Therefore, update the storage handles which may have
      ! changed due a possible resize of rgroupFEMSet. If the edge
      ! list is not shared, then it has been resized by the above call
      ! to subroutine afcstab_resizeStabDirect
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx     = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList        = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag   = ior(rafcstab%iduplicationFlag,&
                                          AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_resizeStabIndGFEM')
        call sys_halt()
      end if
    end if

  end subroutine afcstab_resizeStabIndGFEM

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
      rafcstabDest%bisSymmetricOperator = rafcstabSrc%bisSymmetricOperator
      rafcstabDest%cafcstabType         = rafcstabSrc%cafcstabType
      rafcstabDest%cprelimitingType     = rafcstabSrc%cprelimitingType
      rafcstabDest%climitingType        = rafcstabSrc%climitingType
      rafcstabDest%cdataType            = rafcstabSrc%cdataType
      rafcstabDest%NEQ                  = rafcstabSrc%NEQ
      rafcstabDest%NVAR                 = rafcstabSrc%NVAR
      rafcstabDest%NVARtransformed      = rafcstabSrc%NVARtransformed
      rafcstabDest%NEDGE                = rafcstabSrc%NEDGE
      rafcstabDest%NNVEDGE              = rafcstabSrc%NNVEDGE
      rafcstabDest%istabilisationSpec   = ior(rafcstabDest%istabilisationSpec,&
                                              AFCSTAB_INITIALISED)
    end if

    ! Copy edge structre
    if (check(idupFlag, AFCSTAB_DUP_EDGELIST) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIST)) then
      ! Copy content from source to destination structure
      call storage_copy(rafcstabSrc%h_IedgeListIdx,&
          rafcstabDest%h_IedgeListIdx)
      call storage_copy(rafcstabSrc%h_IedgeList,&
          rafcstabDest%h_IedgeList)
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIST))
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION))
      ! Reset ownership
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          not(AFCSTAB_SHARE_EDGELIST))
    end if

    ! Copy edge values
    if (check(idupFlag, AFCSTAB_DUP_EDGEVALUES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)) then
      ! Copy content from source to destination structure
      call storage_copy(rafcstabSrc%h_CoeffsAtEdge,&
          rafcstabDest%h_CoeffsAtEdge)
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
      rafcstabDest%bisSymmetricOperator = rafcstabSrc%bisSymmetricOperator
      rafcstabDest%cafcstabType         = rafcstabSrc%cafcstabType
      rafcstabDest%cprelimitingType     = rafcstabSrc%cprelimitingType
      rafcstabDest%climitingType        = rafcstabSrc%climitingType
      rafcstabDest%cdataType            = rafcstabSrc%cdataType
      rafcstabDest%NEQ                  = rafcstabSrc%NEQ
      rafcstabDest%NVAR                 = rafcstabSrc%NVAR
      rafcstabDest%NVARtransformed      = rafcstabSrc%NVARtransformed
      rafcstabDest%NEDGE                = rafcstabSrc%NEDGE
      rafcstabDest%NNVEDGE              = rafcstabSrc%NNVEDGE
      rafcstabDest%istabilisationSpec   = ior(rafcstabDest%istabilisationSpec,&
                                              AFCSTAB_INITIALISED)
    end if

    ! Duplicate edge structre
    if (check(idupFlag, AFCSTAB_DUP_EDGELIST) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIST)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_EDGELIST))) then
        call storage_free(rafcstabDest%h_IedgeListIdx)
        call storage_free(rafcstabDest%h_IedgeList)
      end if
      ! Copy handle from source to destination structure
      rafcstabDest%h_IedgeListIdx = rafcstabSrc%h_IedgeListIdx
      rafcstabDest%h_IedgeList = rafcstabSrc%h_IedgeList
      ! Adjust specifier of the destination structure
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGELIST))
      rafcstabDest%istabilisationSpec = ior(rafcstabDest%istabilisationSpec,&
          iand(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION))
      ! Set ownership to shared
      rafcstabDest%iduplicationFlag = iand(rafcstabDest%iduplicationFlag,&
          AFCSTAB_SHARE_EDGELIST)
    end if

    ! Duplicate edge values
    if (check(idupFlag, AFCSTAB_DUP_EDGEVALUES) .and.&
        check(rafcstabSrc%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)) then
      ! Remove existing data owned by the destination structure
      if (.not.(check(rafcstabDest%iduplicationFlag, AFCSTAB_SHARE_EDGEVALUES)))&
          call storage_free(rafcstabDest%h_CoeffsAtEdge)
      ! Copy handle from source to destination structure
      rafcstabDest%h_CoeffsAtEdge = rafcstabSrc%h_CoeffsAtEdge
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
          call storage_free(rafcstabDest%h_BoundsAtEdge)
      ! Copy handle from source to destination structure
      rafcstabDest%h_BoundsAtEdge = rafcstabSrc%h_BoundsAtEdge
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

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_isMatrixCompatibleSc(rafcstab, rmatrix, bcompatible)

!<description>
    ! This subroutine checks if a scalar matrix and a discrete
    ! stabilisation structure are compatible to each other,
    ! i.e. if they share the same structure, size and so on.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Stabilisation structure
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
        call output_line('Matrix/stabilisation structure not compatible, different structure!',&
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
    ! Block matrix
    type(t_matrixBlock), intent(in) :: rmatrixBlock

    ! Stabilisation structure
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
    if (rafcstab%NVAR  .ne. rmatrixBlock%nblocksPerCol         .or.&
        rafcstab%NEQ   .ne. rmatrixBlock%RmatrixBlock(1,1)%NEQ .or.&
        rafcstab%NEDGE .ne. (rmatrixBlock%RmatrixBlock(1,1)%NA-&
                             rmatrixBlock%RmatrixBlock(1,1)%NEQ)/2) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line('Matrix/stabilisation structure not compatible, different structure!',&
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
    ! Scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! Stabilisation structure
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
        call output_line('Vector/stabilisation structure not compatible, different structure!',&
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
        call output_line('Vector/stabuilisation structure not compatible, different structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcstab_isVectorCompatibleBl')
        call sys_halt()
      end if
    end if
  end subroutine afcstab_isVectorCompatibleBl

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_infoStabilisation(rafcstab)

!<description>
  ! This subroutine prints out information about the stabilisation structure
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>
!</subroutine>

    call output_line('AFCstabilisation:')
    call output_line('-----------------')
    call output_line('bisSymmetricOperator:           '//merge('YES','NO ',rafcstab%bisSymmetricOperator))
    call output_line('cafcstabType:                   '//trim(sys_siL(rafcstab%cafcstabType,15)))
    call output_line('cprelimitingType:               '//trim(sys_siL(rafcstab%cprelimitingType,15)))
    call output_line('climitingType:                  '//trim(sys_siL(rafcstab%climitingType,15)))
    call output_line('cdataType:                      '//trim(sys_siL(rafcstab%cdataType,15)))
    call output_line('iduplicationFlag:               '//trim(sys_siL(int(rafcstab%iduplicationFlag),15)))
    call checkAndOutput('AFCSTAB_SHARE_STRUCTURE:        ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_STRUCTURE)
    call checkAndOutput('AFCSTAB_SHARE_EDGELIST:         ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_EDGELIST)
    call checkAndOutput('AFCSTAB_SHARE_EDGEVALUES:       ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_EDGEVALUES)
    call checkAndOutput('AFCSTAB_SHARE_OFFDIAGONALEDGES: ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_OFFDIAGONALEDGES)
    call checkAndOutput('AFCSTAB_SHARE_ADFLUXES:         ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_ADFLUXES)
    call checkAndOutput('AFCSTAB_SHARE_ADINCREMENTS:     ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_ADINCREMENTS)
    call checkAndOutput('AFCSTAB_SHARE_NODEBOUNDS:       ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_NODEBOUNDS)
    call checkAndOutput('AFCSTAB_SHARE_NODELIMITER:      ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_NODELIMITER)
    call checkAndOutput('AFCSTAB_SHARE_EDGELIMITER:      ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_EDGELIMITER)
    call checkAndOutput('AFCSTAB_SHARE_PREDICTOR:        ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_PREDICTOR)
    call checkAndOutput('AFCSTAB_SHARE_NODEVALUES:       ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_NODEVALUES)
    call checkAndOutput('AFCSTAB_SHARE_EDGEBOUNDS:       ',rafcstab%iduplicationFlag,AFCSTAB_SHARE_EDGEBOUNDS)
    call output_line('istabilisationSpec:             '//trim(sys_siL(int(rafcstab%istabilisationSpec),15)))
    call checkAndOutput('AFCSTAB_INITIALISED:            ',rafcstab%istabilisationSpec,AFCSTAB_INITIALISED)
    call checkAndOutput('AFCSTAB_HAS_EDGELIST:           ',rafcstab%istabilisationSpec,AFCSTAB_HAS_EDGELIST)
    call checkAndOutput('AFCSTAB_HAS_EDGEORIENTATION:    ',rafcstab%istabilisationSpec,AFCSTAB_HAS_EDGEORIENTATION)
    call checkAndOutput('AFCSTAB_HAS_EDGEVALUES:         ',rafcstab%istabilisationSpec,AFCSTAB_HAS_EDGEVALUES)
    call checkAndOutput('AFCSTAB_HAS_OFFDIAGONALEDGES:   ',rafcstab%istabilisationSpec,AFCSTAB_HAS_OFFDIAGONALEDGES)
    call checkAndOutput('AFCSTAB_HAS_ADFLUXES:           ',rafcstab%istabilisationSpec,AFCSTAB_HAS_ADFLUXES)
    call checkAndOutput('AFCSTAB_HAS_ADINCREMENTS:       ',rafcstab%istabilisationSpec,AFCSTAB_HAS_ADINCREMENTS)
    call checkAndOutput('AFCSTAB_HAS_NODEBOUNDS:         ',rafcstab%istabilisationSpec,AFCSTAB_HAS_NODEBOUNDS)
    call checkAndOutput('AFCSTAB_HAS_NODELIMITER:        ',rafcstab%istabilisationSpec,AFCSTAB_HAS_NODELIMITER)
    call checkAndOutput('AFCSTAB_HAS_EDGELIMITER:        ',rafcstab%istabilisationSpec,AFCSTAB_HAS_EDGELIMITER)
    call checkAndOutput('AFCSTAB_HAS_PREDICTOR:          ',rafcstab%istabilisationSpec,AFCSTAB_HAS_PREDICTOR)
    call checkAndOutput('AFCSTAB_HAS_NODEVALUES:         ',rafcstab%istabilisationSpec,AFCSTAB_HAS_NODEVALUES)
    call checkAndOutput('AFCSTAB_HAS_EDGEBOUNDS:         ',rafcstab%istabilisationSpec,AFCSTAB_HAS_EDGEBOUNDS)
    call output_line('NEQ:                            '//trim(sys_siL(rafcstab%NEQ,15)))
    call output_line('NEDGE:                          '//trim(sys_siL(rafcstab%NEDGE,15)))
    call output_line('NNVEDGE:                        '//trim(sys_siL(rafcstab%NNVEDGE,15)))
    call output_line('NVAR:                           '//trim(sys_siL(rafcstab%NVAR,15)))
    call output_line('NVARtransformed:                '//trim(sys_siL(rafcstab%NVARtransformed,15)))
    call output_line('ncoeffsAtEdge:                  '//trim(sys_siL(rafcstab%ncoeffsAtEdge,15)))
    call checkAndOutputHandle('IedgeListIdx:                   ', rafcstab%h_IedgeListIdx)
    call checkAndOutputHandle('IedgeList:                      ', rafcstab%h_IedgeList)
    call checkAndOutputHandle('IsuperdiagEdgesIdx:             ', rafcstab%h_IsuperdiagEdgesIdx)
    call checkAndOutputHandle('IsubdiagEdgesIdx:               ', rafcstab%h_IsubdiagEdgesIdx)
    call checkAndOutputHandle('IsubdiagEdges:                  ', rafcstab%h_IsubdiagEdges)
    call checkAndOutputHandle('CoefficientsAtEdge:             ', rafcstab%h_CoeffsAtEdge)
    call checkAndOutputHandle('BoundsAtEdge:                   ', rafcstab%h_BoundsAtEdge)

    if (associated(rafcstab%p_rvectorAlpha)) then
      call output_lbrk()
      call output_line('vectorAlpha')
      call lsyssc_infoVector(rafcstab%p_rvectorAlpha)
    end if

    if (associated(rafcstab%p_rvectorFlux0)) then
      call output_lbrk()
      call output_line('vectorFlux0')
      call lsyssc_infoVector(rafcstab%p_rvectorFlux0)
    end if

    if (associated(rafcstab%p_rvectorFlux)) then
      call output_lbrk()
      call output_line('vectorFlux')
      call lsyssc_infoVector(rafcstab%p_rvectorFlux)
    end if

    if (associated(rafcstab%p_rvectorFluxPrel)) then
      call output_lbrk()
      call output_line('vectorFluxPrel')
      call lsyssc_infoVector(rafcstab%p_rvectorFluxPrel)
    end if

    if (associated(rafcstab%p_rvectorDx)) then
      call output_lbrk()
      call output_line('vectorDx')
      call lsyssc_infoVector(rafcstab%p_rvectorDx)
    end if

    if (associated(rafcstab%p_rvectorPp)) then
      call output_lbrk()
      call output_line('vectorPp')
      call lsyssc_infoVector(rafcstab%p_rvectorPp)
    end if

    if (associated(rafcstab%p_rvectorPm)) then
      call output_lbrk()
      call output_line('vectorPm')
      call lsyssc_infoVector(rafcstab%p_rvectorPm)
    end if

    if (associated(rafcstab%p_rvectorQ)) then
      call output_lbrk()
      call output_line('vectorQ')
      call lsyssc_infoVector(rafcstab%p_rvectorQ)
    end if

    if (associated(rafcstab%p_rvectorQp)) then
      call output_lbrk()
      call output_line('vectorQp')
      call lsyssc_infoVector(rafcstab%p_rvectorQp)
    end if

    if (associated(rafcstab%p_rvectorQm)) then
      call output_lbrk()
      call output_line('vectorQm')
      call lsyssc_infoVector(rafcstab%p_rvectorQm)
    end if

    if (associated(rafcstab%p_rvectorRp)) then
      call output_lbrk()
      call output_line('vectorRp')
      call lsyssc_infoVector(rafcstab%p_rvectorRp)
    end if

    if (associated(rafcstab%p_rvectorRm)) then
      call output_lbrk()
      call output_line('vectorRm')
      call lsyssc_infoVector(rafcstab%p_rvectorRm)
    end if

    if (associated(rafcstab%p_rvectorPredictor)) then
      call output_lbrk()
      call output_line('vectorPredictor')
      call lsysbl_infoVector(rafcstab%p_rvectorPredictor)
    end if

  contains

    ! Here some working routines follow

    !***************************************************************************
    ! Check bitfield and output string

    subroutine checkAndOutput(cstring, iflag, ibitfield)

      ! input parameters
      character(len=*), intent(in) :: cstring
      integer(I32), intent(in) :: iflag, ibitfield

      if (iand(iflag, ibitfield) .eq. ibitfield) then
        call output_line (cstring//'TRUE')
      else
        call output_line (cstring//'FALSE')
      end if

    end subroutine checkAndOutput

    !***************************************************************************
    ! Check handle and output shape of array

    subroutine checkAndOutputHandle(cstring, h_handle)

      ! input parameters
      character(len=*), intent(in) :: cstring
      integer, intent(in) :: h_handle

      ! local variabels
      integer :: isize,idimension
      integer, dimension(2) :: Isize2D
      integer, dimension(3) :: Isize3D

      if (h_handle .ne. ST_NOHANDLE) then
        call storage_getdimension(h_handle, idimension)

        select case(idimension)
        case (1)
          call storage_getsize(h_handle, isize)
          call output_line (cstring//trim(sys_siL(h_handle,15))//&
              ' ('//trim(sys_siL(isize,15))//')')

        case(2)
          call storage_getsize(h_handle, Isize2D)
          call output_line (cstring//trim(sys_siL(h_handle,15))//&
              ' ('//trim(sys_siL(Isize2D(1),15))//','//trim(sys_siL(Isize2D(2),15))//')')

        case (3)
          call storage_getsize(h_handle, Isize3D)
          call output_line (cstring//trim(sys_siL(h_handle,15))//&
              ' ('//trim(sys_siL(Isize3D(1),15))//','//trim(sys_siL(Isize3D(2),15))//&
                    trim(sys_siL(Isize3D(3),15))//')')
        end select
      else
        call output_line (cstring//trim(sys_siL(h_handle,15)))
      end if

    end subroutine checkAndOutputHandle

  end subroutine afcstab_infoStabilisation

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)

!<description>
    ! Returns a pointer to the index pointer for the edge structure
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the edge index structure
    ! NULL() if the stabilisation structure does not provide it.
    integer, dimension(:), pointer :: p_IedgeListIdx
!</output>
!</subroutine>

    ! Do we edge structure at all?
    if (rafcstab%h_IedgeListIdx .eq. ST_NOHANDLE) then
      nullify(p_IedgeListIdx)
      return
    end if

    ! Get the array
    call storage_getbase_int(rafcstab%h_IedgeListIdx,&
        p_IedgeListIdx)

  end subroutine afcstab_getbase_IedgeListIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IedgeList(rafcstab, p_IedgeList)

!<description>
    ! Returns a pointer to the edge structure
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the edge structure
    ! NULL() if the stabilisation structure does not provide it.
    integer, dimension(:,:), pointer :: p_IedgeList
!</output>
!</subroutine>

    ! Do we edge structure at all?
    if ((rafcstab%h_IedgeList .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE             .eq. 0)) then
      nullify(p_IedgeList)
      return
    end if

    ! Get the array
    call storage_getbase_int2D(rafcstab%h_IedgeList,&
        p_IedgeList,rafcstab%NEDGE)

  end subroutine afcstab_getbase_IedgeList

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)

!<description>
    ! Returns a pointer to the index pointer for edge structure
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

  subroutine afcstab_getbase_QcoeffsAtEdge(rafcstab, p_QcoefficientsAtEdge)

!<description>
    ! Returns a pointer to the quad-valued edge data
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the quad-valued edge data
    ! NULL() if the stabilisation structure does not provide it.
    real(QP), dimension(:,:), pointer :: p_QcoefficientsAtEdge
!</output>
!</subroutine>

    ! Do we have double-valued edge data at all?
    if ((rafcstab%h_CoeffsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                .eq. 0)) then
      nullify(p_QcoefficientsAtEdge)
      return
    end if

    ! Get the array
    call storage_getbase_quad2D(rafcstab%h_CoeffsAtEdge,&
        p_QcoefficientsAtEdge, rafcstab%NEDGE)

  end subroutine afcstab_getbase_QcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

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
    if ((rafcstab%h_CoeffsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                .eq. 0)) then
      nullify(p_DcoefficientsAtEdge)
      return
    end if

    ! Get the array
    call storage_getbase_double2D(rafcstab%h_CoeffsAtEdge,&
        p_DcoefficientsAtEdge, rafcstab%NEDGE)

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
    if ((rafcstab%h_CoeffsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE                .eq. 0)) then
      nullify(p_FcoefficientsAtEdge)
      return
    end if

    ! Get the array
    call storage_getbase_single2D(rafcstab%h_CoeffsAtEdge,&
        p_FcoefficientsAtEdge, rafcstab%NEDGE)

  end subroutine afcstab_getbase_FcoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_QboundsAtEdge(rafcstab,p_QboundsAtEdge)

!<description>
    ! Returns a pointer to the quad-valued bounds at edges
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the quad-valued bounds at edges
    ! NULL() if the stabilisation structure does not provide it.
    real(QP), dimension(:,:), pointer :: p_QboundsAtEdge
!</output>
!</subroutine>

    ! Do we have bounds at edges data at all?
    if ((rafcstab%h_BoundsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE          .eq. 0)) then
      nullify(p_QboundsAtEdge)
      return
    end if

    ! Get the array
    call storage_getbase_quad2D(rafcstab%h_BoundsAtEdge,&
        p_QboundsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_QboundsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_DboundsAtEdge(rafcstab,p_DboundsAtEdge)

!<description>
    ! Returns a pointer to the double-valued bounds at edges
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the double-valued bounds at edges
    ! NULL() if the stabilisation structure does not provide it.
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
!</output>
!</subroutine>

    ! Do we have bounds at edges data at all?
    if ((rafcstab%h_BoundsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE          .eq. 0)) then
      nullify(p_DboundsAtEdge)
      return
    end if

    ! Get the array
    call storage_getbase_double2D(rafcstab%h_BoundsAtEdge,&
        p_DboundsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_DboundsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_getbase_FboundsAtEdge(rafcstab,p_FboundsAtEdge)

!<description>
    ! Returns a pointer to the single-valued bounds at edges
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab
!</input>

!<output>
    ! Pointer to the single-valued bounds at edges
    ! NULL() if the stabilisation structure does not provide it.
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
!</output>
!</subroutine>

    ! Do we have bounds at edges data at all?
    if ((rafcstab%h_BoundsAtEdge .eq. ST_NOHANDLE) .or.&
        (rafcstab%NEDGE          .eq. 0)) then
      nullify(p_FboundsAtEdge)
      return
    end if

    ! Get the array
    call storage_getbase_single2D(rafcstab%h_BoundsAtEdge,&
        p_FboundsAtEdge,rafcstab%NEDGE)

  end subroutine afcstab_getbase_FboundsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_genEdgeList(rmatrix, rafcstab)

!<description>
    ! This subroutine generates the list of edges which are
    ! characterised by their two endpoints (i,j) and the absolute
    ! position of matrix entries ij and ji.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
#if defined(USE_OPENMP) || defined(ENABLE_COPROCESSOR_SUPPORT)
    integer, dimension(:,:), pointer :: p_IedgeList
#endif

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_genEdgeList')
      call sys_halt()
    end if

    ! Check if edge structure is owned by the stabilisation structure
    if (iand(rafcstab%iduplicationFlag, AFCSTAB_SHARE_EDGELIST) .eq.&
        AFCSTAB_SHARE_EDGELIST) then
      call output_line('Edge structure is not owned by stabilisation and '//&
          'therefore cannot be generated',&
          OU_CLASS_WARNING,OU_MODE_STD,'afcstab_genEdgeList')
      return
    end if

    ! Generate edge list for a matrix which is structurally symmetric,
    ! i.e. edge (i,j) exists if and only if edge (j,i) exists without
    ! storing the diagonal edges (i,i).
    call lsyssc_genEdgeList(rmatrix, rafcstab%h_IedgeList,&
        LSYSSC_EDGELIST_NODESANDPOS, .true., .true., rafcstab%NEDGE)

    ! Allocate memory
    if (rafcstab%h_IedgeListIdx .eq. ST_NOHANDLE) then
      call storage_new('gfem_genEdgeList', 'IedgeListIdx',&
          2, ST_INT, rafcstab%h_IedgeListIdx, ST_NEWBLOCK_ZERO)
    end if

    ! Set pointer to edge structure
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)

    ! If no OpenMP is used, then all edges belong to the same
    ! group. Otherwise, the edges will be reordered below.
    p_IedgeListIdx    = rafcstab%NEDGE+1
    p_IedgeListIdx(1) = 1

#if defined(USE_OPENMP) || defined(ENABLE_COPROCESSOR_SUPPORT)
    ! OpenMP-Extension: Perform edge-coloring to find groups of
    ! edges which can be processed in parallel, that is, the
    ! nodes of the edges in the group are all distinct
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX1) then
      call lsyssc_regroupEdgeList(rmatrix%NEQ, p_IedgeList,&
          rafcstab%h_IedgeListIdx, 2*(rmatrix%NEQ-1))
    else
      call lsyssc_regroupEdgeList(rmatrix%NEQ, p_IedgeList,&
          rafcstab%h_IedgeListIdx)
    end if
#endif

    ! Set state of stabiliation
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)

  end subroutine afcstab_genEdgeList

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_genOffdiagEdges(rafcstab)

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
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer   :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer   :: p_IsubdiagEdges
    integer, dimension(:), pointer   :: p_IsubdiagEdgesIdx
    integer :: iedge,nedge,istor
    integer :: ieq,jeq,neq
    integer :: isize

    ! Check if edge-based data structure is prepared
    if (iand(rafcstab%istabilisationSpec,AFCSTAB_HAS_EDGELIST) .eq. 0) then
      call output_line('Stabilisation structure does not provide required ' //&
          'edge-based data structure',OU_CLASS_ERROR,OU_MODE_STD,&
          'afcstab_genOffdiagEdge')
      call sys_halt()
    end if

    ! Store dimensions of stabilisation structure
    neq   = rafcstab%NEQ
    nedge = rafcstab%NEDGE

    ! Allocate memory (if required)
    if (rafcstab%h_IsuperdiagEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new('afcstab_genOffdiagEdges','IsuperdiagEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsuperdiagEdgesIdx,ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rafcstab%h_IsuperdiagEdgesIdx,isize)
      if (isize < neq+1) then
        call storage_realloc('afcstab_genOffdiagEdges',neq+1,&
            rafcstab%h_IsuperdiagEdgesIdx,ST_NEWBLOCK_ZERO,.false.)
      else
        call storage_clear(rafcstab%h_IsuperdiagEdgesIdx)
      end if
    end if

    ! Allocate memory (if required)
    if (rafcstab%h_IsubdiagEdgesIdx .eq. ST_NOHANDLE) then
      call storage_new('afcstab_genOffdiagEdges','IsubdiagEdgesIdx',&
          neq+1,ST_INT,rafcstab%h_IsubdiagEdgesIdx,ST_NEWBLOCK_ZERO)
    else
      call storage_getsize(rafcstab%h_IsubdiagEdgesIdx,isize)
      if (isize < neq+1) then
        call storage_realloc('afcstab_genOffdiagEdges',neq+1,&
            rafcstab%h_IsubdiagEdgesIdx,ST_NEWBLOCK_ZERO,.false.)
      else
        call storage_clear(rafcstab%h_IsubdiagEdgesIdx)
      end if
    end if

    ! Allocate memory (if required)
    if (rafcstab%h_IsubdiagEdges .eq. ST_NOHANDLE) then
      call storage_new('afcstab_genOffdiagEdges','IsubdiagEdges',&
          nedge,ST_INT,rafcstab%h_IsubdiagEdges,ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(rafcstab%h_IsubdiagEdges,isize)
      if (isize < nedge) then
        call storage_realloc('afcstab_genOffdiagEdges',nedge,&
            rafcstab%h_IsubdiagEdges,ST_NEWBLOCK_NOINIT,.false.)
      end if
    end if

    ! Set pointers
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab,p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IedgeList(rafcstab,p_IedgeList)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab,p_IsubdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab,p_IsubdiagEdges)

    ! Count the number of edges in each row
    do iedge = 1, nedge

       ! Determine the start-point of the current edge
      ieq = min(p_IedgeList(1,iedge),&
                p_IedgeList(2,iedge))

      ! Determine the end-point of the current edge
      jeq = max(p_IedgeList(1,iedge),&
                p_IedgeList(2,iedge))

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
        jeq = p_IedgeList(1,iedge) + p_IedgeList(2,iedge) - ieq

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

  end subroutine afcstab_genOffdiagEdges

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_genExtSparsity(rmatrix, rmatrixDest)

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

!<inputoutput>
    ! scalar matrix with standard sparsity pattern
    type(t_matrixScalar), intent(inout)    :: rmatrix

    ! OPTIONAL: scalar matrix with extended sparsity pattern
    ! If not present, then the source matrix will be overwritten
    type(t_matrixScalar), intent(inout), optional :: rmatrixDest
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_matrixScalar) :: rmatrixTmp

!!$    ! Clear output matrix
!!$    if (lsyssc_hasMatrixStructure(rmatrixExtended) .or.&
!!$        lsyssc_hasMatrixContent(rmatrixExtended)) then
!!$      call lsyssc_releaseMatrix(rmatrixExtended)
!!$    end if

    ! Compute Z=A*A and let the connectivity graph of Z be the
    ! extended sparsity pattern of the Jacobian matrix
    call lsyssc_multMatMat(rmatrix, rmatrix, rmatrixTmp,&
                           .true., .true., .false.)

    if (present(rmatrixDest)) then
      call lsyssc_moveMatrix(rmatrixTmp, rmatrixDest)
    else
      call lsyssc_moveMatrix(rmatrixTmp, rmatrix)
    end if
    call lsyssc_releaseMatrix(rmatrixTmp)
    
  end subroutine afcstab_genExtSparsity

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_allocInternalData(rafcstab, ballocCoefficients,&
      rblockDiscretisation, rdiscretisation)

!<description>
    ! This subroutine allocates the internal data structures of the
    ! stabilisation structure rafcstab.
    ! Note that this routine does no error checking, and thus, should
    ! not be called from user-defined routine directly.
    ! This routine doe not generate auxiliary data structures such as
    ! the edge-based structure. Depending on the specified type of
    ! stabilisation auxiliary arrays such as the antidiffusive fluxes
    ! or the upper and lower bounds are allocated.
!</description>

!<input>
    ! Flag: if.true. then auxiliary coefficients are allocated
    logical, intent(in) :: ballocCoefficients

    ! OPTIONAL: block discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_blockDiscretisation), intent(in), optional :: rblockDiscretisation

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary 1-block vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>

    ! local variables
    type(t_vectorScalar) :: rvector

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)

      ! No internal data structures are required

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE)

      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      if (ballocCoefficients) call afcstab_allocCoeffsAtEdge(rafcstab,3)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the 3 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux0(rafcstab)
      call afcstab_allocFlux(rafcstab)

      ! We need the edgewise vector for the prelimited fluxes (if any)
      ! or for the semi-implicit version of the FCT algorithm
      if ((rafcstab%cprelimitingType .ne. AFCSTAB_PRELIMITING_NONE) .or.&
          (rafcstab%cafcstabType     .eq. AFCSTAB_NLINFCT_IMPLICIT)) then
        call afcstab_allocFluxPrel(rafcstab)
      end if

      ! We need the nodal vector for the predictor
      allocate(rafcstab%p_rvectorPredictor)

      if (present(rblockDiscretisation)) then
        ! Create block vector by block discretisation
        call lsysbl_createVector(rblockDiscretisation,&
            rafcstab%p_rvectorPredictor, .false., rafcstab%cdataType)
      elseif (present(rdiscretisation)) then
        ! Create scalar vector by spatial discretisation and convert
        ! it into a 1-block vector afterwards
        call lsyssc_createVector(rdiscretisation, rvector, .false.,&
            rafcstab%cdataType)
        call lsysbl_convertVecFromScalar(rvector, rafcstab%p_rvectorPredictor)
        call lsyssc_releaseVector(rvector)
      else
        ! Create 1-block vector directly
        call lsysbl_createVector(rafcstab%p_rvectorPredictor,&
            rafcstab%NEQ, 1, .false., rafcstab%cdataType)
      end if

      !-------------------------------------------------------------------------

    case (AFCSTAB_TVD,&
          AFCSTAB_GP)

      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      if (ballocCoefficients) call afcstab_allocCoeffsAtEdge(rafcstab,3)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the 3 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux0(rafcstab)
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case (AFCSTAB_LINFCT,&
          AFCSTAB_LINFCT_MASS)

      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      if (ballocCoefficients) call afcstab_allocCoeffsAtEdge(rafcstab,3)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)

      ! We need the edgewise vector if raw antidiffusive fluxes should be prelimited
      if (rafcstab%cprelimitingType .ne. AFCSTAB_PRELIMITING_NONE) then
        call afcstab_allocFluxPrel(rafcstab)
      end if

      !-------------------------------------------------------------------------

    case (AFCSTAB_SYMMETRIC)

      ! Handle for DcoefficientsAtEdge: (/d_ij,s_ij/)
      if (ballocCoefficients) call afcstab_allocCoeffsAtEdge(rafcstab,2)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the edgewise vector for the fluxes
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_MASS,&
          AFCSTAB_LINLPT_MASS)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      ! and one extra nodal vector Q
      call afcstab_allocVectorsPQR(rafcstab, ballocCommonQ=.true.)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)

      ! We need the nodal vector for the predictor
      allocate(rafcstab%p_rvectorPredictor)

      ! We need the nodal vector for the predictor
      if (present(rblockDiscretisation)) then
        ! Create block vector by block discretisation
        call lsysbl_createVector(rblockDiscretisation,&
            rafcstab%p_rvectorPredictor, .false., rafcstab%cdataType)
      elseif (present(rdiscretisation)) then
        ! Create scalar vector by spatial discretisation and convert
        ! it into a 1-block vector afterwards
        call lsyssc_createVector(rdiscretisation, rvector, .false.,&
            rafcstab%cdataType)
        call lsysbl_convertVecFromScalar(rvector, rafcstab%p_rvectorPredictor)
        call lsyssc_releaseVector(rvector)
      else
        ! Create 1-block vector directly
        call lsysbl_createVector(rafcstab%p_rvectorPredictor,&
            rafcstab%NEQ, 1, .false., rafcstab%cdataType)
      end if

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_UPWINDBIASED,&
          AFCSTAB_LINLPT_UPWINDBIASED)

      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      if (ballocCoefficients) call afcstab_allocCoeffsAtEdge(rafcstab,3)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      ! and one extra nodal vector Q
      call afcstab_allocVectorsPQR(rafcstab, ballocCommonQ=.true.)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_SYMMETRIC,&
          AFCSTAB_LINLPT_SYMMETRIC)

      ! Handle for DcoefficientsAtEdge: (/d_ij,s_ij/)
      if (ballocCoefficients) call afcstab_allocCoeffsAtEdge(rafcstab,2)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      ! and one extra nodal vector Q
      call afcstab_allocVectorsPQR(rafcstab, ballocCommonQ=.true.)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_allocInternalData')
      call sys_halt()
    end select

  end subroutine afcstab_allocInternalData

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

    ! Handle for IedgeListIdx
    if (rafcstab%h_IedgeListIdx .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IedgeListIdx)
    call storage_new('afcstab_allocEdgeStructure', 'IedgeListIdx',&
        2, ST_INT, rafcstab%h_IedgeListIdx, ST_NEWBLOCK_NOINIT)

    ! Handle for IedgeList: (/i,j,ij,ji/)
    Isize = (/ndata, rafcstab%NEDGE/)
    if (rafcstab%h_IedgeList .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_IedgeList)
    call storage_new('afcstab_allocEdgeStructure', 'IedgeList',&
        Isize, ST_INT, rafcstab%h_IedgeList, ST_NEWBLOCK_NOINIT)

    ! Set ownership
    rafcstab%iduplicationFlag = iand(rafcstab%iduplicationFlag,&
        not(AFCSTAB_SHARE_EDGELIST))

  end subroutine afcstab_allocEdgeStructure

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_allocCoeffsAtEdge(rafcstab, ncoeffsAtEdge, cdataType)

!<description>
    ! This subroutine allocates the coefficients at edge structuree
!</description>

!<input>
    ! Number of coefficients at edge
    integer, intent(in) :: ncoeffsAtEdge

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

    rafcstab%ncoeffsAtEdge = ncoeffsAtEdge
    Isize = (/rafcstab%ncoeffsAtEdge, rafcstab%NEDGE/)
    if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_CoeffsAtEdge)
    call storage_new('afcstab_allocCoeffsAtEdge', 'DcoefficientsAtEdge',&
        Isize, cType, rafcstab%h_CoeffsAtEdge, ST_NEWBLOCK_NOINIT)

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
    if (rafcstab%h_BoundsAtEdge .ne. ST_NOHANDLE)&
        call storage_free(rafcstab%h_BoundsAtEdge)
    call storage_new('afcstab_allocBoundsAtEdge', 'DboundsAtEdge',&
        Isize, cType, rafcstab%h_BoundsAtEdge, ST_NEWBLOCK_NOINIT)

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

  subroutine afcstab_buildBoundsLPT1D(rafcstab, rmatrixCx, rmatrixM,&
      DdofCoords, dscale)

!<description>
    ! This subroutine computes the solution bounds for
    ! linearity-preserving flux correction in 1D
!</description>

!<input>
    ! Coefficient matrix
    type(t_matrixScalar), intent(in) :: rmatrixCx

    ! Diagonal scaling matrix
    type(t_matrixScalar), intent(in) :: rmatrixM

    ! Global coordinates of the degrees of freedom
    real(DP), dimension(:), intent(in) :: DdofCoords

    ! Global scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
    real(DP), dimension(:), pointer :: p_DdataM, p_DdataCx
    real(SP), dimension(:), pointer :: p_FdataM, p_FdataCx
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kdiagonal, p_Kld

    ! Check spatial discretisation structure
    if (.not. rmatrixCx%bidenticalTrialAndTest) then
      call output_line('Test and trial functions must be identical!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_genEdgeList(rmatrixCx, rafcstab)
    end if

    ! Check if stabilisation structure provides array to store bounds;
    ! otherwise allocate new memory using the predefined data type
    if (rafcstab%h_BoundsAtEdge .eq. ST_NOHANDLE)&
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
    if ((rmatrixM%cmatrixFormat .ne. LSYSSC_MATRIXD)  .or.&
        (rmatrixM%cdataType .ne. rmatrixCx%cdataType) .or.&
        (rmatrixM%NEQ .ne. rmatrixCx%NEQ)) then
      call output_line('Coefficient matrices are not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)

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

        call doBoundsMat7DP(p_IedgeList, p_Kld,&
            p_DdataCx, p_DdataM, DdofCoords, dscale, p_DboundsAtEdge)

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call afcstab_getbase_FboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat7SP(p_IedgeList, p_Kld,&
            p_FdataCx, p_FdataM, DdofCoords, dscale, p_FboundsAtEdge)

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

        call doBoundsMat9DP(p_IedgeList, p_Kdiagonal, p_Kld,&
            p_DdataCx, p_DdataM, DdofCoords, dscale, p_DboundsAtEdge)

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call afcstab_getbase_FboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat9SP(p_IedgeList, p_Kdiagonal, p_Kld,&
            p_FdataCx, p_FdataM, DdofCoords, dscale, p_FboundsAtEdge)

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

    subroutine doBoundsMat7DP(IedgeList, Kld,&
        DdataCx, DdataM, DdofCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kld
      real(DP), dimension(:), intent(in) :: DdataCx,DdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DdofCoords(i) - DdofCoords(j)

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

    end subroutine doBoundsMat7DP

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9DP(IedgeList, Kdiagonal, Kld,&
        DdataCx, DdataM, DdofCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DdofCoords(i) - DdofCoords(j)

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

    end subroutine doBoundsMat9DP

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7SP(IedgeList, Kld,&
        FdataCx, FdataM, DdofCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DdofCoords(i) - DdofCoords(j)

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

    end subroutine doBoundsMat7SP

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9SP(IedgeList, Kdiagonal, Kld,&
        FdataCx, FdataM, DdofCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diff
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diff,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diff = DdofCoords(i) - DdofCoords(j)

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

    end subroutine doBoundsMat9SP

  end subroutine afcstab_buildBoundsLPT1D

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_buildBoundsLPT2D(rafcstab, rmatrixCx, rmatrixCy,&
      rmatrixM, DdofCoords, dscale)

!<description>
    ! This subroutine computes the solution bounds for
    ! linearity-preserving flux correction in 2D
!</description>

!<input>
    ! Coefficient matrices
    type(t_matrixScalar), intent(in) :: rmatrixCx, rmatrixCy

    ! Diagonal scaling matrix
    type(t_matrixScalar), intent(in) :: rmatrixM

    ! Global coordinates of the degrees of freedom
    real(DP), dimension(:), intent(in) :: DdofCoords

    ! Global scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
    real(DP), dimension(:), pointer :: p_DdataM, p_DdataCx, p_DdataCy
    real(SP), dimension(:), pointer :: p_FdataM, p_FdataCx, p_FdataCy
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kdiagonal, p_Kld

    ! Check spatial discretisation structure
    if (.not. rmatrixCx%bidenticalTrialAndTest) then
      call output_line('Test and trial functions must be identical!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT2D')
      call sys_halt()
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_genEdgeList(rmatrixCx, rafcstab)
    end if

    ! Check if stabilisation structure provides array to store bounds;
    ! otherwise allocate new memory using the predefined data type
    if (rafcstab%h_BoundsAtEdge .eq. ST_NOHANDLE)&
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
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)

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

        call doBoundsMat7DP(p_IedgeList, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataM, DdofCoords,&
            dscale, p_DboundsAtEdge)

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call afcstab_getbase_FboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat7SP(p_IedgeList, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataM, DdofCoords,&
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

        call doBoundsMat9DP(p_IedgeList, p_Kdiagonal, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataM, DdofCoords,&
            dscale, p_DboundsAtEdge)

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call afcstab_getbase_FboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat9SP(p_IedgeList, p_Kdiagonal, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataM, DdofCoords,&
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

    subroutine doBoundsMat7DP(IedgeList, Kld,&
        DdataCx, DdataCy, DdataM, DdofCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM2D*(i-1)+1) - DdofCoords(NDIM2D*(j-1)+1)
        diffY = DdofCoords(NDIM2D*(i-1)+2) - DdofCoords(NDIM2D*(j-1)+2)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*diffX+&
                            DdataCy(ik)*diffY)
        end do

        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*diffX+&
                            DdataCy(jk)*diffY)
        end do

        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat7DP

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9DP(IedgeList, Kdiagonal, Kld,&
        DdataCx, DdataCy, DdataM, DdofCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM2D*(i-1)+1) - DdofCoords(NDIM2D*(j-1)+1)
        diffY = DdofCoords(NDIM2D*(i-1)+2) - DdofCoords(NDIM2D*(j-1)+2)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(DdataCx(ik)*diffX+&
                            DdataCy(ik)*diffY)
        end do

        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*diffX+&
                            DdataCy(ik)*diffY)
        end do

        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(DdataCx(jk)*diffX+&
                            DdataCy(jk)*diffY)
        end do

        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*diffX+&
                            DdataCy(jk)*diffY)
        end do

        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat9DP

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7SP(IedgeList, Kld,&
        FdataCx, FdataCy, FdataM, DdofCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM2D*(i-1)+1) - DdofCoords(NDIM2D*(j-1)+1)
        diffY = DdofCoords(NDIM2D*(i-1)+2) - DdofCoords(NDIM2D*(j-1)+2)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*diffX+&
                            FdataCy(ik)*diffY)
        end do

        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*diffX+&
                            FdataCy(jk)*diffY)
        end do

        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat7SP

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9SP(IedgeList, Kdiagonal, Kld,&
        FdataCx, FdataCy, FdataM, DdofCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM2D*(i-1)+1) - DdofCoords(NDIM2D*(j-1)+1)
        diffY = DdofCoords(NDIM2D*(i-1)+2) - DdofCoords(NDIM2D*(j-1)+2)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(FdataCx(ik)*diffX+&
                            FdataCy(ik)*diffY)
        end do

        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*diffX+&
                            FdataCy(ik)*diffY)
        end do

        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(FdataCx(jk)*diffX+&
                            FdataCy(jk)*diffY)
        end do

        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*diffX+&
                            FdataCy(jk)*diffY)
        end do

        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat9SP

  end subroutine afcstab_buildBoundsLPT2D

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_buildBoundsLPT3D(rafcstab, rmatrixCx, rmatrixCy,&
      rmatrixCz, rmatrixM, DdofCoords, dscale)

!<description>
    ! This subroutine computes the solution bounds for
    ! linearity-preserving flux correction in 3D
!</description>

!<input>
    ! Coefficient matrices
    type(t_matrixScalar), intent(in) :: rmatrixCx, rmatrixCy, rmatrixCz

    ! Diagonal scaling matrix
    type(t_matrixScalar), intent(in) :: rmatrixM

    ! Global coordinates of the degrees of freedom
    real(DP), dimension(:), intent(in) :: DdofCoords

    ! Global scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(SP), dimension(:,:), pointer :: p_FboundsAtEdge
    real(DP), dimension(:), pointer :: p_DdataM, p_DdataCx, p_DdataCy, p_DdataCz
    real(SP), dimension(:), pointer :: p_FdataM, p_FdataCx, p_FdataCy, p_FdataCz
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kdiagonal, p_Kld

    ! Check spatial discretisation structure
    if (.not. rmatrixCx%bidenticalTrialAndTest) then
      call output_line('Test and trial functions must be identical!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT3D')
      call sys_halt()
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_buildBoundsLPT1D')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_genEdgeList(rmatrixCx, rafcstab)
    end if

    ! Check if stabilisation structure provides array to store bounds;
    ! otherwise allocate new memory using the predefined data type
    if (rafcstab%h_BoundsAtEdge .eq. ST_NOHANDLE)&
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
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)

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

        call doBoundsMat7DP(p_IedgeList, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataCz, p_DdataM,&
            DdofCoords, dscale, p_DboundsAtEdge)

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call lsyssc_getbase_single (rmatrixCz, p_FdataCz)
        call afcstab_getbase_FboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat7SP(p_IedgeList, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataCz, p_FdataM,&
            DdofCoords, dscale, p_FboundsAtEdge)

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

        call doBoundsMat9DP(p_IedgeList, p_Kdiagonal, p_Kld,&
            p_DdataCx, p_DdataCy, p_DdataCz, p_DdataM,&
            DdofCoords, dscale, p_DboundsAtEdge)

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixM, p_FdataM)
        call lsyssc_getbase_single (rmatrixCx, p_FdataCx)
        call lsyssc_getbase_single (rmatrixCy, p_FdataCy)
        call lsyssc_getbase_single (rmatrixCz, p_FdataCz)
        call afcstab_getbase_FboundsAtEdge(rafcstab, p_FboundsAtEdge)

        call doBoundsMat9SP(p_IedgeList, p_Kdiagonal, p_Kld,&
            p_FdataCx, p_FdataCy, p_FdataCz, p_FdataM,&
            DdofCoords, dscale, p_FboundsAtEdge)

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

    subroutine doBoundsMat7DP(IedgeList, Kld,&
        DdataCx, DdataCy, DdataCz, DdataM, DdofCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataCz, DdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY,diffZ
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,diffZ,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM3D*(i-1)+1) - DdofCoords(NDIM3D*(j-1)+1)
        diffY = DdofCoords(NDIM3D*(i-1)+2) - DdofCoords(NDIM3D*(j-1)+2)
        diffZ = DdofCoords(NDIM3D*(i-1)+3) - DdofCoords(NDIM3D*(j-1)+3)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*diffX+&
                            DdataCy(ik)*diffY+&
                            DdataCz(ik)*diffZ)
        end do

        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*diffX+&
                            DdataCy(jk)*diffY+&
                            DdataCz(jk)*diffZ)
        end do

        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat7DP

    !**************************************************************
    ! Compute the bounds in double precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9DP(IedgeList, Kdiagonal, Kld,&
        DdataCx, DdataCy, DdataCz, DdataM, DdofCoords, dscale, DboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(DP), dimension(:), intent(in) :: DdataCx, DdataCy, DdataCz, DdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(DP), dimension(:,:), intent(out) :: DboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY,diffZ
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,diffZ,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM3D*(i-1)+1) - DdofCoords(NDIM3D*(j-1)+1)
        diffY = DdofCoords(NDIM3D*(i-1)+2) - DdofCoords(NDIM3D*(j-1)+2)
        diffZ = DdofCoords(NDIM3D*(i-1)+3) - DdofCoords(NDIM3D*(j-1)+3)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(DdataCx(ik)*diffX+&
                            DdataCy(ik)*diffY+&
                            DdataCz(ik)*diffZ)
        end do

        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(DdataCx(ik)*diffX+&
                            DdataCy(ik)*diffY+&
                            DdataCz(ik)*diffZ)
        end do

        ! Store result into first entry
        DboundsAtEdge(1,iedge) = dscale*daux/DdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(DdataCx(jk)*diffX+&
                            DdataCy(jk)*diffY+&
                            DdataCz(jk)*diffZ)
        end do

        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(DdataCx(jk)*diffX+&
                            DdataCy(jk)*diffY+&
                            DdataCz(jk)*diffZ)
        end do

        ! Store result into first entry
        DboundsAtEdge(2,iedge) = dscale*daux/DdataM(j)
      end do edges

    end subroutine doBoundsMat9DP

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 7.

    subroutine doBoundsMat7SP(IedgeList, Kld,&
        FdataCx, FdataCy, FdataCz, FdataM, DdofCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataCz, FdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY,diffZ
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,diffZ,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM3D*(i-1)+1) - DdofCoords(NDIM3D*(j-1)+1)
        diffY = DdofCoords(NDIM3D*(i-1)+2) - DdofCoords(NDIM3D*(j-1)+2)
        diffZ = DdofCoords(NDIM3D*(i-1)+3) - DdofCoords(NDIM3D*(j-1)+3)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*diffX+&
                            FdataCy(ik)*DiffY+&
                            FdataCz(ik)*DiffZ)
        end do

        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*diffX+&
                            FdataCy(jk)*diffY+&
                            FdataCz(jk)*diffZ)
        end do

        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat7SP

    !**************************************************************
    ! Compute the bounds in single precision.
    ! All matrices are stored in matrix format 9.

    subroutine doBoundsMat9SP(IedgeList, Kdiagonal, Kld,&
        FdataCx, FdataCy, FdataCz, FdataM, DdofCoords, dscale, FboundsAtEdge)

      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal, Kld
      real(SP), dimension(:), intent(in) :: FdataCx, FdataCy, FdataCz, FdataM
      real(DP), dimension(:), intent(in) :: DdofCoords
      real(DP), intent(in) :: dscale

      real(SP), dimension(:,:), intent(out) :: FboundsAtEdge

      ! local variables
      real(DP) :: daux,diffX,diffY,diffZ
      integer :: iedge,i,j,ik,jk

      ! Loop over all edges of the structure
      !$omp parallel do default(shared)&
      !$omp private(daux,diffX,diffY,diffZ,i,ik,j,jk)
      edges: do iedge = 1, size(IedgeList,2)

        ! Get the numbers of the two endpoints
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Compute difference (x_i-x_j)
        diffX = DdofCoords(NDIM3D*(i-1)+1) - DdofCoords(NDIM3D*(j-1)+1)
        diffY = DdofCoords(NDIM3D*(i-1)+2) - DdofCoords(NDIM3D*(j-1)+2)
        diffZ = DdofCoords(NDIM3D*(i-1)+3) - DdofCoords(NDIM3D*(j-1)+3)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{ik}*(x_i-x_j)| $
        do ik = Kld(i), Kdiagonal(i)-1
          daux = daux + abs(FdataCx(ik)*diffX+&
                            FdataCy(ik)*diffY+&
                            FdataCz(ik)*diffZ)
        end do

        do ik = Kdiagonal(i)+1, Kld(i+1)-1
          daux = daux + abs(FdataCx(ik)*diffX+&
                            FdataCy(ik)*diffY+&
                            FdataCz(ik)*diffZ)
        end do

        ! Store result into first entry
        FboundsAtEdge(1,iedge) = real(dscale*daux,SP)/FdataM(i)

        ! Clear temporal variable
        daux = 0.0_DP

        ! Loop over all off-diagonal entries of current row and sum
        ! the absolute values $ |c_{jk}*(x_j-x_i)| $
        do jk = Kld(j), Kdiagonal(j)-1
          daux = daux + abs(FdataCx(jk)*diffX+&
                            FdataCy(jk)*diffY+&
                            FdataCz(jk)*diffZ)
        end do

        do jk = Kdiagonal(j)+1, Kld(j+1)-1
          daux = daux + abs(FdataCx(jk)*diffX+&
                            FdataCy(jk)*diffY+&
                            FdataCz(jk)*diffZ)
        end do

        ! Store result into first entry
        FboundsAtEdge(2,iedge) = real(dscale*daux,SP)/FdataM(j)
      end do edges

    end subroutine doBoundsMat9SP

  end subroutine afcstab_buildBoundsLPT3D

  !*****************************************************************************

!<function>

  elemental function afcstab_limitUnboundedQP(p, q, qval) result(r)

!<description>
    ! This function computes the ratio q/p. If the denominator is
    ! too small, then the default value dval is applied.
!</description>

!<input>
    ! (de)nominator
    real(QP), intent(in) :: p,q

    ! default value
    real(QP), intent(in) :: qval
!</input>

!<result>
    ! limited ratio
    real(QP) :: r
!</result>
!</function>

    if (abs(p) .gt. AFCSTAB_EPSABS) then
      r = q/p
    else
      r = qval
    end if
  end function afcstab_limitUnboundedQP

  !*****************************************************************************

!<function>

  elemental function afcstab_limitUnboundedDP(p, q, dval) result(r)

!<description>
    ! This function computes the ratio q/p. If the denominator is
    ! too small, then the default value dval is applied.
!</description>

!<input>
    ! (de)nominator
    real(DP), intent(in) :: p,q

    ! default value
    real(DP), intent(in) :: dval
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>

    if (abs(p) .gt. AFCSTAB_EPSABS) then
      r = q/p
    else
      r = dval
    end if
  end function afcstab_limitUnboundedDP

  !*****************************************************************************

!<function>

  elemental function afcstab_limitUnboundedSP(p, q, fval) result(r)

!<description>
    ! This function computes the ratio q/p. If the denominator is
    ! too small, then the default value fval is applied.
!</description>

!<input>
    ! (de)nominator
    real(SP), intent(in) :: p,q

    ! default value
    real(SP), intent(in) :: fval
!</input>

!<result>
    ! limited ratio
    real(SP) :: r
!</result>
!</function>

    if (abs(p) .gt. AFCSTAB_EPSABS) then
      r = q/p
    else
      r = fval
    end if
  end function afcstab_limitUnboundedSP

  !*****************************************************************************

!<function>

  elemental function afcstab_limitBoundedQP(p, q, qval, qbound) result(r)

!<description>
    ! This function computes the limited ratio q/p and bounds the
    ! result by the size of dbound. If the denominator is too small
    ! then the default value dval is applied.
!</description>

!<input>
    ! (de)nominator
    real(QP), intent(in) :: p,q

    ! default value
    real(QP), intent(in) :: qval

    ! upper bound
    real(QP), intent(in) :: qbound
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>

    if (abs(p) .gt. AFCSTAB_EPSABS) then
      r = min(q/p, qbound)
    else
      r = qval
    end if
  end function afcstab_limitBoundedQP

  !*****************************************************************************

!<function>

  elemental function afcstab_limitBoundedDP(p, q, dval, dbound) result(r)

!<description>
    ! This function computes the limited ratio q/p and bounds the
    ! result by the size of dbound. If the denominator is too small
    ! then the default value dval is applied.
!</description>

!<input>
    ! (de)nominator
    real(DP), intent(in) :: p,q

    ! default value
    real(DP), intent(in) :: dval

    ! upper bound
    real(DP), intent(in) :: dbound
!</input>

!<result>
    ! limited ratio
    real(DP) :: r
!</result>
!</function>

    if (abs(p) .gt. AFCSTAB_EPSABS) then
      r = min(q/p, dbound)
    else
      r = dval
    end if
  end function afcstab_limitBoundedDP

  !*****************************************************************************

!<function>

  elemental function afcstab_limitBoundedSP(p, q, fval, fbound) result(r)

!<description>
    ! This function computes the limited ratio q/p and bounds the
    ! result by the size of dbound. If the denominator is too small
    ! then the default value dval is applied.
    ! Single valued version
!</description>

!<input>
    ! (de)nominator
    real(SP), intent(in) :: p,q

    ! default value
    real(SP), intent(in) :: fval

    ! upper bound
    real(SP), intent(in) :: fbound
!</input>

!<result>
    ! limited ratio
    real(SP) :: r
!</result>
!</function>

    if (abs(p) .gt. AFCSTAB_EPSABS) then
      r = min(q/p, fbound)
    else
      r = fval
    end if
  end function afcstab_limitBoundedSP

  !*****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine afcstab_combFluxesQP(NEDGE, qscale, Qflux1, Qflux2, Qalpha,&
                                  rperfconfig)

!<description>
    ! This subroutine combines the two fluxes:
    ! Qflux2 := Qflux2 + qscale * Qalpha * Qflux1
!</description>

!<input>
    ! First flux
    real(QP), dimension(:), intent(in) :: Qflux1

    ! Individual scaling factor for each entry of flux1
    real(QP), dimension(:), intent(in), optional :: Qalpha

    ! Global scaling factor for all entries of flux1
    real(QP), intent(in) :: qscale

    ! Number of entries/edges
    integer, intent(in) :: NEDGE

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second flux
    real(QP), dimension(:), intent(inout) :: Qflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    !$ if (present(rperfconfig)) then
    !$  p_rperfconfig => rperfconfig
    !$ else
    !$  p_rperfconfig => afcstab_perfconfig
    !$ end if

    if (present(Qalpha)) then

      if (qscale .eq. 1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(iedge) = Qflux2(iedge)&
                        + Qalpha(iedge) * Qflux1(iedge)
        end do
        !$omp end parallel do

      elseif (qscale .eq. -1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(iedge) = Qflux2(iedge)&
                        - Qalpha(iedge) * Qflux1(iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(iedge) = Qflux2(iedge)&
                        + qscale * Qalpha(iedge) * Qflux1(iedge)
        end do
        !$omp end parallel do

      end if

    else   ! Qalpha not present

      if (qscale .eq. 1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(iedge) = Qflux2(iedge) + Qflux1(iedge)
        end do
        !$omp end parallel do

      elseif (qscale .eq. -1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(iedge) = Qflux2(iedge) - Qflux1(iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(iedge) = Qflux2(iedge) + qscale * Qflux1(iedge)
        end do
        !$omp end parallel do

      end if

    end if

  end subroutine afcstab_combFluxesQP

  !*****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine afcstab_combFluxesDP(NEDGE, dscale, Dflux1, Dflux2, Dalpha,&
                                  rperfconfig)

!<description>
    ! This subroutine combines the two fluxes:
    ! Dflux2 := Dflux2 + dscale * Dalpha * Dflux1
!</description>

!<input>
    ! First flux
    real(DP), dimension(:), intent(in) :: Dflux1

    ! Individual scaling factor for each entry of flux1
    real(DP), dimension(:), intent(in), optional :: Dalpha

    ! Global scaling factor for all entries of flux1
    real(DP), intent(in) :: dscale

    ! Number of entries/edges
    integer, intent(in) :: NEDGE

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second flux
    real(DP), dimension(:), intent(inout) :: Dflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    !$ if (present(rperfconfig)) then
    !$  p_rperfconfig => rperfconfig
    !$ else
    !$  p_rperfconfig => afcstab_perfconfig
    !$ end if

    if (present(Dalpha)) then

      if (dscale .eq. 1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge)&
                        + Dalpha(iedge) * Dflux1(iedge)
        end do
        !$omp end parallel do

      elseif (dscale .eq. -1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge)&
                        - Dalpha(iedge) * Dflux1(iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge)&
                        + dscale * Dalpha(iedge) * Dflux1(iedge)
        end do
        !$omp end parallel do

      end if

    else   ! Dalpha not present

      if (dscale .eq. 1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge) + Dflux1(iedge)
        end do
        !$omp end parallel do

      elseif (dscale .eq. -1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge) - Dflux1(iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge) + dscale * Dflux1(iedge)
        end do
        !$omp end parallel do

      end if

    end if

  end subroutine afcstab_combFluxesDP

  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine afcstab_combFluxesSP(NEDGE, fscale, Fflux1, Fflux2, Falpha,&
                                  rperfconfig)

!<description>
    ! This subroutine combines the two fluxes:
    ! Fflux2 := Fflux2 + fscale * Falpha * Fflux1
!</description>

!<input>
    ! First flux
    real(SP), dimension(:), intent(in) :: Fflux1

    ! Individual scaling factor for each entry of flux1
    real(SP), dimension(:), intent(in), optional :: Falpha

    ! Global scaling factor for all entries of flux1
    real(SP), intent(in) :: fscale

    ! Number of entries/edges
    integer, intent(in) :: NEDGE

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second flux
    real(SP), dimension(:), intent(inout) :: Fflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    !$ if (present(rperfconfig)) then
    !$  p_rperfconfig => rperfconfig
    !$ else
    !$  p_rperfconfig => afcstab_perfconfig
    !$ end if

    if (present(Falpha)) then

      if (fscale .eq. 1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge)&
                        + Falpha(iedge) * Fflux1(iedge)
        end do
        !$omp end parallel do

      elseif (fscale .eq. -1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge)&
                        - Falpha(iedge) * Fflux1(iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge)&
                        + fscale * Falpha(iedge) * Fflux1(iedge)
        end do
        !$omp end parallel do

      end if

    else   ! Falpha not present

      if (fscale .eq. 1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge) + Fflux1(iedge)
        end do
        !$omp end parallel do

      elseif (fscale .eq. -1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge) - Fflux1(iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge) + fscale * Fflux1(iedge)
        end do
        !$omp end parallel do

      end if

    end if

  end subroutine afcstab_combFluxesSP

  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
    pure&
#endif
    subroutine afcstab_combineFluxesQP(NVAR, NEDGE, qscale, Qflux1, Qflux2,&
                                       Qalpha, rperfconfig)

!<description>
    ! This subroutine combines the two fluxes:
    ! Qflux2 := Qflux2 + qscale * Qalpha * Qflux1
!</description>

!<input>
    ! Number of entries/edges
    integer, intent(in) :: NEDGE

    ! Number of variables
    integer, intent(in) :: NVAR

    ! First flux
    real(QP), dimension(NVAR,NEDGE), intent(in) :: Qflux1

    ! Individual scaling factor for each entry of flux1
    real(QP), dimension(:), intent(in), optional :: Qalpha

    ! Global scaling factor for all entries of flux1
    real(QP), intent(in) :: qscale

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second flux
    real(QP), dimension(NVAR,NEDGE), intent(inout) :: Qflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    !$ if (present(rperfconfig)) then
    !$  p_rperfconfig => rperfconfig
    !$ else
    !$  p_rperfconfig => afcstab_perfconfig
    !$ end if

    if (present(Qalpha)) then

      if (qscale .eq. 1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(:,iedge) = Qflux2(:,iedge)&
                          + Qalpha(iedge) * Qflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (qscale .eq. -1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(:,iedge) = Qflux2(:,iedge)&
                          - Qalpha(iedge) * Qflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(:,iedge) = Qflux2(:,iedge)&
                          + qscale * Qalpha(iedge) * Qflux1(:,iedge)
        end do
        !$omp end parallel do

      end if

    else   ! Qalpha not present

      if (qscale .eq. 1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(:,iedge) = Qflux2(:,iedge)&
                          + Qflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (qscale .eq. -1.0_QP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(:,iedge) = Qflux2(:,iedge)&
                          - Qflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Qflux2(:,iedge) = Qflux2(:,iedge)&
                          + qscale * Qflux1(:,iedge)
        end do
        !$omp end parallel do

      end if

    end if

  end subroutine afcstab_combineFluxesQP

  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
    pure&
#endif
    subroutine afcstab_combineFluxesDP(NVAR, NEDGE, dscale, Dflux1, Dflux2,&
                                       Dalpha, rperfconfig)

!<description>
    ! This subroutine combines the two fluxes:
    ! Dflux2 := Dflux2 + dscale * Dalpha * Dflux1
!</description>

!<input>
    ! Number of entries/edges
    integer, intent(in) :: NEDGE

    ! Number of variables
    integer, intent(in) :: NVAR

    ! First flux
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1

    ! Individual scaling factor for each entry of flux1
    real(DP), dimension(:), intent(in), optional :: Dalpha

    ! Global scaling factor for all entries of flux1
    real(DP), intent(in) :: dscale

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second flux
    real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    !$ if (present(rperfconfig)) then
    !$  p_rperfconfig => rperfconfig
    !$ else
    !$  p_rperfconfig => afcstab_perfconfig
    !$ end if

    if (present(Dalpha)) then

      if (dscale .eq. 1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          + Dalpha(iedge) * Dflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (dscale .eq. -1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          - Dalpha(iedge) * Dflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          + dscale * Dalpha(iedge) * Dflux1(:,iedge)
        end do
        !$omp end parallel do

      end if

    else   ! Dalpha not present

      if (dscale .eq. 1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          + Dflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (dscale .eq. -1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          - Dflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          + dscale * Dflux1(:,iedge)
        end do
        !$omp end parallel do

      end if

    end if

  end subroutine afcstab_combineFluxesDP

  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
    pure&
#endif
    subroutine afcstab_combineFluxesSP(NVAR, NEDGE, fscale, Fflux1, Fflux2,&
                                       Falpha, rperfconfig)

!<description>
    ! This subroutine combines the two fluxes:
    ! Fflux2 := Fflux2 + fscale * Falpha * Fflux1
!</description>

!<input>
    ! First flux
    real(SP), dimension(NVAR,NEDGE), intent(in) :: Fflux1

    ! Individual scaling factor for each entry of flux1
    real(SP), dimension(:), intent(in), optional :: Falpha

    ! Global scaling factor for all entries of flux1
    real(SP), intent(in) :: fscale

    ! Number of entries/edges
    integer, intent(in) :: NEDGE

    ! Number of variables
    integer, intent(in) :: NVAR

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second flux
    real(SP), dimension(NVAR,NEDGE), intent(inout) :: Fflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    !$ if (present(rperfconfig)) then
    !$  p_rperfconfig => rperfconfig
    !$ else
    !$  p_rperfconfig => afcstab_perfconfig
    !$ end if

    if (present(Falpha)) then

      if (fscale .eq. 1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          + Falpha(iedge) * Fflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (fscale .eq. -1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          - Falpha(iedge) * Fflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          + fscale * Falpha(iedge) * Fflux1(:,iedge)
        end do
        !$omp end parallel do

      end if

    else   ! Falpha not present

      if (fscale .eq. 1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          + Fflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (fscale .eq. -1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          - Fflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          + fscale * Fflux1(:,iedge)
        end do
        !$omp end parallel do

      end if

    end if

  end subroutine afcstab_combineFluxesSP

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_upwindOrientationDP1(Dcoefficients, IedgeList,&
      DcoefficientsAux, ipos, jpos, rperfconfig)

!<description>
    ! This subroutine orients the edges so that the starting node of
    ! each edge is located upwind. This orientation convention is
    ! introduced in the paper:
    !
    ! D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
    ! conservation laws, In: D. Kuzmin et al. (eds), Flux-Corrected
    ! Transport: Principles, Algorithms, and Applications, Springer,
    ! 2005, 155-206.
!</description>

!<input>
    ! Positions of node I and J
    integer, intent(in) :: ipos,jpos

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Coefficients used to determine the upwind direction
    ! DIMENSION(1:IPOS:JPOS:NPOS,1:NEDGE),
    ! where IPOS and JPOS are the positions where the data
    ! for node I and J is stored, respectively
    real(DP), dimension(:,:), intent(inout) :: Dcoefficients

    ! Additional coefficients at the edges
    ! DIMENSION(1:N2,1:2,1:NEDGE)
    real(DP), dimension(:,:,:), intent(inout) :: DcoefficientsAux

    ! List of edges
    ! DIMENSION(1:N3,1:NEDGE)
    integer, dimension(:,:), intent(inout) :: IedgeList
!</inputoutput>
!</subroutine

    ! local variables
    real(DP) :: daux
    integer :: iedge,iaux,naux,nedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcstab_perfconfig
    end if

    nedge = size(IedgeList,2)
    naux  = size(DcoefficientsAux,1)

    select case(size(IedgeList,1))
    case (2)
      !$omp parallel do default(shared) private(iaux,daux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Dcoefficients(ipos,iedge) .gt. Dcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          daux = Dcoefficients(ipos,iedge)
          Dcoefficients(ipos,iedge) = Dcoefficients(jpos,iedge)
          Dcoefficients(jpos,iedge) = daux

          ! Swap edgewise coefficients Data_ij <-> Data_ji
          do iaux = 1, naux
            daux = DcoefficientsAux(iaux,1,iedge)
            DcoefficientsAux(iaux,1,iedge) = DcoefficientsAux(iaux,2,iedge)
            DcoefficientsAux(iaux,2,iedge) = daux
          end do
        end if
      end do
      !$omp end parallel do

    case (6)
      !$omp parallel do default(shared) private(iaux,daux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Dcoefficients(ipos,iedge) .gt. Dcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edges ij <-> ji
          iaux = IedgeList(3,iedge)
          IedgeList(3,iedge) = IedgeList(4,iedge)
          IedgeList(4,iedge) = iaux

          ! Swap edges ii <-> jj
          iaux = IedgeList(5,iedge)
          IedgeList(5,iedge) = IedgeList(6,iedge)
          IedgeList(6,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          daux = Dcoefficients(ipos,iedge)
          Dcoefficients(ipos,iedge) = Dcoefficients(jpos,iedge)
          Dcoefficients(jpos,iedge) = daux

          ! Swap edgewise coefficients Data_ij <-> Data_ji
          do iaux = 1, naux
            daux = DcoefficientsAux(iaux,1,iedge)
            DcoefficientsAux(iaux,1,iedge) = DcoefficientsAux(iaux,2,iedge)
            DcoefficientsAux(iaux,2,iedge) = daux
          end do
        end if
      end do
      !$omp end parallel do

    case default
      call output_line('SIZE(IedgeList,2) is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_upwindOrientationDP1')
      call sys_halt()
    end select

  end subroutine afcstab_upwindOrientationDP1

!*****************************************************************************

!<subroutine>

  subroutine afcstab_upwindOrientationDP2(Dcoefficients, IedgeList,&
      ipos, jpos, rperfconfig)

!<description>
    ! This subroutine orients the edges so that the starting node of
    ! each edge is located upwind. This orientation convention is
    ! introduced in the paper:
    !
    ! D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
    ! conservation laws, In: D. Kuzmin et al. (eds), Flux-Corrected
    ! Transport: Principles, Algorithms, and Applications, Springer,
    ! 2005, 155-206.
!</description>

!<input>
    ! Positions of node I and J
    integer, intent(in) :: ipos,jpos

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Coefficients used to determine the upwind direction
    ! DIMENSION(1:IPOS:JPOS:NPOS,1:NEDGE),
    ! where IPOS and JPOS are the positions where the data
    ! for node I and J is stored, respectively
    real(DP), dimension(:,:), intent(inout) :: Dcoefficients

    ! List of edges
    ! DIMENSION(1:N3,1:NEDGE)
    integer, dimension(:,:), intent(inout) :: IedgeList
!</inputoutput>
!</subroutine

    ! local variables
    real(DP) :: daux
    integer :: iedge,iaux,nedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcstab_perfconfig
    end if

    nedge = size(IedgeList,2)

    select case(size(IedgeList,1))
    case (2)
      !$omp parallel do default(shared) private(iaux,daux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Dcoefficients(ipos,iedge) .gt. Dcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          daux = Dcoefficients(ipos,iedge)
          Dcoefficients(ipos,iedge) = Dcoefficients(jpos,iedge)
          Dcoefficients(jpos,iedge) = daux
        end if
      end do
      !$omp end parallel do

    case (6)
      !$omp parallel do default(shared) private(iaux,daux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Dcoefficients(ipos,iedge) .gt. Dcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edges ij <-> ji
          iaux = IedgeList(3,iedge)
          IedgeList(3,iedge) = IedgeList(4,iedge)
          IedgeList(4,iedge) = iaux

          ! Swap edges ii <-> jj
          iaux = IedgeList(5,iedge)
          IedgeList(5,iedge) = IedgeList(6,iedge)
          IedgeList(6,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          daux = Dcoefficients(ipos,iedge)
          Dcoefficients(ipos,iedge) = Dcoefficients(jpos,iedge)
          Dcoefficients(jpos,iedge) = daux
        end if
      end do
      !$omp end parallel do

    case default
      call output_line('SIZE(IedgeList,2) is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_upwindOrientationDP2')
      call sys_halt()
    end select

  end subroutine afcstab_upwindOrientationDP2

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_upwindOrientationSP1(Fcoefficients, IedgeList,&
      FcoefficientsAux, ipos, jpos, rperfconfig)

!<description>
    ! This subroutine orients the edges so that the starting node of
    ! each edge is located upwind. This orientation convention is
    ! introduced in the paper:
    !
    ! D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
    ! conservation laws, In: D. Kuzmin et al. (eds), Flux-Corrected
    ! Transport: Principles, Algorithms, and Applications, Springer,
    ! 2005, 155-206.
!</description>

!<input>
    ! Positions of node I and J
    integer, intent(in) :: ipos,jpos

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Coefficients used to determine the upwind direction
    ! DIMENSION(1:IPOS:JPOS:NPOS,1:NEDGE),
    ! where IPOS and JPOS are the positions where the data
    ! for node I and J is stored, respectively
    real(SP), dimension(:,:), intent(inout) :: Fcoefficients

    ! Additional coefficients at the edges
    ! DIMENSION(1:N2,1:2,1:NEDGE)
    real(SP), dimension(:,:,:), intent(inout) :: FcoefficientsAux

    ! List of edges
    ! DIMENSION(1:N3,1:NEDGE)
    integer, dimension(:,:), intent(inout) :: IedgeList
!</inputoutput>
!</subroutine

    ! local variables
    real(SP) :: faux
    integer :: iedge,iaux,naux,nedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcstab_perfconfig
    end if

    nedge = size(IedgeList,2)
    naux  = size(FcoefficientsAux,1)

    select case(size(IedgeList,1))
    case (2)
      !$omp parallel do default(shared) private(iaux,faux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Fcoefficients(ipos,iedge) .gt. Fcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          faux = Fcoefficients(ipos,iedge)
          Fcoefficients(ipos,iedge) = Fcoefficients(jpos,iedge)
          Fcoefficients(jpos,iedge) = faux

          ! Swap edgewise coefficients Data_ij <-> Data_ji
          do iaux = 1, naux
            faux = FcoefficientsAux(iaux,1,iedge)
            FcoefficientsAux(iaux,1,iedge) = FcoefficientsAux(iaux,2,iedge)
            FcoefficientsAux(iaux,2,iedge) = faux
          end do
        end if
      end do
      !$omp end parallel do

    case (6)
      !$omp parallel do default(shared) private(iaux,faux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Fcoefficients(ipos,iedge) .gt. Fcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edges ij <-> ji
          iaux = IedgeList(3,iedge)
          IedgeList(3,iedge) = IedgeList(4,iedge)
          IedgeList(4,iedge) = iaux

          ! Swap edges ii <-> jj
          iaux = IedgeList(5,iedge)
          IedgeList(5,iedge) = IedgeList(6,iedge)
          IedgeList(6,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          faux = Fcoefficients(ipos,iedge)
          Fcoefficients(ipos,iedge) = Fcoefficients(jpos,iedge)
          Fcoefficients(jpos,iedge) = faux

          ! Swap edgewise coefficients Data_ij <-> Data_ji
          do iaux = 1, naux
            faux = FcoefficientsAux(iaux,1,iedge)
            FcoefficientsAux(iaux,1,iedge) = FcoefficientsAux(iaux,2,iedge)
            FcoefficientsAux(iaux,2,iedge) = faux
          end do
        end if
      end do
      !$omp end parallel do

    case default
      call output_line('SIZE(IedgeList,2) is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_upwindOrientationSP1')
      call sys_halt()
    end select

  end subroutine afcstab_upwindOrientationSP1

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_upwindOrientationSP2(Fcoefficients, IedgeList,&
      ipos, jpos, rperfconfig)

!<description>
    ! This subroutine orients the edges so that the starting node of
    ! each edge is located upwind. This orientation convention is
    ! introduced in the paper:
    !
    ! D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
    ! conservation laws, In: D. Kuzmin et al. (eds), Flux-Corrected
    ! Transport: Principles, Algorithms, and Applications, Springer,
    ! 2005, 155-206.
!</description>

!<input>
    ! Positions of node I and J
    integer, intent(in) :: ipos,jpos

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Coefficients used to determine the upwind direction
    ! DIMENSION(1:IPOS:JPOS:NPOS,1:NEDGE),
    ! where IPOS and JPOS are the positions where the data
    ! for node I and J is stored, respectively
    real(SP), dimension(:,:), intent(inout) :: Fcoefficients

    ! List of edges
    ! DIMENSION(1:N3,1:NEDGE)
    integer, dimension(:,:), intent(inout) :: IedgeList
!</inputoutput>
!</subroutine

    ! local variables
    real(SP) :: faux
    integer :: iedge,iaux,nedge

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcstab_perfconfig
    end if

    nedge = size(IedgeList,2)

    select case(size(IedgeList,1))
    case (2)
      !$omp parallel do default(shared) private(iaux,faux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Fcoefficients(ipos,iedge) .gt. Fcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          faux = Fcoefficients(ipos,iedge)
          Fcoefficients(ipos,iedge) = Fcoefficients(jpos,iedge)
          Fcoefficients(jpos,iedge) = faux
        end if
      end do
      !$omp end parallel do

    case (6)
      !$omp parallel do default(shared) private(iaux,faux)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, nedge
        if (Fcoefficients(ipos,iedge) .gt. Fcoefficients(jpos,iedge)) then
          ! Swap nodes i <-> j
          iaux = IedgeList(1,iedge)
          IedgeList(1,iedge) = IedgeList(2,iedge)
          IedgeList(2,iedge) = iaux

          ! Swap edges ij <-> ji
          iaux = IedgeList(3,iedge)
          IedgeList(3,iedge) = IedgeList(4,iedge)
          IedgeList(4,iedge) = iaux

          ! Swap edges ii <-> jj
          iaux = IedgeList(5,iedge)
          IedgeList(5,iedge) = IedgeList(6,iedge)
          IedgeList(6,iedge) = iaux

          ! Swap edge data data_ij <-> data_ji
          faux = Fcoefficients(ipos,iedge)
          Fcoefficients(ipos,iedge) = Fcoefficients(jpos,iedge)
          Fcoefficients(jpos,iedge) = faux
        end if
      end do
      !$omp end parallel do

    case default
      call output_line('SIZE(IedgeList,2) is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcstab_upwindOrientationSP2')
      call sys_halt()
    end select

  end subroutine afcstab_upwindOrientationSP2

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyH2D_IedgeListIdx(rafcstab, btranspose, istream)

!<description>
    ! This subroutine copies the edge index structure from the
    ! host memory to the memory of the coprocessor device. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! OPTIONAL: if true then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>
!</subroutine>

    if (rafcstab%h_IedgeListIdx .ne. ST_NOHANDLE)&
        call storage_syncMemoryHostDevice(rafcstab%h_IedgeListIdx,&
        ST_SYNCBLOCK_COPY_H2D, btranspose, istream)

  end subroutine afcstab_copyH2D_IedgeListIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyD2H_IedgeListIdx(rafcstab, btranspose, istream)

!<description>
    ! This subroutine copies the edge index structure from the
    ! memory of the coprocessor device to the host memory. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! OPTIONAL: if true then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>
!</subroutine>

    if (rafcstab%h_IedgeListIdx .ne. ST_NOHANDLE)&
        call storage_syncMemoryHostDevice(rafcstab%h_IedgeListIdx,&
        ST_SYNCBLOCK_COPY_D2H, btranspose, istream)

  end subroutine afcstab_copyD2H_IedgeListIdx

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyH2D_IedgeList(rafcstab, btranspose, istream)

!<description>
    ! This subroutine copies the edge structure from the
    ! host memory to the memory of the coprocessor device. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! OPTIONAL: if true then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>
!</subroutine>

    if (rafcstab%h_IedgeList .ne. ST_NOHANDLE)&
        call storage_syncMemoryHostDevice(rafcstab%h_IedgeList,&
        ST_SYNCBLOCK_COPY_H2D, btranspose, istream)

  end subroutine afcstab_copyH2D_IedgeList

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyD2H_IedgeList(rafcstab, btranspose, istream)

!<description>
    ! This subroutine copies the edge structure from the
    ! memory of the coprocessor device to the host memory. If no
    ! device is available, then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! OPTIONAL: if true then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>
!</subroutine>

    if (rafcstab%h_IedgeList .ne. ST_NOHANDLE)&
        call storage_syncMemoryHostDevice(rafcstab%h_IedgeList,&
        ST_SYNCBLOCK_COPY_D2H, btranspose, istream)

  end subroutine afcstab_copyD2H_IedgeList

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyH2D_CoeffsAtEdge(rafcstab, btranspose, istream)

!<description>
    ! This subroutine copies the edge data from the host memory to the
    ! memory of the coprocessor device. If no device is available,
    ! then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! OPTIONAL: if true then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>
!</subroutine>

    if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemoryHostDevice(rafcstab%h_CoeffsAtEdge,&
        ST_SYNCBLOCK_COPY_H2D, btranspose, istream)

  end subroutine afcstab_copyH2D_CoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_copyD2H_CoeffsAtEdge(rafcstab, btranspose, istream)

!<description>
    ! This subroutine copies the edge data from the memory of the
    ! coprocessor device to the host memory. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! OPTIONAL: if true then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>
!</subroutine>

    if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_syncMemoryHostDevice(rafcstab%h_CoeffsAtEdge,&
        ST_SYNCBLOCK_COPY_D2H, btranspose, istream)

  end subroutine afcstab_copyD2H_CoeffsAtEdge

  !*****************************************************************************

!<subroutine>

  subroutine afcstab_setDataType(rafcstab, cdataType)

!<description>
    ! This subroutine converts the data arrays of the stabilisation
    ! structure to the given data type cdataType
!</description>

!<input>
    ! Target datatype of the stabilisation structure
    integer, intent(in) :: cdataType
!</input>

!<inputoutput>
    ! Stabilisation structure that should be converts
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if stabilisation structure needs conversion
    if (rafcstab%cdataType .eq. cdataType) return

    ! Set data type
    rafcstab%cdataType = cdataType
    
    ! Convert all available data
    if (rafcstab%h_CoeffsAtEdge .ne. ST_NOHANDLE)&
        call storage_setdatatype (rafcstab%h_CoeffsAtEdge, cdataType)
    if (rafcstab%h_BoundsAtEdge .ne. ST_NOHANDLE)&
        call storage_setdatatype (rafcstab%h_BoundsAtEdge, cdataType)
    if (associated(rafcstab%p_rvectorAlpha))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorAlpha, cdataType)
    if (associated(rafcstab%p_rvectorFlux0))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorFlux0, cdataType)
    if (associated(rafcstab%p_rvectorFlux))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorFlux, cdataType)
    if (associated(rafcstab%p_rvectorFluxPrel))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorFluxPrel, cdataType)
    if (associated(rafcstab%p_rvectorDx))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorDx, cdataType)
    if (associated(rafcstab%p_rvectorPp))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorPp, cdataType)
    if (associated(rafcstab%p_rvectorPm))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorPm, cdataType)
    if (associated(rafcstab%p_rvectorQ))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorQ, cdataType)
    if (associated(rafcstab%p_rvectorQp))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorQp, cdataType)
    if (associated(rafcstab%p_rvectorQm))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorQm, cdataType)
    if (associated(rafcstab%p_rvectorRp))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorRp, cdataType)
    if (associated(rafcstab%p_rvectorRm))&
        call lsyssc_setDataTypeVector (rafcstab%p_rvectorRm, cdataType)
    if (associated(rafcstab%p_rvectorPredictor))&
        call lsysbl_setDataTypeVector (rafcstab%p_rvectorPredictor, cdataType)

  end subroutine afcstab_setDataType

end module afcstabbase

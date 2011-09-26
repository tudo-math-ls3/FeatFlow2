!##############################################################################
!# ****************************************************************************
!# <name> groupfemscalar </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the
!# group-finite element formulation to scalar problems,
!# i.e. conservation laws.  The technique was proposed by
!# C.A.J. Fletcher in:
!#
!#     C.A.J. Fletcher, The group finite element formulation
!#     Computer Methods in Applied Mechanics and Engineering (ISSN
!#     0045-7825), vol. 37, April 1983, p. 225-244.
!#
!# The group finite element formulation uses the same basis functions
!# for the unknown solution and the fluxes. This allows for an
!# efficient matrix assemble, whereby the constant coefficient
!# matrices can be assembled once and for all at the beginning of the
!# simulation and each time the grid is modified.
!#
!# Moreover, this module allows to modifying discrete operators by
!# means of the algebraic flux correction (AFC) methodology proposed
!# by Kuzmin, Moeller and Turek in a series of publications. As a
!# starting point for scalar conservation laws, the reader is referred
!# to the book chapterdiffu
!#
!#     D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
!#     conservation laws, In: D. Kuzmin et al. (eds), Flux-Corrected
!#     Transport: Principles, Algorithms, and Applications, Springer,
!#     2005, 155-206.
!#
!# A more detailed description of the algorithms is given in the
!# comments of the subroutine implementing the corresponding
!# discretisation schemes. All methods are based on the stabilisation
!# structure t_afcstab which is defined in the underlying module
!# afcstabilisation. The initialisation as a scalar stabilisation
!# structure is done by the routine gfsc_initStabilisation.
!# 
!# There are three types of routines. The gfsc_buildOperator routine
!# can be used to assemble discrete diffusion operators resulting from
!# the standard Galerkin finite element discretisation plus some
!# discretely defined artificial diffusion. For convective terms, this
!# technique is termed discrete upwinding, whereas for physical
!# diffusion operators this approach ensures that the discrete maximum
!# principle holds. Consequently, the term DMP is adopted.
!#
!# The second type of routines is given by gfsc_buildVectorXXX. They
!# can be used to update/initialise the rhs/residual applying some
!# sort of algebraic flux correction. Importantly, the family of AFC
!# schemes gives rise to nonlinear algebraic equations that need to be
!# solved iteratively. Thus, it is usefull to build the compensating
!# antidiffusion into the residual vector rather than the right hand
!# side. However, it is still possible to give a negative scaling
!# factor.
!#
!# The third type of routines is used to assemble the Jacobian matrix
!# for Newton`s method. Here, the exact Jacobian matrix is
!# approximated by means of second-order divided differences whereby
!# the perturbation parameter is specified by the user. You should
!# be aware of the fact, that in general the employed flux limiters
!# are not differentiable globally so that the construction of the
!# Jacobian matrix is somehow delicate. Even though the routines
!# will produce some matrix without warnings, this matrix may be
!# singular and/or ill-conditioned.
!#
!# The following routines are available:
!#
!# 1.) gfsc_initStabilisation
!#     -> Initializes the stabilisation structure
!#
!# 2.) gfsc_buildOperator = gfsc_buildOperatorConst /
!#                          gfsc_buildOperatorNodeScalar /
!#                          gfsc_buildOperatorNodeBlock /
!#                          gfsc_buildOperatorEdgeScalar /
!#                          gfsc_buildOperatorEdgeBlock
!#     -> Assembles a discrete operator by the group finite element formulation
!#
!# 4.) gfsc_buildFluxFCT = gfsc_buildFluxFCTScalar /
!#                         gfsc_buildFluxFCTBlock
!#     -> Assembles the raw antidiffusive flux for AFC stabilisation of FCT type
!#
!# 5.) gfsc_buildConvectionVectorFCT = gfsc_buildConvVecFCTScalar /
!#                                     gfsc_buildConvVecFCTBlock /
!#     -> Assembles the convective vector for AFC stabilisation of FCT type
!#
!# 6.) gfsc_buildConvectionVectorTVD = gfsc_buildConvVecTVDScalar /
!#                                     gfsc_buildConvVecTVDBlock
!#     -> Assembles the convective vector for AFC stabilisation of TVD type
!#
!# 7.) gfsc_buildConvectionVectorGP = gfsc_buildConvVecGPScalar /
!#                                    gfsc_buildConvVecGPBlock
!#     -> Assembles the convective vector for AFC stabilisation of general-purpose type
!#
!# 8.) gfsc_buildConvectionVectorSymm = gfsc_buildConvVecSymmScalar /
!#                                      gfsc_buildConvVecSymmBlock
!#     -> Assembles the convective vector for stabilisation by means of
!#        symmetric flux limiting for diffusion operators
!#
!# 9.) gfsc_buildFluxLPT = gfsc_buildFluxLPTScalar /
!#                         gfsc_buildFluxLPTBlock
!#     -> Assembles the raw antidiffusive flux for the
!#        linearity-preserving stabilisation
!#
!# 10.) gfsc_buildVectorLPT = gfsc_buildVecLPTScalar /
!#                            gfsc_buildVecLPTBlock
!#      -> Assembles the vector for stabilisation by means of linearity-
!#         preserving flux correction for mass matrices
!#
!# 11.) gfsc_buildConvectionJacobian = gfsc_buildConvJacobianScalar /
!#                                     gfsc_buildConvJacobianBlock
!#      -> Assembles the Jacobian matrix for the convective part of
!#         the transport operator for a scalar convection equation
!#
!# 12.) gfsc_buildJacobianFCT = gfsc_buildJacLinearFCTScalar /
!#                              gfsc_buildJacLinearFCTBlock /
!#                              gfsc_buildJacobianFCTScalar /
!#                              gfsc_buildJacobianFCTBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of FCT type;
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 13.) gfsc_buildJacobianTVD = gfsc_buildJacLinearTVDScalar /
!#                              gfsc_buildJacLinearTVDBlock /
!#                              gfsc_buildJacobianTVDScalar /
!#                              gfsc_buildJacobianTVDBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of TVD type;
!#         For the first two routines, the velocity is assumed to be linear which
!#         simplifies the evaluation of the Jacobian matrix significantly.
!#         For the second two routines, the velocity can be arbitrary.
!#
!# 14.) gfsc_buildJacobianGP = gfsc_buildJacLinearGPScalar /
!#                             gfsc_buildJacLinearGPBlock /
!#                             gfsc_buildJacobianGPScalar /
!#                             gfsc_buildJacobianGPBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of general 
!#         purpose limiter. For the first two routines, the velocity is assumed
!#         to be linear which simplifies the evaluation of the Jacobian matrix
!#         significantly. For the second two routines, the velocity can be arbitrary.
!#
!# 15.) gfsc_buildJacobianSymm = gfsc_buildJacobianSymmScalar /
!#                               gfsc_buildJacobianSymmBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of symmetric
!#         flux limiting for diffusion operators
!#
!#
!# The following auxiliary routines are available:
!#
!#  1.) gfsc_combineFluxes = gfsc_combineFluxesDble /
!#                           gfsc_combineFluxesSngl
!#      -> Linear combination of fluxes
!#
!#  2.) gfsc_limit = gfsc_limitUnboundedDble /
!#                   gfsc_limitUnboundedSngl /
!#                   gfsc_limitBoundedDble /
!#                   gfsc_limitBoundedSngl
!#      -> Compute the nodal correction factors, i.e., the ratio of
!#         admissible solution increments and raw antidiffusion
!#
!# 3.) gfsc_upwindOrientation = gfsc_upwindOrientationDble /
!#                              gfsc_upwindOrientationSngl
!#     -> Swap edge orientation so that the starting edge is located upwind
!#
!# </purpose>
!##############################################################################

module groupfemscalar

  use afcstabilisation
  use basicgeometry
  use collection
  use fsystem
  use genoutput
  use groupfembase
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use mprimitives
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none

  private
  
  public :: gfsc_initStabilisation
  public :: gfsc_buildOperator
  public :: gfsc_buildFluxFCT
  public :: gfsc_buildFluxLPT
  public :: gfsc_buildConvectionVectorFCT
  public :: gfsc_buildConvectionVectorTVD
  public :: gfsc_buildConvectionVectorGP
  public :: gfsc_buildConvectionVectorSymm
  public :: gfsc_buildConvectionJacobian
  public :: gfsc_buildVectorLPT
  public :: gfsc_buildJacobianFCT
  public :: gfsc_buildJacobianTVD
  public :: gfsc_buildJacobianGP
  public :: gfsc_buildJacobianSymm
 
!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of nodes to handle simultaneously when building matrices
#ifndef GFSC_NEQSIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEQSIM = 128
#else
  integer, public            :: GFSC_NEQSIM = 128
#endif
#endif

  ! Number of edges to handle simultaneously when building matrices
#ifndef GFSC_NEDGESIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEDGESIM = 64
#else
  integer, public            :: GFSC_NEDGESIM = 64
#endif
#endif
  
  ! Minimum number of nodes for OpenMP parallelisation: If the number of
  ! nodes is below this value, then no parallelisation is performed.
#ifndef GFSC_NEQMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEQMIN_OMP = 1000
#else
  integer, public            :: GFSC_NEQMIN_OMP = 1000
#endif
#endif

  ! Minimum number of edges for OpenMP parallelisation: If the number of
  ! edges is below this value, then no parallelisation is performed.
#ifndef GFSC_NEDGEMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSC_NEDGEMIN_OMP = 1000
#else
  integer, public            :: GFSC_NEDGEMIN_OMP = 1000
#endif
#endif
!</constantblock>

!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface gfsc_buildOperator
    module procedure gfsc_buildOperatorConst
    module procedure gfsc_buildOperatorNodeScalar
    module procedure gfsc_buildOperatorNodeBlock
    module procedure gfsc_buildOperatorEdgeScalar
    module procedure gfsc_buildOperatorEdgeBlock
  end interface

  interface gfsc_buildFluxFCT
    module procedure gfsc_buildFluxFCTScalar
    module procedure gfsc_buildFluxFCTBlock
  end interface

  interface gfsc_buildFluxLPT
    module procedure gfsc_buildFluxLPTScalar
    module procedure gfsc_buildFluxLPTBlock
  end interface

  interface gfsc_buildConvectionVectorFCT
    module procedure gfsc_buildConvVecFCTScalar
    module procedure gfsc_buildConvVecFCTBlock
  end interface

  interface gfsc_buildConvectionVectorTVD
    module procedure gfsc_buildConvVecTVDScalar
    module procedure gfsc_buildConvVecTVDBlock
  end interface

  interface gfsc_buildConvectionVectorGP
    module procedure gfsc_buildConvVecGPScalar
    module procedure gfsc_buildConvVecGPBlock
  end interface

  interface gfsc_buildConvectionVectorSymm
    module procedure gfsc_buildConvVecSymmScalar
    module procedure gfsc_buildConvVecSymmBlock
  end interface

  interface gfsc_buildConvectionJacobian
    module procedure gfsc_buildConvJacobianScalar
    module procedure gfsc_buildConvJacobianBlock
  end interface

  interface gfsc_buildVectorLPT
    module procedure gfsc_buildVecLPTScalar
    module procedure gfsc_buildVecLPTBlock
  end interface

  interface gfsc_buildJacobianFCT
    module procedure gfsc_buildJacLinearFCTScalar
    module procedure gfsc_buildJacLinearFCTBlock
    module procedure gfsc_buildJacobianFCTScalar
    module procedure gfsc_buildJacobianFCTBlock
  end interface

  interface gfsc_buildJacobianTVD
    module procedure gfsc_buildJacLinearTVDScalar
    module procedure gfsc_buildJacLinearTVDBlock
    module procedure gfsc_buildJacobianTVDScalar
    module procedure gfsc_buildJacobianTVDBlock
  end interface

  interface gfsc_buildJacobianGP
    module procedure gfsc_buildJacLinearGPScalar
    module procedure gfsc_buildJacLinearGPBlock
    module procedure gfsc_buildJacobianGPScalar
    module procedure gfsc_buildJacobianGPBlock
  end interface

  interface gfsc_buildJacobianSymm
    module procedure gfsc_buildJacobianSymmScalar
    module procedure gfsc_buildJacobianSymmBlock
  end interface

  interface gfsc_combineFluxes
    module procedure gfsc_combineFluxesDble
    module procedure gfsc_combineFluxesSngl
  end interface

  interface gfsc_limit
    module procedure gfsc_limitUnboundedDble
    module procedure gfsc_limitUnboundedSngl
    module procedure gfsc_limitBoundedDble
    module procedure gfsc_limitBoundedSngl
  end interface

  interface gfsc_upwindOrientation
    module procedure gfsc_upwindOrientationDble
    module procedure gfsc_upwindOrientationSngl
  end interface

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

contains

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_initStabilisation(rmatrixTemplate, rafcstab, rblockDiscretisation)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    ! Note that there exists an edge between two nodes i and j of their
    ! basis functions have overlapping support. For linear finite elements
    ! the number of physical edges coincides with the number of edges used
    ! in the stabilisation structure. In general the number of edges equals
    ! (NA-NEQ)/2, that is, it is the number of nonzero matrix entries not
    ! counting the diagonal divided by two.
!</description>

!<input>
    ! The template matrix
    type(t_matrixScalar), intent(in) :: rmatrixTemplate

    ! OPTIONAL: block discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_blockDiscretisation), intent(in), optional :: rblockDiscretisation
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    
    ! Set atomic data
    rafcstab%NVARtransformed = rmatrixTemplate%NVAR
    rafcstab%NVAR            = rmatrixTemplate%NVAR
    rafcstab%NEQ             = rmatrixTemplate%NEQ
    rafcstab%NEDGE           = (rmatrixTemplate%NA-rmatrixTemplate%NEQ)/2
    rafcstab%NNVEDGE         = 0
    
    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)
      
    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      !-------------------------------------------------------------------------
      
    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE)
      
      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)
      
      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      call afcstab_allocCoeffsAtEdge(rafcstab,3)

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
        call lsysbl_createVectorBlock(rblockDiscretisation,&
            rafcstab%p_rvectorPredictor, .false., rafcstab%cdataType)
      else
        call lsysbl_createVectorBlock(rafcstab%p_rvectorPredictor,&
            rafcstab%NEQ, 1, .false., rafcstab%cdataType)
      end if

      !-------------------------------------------------------------------------

    case (AFCSTAB_TVD,&
          AFCSTAB_GP)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      call afcstab_allocCoeffsAtEdge(rafcstab,3)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)
      
      ! We need the 3 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux0(rafcstab)
      call afcstab_allocFlux(rafcstab)
      
      !-------------------------------------------------------------------------

    case (AFCSTAB_LINFCT,&
          AFCSTAB_LINFCT_MASS)
      
      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      call afcstab_allocCoeffsAtEdge(rafcstab,3)

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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! Handle for DcoefficientsAtEdge: (/d_ij,s_ij/)
      call afcstab_allocCoeffsAtEdge(rafcstab,2)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the edgewise vector for the fluxes
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_MASS,&
          AFCSTAB_LINLPT_MASS)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      ! and one extra nodal vector Q
      call afcstab_allocVectorsPQR(rafcstab, ballocCommonQ=.true.)
    
      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)

      ! We need the nodal vector for the predictor
      allocate(rafcstab%p_rvectorPredictor)
      if (present(rblockDiscretisation)) then
        call lsysbl_createVectorBlock(rblockDiscretisation,&
            rafcstab%p_rvectorPredictor, .false., rafcstab%cdataType)
      else
        call lsysbl_createVectorBlock(rafcstab%p_rvectorPredictor,&
            rafcstab%NEQ, 1, .false., rafcstab%cdataType)
      end if

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_UPWINDBIASED,&
          AFCSTAB_LINLPT_UPWINDBIASED)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! Handle for DcoefficientsAtEdge: (/d_ij,k_ij,k_ji/)
      call afcstab_allocCoeffsAtEdge(rafcstab,3)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      ! and one extra nodal vector Q
      call afcstab_allocVectorsPQR(rafcstab, ballocCommonQ=.true.)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)
      
      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_SYMMETRIC,&
          AFCSTAB_LINLPT_SYMMETRIC)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! Handle for DcoefficientsAtEdge: (/d_ij,s_ij/)
      call afcstab_allocCoeffsAtEdge(rafcstab,2)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      ! and one extra nodal vector Q
      call afcstab_allocVectorsPQR(rafcstab, ballocCommonQ=.true.)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_initStabilisation')
      call sys_halt()
    end select
    

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

  end subroutine gfsc_initStabilisation
 
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorConst(rgroupFEMSet, dscale, bclear, rmatrix, rafcstab)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries are taken from
    ! the precomputed values stored in the structure rgroupFEMSet.
    !
    ! This is the simplest of the gfsc_buildOperatorX routines. It can
    ! handle both node-by-node and edge-by-edge matrix assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij=k_ji/)
    !
    ! Otherwise, the standard operator is assembled without stabilisation.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Scaling factor
    real(DP), intent(in) :: dscale
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge

    real(SP), dimension(:), pointer :: p_Fdata
    real(SP), dimension(:,:), pointer :: p_Fcoefficients
    real(SP), dimension(:,:), pointer :: p_FcoeffsAtNode
    real(SP), dimension(:,:,:), pointer :: p_FcoeffsAtEdge

    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer, dimension(:), pointer :: p_Kdiagonal,p_Kld

    logical :: bsymm

    ! Check if matrix and vector have the same data type
    if (rmatrix%cdataType .ne. rgroupFEMSet%cdataType) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Node-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end if

      if (present(rafcstab)) then
        call output_line('Stabilisation is not feasible in node-by-node assembly!',&
            OU_CLASS_WARNING,OU_MODE_STD,'gfsc_buildOperatorConst')
      end if
      
      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)

        ! What data types are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call lsyssc_getbase_double(rmatrix, p_Ddata)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpNodeMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpNodeMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata, p_InodeList)
          end if

        case (ST_SINGLE)
          ! Set pointers
          call gfem_getbase_FcoeffsAtNode(rgroupFEMSet, p_FcoeffsAtNode)
          call lsyssc_getbase_single(rmatrix, p_Fdata)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpNodeMat79Sngl(p_Kld, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpNodeMat79Sngl(p_Kld, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata, p_InodeList)
          end if

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end select

      
    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGESTRUCTURE) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA)      .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA)      .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end if
      
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
      
      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX7) then
          call lsyssc_getbase_Kld(rmatrix, p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        end if
        
        ! What data types are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
          call lsyssc_getbase_double(rmatrix, p_Ddata)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata, p_InodeList)
          end if
          

          ! Do we have to build the stabilisation?
          if (present(rafcstab)) then
            
            ! Check if stabilisation has been prepared
            if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
              call output_line('Stabilisation has not been prepared!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
              call sys_halt()
            end if

            ! Symmetric artificial diffusion?
            bsymm = .not.(rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED)
            
            ! Check if coefficients should be stored in stabilisation
            if (rafcstab%h_CoefficientsAtEdge .ne. ST_NOHANDLE) then
              
              ! Check if stabilisation has the same data type
              if (rafcstab%cdataType .ne. ST_DOUBLE) then
                call output_line('Stabilisation must have double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
                call sys_halt()
              end if
              
              ! Set additional pointers
              call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation and generate coefficients
              !-----------------------------------------------------------------
              call doOpEdgeAFCMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, dscale, bclear, bsymm,&
                  p_Ddata, p_Dcoefficients)
              
              ! Set state of stabilisation
              rafcstab%istabilisationSpec =&
                  ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)
              
              ! Do we need edge orientation?
              if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
                call gfsc_upwindOrientationDble(p_Dcoefficients, p_IedgeList,&
                    p_DcoeffsAtEdge, 2, 3)
                rafcstab%istabilisationSpec =&
                    ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
              else
                rafcstab%istabilisationSpec =&
                    iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
              end if
              
            else
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation but do not generate coeffs
              !-----------------------------------------------------------------
              call doOpEdgeStabMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, dscale, bclear, bsymm, p_Ddata)
            end if

          else   ! no stabilisation structure present

            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOpEdgeMat79Dble(p_IedgeListIdx, p_IedgeList,&
                rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, dscale, bclear, p_Ddata)
          end if

        case (ST_SINGLE)
          ! Set pointers
          call gfem_getbase_FcoeffsAtNode(rgroupFEMSet, p_FcoeffsAtNode)
          call gfem_getbase_FcoeffsAtEdge(rgroupFEMSet, p_FcoeffsAtEdge)
          call lsyssc_getbase_single(rmatrix, p_Fdata)
          
          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpDiagMat79Sngl(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpDiagMat79Sngl(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_FcoeffsAtNode, real(dscale,SP), bclear, p_Fdata, p_InodeList)
          end if

          ! Do we have to build the stabilisation?
          if (present(rafcstab)) then
            
            ! Check if stabilisation has been prepared
            if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
              call output_line('Stabilisation has not been prepared!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
              call sys_halt()
            end if

            ! Symmetric artificial diffusion?
            bsymm = .not.(rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED)
            
            ! Check if coefficients should be stored in stabilisation
            if (rafcstab%h_CoefficientsAtEdge .ne. ST_NOHANDLE) then
              
              ! Check if stabilisation has the same data type
              if (rafcstab%cdataType .ne. ST_SINGLE) then
                call output_line('Stabilisation must have double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
                call sys_halt()
              end if
              
              ! Set additional pointers
              call afcstab_getbase_FcoeffsAtEdge(rafcstab, p_Fcoefficients)
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation
              !-----------------------------------------------------------------
              call doOpEdgeAFCMat79Sngl(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_FcoeffsAtEdge, real(dscale,SP),&
                  bclear, bsymm, p_Fdata, p_Fcoefficients)
              
              ! Set state of stabilisation
              rafcstab%istabilisationSpec =&
                  ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)
              
              ! Do we need edge orientation?
              if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
                call gfsc_upwindOrientationSngl(p_Fcoefficients, p_IedgeList,&
                    p_FcoeffsAtEdge, 2, 3)
                rafcstab%istabilisationSpec =&
                    ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
              else
                rafcstab%istabilisationSpec =&
                    iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
              end if
              
            else
              
              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation but do not generate coeffs
              !-----------------------------------------------------------------
              call doOpEdgeStabMat79Sngl(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_FcoeffsAtEdge, real(dscale,SP),&
                  bclear, bsymm, p_Fdata)
            end if
            
          else   ! no stabilisation structure present

            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOpEdgeMat79Sngl(p_IedgeListIdx, p_IedgeList, rgroupFEMSet%NEDGE,&
                p_FcoeffsAtEdge, real(dscale,SP), bclear, p_Fdata)
          end if
          
        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
          call sys_halt()
        end select
        
      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorConst')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpNodeMat79Dble(Kld, NEQ, DcoeffsAtNode,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      integer :: idx,ieq,ij

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------

      if (bclear) call lalg_clearVector(Ddata)

      if (present(InodeList)) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ij)
        do idx = 1, size(InodeList,2)

          ! Get position of matrix entry
          ij  = InodeList(2,idx)

          ! Update the matrix coefficient
          Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do

      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ij)
        do ieq = 1, NEQ
          do ij = Kld(ieq), Kld(ieq+1)-1

            ! Update the matrix coefficient
            Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtNode(1,ij)
          end do
        end do
        !$omp end parallel do

      end if

    end subroutine doOpNodeMat79Dble

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpNodeMat79Sngl(Kld, NEQ, FcoeffsAtNode,&
        fscale, bclear, Fdata, InodeList)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      integer :: idx,ieq,ij

      !-------------------------------------------------------------------------
      ! Assemble matrix entries
      !-------------------------------------------------------------------------

      if (bclear) call lalg_clearVector(Fdata)

      if (present(InodeList)) then
        
        ! Loop over the subset of equations
        !$omp parallel do default(shared) private(ij)
        do idx = 1, size(InodeList,2)

          ! Get position of matrix entry
          ij  = InodeList(2,idx)

          ! Update the matrix coefficient
          Fdata(ij) = Fdata(ij) + fscale*FcoeffsAtNode(1,idx)
        end do
        !$omp end parallel do

      else
        
        ! Loop over all equations
        !$omp parallel do default(shared) private(ij)
        do ieq = 1, NEQ
          do ij = Kld(ieq), Kld(ieq+1)-1

            ! Update the matrix coefficient
            Fdata(ij) = Fdata(ij) + dscale*FcoeffsAtNode(1,ij)
          end do
        end do
        !$omp end parallel do

      end if

    end subroutine doOpNodeMat79Sngl

    !**************************************************************
    ! Assemble diagonal part of the operator
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpDiagMat79Dble(Kdiagonal, NEQ, DcoeffsAtNode,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! local variables
      integer :: ieq,ii,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and  position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Ddata(ii) = dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do
        end if

      else   ! do not clear matrix

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = Ddata(ii) + dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Ddata(ii) = Ddata(ii) + dscale*DcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        end if

      end if

    end subroutine doOpDiagMat79Dble

    !**************************************************************
    ! Assemble diagonal part of the operator
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpDiagMat79Sngl(Kdiagonal, NEQ, FcoeffsAtNode,&
        fscale, bclear, Fdata, InodeList)

      ! input parameters
      real(SP), dimension(:,:), intent(in) :: FcoeffsAtNode
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! local variables
      integer :: ieq,ii,idx

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and  position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Fdata(ii) = fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Fdata(ii) = fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do
        end if

      else   ! do not clear matrix

        if (present(InodeList)) then

          ! Loop over the subset of equations
          !$omp parallel do default(shared) private(ii,ieq)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do idx = 1, NEQ

            ! Get equation number and position of diagonal entry
            ieq = InodeList(1,idx)
            ii  = InodeList(2,idx)
            
            ! Update the diagonal coefficient
            Fdata(ii) = Fdata(ii) + fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        else

          ! Loop over all equations
          !$omp parallel do default(shared) private(ii)&
          !$omp if (NEQ > GFSC_NEQMIN_OMP)
          do ieq = 1, NEQ
            
            ! Get position of diagonal entry
            ii = Kdiagonal(ieq)
            
            ! Update the diagonal coefficient
            Fdata(ii) = Fdata(ii) + fscale*FcoeffsAtNode(1,ieq)
          end do
          !$omp end parallel do

        end if

      end if

    end subroutine doOpDiagMat79Sngl
    
    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeMat79Dble(IedgeListIdx, IedgeList, NEDGE,&
        DcoeffsAtEdge, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! local variables
      integer :: iedge,igroup,ij,ji     

      !$omp parallel default(shared) private(ij,ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
 
      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
        
        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
        
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Ddata(ij) = dscale*DcoeffsAtEdge(1,1,iedge)
            Ddata(ji) = dscale*DcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      else   ! do not clear matrix

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
            
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Ddata(ij) = Ddata(ij) + dscale*DcoeffsAtEdge(1,1,iedge)
            Ddata(ji) = Ddata(ji) + dscale*DcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      end if
      !$omp end parallel

    end subroutine doOpEdgeMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeMat79Sngl(IedgeListIdx, IedgeList, NEDGE,&
        FcoeffsAtEdge, fscale, bclear, Fdata)
      
      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! local variables
      integer :: iedge,igroup,ij,ji
      
      !$omp parallel default(shared) private(ij,ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
        
        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
        
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Fdata(ij) = fscale*FcoeffsAtEdge(1,1,iedge)
            Fdata(ji) = fscale*FcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      else   ! do not clear matrix

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
            
            ! Get node numbers and matrix positions
            ij = IedgeList(3,iedge)
            ji = IedgeList(4,iedge)
            
            ! Update the global operator
            Fdata(ij) = Fdata(ij) + fscale*FcoeffsAtEdge(1,1,iedge)
            Fdata(ji) = Fdata(ji) + fscale*FcoeffsAtEdge(1,2,iedge)
          end do
          !$omp end do
        end do ! igroup

      end if
      !$omp end parallel

    end subroutine doOpEdgeMat79Sngl

    !**************************************************************
    ! Assemble edge-by-edge operator with stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeStabMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, dscale, bclear, bsymm, Ddata)
      
      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear, bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata

      ! local variables
      real(DP) :: d_ij
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
   
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if
        
      else   ! do not clear matrix
        
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = max(0.0_DP,&
                  -DcoeffsAtEdge(1,1,iedge), -DcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if

      end if
      !$omp end parallel

    end subroutine doOpEdgeStabMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeStabMat79Sngl(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, FcoeffsAtEdge, fscale, bclear, bsymm, Fdata)

      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear, bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata

      ! local variables
      real(SP) :: d_ij
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then
   
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric artificial diffusion coefficient
              d_ij = max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup
          
        end if

      else   ! do not clear matrix
        
        if (bsymm) then

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Symmetric atificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Non-symmetric atificial diffusion coefficient
              d_ij = max(0.0_SP,&
                  -FcoeffsAtEdge(1,1,iedge), -FcoeffsAtEdge(1,2,iedge))
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup
          
        end if

      end if
      !$omp end parallel

    end subroutine doOpEdgeStabMat79Sngl

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeAFCMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, dscale, bclear, bsymm, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear,bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! output parameters
      real(DP), dimension(:,:), intent(out) :: Dcoefficients
      
      ! local variables
      real(DP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (bsymm) then
        
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1       
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Dcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_DP,  DcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Dcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if

      else   ! do not clear matrix
        
        if (bsymm) then
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Dcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_DP, -DcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_DP,  DcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_DP,  DcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Dcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - d_ij
              Ddata(jj) = Ddata(jj) - d_ij
              Ddata(ij) = Ddata(ij) + dscale*(DcoeffsAtEdge(1,1,iedge) + d_ij)
              Ddata(ji) = Ddata(ji) + dscale*(DcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if
        
      end if
      !$omp end parallel

    end subroutine doOpEdgeAFCMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOpEdgeAFCMat79Sngl(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, FcoeffsAtEdge, fscale, bclear, bsymm, Fdata, Fcoefficients)

      ! input parameters
      real(SP), dimension(:,:,:), intent(in) :: FcoeffsAtEdge
      real(SP), intent(in) :: fscale
      logical, intent(in) :: bclear,bsymm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(SP), dimension(:), intent(inout) :: Fdata
      
      ! output parameters
      real(SP), dimension(:,:), intent(out) :: Fcoefficients
      
      ! local variables
      real(SP) :: d_ij,k_ij,k_ji
      integer :: iedge,igroup,ii,ij,ji,jj
      
      !$omp parallel default(shared) private(ii,ij,ji,jj,d_ij,k_ij,k_ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      if (bclear) then

        if (bsymm) then
        
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1       
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Fcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_SP,  FcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Fcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if

      else   ! do not clear matrix
        
        if (bsymm) then
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              
              ! Symmetric AFC w/o edge orientation
              Fcoefficients(1:2,iedge) = (/d_ij, k_ij/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        else

          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Get node numbers and matrix positions
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Artificial diffusion coefficient
              d_ij = max(0.0_SP, -FcoeffsAtEdge(1,1,iedge))
              k_ij = max(0.0_SP,  FcoeffsAtEdge(1,1,iedge))
              k_ji = max(0.0_SP,  FcoeffsAtEdge(1,2,iedge))
              
              ! Non-symmetric AFC w/o edge orientation
              Fcoefficients(1:3,iedge) = (/d_ij, k_ij, k_ji/)
              
              ! Update the global operator
              Fdata(ii) = Fdata(ii) - d_ij
              Fdata(jj) = Fdata(jj) - d_ij
              Fdata(ij) = Fdata(ij) + fscale*(FcoeffsAtEdge(1,1,iedge) + d_ij)
              Fdata(ji) = Fdata(ji) + fscale*(FcoeffsAtEdge(1,2,iedge) + d_ij)
            end do
            !$omp end do
          end do ! igroup

        end if
        
      end if
      !$omp end parallel

    end subroutine doOpEdgeAFCMat79Sngl

  end subroutine gfsc_buildOperatorConst

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeBlock(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, rcollection)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! Note that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called.  Otherwise, an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the matrix entries may depend.
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSc_sim.inc'
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>
    
    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then

      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeBlock')
      call sys_halt()

    else

      call gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, rcollection)

    end if
    
  end subroutine gfsc_buildOperatorNodeBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorNodeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, dscale, bclear, rmatrix, rcollection)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet
    
    ! Vector on which the matrix entries may depend.
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSc_sim.inc'
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata,p_Dx
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:), pointer :: p_Kld

    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
      call sys_halt()
    end if

    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_NODEBASED)
      !-------------------------------------------------------------------------
      ! Node-based assembly
      !-------------------------------------------------------------------------

      ! Check if group finite element set is prepared
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA) .eq. 0) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
        call sys_halt()
      end if

      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)

        ! What data types are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call lsyssc_getbase_double(rmatrix, p_Ddata)

          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOperatorMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOperatorMat79Dble(p_Kld, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, dscale, bclear, p_Ddata, p_InodeList)
          end if

           case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorNodeScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow

    !**************************************************************
    ! Assemble operator node-by-node without stabilisation
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Dble(Kld, NEQ, DcoeffsAtNode,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
    end subroutine doOperatorMat79Dble

  end subroutine gfsc_buildOperatorNodeScalar
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorEdgeBlock(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim,&
      dscale, bclear, rmatrix, rcollection, rafcstab)
    
!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij=k_ji/)
    !
    ! Note that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called.  Otherwise, an error is thrown.
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorBlock), intent(in) :: rx
    
    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSc_sim.inc'
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then

      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeBlock')
      call sys_halt()

    else

      call gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim, dscale,&
          bclear, rmatrix, rcollection, rafcstab)

    end if

  end subroutine gfsc_buildOperatorEdgeBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildOperatorEdgeScalar(rgroupFEMSet, rx,&
      fcb_calcMatrixDiagSc_sim, fcb_calcMatrixSc_sim, dscale,&
      bclear, rmatrix, rcollection, rafcstab)

!<description>
    ! This subroutine assembles a discrete operator by the group
    ! finite element formulation. The matrix entries may depend on the
    ! values of the vector rx (e.g. the solution vector) and they are
    ! computed by calling user-defined callback functions.
    !
    ! This routine supports only node-by-node assembly.
    !
    ! If the optional stabilisation structure rafcstab is present then
    ! for each edge IJ the matrix entries and the artificial diffusion
    ! coefficient is stored according to the following convention:
    !
    ! For upwind-biased stabilisation
    !   Coefficients(1:3, IJ) = (/d_ij, k_ij, k_ji/)
    !
    ! For symmetric stabilisation
    !   Coefficients(1:2, IJ) = (/d_ij, k_ij=k_ji/)
!</description>

!<input>
    ! Group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Vector on which the matrix entries may depend.
    ! Note that this vector is only required for nonlinear
    ! problems which require the evaluation of the velocity
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute matrix entries
    include 'intf_calcMatrixDiagSc_sim.inc'
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! Global operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata,p_Dx
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoeffsAtNode
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:,:), pointer :: p_InodeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer, dimension(:), pointer :: p_Kdiagonal

    ! Check if matrix and vector have the same data type
    if ((rmatrix%cdataType .ne. rx%cdataType) .or.&
        (rmatrix%cdataType .ne. rgroupFEMSet%cdataType)) then
      call output_line('Data types mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
      call sys_halt()
    end if
    
    ! What type of assembly should be performed
    select case(rgroupFEMSet%cassemblyType)

    case (GFEM_EDGEBASED)
      !-------------------------------------------------------------------------
      ! Edge-based assembly
      !-------------------------------------------------------------------------
      
      ! Check if group finite element set is prepared
      if ((iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGESTRUCTURE) .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODEDATA)      .eq. 0) .or.&
          (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGEDATA)      .eq. 0)) then
        call output_line('Group finite element set does not provide required data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
        call sys_halt()
      end if
        
      ! Set pointers
      call gfem_getbase_IedgeListIdx(rgroupFEMSet, p_IedgeListIdx)
      call gfem_getbase_IedgeList(rgroupFEMSet, p_IedgeList)
      
      ! What kind of matrix are we?
      select case(rmatrix%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX7) then
          call lsyssc_getbase_Kld(rmatrix, p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
        end if
        
        ! What data types are we?
        select case(rmatrix%cdataType)
          
        case (ST_DOUBLE)
          ! Set pointers
          call gfem_getbase_DcoeffsAtNode(rgroupFEMSet, p_DcoeffsAtNode)
          call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
          call lsyssc_getbase_double(rmatrix, p_Ddata)
          call lsyssc_getbase_double(rx, p_Dx)
        
          ! Check if only a subset of the matrix is required
          if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_NODESTRUCTURE) .eq. 0) then
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata)
          else
            call gfem_getbase_InodeList(rgroupFEMSet, p_InodeList)
            call doOpDiagMat79Dble(p_Kdiagonal, rgroupFEMSet%NEQ,&
                p_DcoeffsAtNode, p_Dx, dscale, bclear, p_Ddata, p_InodeList)
          end if

          
          ! Do we have to build the stabilisation?
          if (present(rafcstab)) then
            
            ! Check if stabilisation has been prepared
            if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
              call output_line('Stabilisation has not been prepared!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
              call sys_halt()
            end if

            ! Check if coefficients should be stored in stabilisation
            if (rafcstab%h_CoefficientsAtEdge .ne. ST_NOHANDLE) then

              ! Check if stabilisation has the same data type
              if (rafcstab%cdataType .ne. ST_DOUBLE) then
                call output_line('Stabilisation must have double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
                call sys_halt()
              end if
                   
              ! Set additional pointers
              call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)

              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation and generate coefficients
              !-----------------------------------------------------------------
              call doOperatorAFCMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, p_Dx, dscale, bclear,&
                  p_Ddata, p_Dcoefficients)
              
              ! Set state of stabilisation
              rafcstab%istabilisationSpec =&
                  ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)
              
              ! Do we need edge orientation?
              if (rafcstab%climitingType .eq. AFCSTAB_LIMITING_UPWINDBIASED) then
                call gfsc_upwindOrientationDble(p_Dcoefficients, p_IedgeList,&
                    p_DcoeffsAtEdge, 2, 3)
                rafcstab%istabilisationSpec =&
                    ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION)
              else
                rafcstab%istabilisationSpec =&
                    iand(rafcstab%istabilisationSpec, not(AFCSTAB_HAS_EDGEORIENTATION))
              end if
        
            else

              !-----------------------------------------------------------------
              ! Assemble operator with stabilisation but do not generate coeffs
              !-----------------------------------------------------------------
              call doOperatorStabMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                  rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Ddata)
            end if
            
          else   ! no stabilisation structure present

            !-------------------------------------------------------------------
            ! Assemble operator without stabilisation
            !-------------------------------------------------------------------
            call doOperatorMat79Dble(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
                rgroupFEMSet%NEDGE, p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Ddata)
          end if
          
        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
          call sys_halt()
        end select
              
      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported assembly type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildOperatorEdgeScalar')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble diagonal part of the operator
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOpDiagMat79Dble(Kdiagonal, NEQ, DcoeffsAtNode, Dx,&
        dscale, bclear, Ddata, InodeList)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DcoeffsAtNode
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in), optional :: InodeList
      integer, intent(in) :: NEQ

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:), pointer :: DdataAtNode
      real(DP), dimension(:,:), pointer :: DcoefficientsAtNode
      integer, dimension(:,:), pointer  :: IdofsAtNode

      ! local variables
      integer :: idx,IEQset,IEQmax
      integer :: i,iidx,ii

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------
      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtNode,DdataAtNode,IEQmax,&
      !$omp         IdofsAtNode,i,idx,ii)

      ! Allocate temporal memory
      allocate(IdofsAtNode(2,GFSC_NEQSIM))
      allocate(DdataAtNode(GFSC_NEQSIM))
      allocate(DcoefficientsAtNode(1,GFSC_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSC_NEQSIM

        ! We always handle GFSC_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSC_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSC_NEQSIM)
        
        if (present(InodeList)) then

          ! Loop through all equations in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEQmax-IEQset+1
            
            ! Get actual index
            iidx = idx+IEQset-1
            
            ! Get equation number and position of diagonal entry
            i  = InodeList(1,iidx)
            ii = InodeList(2,iidx)
            
            ! Fill auxiliary arrays
            IdofsAtNode(1,idx) = i
            IdofsAtNode(2,idx) = ii
            DdataAtNode(idx)       = Dx(i)
          end do

        else

          ! Loop through all equations in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEQmax-IEQset+1

            ! Get actual equation number
            i = idx+IEQset-1
            
            ! Get position of diagonal entry
            ii = Kdiagonal(i)
            
            ! Fill auxiliary arrays
            IdofsAtNode(1,idx) = i
            IdofsAtNode(2,idx) = ii
            DdataAtNode(idx)       = Dx(i)
          end do

        end if

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSc_sim(&
            DdataAtNode(1:IEQmax-IEQset+1),&
            DcoeffsAtNode(:,IEQset:IEQmax),&
            IdofsAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = DcoefficientsAtNode(1,idx)
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IdofsAtNode(2,idx)
            
            ! Update the diagonal coefficient
            Ddata(ii) = Ddata(ii) + DcoefficientsAtNode(1,idx)
          end do

        end if
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(IdofsAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)
      !$omp end parallel

    end subroutine doOpDiagMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge without stabilisation
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, Dx, dscale, bclear, Ddata)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: Dcoefficients

      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ij,ji
      
      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ij,ji)

      ! Allocate temporal memory
      allocate(DdataAtEdge(2,GFSC_NEDGESIM))
      allocate(Dcoefficients(2,GFSC_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSC_NEDGESIM

          ! We always handle GFSC_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSC_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSC_NEDGESIM)
        
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,1:IEDGEmax-IEDGEset+1), rcollection)
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
                            
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(ij) = Dcoefficients(1,idx)
              Ddata(ji) = Dcoefficients(2,idx)
            end do

          else   ! do not clear matrix

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(ij) = Ddata(ij) + Dcoefficients(1,idx)
              Ddata(ji) = Ddata(ji) + Dcoefficients(2,idx)
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorMat79Dble
    
    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorStabMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, Dx, dscale, bclear, Ddata)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: Dcoefficients
      
      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ii,ij,ji,jj

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(Dcoefficients,DdataAtEdge,IEDGEmax,idx,iedge,ii,ij,ji,jj)
      
      ! Allocate temporal memory
      allocate(DdataAtEdge(2,GFSC_NEDGESIM))
      allocate(Dcoefficients(3,GFSC_NEDGESIM))
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSC_NEDGESIM

          ! We always handle GFSC_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSC_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSC_NEDGESIM)
          
          ! Dcoefficients is not allocated as temporal array
          ! since it is given as global output parameter.
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,1:IEDGEmax-IEDGEset+1), rcollection)
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,idx)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,idx)
              Ddata(ij) = Dcoefficients(2,idx) + Dcoefficients(1,idx) 
              Ddata(ji) = Dcoefficients(3,idx) + Dcoefficients(1,idx)
            end do

          else   ! do not clear matrix

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,idx)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,idx)
              Ddata(ij) = Ddata(ij) + Dcoefficients(2,idx) + Dcoefficients(1,idx) 
              Ddata(ji) = Ddata(ji) + Dcoefficients(3,idx) + Dcoefficients(1,idx)
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(Dcoefficients)
      !$omp end parallel
      
    end subroutine doOperatorStabMat79Dble

    !**************************************************************
    ! Assemble operator edge-by-edge with stabilisation and AFC data.
    ! All matrices are stored in matrix format 7 and 9
    
    subroutine doOperatorAFCMat79Dble(Kdiagonal, IedgeListIdx, IedgeList,&
        NEDGE, DcoeffsAtEdge, Dx, dscale, bclear, Ddata, Dcoefficients)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEDGE

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Ddata
      
      ! output parameters
      real(DP), dimension(:,:), intent(out) :: Dcoefficients
      
      ! auxiliary arras
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      
      ! local variables
      integer :: idx,IEDGEset,IEDGEmax
      integer :: iedge,igroup,ii,ij,ji,jj

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,IEDGEmax,idx,iedge,ii,ij,ji,jj)
      
      ! Allocate temporal memory
      allocate(DdataAtEdge(2,GFSC_NEDGESIM))
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSC_NEDGESIM

          ! We always handle GFSC_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSC_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSC_NEDGESIM)
          
          ! Dcoefficients is not allocated as temporal array
          ! since it is given as global output parameter.
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dcoefficients(:,IEDGEset:IEDGEmax), rcollection)
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global matrix
          if (bclear) then
            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Compute entries of low-order operator
              Dcoefficients(2,iedge) = Dcoefficients(2,iedge) + Dcoefficients(1,iedge)
              Dcoefficients(3,iedge) = Dcoefficients(3,iedge) + Dcoefficients(1,iedge)

              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,iedge)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,iedge)
              Ddata(ij) = Dcoefficients(2,iedge)
              Ddata(ji) = Dcoefficients(3,iedge)
            end do

          else   ! do not clear matrix

            do idx = 1, IEDGEmax-IEDGEset+1
              
              ! Get actual edge number
              iedge = idx+IEDGEset-1
              
              ! Get position of diagonal entries
              ii = Kdiagonal(IedgeList(1,iedge))
              jj = Kdiagonal(IedgeList(2,iedge))
              
              ! Get position of off-diagonal entries
              ij = IedgeList(3,iedge)
              ji = IedgeList(4,iedge)
              
              ! Compute entries of low-order operator
              Dcoefficients(2,iedge) = Dcoefficients(2,iedge) + Dcoefficients(1,iedge)
              Dcoefficients(3,iedge) = Dcoefficients(3,iedge) + Dcoefficients(1,iedge)

              ! Update the global operator
              Ddata(ii) = Ddata(ii) - Dcoefficients(1,iedge)
              Ddata(jj) = Ddata(jj) - Dcoefficients(1,iedge)
              Ddata(ij) = Ddata(ij) + Dcoefficients(2,iedge)
              Ddata(ji) = Ddata(ji) + Dcoefficients(3,iedge)
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      !$omp end parallel
      
    end subroutine doOperatorAFCMat79Dble
        
  end subroutine gfsc_buildOperatorEdgeScalar
 
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildFluxFCTBlock(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, rmatrix, rxTimeDeriv, rxPredictor)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! algebraic flux correction of FCT-type with or without the
    ! contributions of the consistent mass matrix. Note that this
    ! routine serves as a wrapper for block vectors. If there is only
    ! one block, then the corresponding scalar routine is called.
    ! Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep
    
    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Switch for flux assembly
    ! TRUE  : fluxes are not modified externally so that 
    !         quicker assembly procedures may be feasible
    ! FALSE : fluxes are truely assembled even if this
    !         leads to an expensive addition of zeros
    logical, intent(in) :: bquickAssembly

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: Consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorBlock), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorBlock), intent(in), optional :: rxPredictor
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    integer :: nblocks

    ! Check if block vector(s) contains exactly one block
    nblocks = rx%nblocks
    if (present(rxTimeDeriv)) nblocks = max(nblocks, rxTimeDeriv%nblocks)
    if (present(rxPredictor)) nblocks = max(nblocks, rxPredictor%nblocks)

    if (nblocks .ne. 1) then
      call output_line('Vector(s) must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Call subroutine for scalar vectors
    if (present(rxTimeDeriv)) then
      if (present(rxPredictor)) then
        ! ... both approximate time derivative and predictor are present
        call gfsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            rmatrix, rxTimeDeriv=rxTimeDeriv%RvectorBlock(1),&
            rxPredictor=rxPredictor%RvectorBlock(1))
      else
        ! ... only the approximate time derivative is present
        call gfsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            rmatrix, rxTimeDeriv=rxTimeDeriv%RvectorBlock(1))
      end if
    else
      if (present(rxPredictor)) then
        ! ... only the predictor is present
        call gfsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            rmatrix, rxPredictor=rxPredictor%RvectorBlock(1))
      else
        ! ... neither the approximate time derivative nor the predictor is present
        call gfsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            rmatrix)
      end if
    end if
    
  end subroutine gfsc_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildFluxFCTScalar(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, rmatrix, rxTimeDeriv, rxPredictor)
    
!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! algebraic flux correction of FCT-type with or without the
    ! contribution of the consistent mass matrix.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Switch for flux assembly
    ! TRUE  : fluxes are not modified externally so that 
    !         quicker assembly procedures may be feasible
    ! FALSE : fluxes are truely assembled even if this
    !         leads to an expensive addition of zeros
    logical, intent(in) :: bquickAssembly

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorScalar), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorScalar), intent(in), optional :: rxPredictor
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dmatrix,p_Dx
    real(DP), dimension(:), pointer :: p_DxTimeDeriv, p_DxPredictor
    real(DP), dimension(:), pointer :: p_Dflux0,p_Dflux,p_DfluxPrel,p_Dalpha
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
      call sys_halt()
    end if
    
    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge data structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ) .or.&
          (rafcstab%cmatrixFormat .ne. rmatrix%cmatrixFormat)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
        call sys_halt()
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rx, p_Dx)

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)
      
    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE)

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    
      !-------------------------------------------------------------------------
      ! Classical, iterative and semi-implicit nonlinear FEM-FCT algorithm
      !
      ! The raw antidiffusive fluxes for all algorithms can be
      ! assembled essentially in the same way.
      !
      ! $$ f_{ij}^n = -m_{ij}*(u_i^n-u_j^n)+(1-\theta)\Delta t d_{ij}^n(u_i^n-u_j^n) $$
      ! $$ f_{ij}^m = f_{ij}^n + m_{ij}*(u_i^m-u_j^m)+\theta\Delta t d_{ij}^m(u_i^m-u_j^m) $$
      !
      ! The only difference is that the amount of rejected  antidiffusion
      ! is subtracted from the initial fluxes in subsequent iterations if
      ! the iterative FEM-FCT algorithm is applied.
      ! Moreover the initial flux without mass contribution is stored
      ! separately for the semi-implicit FEM-FCT algorithm since it is
      ! used to constrain the raw antidiffusive fluxes in each iteration.
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edge data
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
        call output_line('Stabilisation does not provide edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
        call sys_halt()
      end if
      
      
      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_EXPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble explicit part of raw-antidiffive fluxes
        !-----------------------------------------------------------------------
                
        if (theta .ne. 1.0_DP) then
          ! Assemble the explicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij}^n = (1-\theta)\Delta t d_{ij}^n(u_i^n-u_j^n) $$
          call doFluxesByCoeffDble(p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale*(1.0_DP-theta),&
              bclear, p_Dflux0)
        elseif (.not.bquickAssembly .and. bclear) then
          ! Clear the explicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij}^n = 0 $$
          call lalg_clearVector(p_Dflux0, rafcstab%NEDGE)
          ! if bquickAssembly = TRUE then this step can be skipped
        end if

        !-----------------------------------------------------------------------

        ! Check for special treatment
        if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)
          call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
          
          ! We have to store the raw-antidiffusive fluxes based on the
          ! initial solution without contribution of the consistent
          ! mass matrix and without scaling by the implicitness parameter
          ! $$ f_{ij} = \Delta t d_{ij}^n(u_i^n-u_j^n) $$
          call doFluxesByCoeffDble(p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale, .true., p_DfluxPrel)

        elseif (rafcstab%cprelimitingType .ne. AFCSTAB_PRELIMITING_NONE) then
          
          ! We have to assemble the raw-antidiffusive fluxes for
          ! prelimiting separately based on the low-order predictor
          if (present(rxPredictor)) then
            
            ! Set pointers
            call lsyssc_getbase_double(rxPredictor, p_DxPredictor)
            call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
            
            if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
              ! Compute solution difference for standard prelimiting
              call doDifferencesDble(p_IedgeList, rafcstab%NEDGE,&
                  p_DxPredictor, p_DfluxPrel)
            elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
              ! Compute fluxes for minmod prelimiting
              call doFluxesByCoeffDble(p_IedgeList, rafcstab%NEDGE,&
                  p_DcoefficientsAtEdge, p_DxPredictor, dscale, .true., p_DfluxPrel)
            else
              call output_line('Invalid type of prelimiting!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
              call sys_halt()
            end if
          else
            call output_line('Fluxes for prelimiting cannot be assembled without predictor!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
            call sys_halt()
          end if
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the explicit part of the mass-antidiffusive fluxes
          ! $$ f_{ij}^n := f_{ij}^n - m_{ij}(u_i^n-u_j^n) $$
          call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
              p_Dmatrix, p_Dx, -dscale/tstep, .false., p_Dflux0)
        end if

      end if


      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_IMPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble implicit part of raw-antidiffusive fluxes
        !-----------------------------------------------------------------------
        
        if ((rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_ITERATIVE) .and.&
            iand(ioperationSpec, AFCSTAB_FCTFLUX_REJECTED) .ne. 0) then
          !---------------------------------------------------------------------
          ! Apply the rejected antidiffusive fluxes from the previous limiting
          ! step to the implicit part of the raw-antidiffusive fluxes
          ! --------------------------------------------------------------------
          
          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
            call output_line('Stabilisation does not provide correction factors!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
            call sys_halt()
          end if
          
          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
            call sys_halt()
          end if
          
          ! Set pointer
          call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
          
          ! Subtract amount of rejected antidiffusion
          call gfsc_combineFluxes(rafcstab%NEDGE, -1.0_DP, p_Dflux, p_Dflux0, p_Dalpha)
        end if
        
        !-----------------------------------------------------------------------
        
        if (theta .ne. 0.0_DP) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij} = \theta\Delta t d_{ij}(u_i-u_j) $$
          call doFluxesByCoeffDble(p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale*theta, bclear, p_Dflux)
        end if
        
        if (bquickAssembly) then
          ! We may check of either the implicit or explicit part are
          ! missing so that some redundant computations may be skipped
          if (theta .ne. 1.0_DP) then
            ! The explicit part of the raw-antidiffusive fluxes exists
            if (theta .ne. 0.0_DP) then
              ! The implicit part of the raw-antidiffusive fluxes
              ! exists; so combine them both into common fluxes
              call gfsc_combineFluxes(rafcstab%NEDGE, 1.0_DP, p_Dflux0, p_Dflux)
            else
              ! The implicit part of the raw-antidiffusive fluxes does
              ! not exists; the fluxes should be cleared so just
              ! overwrite them by the explicit part 
              call lalg_copyVector(p_Dflux0, p_Dflux)
            end if
            ! if theta = 1 then the explicit part does not exist
          end if
        else
          ! Truely combine both parts of the raw-antidiffusive fluxes
          call gfsc_combineFluxes(rafcstab%NEDGE, 1.0_DP, p_Dflux0, p_Dflux)
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the implicit part of the mass-antidiffusive fluxes
          ! $$ f_{ij}^m := f_{ij}^m + m_{ij}(u_i^m-u_j^m) $$
          call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
              p_Dmatrix, p_Dx, dscale/tstep, .false., p_Dflux)
    
        end if

      end if

      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
      

    case (AFCSTAB_LINFCT)
      
      !-------------------------------------------------------------------------
      ! Linearised FEM-FCT algorithm
      !-------------------------------------------------------------------------
      
      ! Check if stabilisation provides edge data
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
        call output_line('Stabilisation does not provide edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
        call sys_halt()
      end if

      !-------------------------------------------------------------------------
      
      ! Set pointer
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      
      ! Assemble spatial part of raw-antidiffusive fluxes
      call doFluxesByCoeffDble(p_IedgeList, rafcstab%NEDGE,&
          p_DcoefficientsAtEdge, p_Dx, dscale, bclear, p_Dflux)
      
      !-------------------------------------------------------------------------
      
      if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
        ! Compute fluxes for standard prelimiting based on the
        ! low-order solution which serves as predictor
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        call doDifferencesDble(p_IedgeList, rafcstab%NEDGE,&
            p_Dx, p_DfluxPrel)
        
      elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
        ! Make a backup of the spatial part of the raw-antidiffusive
        ! fluxes which are used for minmod prelimiting
        call lsyssc_copyVector(rafcstab%p_rvectorFlux,&
            rafcstab%p_rvectorFluxPrel)
      end if

      !-------------------------------------------------------------------------
      
      ! Do we have to include mass antidiffusion?
      if (present(rmatrix) .and. present(rxTimeDeriv)) then

        ! Set pointer
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        call lsyssc_getbase_double(rxTimeDeriv, p_DxTimeDeriv)
      
        ! Apply mass antidiffusion to antidiffusive fluxes based on
        ! the approximation to the time derivative
        call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
            p_Dmatrix, p_DxTimeDeriv, dscale, .false., p_Dflux)
      end if
      
      ! Set specifiers for raw antidiffusive fluxes
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
      

    case (AFCSTAB_LINFCT_MASS)
      
      !-------------------------------------------------------------------------
      ! FEM-FCT algorithm for mass antidiffusion
      !-------------------------------------------------------------------------
      
      if (present(rmatrix)) then
        
        ! Set pointers
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        
        ! Clear vector and assemble antidiffusive fluxes
        call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
            p_Dmatrix, p_Dx, dscale, .true., p_Dflux)

        ! Set specifiers for raw antidiffusive fluxes
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        
      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
        call sys_halt()
      end if
      
      
    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxFCTScalar')
      call sys_halt()
    end select

  contains
    
    ! Here, the working routines follow

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the edge-by-edge array DcoefficientsAtEdge
    
    subroutine doFluxesByCoeffDble(IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear
      
      real(DP), dimension(:), intent(inout) :: Dflux
      
      ! local variables
      integer :: iedge,i,j
      
      if (dscale .eq. 0.0_DP) then
        
        if (bclear) call lalg_clearVector(Dflux, NEDGE)
        
      elseif (dscale .eq. 1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if
        
      elseif (dscale .eq. -1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        end if
        
      else
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if
        
      end if
      
    end subroutine doFluxesByCoeffDble

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix Dmatrix
    
    subroutine doFluxesByMatrixDble(IedgeList, NEDGE,&
        Dmatrix, Dx, dscale, bclear, Dflux)
      
      real(DP), dimension(:), intent(in) :: Dmatrix, Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear
      
      real(DP), dimension(:), intent(inout) :: Dflux
    
      ! local variables
      integer :: iedge,i,j,ij
      
      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      elseif (dscale .eq. -1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        end if

      else

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesByMatrixDble
    
    !**************************************************************
    ! Assemble solution difference used for classical prelimiting.

    subroutine doDifferencesDble(IedgeList, NEDGE, Dx, Dflux)
      
      real(DP), dimension(:), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(out) :: Dflux

      ! local variables
      integer :: iedge,i,j

      !$omp parallel do default(shared) private(i,j)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Determine indices
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)
        
        ! Compute solution difference; in contrast to the literature,
        ! we compute the solution difference $u_i-u_j$ and check if
        ! $f_{ij}(u_i-u_j)<0$ in the prelimiting step.
        Dflux(iedge) = Dx(i)-Dx(j)
      end do
      !$omp end parallel do

    end subroutine doDifferencesDble
       
  end subroutine gfsc_buildFluxFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecFCTBlock(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-FCT type. Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called. Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx
    
    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
    
    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
    
    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry
    !</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks .ne. 1 .or. ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvVecFCTScalar(rafcstab, rmatrix, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1))
      
    end if
    
  end subroutine gfsc_buildConvVecFCTBlock

  ! *****************************************************************************
  
!<subroutine>
  
  subroutine gfsc_buildConvVecFCTScalar(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-FCT type. The idea of flux corrected
    ! transport can be traced back to the early SHASTA algorithm by
    ! Boris and Bock in the early 1970s. Zalesak suggested a fully
    ! multi-dimensional generalisation of this approach and paved the
    ! way for a large family of FCT algorithms.
    !
    ! This subroutine provides different algorithms:
    !
    ! Nonlinear FEM-FCT
    ! ~~~~~~~~~~~~~~~~~
    ! 1. Semi-explicit FEM-FCT algorithm
    !
    !    This is the classical algorithm which makes use of Zalesak`s
    !    flux limiter and recomputes and auxiliary positivity-
    !    preserving solution in each iteration step. 
    !    The details of this method can be found in:
    !
    !    D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
    !    conservation laws, Ergebnisberichte Angew. Math. 249,
    !    University of Dortmund, 2004.
    !
    ! 2. Iterative FEM-FCT algorithm
    !
    !    This is an extension of the classical algorithm which makes
    !    use of Zalesak`s flux limiter and tries to include the
    !    amount of rejected antidiffusion in subsequent iteration
    !    steps. The details of this method can be found in:
    !
    !    D. Kuzmin and M. Moeller, Algebraic flux correction I. Scalar
    !    conservation laws, Ergebnisberichte Angew. Math. 249,
    !    University of Dortmund, 2004.
    !
    ! 3. Semi-implicit FEM-FCT algorithm
    !
    !    This is the FCT algorithm that should be used by default. It
    !    is quite efficient since the nodal correction factors are
    !    only computed in the first iteration and used to limit the
    !    antidiffusive flux from the first iteration. This explicit
    !    predictor is used in all subsequent iterations to constrain
    !    the actual target flux.
    !    The details of this method can be found in:
    !
    !    D. Kuzmin and D. Kourounis, A semi-implicit FEM-FCT
    !    algorithm for efficient treatment of time-dependent
    !    problems, Ergebnisberichte Angew. Math. 302, University of
    !    Dortmund, 2005.
    !
    ! Linearised FEM-FCT
    ! ~~~~~~~~~~~~~~~~~~
    !
    ! 4. Linearised FEM-FCT algorithm
    !
    !    A new trend in the development of FCT algorithms is to
    !    linearise the raw antidiffusive fluxes about an intermediate
    !    solution computed by a positivity-preserving low-order
    !    scheme. By virtue of this linearisation, the costly
    !    evaluation of correction factors needs to be performed just
    !    once per time step. Furthermore, no questionable
    !    `prelimiting` of antidiffusive fluxes is required, which
    !    eliminates the danger of artificial steepening.
    !    The details of this method can be found in:
    !
    !    D. Kuzmin, Explicit and implicit FEM-FCT algorithms with
    !    flux linearization, Ergebnisberichte Angew. Math. 358,
    !    University of Dortmund, 2008.
!</description>

!<input>
    ! lumped mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
    
    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_DfluxPrel
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Check if stabilisation is prepeared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
      call sys_halt()
    end if

    ! Clear convective vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)

    !---------------------------------------------------------------------------
    ! The FEM-FCT algorithm is split into the following steps which
    ! can be skipped and performed externally by the user:
    !
    ! 1) Initialise the edgewise correction factors (alpha).
    !
    ! 2) Prelimit the antidiffusive fluxes (alpha).
    !
    ! 3) Compute the antidiffusive increments (Pp, Pm)
    !
    ! 4) Compute the local solution bounds (Qp, Qm).
    !
    ! 5) Compute the nodal correction factors (Rp, Rm).
    !
    ! 6) Compute edgewise correction factors (Alpha).
    !
    ! 7) Apply the limited antiddifusive fluxes.
    !-------------------------------------------------------------------------
    
    if (iand(ioperationSpec, AFCSTAB_FCTALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 1) Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------
            
      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if

    
    if (iand(ioperationSpec, AFCSTAB_FCTALGO_PRELIMIT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 2) Prelimit the raw antidiffusive fluxes (if required)
      !-------------------------------------------------------------------------
      if (rafcstab%cprelimitingType .ne. AFCSTAB_PRELIMITING_NONE) then
        
        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
          call sys_halt()
        end if
        
        ! Check if stabilisation provides edge-based structure
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide edge structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
          call sys_halt()
        end if
        
        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        
        if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
          ! Perform standard prelimiting
          call doStdPrelimitDble(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, p_Dflux, p_DfluxPrel, p_Dalpha)
        elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
          ! Perform minmod prelimiting
          call doMinModPrelimitDble(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, p_Dflux, p_DfluxPrel, p_Dalpha)
        else
          call output_line('Invalid type of prelimiting!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
          call sys_halt()
        end if
      end if
    end if
    

    if (iand(ioperationSpec, AFCSTAB_FCTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 3) Compute sums of antidiffusive increments 
      !-------------------------------------------------------------------------
      
      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if
      
      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        
        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        
        ! Compute sums of antidiffusive increments
        ! based on the prelimiting fluxes
        call doADIncrementsDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
      else
        ! Compute sums of antidiffusive increments
        ! based on the raw-antidiffusive fluxes
        call doADIncrementsDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
      end if

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 4) Compute local solution bounds
      !-------------------------------------------------------------------------
      
      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if
      
      ! Compute bounds
      call doBoundsDble(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEDGE, p_Dx, p_Dqp, p_Dqm)

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 5) Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)   .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Set additional pointers
      call lsyssc_getbase_double(rmatrix, p_ML)

      ! Compute nodal correction factors
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        call doLimitNodalDble(rafcstab%NEQ, dscale,&
            p_ML, p_Dx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        call doLimitNodalConstrainedDble(rafcstab%NEQ, dscale,&
            p_ML, p_Dx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      end if
      
      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 6) Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides nodal correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provide nodal correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Compute edgewise correction factors
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        call doLimitEdgewiseConstrainedDble(p_IedgeList,&
            rafcstab%NEDGE, p_DfluxPrel, p_Dflux, p_Drp, p_Drm, p_Dalpha)
      else
        call doLimitEdgewiseDble(p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dpp, p_Dpm, p_Drp, p_Drm, p_Dalpha)
      end if

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_CORRECT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 7) Correct antidiffusive fluxes and apply them
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edgewise correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provide edgewise correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then
        call lsyssc_getbase_double(rmatrix, p_ML)
        call doCorrectScaleByMassDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
      else
        call doCorrectDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, dscale, p_Dalpha, p_Dflux, p_Dy)
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes the standard way, as
    ! suggested by Boris and Book in their first FCT algorithm

    subroutine doStdPrelimitDble(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, DfluxPrel, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux,DfluxPrel
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      integer :: iedge

      ! Loop over all edges
      !$omp parallel do default(shared)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Check if the antidiffusive flux is directed down the gradient
        !   $f_ij*(u_i-u_j) < 0$
        ! and if its magnitude is larger than an absolute tolerance
        !  $ |f_ij| > tol$
        ! In this case, cancel the flux completely.
        if ((Dflux(iedge)*DfluxPrel(iedge) .lt. 0.0_DP) .and.&
            abs(Dflux(iedge)) .gt. AFCSTAB_PRELIMABS)&
            Dalpha(iedge) = 0.0_DP
      end do
      !$omp end parallel do

    end subroutine doStdPrelimitDble
    
    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes using minmod limiter

    subroutine doMinModPrelimitDble(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, DfluxPrel, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux,DfluxPrel
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      integer :: iedge

      ! Loop over all edges
      !$omp parallel do default(shared)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Check if the magnitude of the antidiffusive flux is larger
        ! than an absolute tolerance; otherwise no prelimiting is done
        if (abs(Dflux(iedge)) .gt. AFCSTAB_PRELIMABS) then
          ! Check if the antidiffusive flux is directed down the gradient
          !   $f_ij*fp_ij < 0$
          if (Dflux(iedge)*DfluxPrel(iedge) .lt. 0.0_DP) then
            ! Then, cancel the antidiffusive flux completely
            Dalpha(iedge) = 0.0_DP
          elseif (abs(Dflux(iedge)) .gt. abs(DfluxPrel(iedge))) then
            ! Check if the magnitude of the raw antidiffusive flux
            ! exceeds the magnitude of the prelimiting flux
            !   $|f_ij| > |fp_ij|$
            ! then set the correction factor as follows
            Dalpha(iedge) = min(Dalpha(iedge),DfluxPrel(iedge)/Dflux(iedge))
          end if
        end if
      end do
      !$omp end parallel do
      
    end subroutine doMinModPrelimitDble
    
    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes without prelimiting
    
    subroutine doADIncrementsDble(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, Dalpha, Dpp, Dpm)
      
      real(DP), dimension(:), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! The sums of positive/negative antidiffusive increments
      real(DP), dimension(:), intent(out) :: Dpp,Dpm
      
      ! local variables
      real(DP) :: f_ij,fp_ij,fm_ij
      integer :: i,iedge,igroup,j
      
      !$omp parallel default(shared) private(i,j,f_ij,fp_ij,fm_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          
          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)
      
          ! Separate fluxes into positive/negative contributions
          fp_ij = max(0.0_DP,f_ij)
          fm_ij = min(0.0_DP,f_ij)
          
          ! Compute the sums of antidiffusive increments
          Dpp(i) = Dpp(i) + fp_ij   ! += max(0.0_DP, f_ij)
          Dpp(j) = Dpp(j) - fm_ij   ! += max(0.0_DP,-f_ij)
          Dpm(i) = Dpm(i) + fm_ij   ! += min(0.0_DP, f_ij)
          Dpm(j) = Dpm(j) - fp_ij   ! += min(0.0_DP,-f_ij)
        end do
        !$omp end do
        
      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsDble

    !**************************************************************
    ! Assemble the local bounds from the predicted solution
    
    subroutine doBoundsDble(IedgeListIdx, IedgeList, NEDGE, Dx, Dqp, Dqm)
      
      real(DP), dimension(:), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      ! The local upper/lower bounds computed from Dx
      real(DP), dimension(:), intent(out) :: Dqp,Dqm
      
      ! local variables
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      ! Initialise Q`s by solution
      !$omp sections
      !$omp section
      call lalg_copyVector(Dx, Dqp)
      !$omp section
      call lalg_copyVector(Dx, Dqm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
        
        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          
          ! Compute local upper and lower bounds
          Dqp(i) = max(Dqp(i), Dx(j))
          Dqm(i) = min(Dqm(i), Dx(j))
          Dqp(j) = max(Dqp(j), Dx(i))
          Dqm(j) = min(Dqm(j), Dx(i))
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doBoundsDble

    !**************************************************************
    ! Compute the nodal correction factors without constraints
    
    subroutine doLimitNodalDble(NEQ, dscale,&
        ML, Dx, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      real(DP), dimension(:), intent(in) :: ML,Dx
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ
      
      ! The nodal correction factors for positive/negative fluxes
      real(DP), dimension(:), intent(inout) :: Drp,Drm
      
      ! local variables
      real(DP) :: diff
      integer :: ieq
      
      !$omp parallel sections default(shared) private(ieq,diff)
      
      !$omp section

      !$omp parallel do default(shared) private(ieq,diff)
      do ieq = 1, NEQ
        diff = Dqp(ieq)-Dx(ieq)
        if (dscale*Dpp(ieq) .gt. AFCSTAB_EPSABS) then
          Drp(ieq) = ML(ieq)*diff/(dscale*Dpp(ieq))
        else
          Drp(ieq) = 1.0_DP
        end if
      end do
      !$omp end parallel do

      !$omp section

      !$omp parallel do default(shared) private(ieq,diff)
      do ieq = 1, NEQ
        diff = Dqm(ieq)-Dx(ieq)
        if (dscale*Dpm(ieq) .lt. -AFCSTAB_EPSABS) then
          Drm(ieq) = ML(ieq)*diff/(dscale*Dpm(ieq))
        else
          Drm(ieq) = 1.0_DP
        end if
      end do
      !$omp end parallel do

      !$omp end parallel sections

    end subroutine doLimitNodalDble

    !**************************************************************
    ! Compute nodal correction factors with constraints
    
    subroutine doLimitNodalConstrainedDble(NEQ, dscale,&
        ML, Dx, Dpp, Dpm, Dqp, Dqm, Drp, Drm)
      
      real(DP), dimension(:), intent(in) :: ML,Dx
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ
      
      ! The nodal correction factors for positive/negative fluxes
      real(DP), dimension(:), intent(out) :: Drp,Drm
      
      ! local variables
      real(DP) :: diff
      integer :: ieq
      
      !$omp parallel sections default(shared) private(ieq,diff)
      
      !$omp section

      !$omp parallel do default(shared) private(ieq,diff)
      do ieq = 1, NEQ
        diff = Dqp(ieq)-Dx(ieq)
        if (dscale*Dpp(ieq) .gt. ML(ieq)*diff) then
          Drp(ieq) = ML(ieq)*diff/(dscale*Dpp(ieq))
        else
          Drp(ieq) = 1.0_DP
        end if
      end do
      !$omp end parallel do

      !$omp section

      !$omp parallel do default(shared) private(ieq,diff)
      do ieq = 1, NEQ
        diff = Dqm(ieq)-Dx(ieq)
        if (dscale*Dpm(ieq) .lt. ML(ieq)*diff) then
          Drm(ieq) = ML(ieq)*diff/(dscale*Dpm(ieq))
        else
          Drm(ieq) = 1.0_DP
        end if
      end do
      !$omp end parallel do

      !$omp end parallel sections

    end subroutine doLimitNodalConstrainedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    
    subroutine doLimitEdgewiseDble(IedgeList,&
        NEDGE, Dflux, Dpp, Dpm, Drp, Drm, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dpp,Dpm
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      
      ! On input: the edge-wise correction factors from previous
      !           multiplicative correction steps
      ! On exit: the edge-wise correction factors resulting from
      !          the nodal correction factors Rp and Rm
      real(DP), dimension(:), intent(inout) :: Dalpha
      
      ! local variables
      real(DP) :: f_ij,r_ij
      integer :: iedge,i,j
      
      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,f_ij,r_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)
        
        ! Get precomputed raw antidiffusive fluxes
        f_ij = Dflux(iedge)
 
        ! Compute nodal correction factors
        if (f_ij .gt. AFCSTAB_EPSABS) then
          r_ij = min(Drp(i),Drm(j))
        elseif (f_ij .lt. -AFCSTAB_EPSABS) then
          r_ij = min(Drp(j),Drm(i))
        else
          r_ij = 1.0_DP
        end if
        
        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes
    
    subroutine doLimitEdgewiseConstrainedDble(IedgeList,&
        NEDGE, Dflux1, Dflux2, Drp, Drm, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP) :: f1_ij,f2_ij,r_ij
      integer :: iedge,i,j
      
      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,f1_ij,f2_ij,r_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)
        
        ! Get precomputed raw antidiffusive fluxes
        f1_ij = Dflux1(iedge)
        f2_ij = Dflux2(iedge)
 
        ! Compute nodal correction factors
        if (f1_ij*f2_ij .le. 0.0_DP) then
          r_ij = 0.0_DP
        else
          if (f1_ij .ge. 0.0_DP) then
            r_ij = min(1.0_DP, f1_ij/f2_ij*min(Drp(i),Drm(j)))
          else
            r_ij = min(1.0_DP, f1_ij/f2_ij*min(Drp(j),Drm(i)))
          end if
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseConstrainedDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    
    subroutine doCorrectDble(IedgeListIdx, IedgeList,&
        NEDGE, dscale, Dalpha, Dflux, Dy)
      
      real(DP), dimension(:), intent(in) :: Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dy
      
      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          
          ! Correct antidiffusive flux
          f_ij = dscale * Dalpha(iedge) * Dflux(iedge)
          
          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij
          Dy(j) = Dy(j) - f_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel
      
    end subroutine doCorrectDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix
    
    subroutine doCorrectScaleByMassDble(IedgeListIdx,&
        IedgeList, NEDGE, dscale, ML, Dalpha, Dflux, Dy)
      
      real(DP), dimension(:), intent(in) :: ML,Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dy
      
      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup,j


      !$omp parallel default(shared) private(i,j,f_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
        
          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          
          ! Correct antidiffusive flux
          f_ij = dscale * Dalpha(iedge) * Dflux(iedge)
          
          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij/ML(i)
          Dy(j) = Dy(j) - f_ij/ML(j)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectScaleByMassDble

  end subroutine gfsc_buildConvVecFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecTVDBlock(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-TVD type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called. Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
    
    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_TVDALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry    
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks   .ne. 1 .or.&
        ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecTVDBlock')
      call sys_halt()

    else

      call gfsc_buildConvVecTVDScalar(rafcstab, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1))
      
    end if

  end subroutine gfsc_buildConvVecTVDBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecTVDScalar(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-TVD type.
    !
    ! A detailed description of the FEM-TVD limiter in general is given in:
    !
    !     D. Kuzmin and S. Turek, Multidimensional FEM-TVD paradigm
    !     for convection-dominated flows In:  Proceedings of the 
    !     IV European Congress on Computational Methods in Applied Sciences
    !     and Engineering (ECCOMAS 2004). Vol. II, ISBN 951-39-1869-6.
    !
    ! The method actually implemented in this routine is described in:
    !
    !     D. Kuzmin, Algebraic flux correction for finite element
    !     discretizations of coupled systems In: E. Onate,
    !     M. Papadrakakis and B. Schrefler (eds.) Computational
    !     Methods for Coupled Problems in Science and Engineering II,
    !     CIMNE, Barcelona, 2007, 653-656.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
    
    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_TVDALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm
    real(DP), dimension(:), pointer :: p_Dqp,p_Dqm
    real(DP), dimension(:), pointer :: p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dx,p_Dy,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq.0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq.0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq.0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecTVDScalar')
      call sys_halt()
    end if

    ! Clear destination vector?
    if (bclear) call lsyssc_clearVector(ry)
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)
    
    ! Perform flux limiting of TVD-type
    call doLimitDble(p_IedgeListIdx, p_IedgeList,&
        p_DcoefficientsAtEdge, p_Dx, dscale, rafcstab%NEDGE,&
        p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dflux, p_Dy)
    
  contains

    ! Here, the working routine follows
    
    !**************************************************************
    ! The upwind-biased FEM-TVD limiting procedure
    
    subroutine doLimitDble(IedgeListIdx, IedgeList,&
        DcoefficientsAtEdge, Dx, dscale, NEDGE,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dflux, Dy)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(:), intent(inout) :: Dflux,Dy
      
      ! local variables
      real(DP) :: d_ij,diff,f_ij,fm_ij,fp_ij,l_ji
      integer :: i,iedge,igroup,j
      
      
      if ((iand(ioperationSpec, AFCSTAB_TVDALGO_ADINCREMENTS) .ne. 0) .and.&
          (iand(ioperationSpec, AFCSTAB_TVDALGO_BOUNDS)       .ne. 0)) then
        
        ! Clear nodal vectors
        call lalg_clearVectorDble(Dpp)
        call lalg_clearVectorDble(Dpm)
        call lalg_clearVectorDble(Dqp)
        call lalg_clearVectorDble(Dqm)

        if (iand(ioperationSpec, AFCSTAB_TVDALGO_ADFLUXES) .ne. 0) then

          ! Assemble antidiffusive fluxes, sums of antidiffusive
          ! increments and the upper and lower bounds simultaneously
          ! to achieve higher efficiency.
        
          !$omp parallel default(shared)&
          !$omp private(i,j,d_ij,l_ji,diff,f_ij,fm_ij,fp_ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            
            ! Do nothing for empty groups
            if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
            
            ! Loop over the edges
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Determine indices
              i = IedgeList(1,iedge)
              j = IedgeList(2,iedge)
              
              ! Determine coefficients
              d_ij = DcoefficientsAtEdge(1,iedge)
              l_ji = DcoefficientsAtEdge(3,iedge)
              
              ! Determine solution difference
              diff = Dx(i)-Dx(j)
              
              ! Prelimit the antidiffusive flux
              ! F`_IJ=MIN(-P_IJ,L_JI)(DX_I-DX_J)
              f_ij = dscale*min(d_ij,l_ji)*diff
              
              ! And store it
              Dflux(iedge) = f_ij

              ! Separate fluxes into positive/negative contributions
              fp_ij = max(0.0_DP,f_ij)
              fm_ij = min(0.0_DP,f_ij)
              
              ! Assemble P`s accordingly
              Dpp(i) = Dpp(i) + fp_ij   ! += max(0.0_DP, f_ij)
              Dpm(i) = Dpm(i) + fm_ij   ! += min(0.0_DP, f_ij)
              
              ! Assemble Q`s
              Dqp(i) = Dqp(i) - fm_ij   ! += max(0.0_DP,-f_ij)
              Dqp(j) = Dqp(j) + fp_ij   ! += max(0.0_DP, f_ij)
              Dqm(i) = Dqm(i) - fp_ij   ! += min(0.0_DP,-f_ij)
              Dqm(j) = Dqm(j) + fm_ij   ! += min(0.0_DP, f_ij)
            end do
            !$omp end do
          end do ! igroup
          !$omp end parallel

          ! Set specifier
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
          
        else

          ! Check if stabilisation is prepared
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq.0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecTVDScalar')
            call sys_halt()
          end if
          
          ! Assemble antidiffusive fluxes, sums of antidiffusive
          ! increments and the upper and lower bounds simultaneously
          
          !$omp parallel default(shared)&
          !$omp private(i,j,f_ij,fm_ij,fp_ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          
          ! Loop over the edge groups and process all edges of one group
          ! in parallel without the need to synchronize memory access
          do igroup = 1, size(IedgeListIdx)-1
            
            ! Do nothing for empty groups
            if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
            
            ! Loop over the edges
            !$omp do
            do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
              
              ! Determine indices
              i = IedgeList(1,iedge)
              j = IedgeList(2,iedge)
              
              ! Get antidiffusive flux
              f_ij = Dflux(iedge)
              
              ! Assemble P`s accordingly
              Dpp(i) = Dpp(i) + fp_ij   ! += max(0.0_DP, f_ij)
              Dpm(i) = Dpm(i) + fm_ij   ! += min(0.0_DP, f_ij)
              
              ! Assemble Q`s
              Dqp(i) = Dqp(i) - fm_ij   ! += max(0.0_DP,-f_ij)
              Dqp(j) = Dqp(j) + fp_ij   ! += max(0.0_DP, f_ij)
              Dqm(i) = Dqm(i) - fp_ij   ! += min(0.0_DP,-f_ij)
              Dqm(j) = Dqm(j) + fm_ij   ! += min(0.0_DP, f_ij)
            end do
            !$omp end do
          end do ! igroup
          !$omp end parallel
          
          ! Set specifier
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)

        end if

      elseif (iand(ioperationSpec, AFCSTAB_TVDALGO_ADFLUXES) .ne. 0) then

        ! Assemble antidiffusive fluxes
        
        !$omp parallel default(shared)&
        !$omp private(i,j,d_ij,l_ji,diff,f_ij)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        
        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          
          ! Do nothing for empty groups
          if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
          
          ! Loop over the edges
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)
            
            ! Determine solution difference
            diff = dscale*(Dx(i)-Dx(j))
            
            ! Prelimit the antidiffusive flux
            ! F`_IJ=MIN(-P_IJ,L_JI)(DX_I-DX_J)
            f_ij = min(d_ij,l_ji)*diff
            
            ! And store it
            Dflux(iedge) = f_ij
          end do
          !$omp end do
        end do ! igroup
        !$omp end parallel
          
        ! Set specifier
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)

      end if

      !-------------------------------------------------------------------------

      if (iand(ioperationSpec, AFCSTAB_TVDALGO_LIMIT) .ne. 0) then
        
        ! Check if stabilisation is prepared
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)  .eq.0) .or.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS).eq.0)) then
          call output_line('Stabilisation does not provide bounds or increments!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecTVDScalar')
          call sys_halt()
        end if

        ! Apply the nodal limiter
        
        !$omp parallel sections
        !$omp section
        Drp = gfsc_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)

        !$omp section
        Drm = gfsc_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)
        !$omp end parallel sections

        ! Set specifier
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
        
      end if

      !-------------------------------------------------------------------------

      if (iand(ioperationSpec, AFCSTAB_TVDALGO_CORRECT) .ne. 0) then
        
        ! Check if stabilisation is prepared
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq.0) .or.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)    .eq.0)) then
          call output_line('Stabilisation does not provide fluxes or limiting factors!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecTVDScalar')
          call sys_halt()
        end if

        !$omp parallel default(shared)&
        !$omp private(i,j,d_ij,l_ji,diff,f_ij)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        
        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronize memory access
        do igroup = 1, size(IedgeListIdx)-1
          
          ! Do nothing for empty groups
          if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
          
          ! Loop over the edges
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Get precomputed raw antidiffusive flux
            f_ij = Dflux(iedge)
            
            ! Apply correction factor and store limited flux
            if (f_ij .gt. 0.0_DP) then
              f_ij = Drp(i)*f_ij
            else
              f_ij = Drm(i)*f_ij
            end if
            
            ! Update the vector
            Dy(i) = Dy(i)+f_ij
            Dy(j) = Dy(j)-f_ij
          end do
          !$omp end do
          
        end do ! igroup
        !$omp end parallel
        
      end if

    end subroutine doLimitDble
    
  end subroutine gfsc_buildConvVecTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecGPBlock(rafcstab, rmatrix,&
      rx, rx0, theta, dscale, bclear, ioperationSpec, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation of FEM-GP type.  Note that this routine serves as
    ! a wrapper for block vectors. If there is only one block, then
    ! the corresponding scalar routine is called.  Otherwise, an error
    ! is thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
    
    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_GPALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry   
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if ((rx%nblocks   .ne. 1) .or.&
        (rx0%nblocks  .ne. 1) .or.&
        (ry%nblocks .ne. 1)) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecGPBlock')
      call sys_halt()

    else

      call gfsc_buildConvVecGPScalar(rafcstab, rmatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), theta, dscale, bclear, ioperationSpec,&
          ry%RvectorBlock(1))
      
    end if
  end subroutine gfsc_buildConvVecGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecGPScalar(rafcstab, rmatrix,&
      rx, rx0, theta, dscale, bclear, ioperationSpec, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation using the general purpose limiter.
    !
    ! A detailed description of the FEM-GP limiter in general is given in:
    !
    !     D. Kuzmin, On the design of general-purpose flux 
    !     limiters for implicit FEM with a consistent mass matrix.
    !     I. Scalar convection.
    !     J. Comput. Phys.  219  (2006) 513-531.
    !
    ! Note however, that is is quite expensive and not recommended as
    ! a standard limiter. In fact, it is only implemented to
    ! demonstrate that the construction of general-purpose flux
    ! limiters is possible. If you want to recover the consistent mass
    ! matrix for time-dependent problems, you should apply flux
    ! correction of FCT type.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! initial solution vector
    type(t_vectorScalar), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
    
    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_TVDALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_MC,p_Dx,p_Dx0,p_Dy,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecGPScalar')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rmatrix, p_MC)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(rx0, p_Dx0)
    call lsyssc_getbase_double(ry, p_Dy)

    ! Perform flux limiting by the general purpose limiter
    call doLimitDble(p_IedgeListIdx, p_IedgeList,&
        p_DcoefficientsAtEdge, p_MC, p_Dx, p_Dx0, theta, dscale,&
        rafcstab%NEDGE, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm,&
        p_Dflux, p_Dflux0, p_Dy)
    
    ! Set specifier
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
    
  contains

    ! Here, the working routine follows
    
    !**************************************************************
    ! The FEM-GP limiting procedure
    
    subroutine doLimitDble(IedgeListIdx, IedgeList,&
        DcoefficientsAtEdge, MC, Dx, Dx0, theta, dscale, NEDGE,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dflux, Dflux0, Dy)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0
      real(DP), intent(in) :: theta,dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(:), intent(inout) :: Dflux,Dflux0,Dy
      
      ! local variables
      real(DP) :: d_ij,df_ij,f_ij,l_ij,l_ji,m_ij,p_ij,pf_ij,q_ij,q_ji
      real(DP) :: diff,diff0,diff1
      integer :: i,iedge,igroup,ij,j
      

      ! Clear nodal vectors
      call lalg_clearVectorDble(Dpp)
      call lalg_clearVectorDble(Dpm)
      call lalg_clearVectorDble(Dqp)
      call lalg_clearVectorDble(Dqm)
      
      !$omp parallel default(shared)&
      !$omp private(d_ij,df_ij,diff,diff0,diff1,&
      !$omp         f_ij,i,ij,j,l_ij,l_ji,m_ij,p_ij,pf_ij,q_ij,q_ji)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Determine indices
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          ij = IedgeList(3,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          l_ij = DcoefficientsAtEdge(2,iedge)
          l_ji = DcoefficientsAtEdge(3,iedge)
          m_ij = MC(ij)
        
          ! Compute: diff1 = dt*theta*(Dx_i-Dx_j) + dt*(1-theta)*(Dx0_i-Dx0_j)
          diff1 = Dx(i)-Dx(j); diff0 = Dx0(i)-Dx0(j)
          diff  = dscale*(theta*diff1+(1.0_DP-theta)*diff0)
          
          ! Compute antidiffusive flux f_ij=min(0,p_ij)*(Dx_j-Dx_i)
          if (abs(diff) < AFCSTAB_EPSABS) then
            p_ij = 0
            f_ij = 0
          else
            p_ij = max(0.0_DP,m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = p_ij*diff
          end if
          
          ! Prelimit the antidiffusive flux F`_IJ=MIN(-P_IJ,L_JI)(DX_I-DX_J)
          pf_ij = min(p_ij,l_ji)*diff; Dflux0(iedge) = pf_ij
          
          ! Compute the remaining flux dF_IJ=F_IJ-F`_IJ
          df_ij = f_ij-pf_ij; Dflux(iedge) = df_ij
          
          ! Assemble P`s accordingly
          Dpp(i) = Dpp(i)+max(0.0_DP,  f_ij)
          Dpm(i) = Dpm(i)+min(0.0_DP,  f_ij)
          Dpp(j) = Dpp(j)+max(0.0_DP,-df_ij)
          Dpm(j) = Dpm(j)+min(0.0_DP,-df_ij)
          
          q_ij = m_ij/dscale+l_ij
          q_ji = m_ij/dscale+l_ji
          
          ! Assemble Q`s
          Dqp(i) = Dqp(i)+q_ij*max(0.0_DP,-diff)
          Dqm(i) = Dqm(i)+q_ij*min(0.0_DP,-diff)
          Dqp(j) = Dqp(j)+q_ji*max(0.0_DP, diff)
          Dqm(j) = Dqm(j)+q_ji*min(0.0_DP, diff)
        end do
        !$omp end do
      end do ! igroup

      ! Apply nodal limiter
      !$omp single
      Drp = gfsc_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)
      Drm = gfsc_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
        
        ! Loop over the edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
          ! Get precomputed fluxes
          pf_ij = Dflux0(iedge); df_ij = Dflux(iedge)
          
          ! Limit upwind contribution
          if (pf_ij > 0.0_DP) then
            pf_ij = Drp(i)*pf_ij
          else
            pf_ij = Drm(i)*pf_ij
          end if
          
          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = min(Drp(i), Drm(j))*df_ij
          else
            df_ij = min(Drm(i), Drp(j))*df_ij
          end if
          
          f_ij = pf_ij+df_ij
          
          ! Update the vector
          Dy(i) = Dy(i)+f_ij
          Dy(j) = Dy(j)-f_ij
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doLimitDble
    
  end subroutine gfsc_buildConvVecGPScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecSymmBlock(rafcstab, rx, dscale, ry)

!<description>
    ! This subroutine assembles the convective vector and applies
    ! stabilisation by means of symmetric flux limiting for diffusion
    ! operators.  Note that this routine serves as a wrapper for block
    ! vectors. If there is only one block, then the corresponding
    ! scalar routine is called.  Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks .ne. 1 .or. ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecSymmBlock')
      call sys_halt()

    else

      call gfsc_buildConvVecSymmScalar(rafcstab,&
          rx%RvectorBlock(1), dscale, ry%RvectorBlock(1))

    end if
  end subroutine gfsc_buildConvVecSymmBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvVecSymmScalar(rafcstab, rx, dscale, ry)

!<description>
    ! This subroutine assembles the convective vector and applies stabilisation
    ! by means of symmetric flux limiting for diffusion operators.
    !
    ! Yet, there is no publication available. This routine is based on
    ! private communication with D. Kuzmin.
    !
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dx,p_Dy,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    
    
    ! Check if stabilisation is prepared
    if ((rafcstab%cafcstabType .ne. AFCSTAB_SYMMETRIC) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0)    .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecSymmScalar')
      call sys_halt()
    end if
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)
    
    ! Perform symmetric flux limiting
    call doLimitDble(p_IedgeListIdx, p_IedgeList,&
        p_DcoefficientsAtEdge, p_Dx, dscale, rafcstab%NEDGE,&
        p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dflux, p_Dy)

    ! Set specifier
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
    
  contains
    
    ! Here, the working routine follows
    
    !**************************************************************
    ! Perform symmetric flux limiting
    
    subroutine doLimitDble(IedgeListIdx, IedgeList,&
        DcoefficientsAtEdge, Dx, dscale, NEDGE,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dflux, Dy)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      real(DP), dimension(:), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(:), intent(inout) :: Dflux,Dy

      ! local variables
      real(DP) :: d_ij,f_ij,s_ij,diff
      integer :: i,iedge,igroup,j
      
      
      ! Clear nodal vectors
      call lalg_clearVectorDble(Dpp)
      call lalg_clearVectorDble(Dpm)
      call lalg_clearVectorDble(Dqp)
      call lalg_clearVectorDble(Dqm)

      !$omp parallel default(shared)&
      !$omp private(d_ij,diff,f_ij,i,j,s_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
          ! Determine coefficients
          d_ij = DcoefficientsAtEdge(1,iedge)
          s_ij = DcoefficientsAtEdge(2,iedge)
          
          ! Determine fluxes
          diff = Dx(i)-Dx(j); f_ij = d_ij*diff
          Dflux(iedge) = f_ij
          
          ! Sums of raw positive/negative fluxes
          Dpp(i) = Dpp(i)+max(0.0_DP, f_ij)
          Dpp(j) = Dpp(j)+max(0.0_DP,-f_ij)
          Dpm(i) = Dpm(i)+min(0.0_DP, f_ij)
          Dpm(j) = Dpm(j)+min(0.0_DP,-f_ij)
          
          ! Upper/lower bounds
          f_ij = -s_ij*diff
          Dqp(i) = Dqp(i)+max(0.0_DP, f_ij)
          Dqp(j) = Dqp(j)+max(0.0_DP,-f_ij)
          Dqm(i) = Dqm(i)+min(0.0_DP, f_ij)
          Dqm(j) = Dqm(j)+min(0.0_DP,-f_ij)
        end do
        !$omp end do
      end do ! igroup
      
      ! Apply the nodal limiter
      !$omp single
      Drp = gfsc_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)
      Drm = gfsc_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
        
        ! Loop over the edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
          ! Get precomputed raw antidiffusive flux
          f_ij = Dflux(iedge)
          
          if (f_ij > 0.0_DP) then
            f_ij = dscale*min(Drp(i), Drm(j))*f_ij
          else
            f_ij = dscale*min(Drm(i), Drp(j))*f_ij
          end if
          
          ! Update the vector
          Dy(i) = Dy(i)+f_ij
          Dy(j) = Dy(j)-f_ij
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doLimitDble

  end subroutine gfsc_buildConvVecSymmScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildFluxLPTBlock(rafcstab, rx, dscale, bclear,&
      ioperationSpec, rmatrix)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! linearity-preserving flux correction, whereby the coefficients
    ! are determined from the off-diagonal entries of the matrix. Note
    ! that this routine serves as a wrapper for block vectors. If
    ! there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, an error is thrown. 
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: coefficient matrix which is used instead of the
    ! coefficients at edges provided by the stabilisation structure
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    if (rx%nblocks .ne. 1) then
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTBlock')
      call sys_halt()
    end if

    ! Call subroutine for scalar vectors
    call gfsc_buildFluxLPTScalar(rafcstab, rx%RvectorBlock(1),&
        dscale, bclear, ioperationSpec, rmatrix)
    
  end subroutine gfsc_buildFluxLPTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildFluxLPTScalar(rafcstab, rx, dscale, bclear,&
      ioperationSpec, rmatrix)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! linearity-preserving flux correction, whereby the coefficients
    ! are determined from the off-diagonal entries of the matrix.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: coefficient matrix which is used instead of the
    ! coefficients at edges provided by the stabilisation structure
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dmatrix,p_Dx
    real(DP), dimension(:), pointer :: p_Dflux,p_Dq
    real(DP), dimension(:,:), pointer :: p_DboundsAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer, dimension(:,:), pointer :: p_IedgeList

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge data structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ) .or.&
          (rafcstab%cmatrixFormat .ne. rmatrix%cmatrixFormat)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
        call sys_halt()
      end if
    end if

    
    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case(AFCSTAB_NLINLPT_MASS,&
         AFCSTAB_LINLPT_MASS)
      
      if (present(rmatrix)) then

        ! Set pointer
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX) .ne. 0) then

          ! Set pointers
          call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
          call lsyssc_getbase_double(rx, p_Dx)

          ! Assemble raw antidiffusive fluxes for mass antidiffusion
          call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
              p_Dmatrix, p_Dx, dscale, bclear, p_Dflux)

          ! Set specifiers for raw antidiffusive fluxes
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        end if

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX_BOUNDS) .ne. 0) then

          ! Check if stabilisation provides edgewise bounds
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS) .ne. 0) then

            ! Set pointers
            call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
            call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
            call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
            call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)
            
            ! Assemble bounds based on matrix coefficients
            call doBoundsByMatrixDble(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_Dmatrix, p_DboundsAtEdge, p_Dq)
            
            ! Set specifiers for nodal bounds
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)

          else

            call output_line('Stabilisation does not provide bounds at edges!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if

      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
        call sys_halt()
      end if

      !-------------------------------------------------------------------------
      
    case(AFCSTAB_NLINLPT_UPWINDBIASED)

      ! Check if stabilisation provides edge data
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .ne. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .ne. 0)) then

        ! Set pointer
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX) .ne. 0) then
          
          ! Set pointers
          call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
          call lsyssc_getbase_double(rx, p_Dx)

          ! Assemble upwind-biased raw antidiffusive fluxes
          call doFluxesUpwindBiasedDble(p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale, bclear, p_Dflux)

          ! Set specifiers for raw antidiffusive fluxes
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        end if

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX_BOUNDS) .ne. 0) then
        
          ! Check if stabilisation provides edgewise bounds
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS) .ne. 0) then

            ! Set pointers
            call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
            call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
            call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
            call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)

            ! Assemble bounds based on edge coefficients
            call doBoundsByCoeffDble(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_DcoefficientsAtEdge, p_DboundsAtEdge, p_Dq)

            ! Set specifiers for nodal bounds
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)

          else
            
            call output_line('Stabilisation does not provide bounds at edges!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if
        
      else
        call output_line('Stabilisation does not provide (oriented) edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
        call sys_halt()
      end if
      
      !-------------------------------------------------------------------------

    case(AFCSTAB_NLINLPT_SYMMETRIC,&
         AFCSTAB_LINLPT_SYMMETRIC,&
         AFCSTAB_LINLPT_UPWINDBIASED)

      ! Remark: in the linearised version, all fluxes must be limited
      ! in a symmetric fashion since the flux into the downwind node
      ! is no longer balanced by the nonoscillatory part of the
      ! Galerkin operator; cf. comment by D. Kuzmin, 2011

      ! Check if stabilisation provides edge data
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .ne. 0) then

        ! Set pointer
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)

        !-----------------------------------------------------------------------

        if (iand(ioperationSpec, AFCSTAB_LPTFLUX) .ne. 0) then

          ! Set pointers
          call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
          call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
          call lsyssc_getbase_double(rx, p_Dx)

          ! Assemble upwind-biased raw antidiffusive fluxes
          call doFluxesSymmetricDble(p_IedgeList, rafcstab%NEDGE,&
              p_DcoefficientsAtEdge, p_Dx, dscale, bclear, p_Dflux)

          ! Set specifiers for raw antidiffusive fluxes
          rafcstab%istabilisationSpec =&
              ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        end if

        !-----------------------------------------------------------------------
        
        if (iand(ioperationSpec, AFCSTAB_LPTFLUX_BOUNDS) .ne. 0) then
          
          ! Check if stabilisation provides edgewise bounds
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS) .ne. 0) then
            
            ! Set pointers
            call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
            call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
            call afcstab_getbase_DboundsAtEdge(rafcstab, p_DboundsAtEdge)
            call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)
            
            ! Assemble bounds based on edge coefficients
            call doBoundsByCoeffDble(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, p_DcoefficientsAtEdge, p_DboundsAtEdge, p_Dq)

            ! Set specifiers for nodal bounds
            rafcstab%istabilisationSpec =&
                ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
            
          else
            
            call output_line('Stabilisation does not provide bounds at edges!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if

      else
        call output_line('Stabilisation does not provide edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
        call sys_halt()
      end if
      
    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildFluxLPTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble symmetric raw antidiffusive fluxes using the
    ! coefficients supplied by the CSR-matrix Dmatrix
    
    subroutine doFluxesByMatrixDble(IedgeList, NEDGE,&
        Dmatrix, Dx, dscale, bclear, Dflux)
      
      real(DP), dimension(:), intent(in) :: Dmatrix, Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear
      
      real(DP), dimension(:), intent(inout) :: Dflux

      ! local variables
      integer :: iedge,i,j,ij
      
      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      elseif (dscale .eq. -1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = Dmatrix(ij) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + Dmatrix(ij) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        end if

      else

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * Dmatrix(ij) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesByMatrixDble

    !**************************************************************
    ! Assemble upwind-biased raw antidiffusive fluxes using the
    ! coefficients supplied by the edge-by-edge array DcoefficientsAtEdge
    
    subroutine doFluxesUpwindBiasedDble(IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear
      
      real(DP), dimension(:), intent(inout) :: Dflux

      ! local variables
      real(DP) :: f_ij,d_ij,l_ji,g_ij
      integer :: iedge,i,j
      
      if (dscale .eq. 0.0_DP) then
        
        if (bclear) call lalg_clearVector(Dflux, NEDGE)
        
      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))
            
            if (l_ji .lt. d_ij) then
              ! f_ij := minmod(f_ij,g_ij)
              if (f_ij*g_ij .le. 0.0_DP) then
                f_ij = 0.0_DP
              elseif (abs(g_ij) .lt. abs(f_ij)) then
                f_ij = g_ij
              end if
            end if
            
            ! Store raw antidiffusive fluxes
            Dflux(iedge) = f_ij
          end do
          !$omp end parallel do

        else
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))
            
            if (l_ji .lt. d_ij) then
              ! f_ij := minmod(f_ij,g_ij)
              if (f_ij*g_ij .le. 0.0_DP) then
                f_ij = 0.0_DP
              elseif (abs(g_ij) .lt. abs(f_ij)) then
                f_ij = g_ij
              end if
            end if

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + f_ij
          end do
          !$omp end parallel do
        end if

      elseif (dscale .eq. -1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))
            
            if (l_ji .lt. d_ij) then
              ! f_ij := minmod(f_ij,g_ij)
              if (f_ij*g_ij .le. 0.0_DP) then
                f_ij = 0.0_DP
              elseif (abs(g_ij) .lt. abs(f_ij)) then
                f_ij = g_ij
              end if
            end if

            ! Store raw antidiffusive flux
            Dflux(iedge) = -f_ij
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))
            
            if (l_ji .lt. d_ij) then
              ! f_ij := minmod(f_ij,g_ij)
              if (f_ij*g_ij .le. 0.0_DP) then
                f_ij = 0.0_DP
              elseif (abs(g_ij) .lt. abs(f_ij)) then
                f_ij = g_ij
              end if
            end if

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) - f_ij
          end do
          !$omp end parallel do
        end if
        
      else

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))
            
            if (l_ji .lt. d_ij) then
              ! f_ij := minmod(f_ij,g_ij)
              if (f_ij*g_ij .le. 0.0_DP) then
                f_ij = 0.0_DP
              elseif (abs(g_ij) .lt. abs(f_ij)) then
                f_ij = g_ij
              end if
            end if

            ! Store raw antidiffusive flux
            Dflux(iedge) = dscale * f_ij
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j,d_ij,f_ij,g_ij,l_ji)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)

            ! Determine coefficients
            d_ij = DcoefficientsAtEdge(1,iedge)
            l_ji = DcoefficientsAtEdge(3,iedge)

            ! Determine fluxes
            f_ij = d_ij * (Dx(i)-Dx(j))
            g_ij = l_ji * (Dx(i)-Dx(j))
            
            if (l_ji .lt. d_ij) then
              ! f_ij := minmod(f_ij,g_ij)
              if (f_ij*g_ij .le. 0.0_DP) then
                f_ij = 0.0_DP
              elseif (abs(g_ij) .lt. abs(f_ij)) then
                f_ij = g_ij
              end if
            end if

            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge) + dscale*f_ij
          end do
          !$omp end parallel do
        end if
        
      end if

    end subroutine doFluxesUpwindBiasedDble
    
    !**************************************************************
    ! Assemble symmetric raw antidiffusive fluxes using the
    ! coefficients supplied by the edge-by-edge array DcoefficientsAtEdge
    
    subroutine doFluxesSymmetricDble(IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear
      
      real(DP), dimension(:), intent(inout) :: Dflux

            ! local variables
      integer :: iedge,i,j
      
      if (dscale .eq. 0.0_DP) then
        
        if (bclear) call lalg_clearVector(Dflux, NEDGE)
        
      elseif (dscale .eq. 1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if
        
      elseif (dscale .eq. -1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        end if
        
      else
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if
        
      end if

    end subroutine doFluxesSymmetricDble
    
    !**************************************************************
    ! Assemble nodal bounds using the coefficients supplied by the
    ! CSR-matrix Dmatrix and the precomputed bounds at edges

    subroutine doBoundsByMatrixDble(IedgeListIdx, IedgeList,&
        NEDGE, Dmatrix, DboundsAtEdge, Dq)
      
      real(DP), dimension(:,:), intent(in) :: DboundsAtEdge
      real(DP), dimension(:), intent(in) :: Dmatrix
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(out) :: Dq

      ! local variables
      integer :: i,iedge,igroup,ij,j,ji

      ! Clear bounds`s
      call lalg_clearVector(Dq)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      !$omp parallel default(shared) private(iedge,i,j,ij,ji)
      do igroup = 1, size(IedgeListIdx)-1
        
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers and matrix positions
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Update bounds
          Dq(i) = Dq(i) + DboundsAtEdge(1,iedge)*Dmatrix(ij)
          Dq(j) = Dq(j) + DboundsAtEdge(2,iedge)*Dmatrix(ji)
        end do
        !$omp end do
      end do
      !$omp end parallel
        
    end subroutine doBoundsByMatrixDble

    !**************************************************************
    ! Assemble nodal bounds using the coefficients supplied by the
    !  edge-by-edge array DcoefficientsAtEdge

    subroutine doBoundsByCoeffDble(IedgeListIdx, IedgeList,&
        NEDGE, DcoefficientsAtEdge, DboundsAtEdge, Dq)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:,:), intent(in) :: DboundsAtEdge
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(out) :: Dq

      ! local variables
      integer :: i,iedge,igroup,j

      ! Clear bounds`s
      call lalg_clearVector(Dq)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      !$omp parallel default(shared) private(iedge,i,j)
      do igroup = 1, size(IedgeListIdx)-1
        
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)

          ! Update bounds
          Dq(i) = Dq(i) + DboundsAtEdge(1,iedge)*DcoefficientsAtEdge(1,iedge)
          Dq(j) = Dq(j) + DboundsAtEdge(2,iedge)*DcoefficientsAtEdge(1,iedge)
        end do
        !$omp end do
      end do
      !$omp end parallel
      
    end subroutine doBoundsByCoeffDble

  end subroutine gfsc_buildFluxLPTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVecLPTBlock(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rmatrix)

!<description>
    ! This subroutine assembles the vector resulting from the
    ! application of linearyity-preserving flux correction. Note that
    ! this routine serves as a wrapper for block vectors. If there is
    ! only one block, then the corresponding scalar routine is
    ! called. Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
    
    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: lumped mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorBlock), intent(inout) :: ry    
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks .ne. 1 .or.&
        ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTBlock')
      call sys_halt()

    else

      call gfsc_buildVecLPTScalar(rafcstab, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1), rmatrix)
      
    end if

  end subroutine gfsc_buildVecLPTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildVecLPTScalar(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rmatrix)

!<description>
    ! This subroutine assembles the vector resulting from the
    ! application of linearyity-preserving flux correction.
    !
    ! This subroutine provides different algorithms:
    !
    ! Nonlinear/Linearised FEM-LPT
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Linearity-perserving flux correction for upwind-biased and
    ! symmetric antidiffusive fluxes.
    !
    ! The details of this method can be found in:
    !
    ! D. Kuzmin, Linearity-preserving flux correction and convergence
    ! acceleration for constrained Galerkin schemes, Ergebnisberichte
    ! Angew. Math. 430, Technische Universität Dortmund, 2011.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear
    
    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_LPTALGO_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: lumped mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! convective vector
    type(t_vectorScalar), intent(inout) :: ry    
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy,p_Dq
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Check if stabilisation is prepeared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
      call sys_halt()
    end if

    ! Clear convective vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)

    !---------------------------------------------------------------------------
    ! The FEM-LPT algorithm is split into the following steps which
    ! can be skipped and performed externally by the user:
    !
    ! 1) Initialise the edgewise correction factors (alpha).
    !
    ! 2) Compute the antidiffusive increments (Pp, Pm)
    !
    ! 3) Compute the local solution bounds (Qp, Qm).
    !
    ! 4) Compute the nodal correction factors (Rp, Rm).
    !
    ! 5) Compute edgewise correction factors (Alpha).
    !
    ! 6) Apply the limited antiddifusive fluxes
    !-------------------------------------------------------------------------
    
    if (iand(ioperationSpec, AFCSTAB_LPTALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 1) Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------

      ! Initialise alpha by unity
      call lalg_setVector(p_Dalpha, 1.0_DP)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_ADINCREMENTS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 2) Compute sums of antidiffusive increments 
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Compute sums of antidiffusive increments based on the
      ! raw-antidiffusive fluxes ...
      select case(rafcstab%climitingType)
      case (AFCSTAB_LIMITING_SYMMETRIC)
        ! ... in a symmetric fashion
        call doADIncrementsSymmDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
        
      case (AFCSTAB_LIMITING_UPWINDBIASED)
        ! ... in an upwind-biased fashion
        call doADIncrementsUpwDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)

      case default
        call output_line('Unsupported type of flux limiting!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end select

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 3) Compute local solution bounds
      !-------------------------------------------------------------------------
      
      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if
      
      ! Compute bounds
      call doBoundsDble(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEDGE, p_Dx, p_Dqp, p_Dqm)
      
      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 4) Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)   .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvVecFCTScalar')
        call sys_halt()
      end if

      ! Set additional pointer
      call lsyssc_getbase_double(rafcstab%p_rvectorQ, p_Dq)

      ! Compute nodal correction factors
      call doLimitNodalConstrainedDble(rafcstab%NEQ,&
          p_Dq, p_Dx, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      
      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_LPTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 5) Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides nodal correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provide nodal correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Compute edgewise correction factors ...
      select case(rafcstab%climitingType)
      case (AFCSTAB_LIMITING_SYMMETRIC)
        ! ... in a symmetric fashion
        call doLimitEdgewiseSymmDble(p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dpp, p_Dpm, p_Drp, p_Drm, p_Dalpha)
        
      case (AFCSTAB_LIMITING_UPWINDBIASED)
        ! ... in an upwind-biased fashion
        call doLimitEdgewiseUpwDble(p_IedgeList,&
            rafcstab%NEDGE, p_Dflux, p_Dpp, p_Dpm, p_Drp, p_Drm, p_Dalpha)
        
      case default
        call output_line('Unsupported type of flux limiting!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end select
      
      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
    end if
    

    if (iand(ioperationSpec, AFCSTAB_LPTALGO_CORRECT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 6) Correct antidiffusive fluxes and apply them
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edgewise correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provide edgewise correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_LPTALGO_SCALEBYMASS) .ne. 0) then
        if (present(rmatrix)) then
          call lsyssc_getbase_double(rmatrix, p_ML)

          call doCorrectScaleByMassDble(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, dscale, p_ML, p_Dalpha, p_Dflux, p_Dy)
        else
          call output_line('Lumped mass matrix is not provided!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildVecLPTScalar')
          call sys_halt()
        end if
      else

        call doCorrectDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, dscale, p_Dalpha, p_Dflux, p_Dy)
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes in a symmetric fashion
    
    subroutine doADIncrementsSymmDble(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, Dalpha, Dpp, Dpm)
      
      real(DP), dimension(:), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! The sums of positive/negative antidiffusive increments
      real(DP), dimension(:), intent(out) :: Dpp,Dpm
      
      ! local variables
      real(DP) :: f_ij,fp_ij,fm_ij
      integer :: i,iedge,igroup,j
      
      !$omp parallel default(shared) private(i,j,f_ij,fp_ij,fm_ij)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)
          
          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)
      
          ! Separate fluxes into positive/negative contributions
          fp_ij = max(0.0_DP,f_ij)
          fm_ij = min(0.0_DP,f_ij)
          
          ! Compute the sums of antidiffusive increments
          Dpp(i) = Dpp(i) + fp_ij   ! += max(0.0_DP, f_ij)
          Dpp(j) = Dpp(j) - fm_ij   ! += max(0.0_DP,-f_ij)
          Dpm(i) = Dpm(i) + fm_ij   ! += min(0.0_DP, f_ij)
          Dpm(j) = Dpm(j) - fp_ij   ! += min(0.0_DP,-f_ij)
        end do
        !$omp end do
      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsSymmDble
    
    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes in an upwind-biased fashion
    
    subroutine doADIncrementsUpwDble(IedgeListIdx, IedgeList,&
        NEDGE, Dflux, Dalpha, Dpp, Dpm)
      
      real(DP), dimension(:), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! The sums of positive/negative antidiffusive increments
      real(DP), dimension(:), intent(out) :: Dpp,Dpm
      
      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup
      
      !$omp parallel default(shared) private(i,f_ij)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node number of the upwind nodeee
          i = IedgeList(1,iedge)
          
          ! Apply multiplicative correction factor
          f_ij = Dalpha(iedge) * Dflux(iedge)
      
          ! Compute the sums of antidiffusive increments
          Dpp(i) = Dpp(i) + max(0.0_DP,f_ij)
          Dpm(i) = Dpm(i) + min(0.0_DP,f_ij)
        end do
        !$omp end do
        
      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsUpwDble

    !**************************************************************
    ! Assemble the local bounds from the predicted solution
    
    subroutine doBoundsDble(IedgeListIdx, IedgeList, NEDGE, Dx, Dqp, Dqm)
      
      real(DP), dimension(:), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      ! The local upper/lower bounds computed from Dx
      real(DP), dimension(:), intent(out) :: Dqp,Dqm
      
      ! local variables
      integer :: i,iedge,igroup,j

      
      !$omp parallel default(shared) private(i,j)

      ! Initialise Q`s by solution
      !$omp sections
      !$omp section
      call lalg_copyVector(Dx, Dqp)
      !$omp section
      call lalg_copyVector(Dx, Dqm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle
        
        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node numbers
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
          ! Compute local upper and lower bounds
          Dqp(i) = max(Dqp(i), Dx(j))
          Dqm(i) = min(Dqm(i), Dx(j))
          Dqp(j) = max(Dqp(j), Dx(i))
          Dqm(j) = min(Dqm(j), Dx(i))
        end do
        !$omp end do

      end do ! igroup     
      !$omp end parallel

    end subroutine doBoundsDble

    !**************************************************************
    ! Compute nodal correction factors with constraints
    
    subroutine doLimitNodalConstrainedDble(NEQ, Dq, Dx,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm)
      
      real(DP), dimension(:), intent(in) :: Dq,Dx
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm
      integer, intent(in) :: NEQ
      
      ! The nodal correction factors for positive/negative fluxes
      real(DP), dimension(:), intent(out) :: Drp,Drm
      
      ! local variables
      real(DP) :: diff
      integer :: ieq
      
      !$omp parallel sections default(shared) private(ieq,diff)
      
      !$omp section

      !$omp parallel do default(shared) private(ieq,diff)
      do ieq = 1, NEQ
        diff = Dq(ieq)*(Dqp(ieq)-Dx(ieq))
        if (Dpp(ieq) .gt. diff) then
          Drp(ieq) = diff/Dpp(ieq)
        else
          Drp(ieq) = 1.0_DP
        end if
      end do
      !$omp end parallel do

      !$omp section

      !$omp parallel do default(shared) private(ieq,diff)
      do ieq = 1, NEQ
        diff = Dq(ieq)*(Dqm(ieq)-Dx(ieq))
        if (Dpm(ieq) .lt. diff) then
          Drm(ieq) = diff/Dpm(ieq)
        else
          Drm(ieq) = 1.0_DP
        end if
      end do
      !$omp end parallel do

      !$omp end parallel sections

    end subroutine doLimitNodalConstrainedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes in
    ! a symmetric fashion, e.g., consider the nodal correction factors
    ! at both endpoints of the edge
    
    subroutine doLimitEdgewiseSymmDble(IedgeList,&
        NEDGE, Dflux, Dpp, Dpm, Drp, Drm, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dpp,Dpm
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      
      ! On input: the edge-wise correction factors from previous
      !           multiplicative correction steps
      ! On exit: the edge-wise correction factors resulting from
      !          the nodal correction factors Rp and Rm
      real(DP), dimension(:), intent(inout) :: Dalpha
      
      ! local variables
      real(DP) :: f_ij,r_ij
      integer :: iedge,i,j
      
      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,f_ij,r_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Get node numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Get precomputed raw antidiffusive fluxes
        f_ij = Dflux(iedge)

        ! Compute nodal correction factors
        if (f_ij .gt. AFCSTAB_EPSABS) then
          r_ij = min(Drp(i),Drm(j))
        elseif (f_ij .lt. -AFCSTAB_EPSABS) then
          r_ij = min(Drp(j),Drm(i))
        else
          r_ij = 1.0_DP
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do
      
    end subroutine doLimitEdgewiseSymmDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes in
    ! an upwind-biased fashion, e.g., consider the nodal correction
    ! factors only at the endpoint of the edge located upwind
    
    subroutine doLimitEdgewiseUpwDble(IedgeList,&
        NEDGE, Dflux, Dpp, Dpm, Drp, Drm, Dalpha)
      
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dpp,Dpm
      real(DP), dimension(:), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      
      ! On input: the edge-wise correction factors from previous
      !           multiplicative correction steps
      ! On exit: the edge-wise correction factors resulting from
      !          the nodal correction factors Rp and Rm
      real(DP), dimension(:), intent(inout) :: Dalpha
      
      ! local variables
      real(DP) :: f_ij,r_ij
      integer :: iedge,i
      
      ! Loop over all edges
      !$omp parallel do default(shared) private(i,f_ij,r_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Get node number of the upwind node
        i = IedgeList(1,iedge)
        
        ! Get precomputed raw antidiffusive fluxes
        f_ij = Dflux(iedge)

        ! Compute nodal correction factors
        if (f_ij .gt. AFCSTAB_EPSABS) then
          r_ij = Drp(i)
        elseif (f_ij .lt. -AFCSTAB_EPSABS) then
          r_ij = Drm(i)
        else
          r_ij = 1.0_DP
        end if

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * r_ij
      end do
      !$omp end parallel do
      
    end subroutine doLimitEdgewiseUpwDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    
    subroutine doCorrectDble(IedgeListIdx, IedgeList,&
        NEDGE, dscale, Dalpha, Dflux, Dy)
      
      real(DP), dimension(:), intent(in) :: Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dy
      
      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node numbers
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
          ! Correct antidiffusive flux
          f_ij = dscale * Dalpha(iedge) * Dflux(iedge)
          
          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij
          Dy(j) = Dy(j) - f_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel
      
    end subroutine doCorrectDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    
    subroutine doCorrectScaleByMassDble(IedgeListIdx, IedgeList,&
        NEDGE, dscale, ML, Dalpha, Dflux, Dy)
      
      real(DP), dimension(:), intent(in) :: ML,Dalpha,Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Dy
      
      ! local variables
      real(DP) :: f_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,f_ij)&
      !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1
          
          ! Get node numbers
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
          ! Correct antidiffusive flux
          f_ij = dscale * Dalpha(iedge) * Dflux(iedge)
          
          ! Apply limited antidiffusive fluxes
          Dy(i) = Dy(i) + f_ij/ML(i)
          Dy(j) = Dy(j) - f_ij/ML(j)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel
      
    end subroutine doCorrectScaleByMassDble

  end subroutine gfsc_buildVecLPTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvJacobianBlock(RcoeffMatrices, rx, fcb_calcMatrixSc_sim,&
      hstep, dscale, bbuildStabilisation, bclear, rjacobian, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective
    ! part of the discrete transport operator for a scalar convection
    ! equation.  Note that this routine serves as a wrapper for block
    ! vectors. If there is only one block, then the corresponding
    ! scalar routine is called.  Otherwise, an error is thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx
    
    ! perturbation parameter
    real(DP), intent(in) :: hstep

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! Check if block vector contains exactly one block
    if (rx%nblocks .ne. 1) then
      
      call output_line('Solution vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianBlock')
      call sys_halt()

    else
      
      call gfsc_buildConvJacobianScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrixSc_sim, hstep,&
          dscale, bbuildStabilisation, bclear, rjacobian, rcollection)
      
    end if

  end subroutine gfsc_buildConvJacobianBlock

   !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildConvJacobianScalar(RcoeffMatrices, rx, fcb_calcMatrixSc_sim,&
      hstep, dscale, bbuildStabilisation, bclear, rjacobian, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the convective part
    ! of the discrete transport operator for a scalar convection equation.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx
    
    ! perturbation parameter
    real(DP), intent(in) :: hstep

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for stabilisation
    ! TRUE  : perform stabilisation
    ! FALSE : perform no stabilisation
    logical, intent(in) :: bbuildStabilisation

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Jac,p_Dx
    integer :: h_Ksep,ndim
    
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)
    
    ! Set pointers
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      
    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_DcoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianScalar')
      call sys_halt()
    end select
    
    
    ! What kind of matrix are we?
    select case(rjacobian%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      
      ! Do we have to build the upwind Jacobian?
      if (bbuildStabilisation) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doUpwindMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doUpwindMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select

      else   ! bbuildStabilisation

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat7_1D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doGalerkinMat7_2D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doGalerkinMat7_3D(p_Kld, p_Kcol, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select

      end if   ! bbuildStabilisation

      ! Release diagonal separator
      call storage_free(h_Ksep)
      
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      
      ! Do we have to build the upwind Jacobian?
      if (bbuildStabilisation) then
        
        select case(ndim)
        case (NDIM1D)
          call doUpwindMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doUpwindMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doUpwindMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select
      
      else   ! bbuildStabilisation

        select case(ndim)
        case (NDIM1D)
          call doGalerkinMat9_1D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_Dx, p_Jac)
        case (NDIM2D)
          call doGalerkinMat9_2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_Dx, p_Jac)
        case (NDIM3D)
          call doGalerkinMat9_3D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              rjacobian%NEQ, p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Jac)
        end select

      end if   ! bbuildStabilisation

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildConvJacobianScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow
    
    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doGalerkinMat7_1D(Kld, Kcol, Ksep, NEQ, DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients k_ij and k_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij ,k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat7_1D


    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doGalerkinMat7_2D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients k_ij and k_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij ,k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat7_2D

    
    !**************************************************************
    ! Assemble standard Jacobian matrix for convective 
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doGalerkinMat7_3D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables     
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients k_ij and k_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat7_3D

    
    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doGalerkinMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      
      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat9_1D


    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doGalerkinMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat9_2D


    !**************************************************************
    ! Assemble standard Jacobian matrix for convective
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doGalerkinMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: k_ij,k_ji,a_ij,a_ji,b_ij,b_ji,d_ij,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+k_ji)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = k_ij; a_ji = k_ji
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+k_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-k_ij)/hstep
          a_ji = 0.5_DP*(a_ji-k_ji)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doGalerkinMat9_3D


    !**************************************************************
    ! Assemble upwind Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_1D(Kld, Kcol, Ksep, NEQ, DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      
      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat7_1D
    

    !**************************************************************
    ! Assemble upwind Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_2D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
                    
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat7_2D


    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 7

    subroutine doUpwindMat7_3D(Kld, Kcol, Ksep, NEQ, DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables     
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j

      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kld(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kld(j); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat7_3D
    
    
    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 1D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_1D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j

      
      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat9_1D


    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 2D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_2D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij = (a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat9_2D


    !**************************************************************
    ! Assemble Jacobian matrix for convective
    ! operator in 3D and assume zero row-sums.
    ! All matrices are stored in matrix format 9

    subroutine doUpwindMat9_3D(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        DcoeffX, DcoeffY, DcoeffZ, Dx, Jac)

      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: d_ij,l_ij,l_ji,a_ij,a_ji,b_ij,b_ji,diff
      integer :: ii,ij,ji,jj,i,j
      

      ! Loop over all rows I of Jacobian matrix
      do i = 1, NEQ
        
        ! Get position of diagonal entry II
        ii = Kdiagonal(i)
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get row number J, the corresponding matrix position JI, and
          ! let the separator point to the next entry
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          ! Now, we have the global position of the matrix entries IJ
          ! and JI (!!!) as well as the numbers I and J for which I < J.
          ! Next, we need to consider all matrix entries of rows I and
          ! J and perturb the matrix coefficients l_ij(u) and l_ji(u)
          ! by +/-h*e_k, whebery K stands for the column number of the
          ! current matrix entry. For the computation of the low-order
          ! coefficient l_ij we have to consider k_ij and k_ji in order
          ! to determine the artificial diffusion coefficient d_ij.
          ! However, the coefficients a_ij and a_ji resulting from a
          ! divided difference approximation of the derivatives of the
          ! transport operator need to be handled separately.
          
          ! Due to the fact, that we need the unperturbed quantities
          ! quite frequently, we store them in local auxiliary variables
          
          ! Compute solution difference Dx_j-Dx_i
          diff = Dx(j)-Dx(i)

          ! Compute coefficients
          C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
          C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
          C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

          ! We have to loop over all columns K of the I-th and J-th row
          ! of the Jacobian matrix and update th positions IK and JK,
          ! respectively. The perturbation +/-h*e_k only influences the
          ! coefficients a_ij^k which s defined as
          !   a_ji^k:=\frac{l_ij(u+h*e_k)-l_ij(u-h*e_k)}{2h}
          ! if either I=K or J=K. In short the loop over the I-th and 
          ! J-th row only affects the matrix position II, IJ, JI, and JJ
          ! which are known a priori(!!)

          ! (1) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_I
          
!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_I" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients
          b_ji = (a_ji+l_ji+d_ij)/2._DP
          Jac(ji) = Jac(ji)+b_ji
          Jac(jj) = Jac(jj)-b_ji
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the I-th column of the I-th and J-th row
          Jac(ii) = Jac(ii)+a_ij*diff
          Jac(ji) = Jac(ji)-a_ji*diff

          
          ! (2) Update Jac(II,IJ,JI,JJ) for perturbation +/-h*e_J

!!$          ! Compute perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply perturbed coefficient to a_ij and a_ji
          a_ij = l_ij+d_ij; a_ji = l_ji+d_ij
          
!!$          ! Compute "-h*e_J" perturbed coefficients l_ij and l_ji
!!$          call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$              C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
          
          ! Apply the average of the perturbed coefficients for J=K
          b_ij =(a_ij+l_ij+d_ij)/2._DP
          Jac(ij) = Jac(ij)+b_ij
          Jac(ii) = Jac(ii)-b_ij
          
          ! Compute final coefficients a_ij and a_ji as the second
          ! order divided differences of the low-order coefficients
          a_ij = 0.5_DP*(a_ij-l_ij-d_ij)/hstep
          a_ji = 0.5_DP*(a_ji-l_ji-d_ij)/hstep
          
          ! Update the K-th column of the I-th row, that is, the
          ! entriy IK of the Jacobian matrix
          Jac(ij) = Jac(ij)+a_ij*diff
          Jac(jj) = Jac(jj)-a_ji*diff
        end do
      end do
    end subroutine doUpwindMat9_3D

  end subroutine gfsc_buildConvJacobianScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearFCTBlock(rx, theta, tstep, hstep,&
      bclear, rafcstab, rjacobian, rmatrix)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  Note that the velocity is assumed
    ! to be linear.  Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the
    ! corresponding scalar routine is called.  Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearFCTBlock')
      call sys_halt()

    else

      call gfsc_buildJacLinearFCTScalar(&
          rx%RvectorBlock(1), theta, tstep, hstep, bclear,&
          rafcstab, rjacobian, rmatrix)

    end if
  end subroutine gfsc_buildJacLinearFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearFCTScalar(rx, theta, tstep, hstep,&
      bclear, rafcstab, rjacobian, rmatrix)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  Note that the velocity is assumed
    ! to be linear.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dflux,p_Dflux0
    real(DP), dimension(:), pointer :: p_MC,p_Jac,p_Dx
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearFCTScalar')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)
      
    case (AFCSTAB_NLINFCT_IMPLICIT)
      
      ! What kind of matrix are we?
      select case(rjacobian%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        !-------------------------------------------------------------------------
        ! Matrix format 7
        !-------------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobian, p_Kld)
        
        if (present(rmatrix)) then
          call lsyssc_getbase_double(rmatrix, p_MC)
          call doJacobian_implFCTconsMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kld, p_MC, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kld, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        end if
        
      case(LSYSSC_MATRIX9)
        !-------------------------------------------------------------------------
        ! Matrix format 9
        !-------------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobian, p_Kld)
        call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
        
        if (present(rmatrix)) then
          call lsyssc_getbase_double(rmatrix, p_MC)
          call doJacobian_implFCTconsMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal, p_MC, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        else
          call doJacobian_implFCTnoMass(&
              p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal, p_Dx,&
              p_Dflux, p_Dflux0, theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
        end if

        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearFCTScalar')
        call sys_halt()
      end select
      
    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearFCTScalar')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow

    !**************************************************************
    ! Assemble the Jacobian matrix for semi-implicit FEM-FCT,
    ! whereby no mass antidiffusion is built into the matrix

    subroutine doJacobian_implFCTnoMass(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Jac
      
      ! local variables
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        a_ij = theta*d_ij
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Compute perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        ! Compute limited antidiffusive flux f(Dx_ij+h*e_i) 
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        ! Compute limited antidiffusive flux f(Dx_ij-h*e_j)
        f_j = a_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)+f_ij
        Jac(jj) = Jac(jj)-f_ij
      end do

    end subroutine doJacobian_implFCTnoMass
    
    
    !**************************************************************
    ! Assemble the Jacobian matrix for semi-implicit FEM-FCT,
    ! whereby consistent mass antidiffusion is built into the matrix

    subroutine doJacobian_implFCTconsMass(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, MC, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dflux,Dflux0
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Jac
      
      ! local variables
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j
      
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Determine coefficients
        d_ij = DcoefficientsAtEdge(1,iedge)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Compute perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        ! Compute limited antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        ! Compute limited antidiffusive flux f(Dx_ij-h*e_j)
        f_j = a_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)+f_ij
        Jac(jj) = Jac(jj)-f_ij
      end do

    end subroutine doJacobian_implFCTconsMass

  end subroutine gfsc_buildJacLinearFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearTVDBlock(rx, tstep, hstep,&
      bclear, rafcstab, rjacobian, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  Note that the velocity is assumed
    ! to be linear.  Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the
    ! corresponding scalar routine is called.  Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearTVDBlock')
      call sys_halt()

    else

      call gfsc_buildJacLinearTVDScalar(rx%RvectorBlock(1), tstep,&
          hstep, bclear, rafcstab, rjacobian, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacLinearTVDBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearTVDScalar(rx, tstep, hstep,&
      bclear, rafcstab, rjacobian, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm
    real(DP), dimension(:), pointer :: p_Jac,p_Dx,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    
    integer :: h_Ksep
    logical :: bisExtended


    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearTVDScalar')
      call sys_halt()
    end if
      
    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_genOffdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if


    ! What kind of matrix format are we?
    select case(rjacobian%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      
      ! Create diagonal separator and increase it by one
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_TVD(&
          p_IsuperdiagEdgesIdx, p_IedgeList,&
          p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
          p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
          p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
          tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
          rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      
      ! Free storage
      call storage_free(h_Ksep)
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      
      call doJacobianMat79_TVD(&
          p_IsuperdiagEdgesIdx, p_IedgeList,&
          p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
          p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
          p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
          tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
          rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
        
      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearTVDScalar')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    pure subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,isep,l,iloc
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k)+1, Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat7
    
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    pure subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k
      
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,isep,l,iloc
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k), Kld(k+1)-1

        ! Get the column number
        l = Kcol(ild)

        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat9


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
        Dx, Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Drploc,Drmloc,Dfluxloc
      integer, dimension(NNVEDGE) :: Kloc
      integer :: k,l,ild,iedge,iloc,nloc
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ

        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)

          ! Increase local counter
          iloc = iloc+1
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              IedgeList, DcoefficientsAtEdge, Dx,&
              Dpp, Dpm, Dqp, Dqm, tstep, hstep, iedge, iloc, k,&
              Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              IedgeList, DcoefficientsAtEdge, Dx,&
              Dpp, Dpm, Dqp, Dqm, tstep, hstep, iedge, iloc, k,&
              Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(iloc)

          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_TVD


    !**************************************************************
    ! Update the local coefficients for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine updateJacobianMat79_TVD(IedgeList,&
        DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        iedge, iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: iedge,k,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(inout) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Dfluxloc
      integer, dimension(:), intent(inout) :: Kloc

      ! local variables
      real(DP) :: d_ij,f_ij,l_ij,l_ji,diff,dsign
      integer :: i,j,iperturb
      
      
      ! Determine indices. Obviously, either i or j must be equal to k.
      ! Otherwise, the edge ij would not be present in the list of 
      ! incident edges for node k.
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)

      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Determine prelimited antidiffusive flux
      diff = tstep*(Dx(i)-Dx(j))
      f_ij = min(d_ij, l_ji)*diff

      !-------------------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     eliminate the contribution of the edge IJ for the
      !     unperturbed solution values Dx_i and Dx_j.
      !
      ! (2) perturbed values: The local Ps and Qs require the 
      !     contribution of the perturbed solution values u +/- h*e_k,
      !     whereby e_k denotes the k-th unit vector and h stands
      !     for the  perturbation step length.
      !-------------------------------------------------------------------------
      
      ! Which is the upwind node?
      if (i .eq. k) then

        ! Store global node number of the opposite node
        Kloc(iloc) = j

        ! Update nodal coefficients for vertex j (!) which is the downwind node
        Dpploc(:,iloc) = Dpp(j)
        Dpmloc(:,iloc) = Dpm(j)
        Dqploc(:,iloc) = Dqp(j)-max(0.0_DP, f_ij)
        Dqmloc(:,iloc) = Dqm(j)-min(0.0_DP, f_ij)

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Compute perturbed antidiffusive flux
          f_ij = min(d_ij,l_ji)*(diff+tstep*dsign*hstep)
          Dfluxloc(iperturb,iloc) = f_ij

          ! For node k which is the upwind node
          Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP, f_ij)
          Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP, f_ij)
          Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP,-f_ij)
          Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP,-f_ij)
          
          ! For node l opposite to k which is the downwind node
          Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, f_ij)
          Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, f_ij)
        end do

      else

        ! Store global node number of the opposite node
        Kloc(iloc) = i
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        Dpploc(:,iloc) = Dpp(i)-max(0.0_DP, f_ij)
        Dpmloc(:,iloc) = Dpm(i)-min(0.0_DP, f_ij)
        Dqploc(:,iloc) = Dqp(i)-max(0.0_DP,-f_ij)
        Dqmloc(:,iloc) = Dqm(i)-min(0.0_DP,-f_ij)

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Compute perturbed antidiffusive flux
          f_ij = min(d_ij,l_ji)*(diff-tstep*dsign*hstep)
          Dfluxloc(iperturb,iloc) = f_ij
          
          ! For node k which is the downwind node
          Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP, f_ij)
          Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP, f_ij)
          
          ! For node l opposite to k
          Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, f_ij)
          Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP,-f_ij)
          Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP,-f_ij)
        end do
      end if
    end subroutine updateJacobianMat79_TVD


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine assembleJacobianMat79_TVD(IedgeList, Kdiagonal,&
        Dflux, Kloc, Drploc, Drmloc, Dfluxloc, hstep, iedge, iloc, k, l,&
        bisExtended, Ksep, Jac)

      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal,Kloc
      integer, intent(in) :: iedge,k,l,iloc
      logical, intent(in) :: bisExtended

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      

      ! local variables
      real(DP) :: f_ij,df_ij
      integer :: ik,jk,i,j,m,iperturb
      
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      if (m .eq. k) then

        !-----------------------------------------------------------------------
        ! 1. Case: primary edge
        !-----------------------------------------------------------------------
        ! The current edge connects the perturbed node k with its direct
        ! neighbor l. Hence, all required information can be extracted from 
        ! the local arrays and no global data retrieval has to be performed.

        ! Initilaize flux difference
        df_ij = 0.0_DP

        ! Which node is located upwind?
        if (i .eq. k) then
          
          do iperturb = 1, 2
          
            ! Retrieve precomputed flux
            f_ij = Dfluxloc(iperturb,iloc)

            ! Limit flux 
            if (f_ij > 0.0_DP) then
              f_ij = Drploc(iperturb,0)*f_ij
            else
              f_ij = Drmloc(iperturb,0)*f_ij
            end if

            ! Adopt sign for perturbation direction
            df_ij = df_ij-(iperturb-1.5_DP)*f_ij/hstep
          end do
          
          ! Get corresponding matrix indices
          ik = Kdiagonal(i); jk = Ksep(j)
        else

          do iperturb = 1, 2
            
            ! Retrieve precomputed flux
            f_ij = Dfluxloc(iperturb,iloc)

            ! Limit flux
            if (f_ij > 0.0_DP) then
              f_ij = Drploc(iperturb,iloc)*f_ij
            else
              f_ij = Drmloc(iperturb,iloc)*f_ij
            end if
            
            ! Adopt sign for perturbation direction
            df_ij = df_ij-(iperturb-1.5_DP)*f_ij/hstep
          end do

          ! Get corresponding matrix indices
          jk = Kdiagonal(j); ik = Ksep(i)
        end if

        ! Apply perturbed antidiffusive contribution
        Jac(ik) = Jac(ik)-df_ij
        Jac(jk) = Jac(jk)+df_ij
                
      elseif (bisExtended) then
       
        !-----------------------------------------------------------------------
        ! 2. Case: secondary edge
        !-----------------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are not equal
        ! to the perturbed vertex k. Thus, the influence of the solution
        ! perturbation can only be due to a change in the correction factors
        ! alpha_ij. Moreover, for upwind-biased flux limiting techniques only
        ! the nodal correction factors for the upwind node i is used. Hence, it
        ! suffices to check if node i corresponds to the modified vertex l.
        ! Interestingly enough, some edge LM which connects two direct neighbors
        ! of the perturbed vertex k is only processed once due to the fact that
        ! either l or (!) m corresponds to the upwind node.

        if (i .eq. l) then

          if (Dflux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*(Drploc(1,iloc)-Drploc(2,iloc))*Dflux(iedge)/hstep
          else
            f_ij = 0.5_DP*(Drmloc(1,iloc)-Drmloc(2,iloc))*Dflux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end if
      end if
    end subroutine assembleJacobianMat79_TVD

  end subroutine gfsc_buildJacLinearTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearGPBlock(rmatrix, rx,&
      rx0, theta, tstep, hstep, bclear, rafcstab, rjacobian,&
      bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  Note that the velocity is assumed
    ! to be linear.  Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the
    ! corresponding scalar routine is called.  Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1 .or.&
        rx0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearGPBlock')
      call sys_halt()

    else

      call gfsc_buildJacLinearGPScalar(&
          rmatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), theta, tstep, hstep,&
          bclear, rafcstab, rjacobian, bextendedSparsity)
      
    end if
  end subroutine gfsc_buildJacLinearGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacLinearGPScalar(rmatrix, rx,&
      rx0, theta, tstep, hstep, bclear, rafcstab, rjacobian,&
      bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation part
    ! of the discrete transport operator for a scalar convection equation.
    ! Note that the velocity is assumed to be linear.
!</description>

!<input>
    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! initial solution vector
    type(t_vectorScalar), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_MC,p_Jac,p_Dx,p_Dx0,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep
    logical :: bisExtended

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearGPScalar')
      call sys_halt()
    end if
    
    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_genOffdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rmatrix, p_MC)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(rx0, p_Dx0)

    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if


    ! What kind of matrix format are we?
    select case(rjacobian%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_GP(&
          p_IsuperdiagEdgesIdx, p_IedgeList, p_IsubdiagEdgesIdx,&
          p_IsubdiagEdges, p_DcoefficientsAtEdge, p_Kld, p_Kcol,&
          p_Kld, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta, tstep, hstep,&
          rafcstab%NEQ, rafcstab%NEDGE, rafcstab%NNVEDGE,&
          bisExtended, .true., p_Ksep, p_Jac)
      
      ! Free storage
      call storage_free(h_Ksep)
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      
      call doJacobianMat79_GP(&
          p_IsuperdiagEdgesIdx, p_IedgeList, p_IsubdiagEdgesIdx,&
          p_IsubdiagEdges, p_DcoefficientsAtEdge, p_Kld, p_Kcol,&
          p_Kdiagonal, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm,&
          theta, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
          rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      
      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacLinearGPScalar')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    pure subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,isep,l,iloc
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k)+1, Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat7

    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    pure subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k
      
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,isep,l,iloc
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k), Kld(k+1)-1

        ! Get the column number
        l = Kcol(ild)

        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat9
    
    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_GP(IsuperdiagEdgesIdx, IedgeList,&
        IsubdiagEdgesIdx, IsubdiagEdges, DcoefficientsAtEdge, Kld,&
        Kcol, Kdiagonal, MC, Dx, Dx0, Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp,&
        Drm, theta, tstep, hstep, NEQ, NEDGE, NNVEDGE, bisExtended,&
        bisMat7, Ksep, Jac)
    
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0,Dflux,Dflux0,Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc
      real(DP), dimension(2,0:NNVEDGE) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      integer, dimension(NNVEDGE) :: Kloc
      integer :: k,l,ild,iedge,iloc,nloc
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              IedgeList, DcoefficientsAtEdge, MC, Dx, Dx0,&
              Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, theta, tstep,&
              hstep, iedge, iloc, k, Dpploc, Dpmloc, Dqploc,&
              Dqmloc, Dfluxloc, Dfluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              IedgeList, DcoefficientsAtEdge, MC, Dx, Dx0,&
              Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, theta, tstep,&
              hstep, iedge, iloc, k, Dpploc, Dpmloc, Dqploc,&
              Dqmloc, Dfluxloc, Dfluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux, Dflux0, Drp, Drm,&
                Kloc, Drploc, Drmloc, Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux, Dflux0, Drp, Drm,&
                Kloc, Drploc, Drmloc, Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_GP

    
    !**************************************************************
    ! Update the local coefficients for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine updateJacobianMat79_GP(IedgeList,&
        DcoefficientsAtEdge, MC, Dx, Dx0, Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm,&
        theta, tstep, hstep, iedge, iloc, k, Dpploc, Dpmloc,&
        Dqploc, Dqmloc, Dfluxloc, Dfluxloc0, Kloc)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0,Dflux,Dflux0,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: iedge,k,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(inout) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Dfluxloc,Dfluxloc0
      integer, dimension(:), intent(inout)    :: Kloc

      ! local variables
      real(DP) :: m_ij,d_ij,df_ij,f_ij,l_ij,l_ji,p_ij,pf_ij,q_ij,q_ji,diff,diff1,diff0,dsign
      integer :: i,j,ij,iperturb
      
      
      ! Determine indices. Obviously, either i or j must be equal
      ! to k. Otherwise, the edge ij would not be present in the
      ! list of incident edges for node k.
      i  = IedgeList(1,iedge)
      j  = IedgeList(2,iedge)
      ij = IedgeList(3,iedge)

      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Include consistent mass matrix
      m_ij = MC(ij)
      q_ij = m_ij/tstep+l_ij
      q_ji = m_ij/tstep+l_ji

      ! Determine solution differences
      diff1 = Dx(i)-Dx(j)
      diff0 = Dx0(i)-Dx0(j)

      ! Determine total solution difference
      diff = tstep*(theta*diff1+(1.0_DP-theta)*diff0)

      ! Compute antidiffusive flux 
      if (abs(diff) < AFCSTAB_EPSABS) then
        p_ij = 0.0_DP
        f_ij = 0.0_DP
      else
        p_ij = max(0.0_DP, m_ij*(diff1-diff0)/diff+d_ij)
        f_ij = p_ij*diff
      end if

      ! Prelimit the antidiffusive flux
      pf_ij = min(p_ij, l_ji)*diff
      
      ! Compute the remaining flux
      df_ij = f_ij-pf_ij

      !-------------------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     eliminate the contribution of the edge IJ for the
      !     unperturbed solution values Dx_i and Dx_j.
      !
      ! (2) perturbed values: The local Ps and Qs require the 
      !     contribution of the perturbed solution values u +/- h*e_k,
      !     whereby e_k denotes the k-th unit vector and h stands
      !     for the  perturbation step length.
      !-------------------------------------------------------------------------

      ! Which is the upwind node?
      if (i .eq. k) then

        ! Store global node number of the opposite node
        Kloc(iloc) = j

        ! Update nodal coefficients for vertex j (!) which is the downwind node
        Dpploc(:,iloc) = Dpp(j)-max(0.0_DP,-df_ij)
        Dpmloc(:,iloc) = Dpm(j)-min(0.0_DP,-df_ij)
        Dqploc(:,iloc) = Dqp(j)-max(0.0_DP, diff)*q_ji
        Dqmloc(:,iloc) = Dqm(j)-min(0.0_DP, diff)*q_ji

        do iperturb = 1, 2
          
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Update solution difference
          diff1 = Dx(i)-Dx(j)+dsign*hstep
          
          ! Update total solution difference
          diff = tstep*(theta*diff1+(1.0_DP-theta)*diff0)
          
          ! Compute antidiffusive flux
          if (abs(diff) < AFCSTAB_EPSABS) then
            p_ij = 0.0_DP
            f_ij = 0.0_DP
          else
            p_ij = max(0.0_DP, m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = p_ij*diff
          end if
          
          ! Prelimit the antidiffusive flux
          pf_ij = min(p_ij, l_ji)*diff
          Dfluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          Dfluxloc(iperturb,iloc) = df_ij
          
          ! For node k which is the upwind node
          Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP, f_ij)
          Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP, f_ij)
          Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP,-diff)*q_ij
          Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP,-diff)*q_ij
          
          ! For node l opposite to k which is the downwind node
          Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP,-df_ij)
          Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP,-df_ij)
          Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, diff)*q_ji
          Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, diff)*q_ji
        end do

      else

        ! Store global node number of the opposite node
        Kloc(iloc) = i

        ! Update nodal coefficients for vertex i (!) which is the upwind node
        Dpploc(:,iloc) = Dpp(i)-max(0.0_DP, f_ij)
        Dpmloc(:,iloc) = Dpm(i)-min(0.0_DP, f_ij)
        Dqploc(:,iloc) = Dqp(i)-max(0.0_DP,-diff)*q_ij
        Dqmloc(:,iloc) = Dqm(i)-min(0.0_DP,-diff)*q_ij

        do iperturb = 1, 2
        
          ! Compute correct sign of perturbation
          dsign = 3-2*iperturb

          ! Update solution difference
          diff1 = Dx(i)-Dx(j)-dsign*hstep
          
          ! Update total solution difference
          diff = tstep*(theta*diff1+(1.0_DP-theta)*diff0)
          
          ! Compute antidiffusive flux
          if (abs(diff) < AFCSTAB_EPSABS) then
            p_ij = 0.0_DP
            f_ij = 0.0_DP
          else
            p_ij = max(0.0_DP, m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = p_ij*diff
          end if
          
          ! Prelimit the antidiffusive flux
          pf_ij = min(p_ij, l_ji)*diff
          Dfluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          Dfluxloc(iperturb,iloc) = df_ij

          ! For node k which is the downwind node
          Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP,-df_ij)
          Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP,-df_ij)
          Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP, diff)*q_ji
          Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP, diff)*q_ji
          
          ! For node l opposite to k
          Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, f_ij)
          Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP,-diff)*q_ij
          Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP,-diff)*q_ij
        end do
      end if
      
    end subroutine updateJacobianMat79_GP
    

    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine assembleJacobianMat79_GP(IedgeList, Kdiagonal,&
        Dflux, Dflux0, Drp, Drm, Kloc, Drploc, Drmloc, Dfluxloc, Dfluxloc0,&
        hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      real(DP), dimension(:), intent(in) :: Drp,Drm,Dflux,Dflux0
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal,Kloc
      integer, intent(in) :: iedge,iloc,k,l
      logical, intent(in) :: bisExtended

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP) :: f_ij,pf_ij,df_ij
      integer :: ik,jk,i,j,m,iperturb
      
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i=IedgeList(1,iedge)
      j=IedgeList(2,iedge)
      m=(i+j)-l
      
      ! We need to find out, which kind of edge is processed
      if (m .eq. k) then

        !-----------------------------------------------------------------------
        ! 1. Case: primary edge
        !-----------------------------------------------------------------------
        ! The current edge connects the perturbed node k with its direct
        ! neighbor l. Hence, all required information can be extracted from 
        ! the local arrays and no global data retrieval has to be performed.
                
        do iperturb = 1, 2
          
          ! Retrieve precomputed fluxes
          df_ij = Dfluxloc(iperturb,iloc)
          pf_ij = Dfluxloc0(iperturb,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = Drploc(iperturb,0)*pf_ij
            else
              pf_ij = Drmloc(iperturb,0)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(Drploc(iperturb,0), Drmloc(iperturb,iloc))*df_ij
            else
              df_ij = min(Drmloc(iperturb,0), Drploc(iperturb,iloc))*df_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = Drploc(iperturb,iloc)*pf_ij
            else
              pf_ij = Drmloc(iperturb,iloc)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(Drmloc(iperturb,0), Drploc(iperturb,iloc))*df_ij
            else
              df_ij = min(Drploc(iperturb,0), Drmloc(iperturb,iloc))*df_ij
            end if
            
          end if
          
          ! Combine both contributions and 
          ! adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*(pf_ij+df_ij)/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end do
        
      elseif (bisExtended) then
        
        !-----------------------------------------------------------------------
        ! 2. Case: secondary edge
        !-----------------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are not equal
        ! to the perturbed vertex k. Thus, the influence of the solution
        ! perturbation can only be due to a change in the correction factors
        ! alpha_ij. Moreover, for upwind-biased flux limiting techniques only
        ! the nodal correction factors for the upwind node i is used. Hence, it
        ! suffices to check if node i corresponds to the modified vertex l.
        ! Interestingly enough, some edge LM which connects two direct neighbors
        ! of the perturbed vertex k is only processed once due to the fact that
        ! either l or (!) m corresponds to the upwind node.

        if (i .eq. l) then

          ! Get precomputed fluxes
          pf_ij = Dflux0(iedge)
          df_ij = Dflux(iedge)

          ! Limit upwind contribution
          if (pf_ij > 0.0_DP) then
            pf_ij = (Drploc(1,iloc)-Drploc(2,iloc))*pf_ij
          else
            pf_ij = (Drmloc(1,iloc)-Drmloc(2,iloc))*pf_ij
          end if

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(Drploc(1,iloc), Drm(j))-&
                     min(Drploc(2,iloc), Drm(j)))*df_ij
          else
            df_ij = (min(Drmloc(1,iloc), Drp(j))-&
                     min(Drmloc(2,iloc), Drp(j)))*df_ij
          end if

          ! Combine both contributions
          f_ij = 0.5_DP*(pf_ij+df_ij)/hstep

          ! Get corresponding matrix indices
          ik=Ksep(i); jk=Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        else

          ! Get precomputed flux (only symmetric part)
          df_ij = Dflux(iedge)

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(Drp(i), Drmloc(1,iloc))-&
                     min(Drp(i), Drmloc(2,iloc)))*df_ij
          else
            df_ij = (min(Drm(i), Drploc(1,iloc))-&
                     min(Drm(i), Drploc(2,iloc)))*df_ij
          end if

          ! Compute divided difference
          f_ij = 0.5_DP*df_ij/hstep

          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)
          
          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
          
        end if
      end if
    end subroutine assembleJacobianMat79_GP
  end subroutine gfsc_buildJacLinearGPScalar

  !*****************************************************************************
  
!<subroutine>

  subroutine gfsc_buildJacobianFCTBlock(RcoeffMatrices, rx,&
      fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, rmatrix, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  The velocity is assumed to be
    ! nonlinear/arbitrary.  Note that this routine serves as a wrapper
    ! for block vectors. If there is only one block, then the
    ! corresponding scalar routine is called.  Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices
    
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTBlock')
      call sys_halt()

    else

      call gfsc_buildJacobianFCTScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrixSc_sim,&
          theta, tstep, hstep, bclear, rafcstab, rjacobian,&
          rmatrix, rcollection)

    end if
  end subroutine gfsc_buildJacobianFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianFCTScalar(RcoeffMatrices, rx,&
      fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, rmatrix, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  The velocity is assumed to be
    ! nonlinear/arbitrary.  This routine will also work for linear
    ! velocities but then it is inefficient since the solution
    ! perturbation does not affect the velocity.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear
    
    ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'
    
    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dflux,p_Dflux0,p_Dx,p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_MC,p_Jac
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer :: ndim
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTScalar')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_DcoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTScalar')
      call sys_halt()
    end select
    
    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)
      
    case (AFCSTAB_NLINFCT_IMPLICIT)
      
      ! What kind of matrix format are we?
      select case(rjacobian%cmatrixFormat)
      case(LSYSSC_MATRIX7)
        !-----------------------------------------------------------------------
        ! Matrix format 7
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kld(rjacobian, p_Kld)
        
        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM2D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE,  p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM3D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kld,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
        end select
        
        
      case(LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 9
        !-----------------------------------------------------------------------

        ! Set pointers
        call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
        
        ! How many dimensions do we have?
        select case(ndim)
        case (NDIM1D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_1D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
            
        case (NDIM2D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_2D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
          
        case (NDIM3D)
          if (present(rmatrix)) then
            call doJacobian_implFCTconsMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_MC, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          else
            call doJacobian_implFCTnoMass_3D(&
                p_IedgeList, p_DcoefficientsAtEdge, p_Kdiagonal,&
                p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dflux0,&
                theta, tstep, hstep, rafcstab%NEDGE, p_Jac)
          end if
        end select
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTScalar')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianFCTScalar')
      call sys_halt()
    end select

  contains
    
    ! Here, the working routine follow

    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 1D,
    ! whereby no mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTnoMass_1D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j
      
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_i+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTnoMass_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 1D,
    ! whereby consistent mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTconsMass_1D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, MC, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM1D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j
      
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j) 
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_i+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTconsMass_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 2D,
    ! whereby no mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTnoMass_2D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in)   :: Kdiagonal
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j
      
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_i+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTnoMass_2D

    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 2D,
    ! whereby consistent mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTconsMass_2D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, MC, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in)   :: Kdiagonal
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j
      
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j) 
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_i+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do

    end subroutine doJacobian_implFCTconsMass_2D
    
    
    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 3D
    ! whereby no mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTnoMass_3D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, DcoeffZ, Dx, Dflux, Dflux0,&
        theta, tstep, hstep, NEDGE, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE

      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j

      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)

        ! Compute solution difference
        diff = Dx(i)-Dx(j)

        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep


        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if


!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if


        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep

        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij

        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------

!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = theta*d_ij

        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j)
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_i+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
      
    end subroutine doJacobian_implFCTnoMass_3D
    

    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-FCT in 3D
    ! whereby consistent mass antidiffusion is built into the Jacobian.
    ! All matrices can be stored in matrix format 7 or 9
    subroutine doJacobian_implFCTconsMass_3D(IedgeList,&
        DcoefficientsAtEdge, Kdiagonal, DcoeffX, DcoeffY, DcoeffZ, MC, Dx, Dflux,&
        Dflux0, theta, tstep, hstep, NEDGE, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,MC,Dx,Dflux,Dflux0
      real(DP), intent(in) :: theta,tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: NEDGE
      
      real(DP), dimension(:), intent(inout) :: Jac

      ! local variables
      real(DP), dimension(NDIM3D) :: C_ij,C_ji
      real(DP) :: f_i,f_j,f_ij,d_ij,a_ij,b_ij,l_ij,l_ji,diff,diff_i,diff_j
      integer :: iedge,ij,ji,ii,jj,i,j
      
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Determine vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Determine matrix indices
        ij = IedgeList(3,iedge)
        ji = IedgeList(4,iedge)
        
        ! Determine diagonal indices
        ii = Kdiagonal(i); jj = Kdiagonal(j)
        
        ! Compute coefficients
        C_ij(1) = DcoeffX(ij); C_ji(1) = DcoeffX(ji)
        C_ij(2) = DcoeffY(ij); C_ji(2) = DcoeffY(ji)
        C_ij(3) = DcoeffZ(ij); C_ji(3) = DcoeffZ(ji)
        
        ! Compute solution difference
        diff = Dx(i)-Dx(j)
        
        ! Determine perturbed solution differences
        diff_i = diff+hstep
        diff_j = diff-hstep
        
        !------------------------------------------------------------
        ! Compute flux for i-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_i)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_i)
        f_i = a_ij*diff_i+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)-hstep, Dx(j),&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_i)
        b_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_j+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply i-th column
        Jac(ii) = Jac(ii)-f_ij
        Jac(ji) = Jac(ji)+f_ij
        
        !------------------------------------------------------------
        ! Compute flux for j-th column
        !------------------------------------------------------------
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)+hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient a_ij(u+hstep*e_j)
        a_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij+h*e_j) 
        f_i = a_ij*diff_j+Dflux0(iedge)
        if (f_i > 0.0_DP) then
          f_i = min(f_i, max(Dflux(iedge), 0.0_DP))
        else
          f_i = max(f_i, min(Dflux(iedge), 0.0_DP))
        end if
        
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i), Dx(j)-hstep,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Compute perturbed coefficient b_ij(u-hstep*e_j)
        b_ij = MC(ij)/tstep+theta*d_ij
        
        ! Compute and limit raw antidiffusive flux f(Dx_ij-h*e_j)
        f_j = b_ij*diff_i+Dflux0(iedge)
        if (f_j > 0.0_DP) then
          f_j = min(f_j, max(Dflux(iedge), 0.0_DP))
        else
          f_j = max(f_j, min(Dflux(iedge), 0.0_DP))
        end if
        
        
        ! Compute divided differences of fluxes
        f_ij = 0.5_DP*tstep*(f_i-f_j)/hstep
        
        ! Apply j-th column
        Jac(ij) = Jac(ij)-f_ij
        Jac(jj) = Jac(jj)+f_ij
      end do
    
    end subroutine doJacobian_implFCTconsMass_3D

  end subroutine gfsc_buildJacobianFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianTVDBlock(RcoeffMatrices, rx,&
      fcb_calcMatrixSc_sim, tstep, hstep, bclear, rafcstab,&
      rjacobian, bextendedSparsity, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  The velocity is assumed to be
    ! nonlinear/arbitrary.  Note that this routine serves as a wrapper
    ! for block vectors. If there is only one block, then the
    ! corresponding scalar routine is called.  Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
    
    ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1) then
      
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianTVDBlock')
      call sys_halt()

    else

      call gfsc_buildJacobianTVDScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrixSc_sim,&
          tstep, hstep, bclear, rafcstab, rjacobian,&
          bextendedSparsity, rcollection)

    end if
  end subroutine gfsc_buildJacobianTVDBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianTVDScalar(RcoeffMatrices, rx,&
      fcb_calcMatrixSc_sim, tstep, hstep, bclear, rafcstab,&
      rjacobian, bextendedSparsity, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary. 
    ! This routine will also work for linear velocities but then it is inefficient
    ! since the solution perturbation does not affect the velocity.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear
    
    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
    
    ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Dflux
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_Jac,p_Dx
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    logical :: bisExtended
    

    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianTVDScalar')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)

    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_DcoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianTVDScalar')
      call sys_halt()
    end select

    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_genOffdiagEdges(rafcstab)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if

    
    ! What kind of matrix format are we?
    select case(rjacobian%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_TVD_1D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dpp, p_Dpm,&
            p_Dqp, p_Dqm, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      end select
      
      ! Free storage
      call storage_free(h_Ksep)
      
        
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------

      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian,   p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_TVD_1D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_TVD_2D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_DcoeffY, p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_TVD_3D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_Dx, p_Dflux, p_Dpp, p_Dpm,&
            p_Dqp, p_Dqm, tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      end select
        
      ! Free storage
      call storage_free(h_Ksep)

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianTVDScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,l
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k)+1, Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat7

    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,l
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k), Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat9


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD in 1D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD_1D(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal,&
        DcoeffX, Dx, Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Drploc,Drmloc,Dfluxloc
      real(DP), dimension(NDIM1D) :: c_ij, c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc

      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = DcoeffX(ij)
          c_ji = DcoeffX(ji)
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = DcoeffX(ij)
          c_ji = DcoeffX(ji)

          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux,&
                Kloc, Drploc, Drmloc, Dfluxloc,&
                hstep, iedge, iloc, k, l,&
                bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_TVD_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD in 2D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD_2D(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, DcoeffX, DcoeffY, Dx, Dflux,&
        Dpp, Dpm, Dqp, Dqm, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Drploc,Drmloc,Dfluxloc
      real(DP), dimension(NDIM2D) :: c_ij, c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc      
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji)/)

          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm,&
              c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_TVD_2D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-TVD in 3D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_TVD_3D(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, DcoeffX, DcoeffY, DcoeffZ, Dx,&
        Dflux, Dpp, Dpm, Dqp, Dqm, tstep, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,Dx,Dflux,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: tstep,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Drploc,Drmloc,Dfluxloc
      real(DP), dimension(NDIM3D) :: c_ij, c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc


      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij),DcoeffZ(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji),DcoeffZ(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              tstep, hstep, iedge, i, j, ij, ji, iloc, k, Dpploc,&
              Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                    
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij),DcoeffZ(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji),DcoeffZ(ji)/)

          ! Update local coefficients
          call updateJacobianMat79_TVD(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              tstep, hstep, iedge, i, j, ij, ji, iloc, k, Dpploc,&
              Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc
        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_TVD(&
                IedgeList, Kdiagonal, Dflux, Kloc, Drploc, Drmloc,&
                Dfluxloc, hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_TVD_3D


    !**************************************************************
    ! Update the local coefficients for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.    
    subroutine updateJacobianMat79_TVD(DcoefficientsAtEdge, Dx,&
        Dpp, Dpm, Dqp, Dqm, c_ij, c_ji, tstep, hstep, iedge, i, j, ij, ji,&
        iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dpp,Dpm,Dqp,Dqm,C_ij,C_ji
      real(DP), intent(in) :: tstep,hstep
      integer, intent(in) :: iedge,i,j,k,ij,ji,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(inout) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Dfluxloc
      integer, dimension(:,:), intent(inout) :: Kloc

      ! local variables
      real(DP) :: d_ij,f_ij,l_ij,l_ji,diff,hstep_ik,hstep_jk,dsign
      integer  :: iperturb
      

      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     eliminate the contribution of the edge IJ for the
      !     unperturbed solution values Dx_i and Dx_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Determine prelimited antidiffusive flux
      diff = tstep*(Dx(i)-Dx(j))
      f_ij = min(d_ij,l_ji)*diff
      
      if (i .eq. k) then
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0.0_DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        Dpploc(:,iloc) = Dpp(j)
        Dpmloc(:,iloc) = Dpm(j)
        Dqploc(:,iloc) = Dqp(j)-max(0.0_DP, f_ij)
        Dqmloc(:,iloc) = Dqm(j)-min(0.0_DP, f_ij)

      else

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0.0_DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        Dpploc(:,iloc) = Dpp(i)-max(0.0_DP, f_ij)
        Dpmloc(:,iloc) = Dpm(i)-min(0.0_DP, f_ij)
        Dqploc(:,iloc) = Dqp(i)-max(0.0_DP,-f_ij)
        Dqmloc(:,iloc) = Dqm(i)-min(0.0_DP,-f_ij)
      end if

      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------

      do iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = 3-2*iperturb
        
!!$        ! Compute perturbed coefficients k_ij and k_ji
!!$        call fcb_calcMatrix(Dx(i)+dsign*hstep_ik, Dx(j)+dsign*hstep_jk,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)
        
        ! Apply discrete upwinding
        l_ij = l_ij+d_ij
        l_ji = l_ji+d_ij
        
        ! Due to the (possible) nonlinearity of the velocity vector
        ! the orientation convention for the edge ij may be violated,
        ! that is, the condition 0=l_ij < l_ji may not be valid. In this
        ! case the node number i and j must be swapped logically
        if (l_ij .le. l_ji) then
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
          
          ! In this case the orientation of edge ij remains unchanged
          f_ij = min(d_ij,l_ji)*(diff+tstep*dsign*(hstep_ik-hstep_jk))
          Dfluxloc(iperturb,iloc) = f_ij
        
          if (i .eq. k) then

            ! For node k which is the upwind node
            Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP, f_ij)
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP,-f_ij)
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP,-f_ij)
            
            ! For node l opposite to k which is the downwind node
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, f_ij)
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, f_ij)

          else

            ! For node k which is the downwind node
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP, f_ij)
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP, f_ij)
            
            ! For node l opposite to k
            Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP,-f_ij)
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP,-f_ij)

          end if
          
        else
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/j,i/)
          
          ! In this case the orientation of edge ij needs to be
          ! reverted so as to let i denote the 'upwind' node
          f_ij = -min(d_ij,l_ij)*(diff+tstep*dsign*(hstep_ik-hstep_jk))
          Dfluxloc(iperturb,iloc) = f_ij
          
          if (j .eq. k) then

            ! For node k which is the upwind node
            Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP, f_ij)
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP,-f_ij)
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP,-f_ij)
            
            ! For node l opposite to k which is the downwind node
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, f_ij)
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, f_ij)

          else

            ! For node k which is the downwind node
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP, f_ij)
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP, f_ij)
            
            ! For node l opposite to k which is the upwind node
            Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP,-f_ij)
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP,-f_ij)
          end if
        end if
      end do
    end subroutine updateJacobianMat79_TVD


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-TVD,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine assembleJacobianMat79_TVD(IedgeList, Kdiagonal,&
        Dflux, Kloc, Drploc, Drmloc, Dfluxloc, hstep, iedge, iloc, k, l,&
        bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), dimension(:), intent(in) :: Dflux
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IedgeList,Kloc
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: iedge,iloc,k,l
      logical, intent(in) :: bisExtended

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP) :: f_ij
      integer :: ik,jk,i,j,m,iperturb
      
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      if (m .eq. k) then

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        do iperturb = 1, 2
          
          ! Retrieve precomputed flux
          f_ij = Dfluxloc(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit flux 
            if (f_ij > 0.0_DP) then
              f_ij = Drploc(iperturb,0)*f_ij
            else
              f_ij = Drmloc(iperturb,0)*f_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit flux
            if (f_ij > 0.0_DP) then
              f_ij = Drploc(iperturb,iloc)*f_ij
            else
              f_ij = Drmloc(iperturb,iloc)*f_ij
            end if
            
          end if
          
          ! Adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*f_ij/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end do
        
      elseif (bisExtended) then
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! Interestingly enough, some edge LM which connects two
        ! direct neighbors of the perturbed vertex k is only
        ! processed once due to the fact that either l or (!) m
        ! corresponds to the upwind node.

        if (i .eq. l) then

          if (Dflux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*(Drploc(1,iloc)-Drploc(2,iloc))*Dflux(iedge)/hstep
          else
            f_ij = 0.5_DP*(Drmloc(1,iloc)-Drmloc(2,iloc))*Dflux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end if
      end if
    end subroutine assembleJacobianMat79_TVD
    
  end subroutine gfsc_buildJacobianTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianGPBlock(RcoeffMatrices, rmatrix,&
      rx, rx0, fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, bextendedSparsity, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete transport operator for a
    ! scalar convection equation.  The velocity is assumed to be
    ! nonlinear/arbitrary.  Note that this routine serves as a wrapper
    ! for block vectors. If there is only one block, then the
    ! corresponding scalar routine is called.  Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity

     ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    if (rx%nblocks  .ne. 1 .or.&
        rx0%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianGPBlock')
      call sys_halt()
      
    else
      
      call gfsc_buildJacobianGPScalar(&
          RcoeffMatrices, rmatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), fcb_calcMatrixSc_sim, theta, tstep, hstep,&
          bclear, rafcstab, rjacobian,bextendedSparsity, rcollection)

    end if
  end subroutine gfsc_buildJacobianGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianGPScalar(RcoeffMatrices, rmatrix,&
      rx, rx0, fcb_calcMatrixSc_sim, theta, tstep, hstep, bclear, rafcstab,&
      rjacobian, bextendedSparsity, rcollection)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete transport operator for a scalar convection equation.
    ! The velocity is assumed to be nonlinear/arbitrary. 
    ! This routine will also work for linear velocities but then it is inefficient
    ! since the solution perturbation does not affect the velocity.
!</description>

!<input>
    ! array of coefficient matrices C = (phi_i,D phi_j)
    type(t_matrixScalar), dimension(:), intent(in) :: RcoeffMatrices

    ! consistent mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! initial solution vector
    type(t_vectorScalar), intent(in) :: rx0

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity

    ! callback functions to compute velocity
    include 'intf_calcMatrixSc_sim.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm,p_Dflux,p_Dflux0
    real(DP), dimension(:), pointer :: p_DcoeffX,p_DcoeffY,p_DcoeffZ,p_MC,p_Jac,p_Dx,p_Dx0
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep,ndim
    logical :: bisExtended
    
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianGPScalar')
      call sys_halt()
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    call lsyssc_getbase_double(rmatrix, p_MC)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(rx0, p_Dx0)
    
    ! How many dimensions do we have?
    ndim = size(RcoeffMatrices,1)
    select case(ndim)
    case (NDIM1D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)

    case (NDIM2D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)

    case (NDIM3D)
      call lsyssc_getbase_double(RcoeffMatrices(1), p_DcoeffX)
      call lsyssc_getbase_double(RcoeffMatrices(2), p_DcoeffY)
      call lsyssc_getbase_double(RcoeffMatrices(3), p_DcoeffZ)

    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianGPScalar')
      call sys_halt()
    end select

    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_genOffdiagEdges(rafcstab)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if

    
    ! What kind of matrix format are we?
    select case(rjacobian%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_GP_1D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_GP_2D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_DcoeffY, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_GP_3D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)
      end select
      
      ! Free storage
      call storage_free(h_Ksep)
      
    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
            
      ! How many dimensions do we have?
      select case(ndim)
      case (NDIM1D)
        call doJacobianMat79_GP_1D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM2D)
        call doJacobianMat79_GP_2D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_DcoeffY, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      case (NDIM3D)
        call doJacobianMat79_GP_3D(&
            p_IsuperdiagEdgesIdx, p_IedgeList,&
            p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
            p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
            p_DcoeffX, p_DcoeffY, p_DcoeffZ, p_MC, p_Dx, p_Dx0, p_Dflux, p_Dflux0,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, theta,&
            tstep, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
            rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)
      end select

      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianGPScalar')
      call sys_halt()
    end select
    
  contains

    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,l
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k)+1, Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat7

    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,l
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k), Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat9

    
    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP in 1D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_GP_1D(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, DcoeffX, MC, Dx, Dx0,&
        Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp, Drm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,MC,Dx,Dx0,Dflux,Dflux0
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc
      real(DP), dimension(2,0:NNVEDGE) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      real(DP), dimension(NDIM1D) :: c_ij,c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc
      
      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = DcoeffX(ij)
          c_ji = DcoeffX(ji)
          
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              DcoefficientsAtEdge, MC, Dx, Dx0, Dflux,&
              Dflux0, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              theta, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k,  Dpploc, Dpmloc, Dqploc, Dqmloc,&
              Dfluxloc, Dfluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = DcoeffX(ij)
          c_ji = DcoeffX(ji)
   
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              DcoefficientsAtEdge, MC, Dx, Dx0, Dflux,&
              Dflux0, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              theta, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc,&
              Dfluxloc, Dfluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_GP_1D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP in 2D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_GP_2D(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, DcoeffX, DcoeffY, MC, Dx, Dx0,&
        Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp, Drm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,MC,Dx,Dx0,Dflux,Dflux0
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc
      real(DP), dimension(2,0:NNVEDGE) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      real(DP), dimension(NDIM2D) :: c_ij,c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc

      
      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              DcoefficientsAtEdge, MC, Dx, Dx0, Dflux,&
              Dflux0, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              theta, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc,&
              Dfluxloc, Dfluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji)/)
   
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              DcoefficientsAtEdge, MC, Dx, Dx0, Dflux,&
              Dflux0, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              theta, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc,&
              Dfluxloc, Dfluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_GP_2D


    !**************************************************************
    ! Assemble the Jacobian matrix for FEM-GP in 3D,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine doJacobianMat79_GP_3D(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, DcoeffX, DcoeffY, DcoeffZ, MC, Dx,&
        Dx0, Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, Drp, Drm, theta, tstep, hstep,&
        NEQ, NEDGE, NNVEDGE, bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: DcoeffX,DcoeffY,DcoeffZ,MC,Dx,Dx0,Dflux,Dflux0
      real(DP), dimension(:), intent(in) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: theta,tstep,hstep  
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7
      
      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc
      real(DP), dimension(2,0:NNVEDGE) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      real(DP), dimension(NDIM3D) :: c_ij,c_ji
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ij,ji,ild,iedge,i,j,k,l,iloc,nloc


      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1
          
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij),DcoeffZ(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji),DcoeffZ(ji)/)
          
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              DcoefficientsAtEdge, MC, Dx, Dx0, Dflux,&
              Dflux0, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              theta, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc,&
              Dfluxloc, Dfluxloc0, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
                 
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Determine matrix indices
          ij = IedgeList(3,iedge)
          ji = IedgeList(4,iedge)

          ! Determine matrix coefficients
          c_ij = (/DcoeffX(ij),DcoeffY(ij),DcoeffZ(ij)/)
          c_ji = (/DcoeffX(ji),DcoeffY(ji),DcoeffZ(ji)/)
   
          ! Update local coefficients
          call updateJacobianMat79_GP(&
              DcoefficientsAtEdge, MC, Dx, Dx0, Dflux,&
              Dflux0, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji,&
              theta, tstep, hstep, iedge, i, j, ij, ji,&
              iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc,&
              Dfluxloc, Dfluxloc0, Kloc)
        end do
        
        ! Save total number of local neighbors
        nloc = iloc

        
        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1

            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do

          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_GP(&
                IedgeList, Kdiagonal, Dflux,&
                Dflux0, Drp, Drm, Kloc, Drploc, Drmloc,&
                Dfluxloc, Dfluxloc0, hstep, iedge,&
                iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do

        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_GP_3D

    
    !**************************************************************
    ! Update the local coefficients for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine updateJacobianMat79_GP(DcoefficientsAtEdge, MC, Dx, Dx0,&
        Dflux, Dflux0, Dpp, Dpm, Dqp, Dqm, c_ij, c_ji, theta, tstep, hstep,&
        iedge, i, j, ij, ji, iloc, k, Dpploc, Dpmloc, Dqploc, Dqmloc,&
        Dfluxloc, Dfluxloc0, Kloc)
      
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: MC,Dx,Dx0,Dflux,Dflux0,Dpp,Dpm,Dqp,Dqm,C_ij,C_ji
      real(DP), intent(in) :: theta,tstep,hstep
      integer, intent(in) :: iedge,i,j,k,ij,ji,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(inout) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Dfluxloc,Dfluxloc0
      integer, dimension(:,:), intent(inout)  :: Kloc

      ! local variables
      real(DP) :: m_ij,d_ij,df_ij,f_ij,l_ij,l_ji,p_ij,pf_ij,q_ij,q_ji
      real(DP) :: diff,diff1,diff0,hstep_ik,hstep_jk,dsign
      integer :: iperturb

      
      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     eliminate the contribution of the edge IJ for the
      !     unperturbed solution values Dx_i and Dx_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      l_ij = DcoefficientsAtEdge(2,iedge)
      l_ji = DcoefficientsAtEdge(3,iedge)
      
      ! Include consistent mass matrix
      m_ij = MC(ij)
      q_ij = m_ij/tstep+l_ij
      q_ji = m_ij/tstep+l_ji

      ! Determine solution differences
      diff1 = Dx(i)-Dx(j)
      diff0 = Dx0(i)-Dx0(j)

      ! Determine total solution difference
      diff = tstep*(theta*diff1+(1.0_DP-theta)*diff0)

      ! Compute antidiffusive flux 
      if (abs(diff) < AFCSTAB_EPSABS) then
        p_ij = 0
        f_ij = 0
      else
        p_ij = max(0.0_DP, m_ij*(diff1-diff0)/diff+d_ij)
        f_ij = p_ij*diff
      end if

      ! Prelimit the antidiffusive flux
      pf_ij = min(p_ij, l_ji)*diff
      
      ! Compute the remaining flux
      df_ij = f_ij-pf_ij

      if (i .eq. k) then
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0.0_DP
        
        ! Update nodal coefficients for vertex j (!) which is the downwind node
        Dpploc(:,iloc) = Dpp(j)-max(0.0_DP,-df_ij)
        Dpmloc(:,iloc) = Dpm(j)-min(0.0_DP,-df_ij)
        Dqploc(:,iloc) = Dqp(j)-max(0.0_DP, diff)*q_ji
        Dqmloc(:,iloc) = Dqm(j)-min(0.0_DP, diff)*q_ji

      else

        ! Store global node number of the opposite node
        Kloc(1,iloc) = i

        ! Compute signed perturbation parameters
        hstep_ik = 0.0_DP; hstep_jk = hstep
        
        ! Update nodal coefficients for vertex i (!) which is the upwind node
        Dpploc(:,iloc) = Dpp(i)-max(0.0_DP, f_ij)
        Dpmloc(:,iloc) = Dpm(i)-min(0.0_DP, f_ij)
        Dqploc(:,iloc) = Dqp(i)-max(0.0_DP,-diff)*q_ij
        Dqmloc(:,iloc) = Dqm(i)-min(0.0_DP,-diff)*q_ij
      end if

      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------
      
      do iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = -2*iperturb+3
        
!!$        ! Compute perturbed velocity
!!$        call fcb_calcMatrix(Dx(i)+dsign*hstep_ik, Dx(j)+dsign*hstep_jk,&
!!$            C_ij, C_ji, i, j, l_ij, l_ji, d_ij)

        ! Perform discrete upwinding
        l_ij = l_ij+d_ij
        l_ji = l_ji+d_ij

        q_ij = m_ij/tstep+l_ij
        q_ji = m_ij/tstep+l_ji

        ! Due to the (possible) nonlinearity of the velocity vector
        ! the orientation convention for the edge ij may be violated,
        ! that is, the condition 0=l_ij < l_ji may not be valid. In this
        ! case the node number i and j must be swapped logically
        if (l_ij .le. l_ji) then
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
          
          ! Update solution difference
          diff1 = Dx(i)-Dx(j)+dsign*(hstep_ik-hstep_jk)

          ! Update total solution difference
          diff = tstep*(theta*diff1+(1.0_DP-theta)*diff0)

          ! Compute antidiffusive flux
          if (abs(diff) < AFCSTAB_EPSABS) then
            p_ij = 0
            f_ij = 0
          else
            p_ij = max(0.0_DP,m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = p_ij*diff
          end if
          
          ! Prelimit the antidiffusive flux
          pf_ij = min(p_ij,l_ji)*diff
          Dfluxloc0(iperturb,iloc) = pf_ij
          
          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          Dfluxloc(iperturb,iloc) = df_ij
        
          if (i .eq. k) then

            ! For node k which is the upwind node
            Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP, f_ij)
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP,-diff)*q_ij
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP,-diff)*q_ij
            
            ! For node l opposite to k which is the downwind node
            Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP,-df_ij)
            Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP,-df_ij)
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, diff)*q_ji
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, diff)*q_ji

          else

            ! For node k which is the downwind node
            Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP,-df_ij)
            Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP,-df_ij)
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP, diff)*q_ji
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP, diff)*q_ji
            
            ! For node l opposite to k
            Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP,-diff)*q_ij
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP,-diff)*q_ij

          end if
          
        else
          
          ! Save oriented node numbers
          Kloc(2*iperturb:2*iperturb+1,iloc) = (/j,i/)
          
          ! Update solution difference
          diff1 = Dx(i)-Dx(j)+dsign*(hstep_ik-hstep_jk)
          
          ! Update total solution difference
          diff = tstep*(theta*diff1+(1.0_DP-theta)*diff0)

          ! Compute antidiffusive flux
          if (abs(diff) < AFCSTAB_EPSABS) then
            p_ij = 0
            f_ij = 0
          else
            p_ij = max(0.0_DP, m_ij*(diff1-diff0)/diff+d_ij)
            f_ij = -p_ij*diff
          end if

          ! Prelimit the antidiffusive flux
          pf_ij = -min(p_ij,l_ij)*diff
          Dfluxloc0(iperturb,iloc) = pf_ij

          ! Compute the remaining flux
          df_ij = f_ij-pf_ij
          Dfluxloc(iperturb,iloc) = df_ij
          
          if (j .eq. k) then
            
            ! For node k which is the upwind node
            Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP, f_ij)
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP, diff)*q_ij
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP, diff)*q_ij
                       
            ! For node l opposite to k which is the downwind node
            Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP,-df_ij)
            Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP,-df_ij)
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP,-diff)*q_ji
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP,-diff)*q_ji

          else

            ! For node k which is the downwind node
            Dpploc(iperturb,0) = Dpploc(iperturb,0)+max(0.0_DP,-df_ij)
            Dpmloc(iperturb,0) = Dpmloc(iperturb,0)+min(0.0_DP,-df_ij)
            Dqploc(iperturb,0) = Dqploc(iperturb,0)+max(0.0_DP,-diff)*q_ji
            Dqmloc(iperturb,0) = Dqmloc(iperturb,0)+min(0.0_DP,-diff)*q_ji
            
            ! For node l opposite to k which is the upwind node
            Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, f_ij)
            Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, f_ij)
            Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, diff)*q_ij
            Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, diff)*q_ij

          end if
        end if
      end do
    end subroutine updateJacobianMat79_GP


    !**************************************************************
    ! Assemble the given column of the Jacobian for FEM-GP,
    ! whereby the matrix can be stored in format 7 or 9.
    subroutine assembleJacobianMat79_GP(IedgeList, Kdiagonal,&
        Dflux, Dflux0, Drp, Drm, Kloc, Drploc, Drmloc, Dfluxloc, Dfluxloc0,&
        hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc,Dfluxloc0
      real(DP), dimension(:), intent(in) :: Dflux,Dflux0,Drp,Drm
      real(DP), intent(in) :: hstep
      integer, dimension(:,:), intent(in) :: IedgeList,Kloc
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: iedge,iloc,k,l
      logical, intent(in) :: bisExtended

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP) :: f_ij,pf_ij,df_ij
      integer :: ik,jk,i,j,m,iperturb
      
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      if (m .eq. k) then

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        do iperturb = 1, 2
          
          ! Retrieve precomputed fluxes
          df_ij = Dfluxloc(iperturb,iloc)
          pf_ij = Dfluxloc0(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = Drploc(iperturb,0)*pf_ij
            else
              pf_ij = Drmloc(iperturb,0)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(Drploc(iperturb,0), Drmloc(iperturb,iloc))*df_ij
            else
              df_ij = min(Drmloc(iperturb,0), Drploc(iperturb,iloc))*df_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit upwind contribution
            if (pf_ij > 0.0_DP) then
              pf_ij = Drploc(iperturb,iloc)*pf_ij
            else
              pf_ij = Drmloc(iperturb,iloc)*pf_ij
            end if

            ! Limit symmetric contribution
            if (df_ij > 0.0_DP) then
              df_ij = min(Drmloc(iperturb,0), Drploc(iperturb,iloc))*df_ij
            else
              df_ij = min(Drploc(iperturb,0), Drmloc(iperturb,iloc))*df_ij
            end if
            
          end if
          
          ! Combine both contributions and 
          ! adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*(pf_ij+df_ij)/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end do
        
      elseif (bisExtended) then
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! For the symmetric part, we must also check if node j
        ! correspond to the modified vertex l.
        
        if (i .eq. l) then

          ! Get precomputed fluxes
          pf_ij = Dflux0(iedge)
          df_ij = Dflux(iedge)

          ! Limit upwind contribution
          if (pf_ij > 0.0_DP) then
            pf_ij = (Drploc(1,iloc)-Drploc(2,iloc))*pf_ij
          else
            pf_ij = (Drmloc(1,iloc)-Drmloc(2,iloc))*pf_ij
          end if

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(Drploc(1,iloc), Drm(j))-&
                     min(Drploc(2,iloc), Drm(j)))*df_ij
          else
            df_ij = (min(Drmloc(1,iloc), Drp(j))-&
                     min(Drmloc(2,iloc), Drp(j)))*df_ij
          end if

          ! Combine both contributions
          f_ij = 0.5_DP*(pf_ij+df_ij)/hstep

          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        else

          ! Get precomputed flux (only symmetric part)
          df_ij = Dflux(iedge)

          ! Limit symmetric contribution
          if (df_ij > 0.0_DP) then
            df_ij = (min(Drp(i), Drmloc(1,iloc))-&
                     min(Drp(i), Drmloc(2,iloc)))*df_ij
          else
            df_ij = (min(Drm(i), Drploc(1,iloc))-&
                     min(Drm(i), Drploc(2,iloc)))*df_ij
          end if

          ! Compute divided difference
          f_ij = 0.5_DP*df_ij/hstep

          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)
          
          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
          
        end if
      end if
    end subroutine assembleJacobianMat79_GP
    
  end subroutine gfsc_buildJacobianGPScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianSymmBlock(rx, dscale, hstep, bclear,&
      rafcstab, rjacobian, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the
    ! stabilisation part of the discrete diffusion operator for a
    ! scalar convection equation.  Note that this routine serves as a
    ! wrapper for block vectors. If there is only one block, then the
    ! corresponding scalar routine is called.  Otherwise, an error is
    ! thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   
!</inputoutput>
!</subroutine>

    if (rx%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianSymmBlock')
      call sys_halt()

    else

      call gfsc_buildJacobianSymmScalar(rx%RvectorBlock(1), dscale,&
          hstep, bclear, rafcstab, rjacobian, bextendedSparsity)

    end if
  end subroutine gfsc_buildJacobianSymmBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_buildJacobianSymmScalar(rx, dscale, hstep, bclear,&
      rafcstab, rjacobian, bextendedSparsity)

!<description>
    ! This subroutine assembles the Jacobian matrix for the stabilisation
    ! part of the discrete diffusion operator for a scalar convection equation.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! perturbation parameter
    real(DP), intent(in) :: hstep
    
    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FALSE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! OPTIONAL: Switch for matrix assembly
    ! TRUE  : assemble the Jacobian matrix with extended sparsity pattern (default)
    ! FALSE : assemble the Jacobian matrix with standard sparsity pattern
    logical, intent(in), optional :: bextendedSparsity
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! Jacobian matrix
    type(t_matrixScalar), intent(inout) :: rjacobian   
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dflux,p_Dx,p_Jac
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IsuperdiagEdgesIdx
    integer, dimension(:), pointer :: p_IsubdiagEdges
    integer, dimension(:), pointer :: p_IsubdiagEdgesIdx
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer :: h_Ksep
    logical :: bisExtended

    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)  .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianSymmScalar')
      call sys_halt()
    end if
    
    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rjacobian)
    
    ! Check if off-diagonal edges need to be generated
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_OFFDIAGONALEDGES) .eq. 0)&
        call afcstab_genOffdiagEdges(rafcstab)
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IsupdiagEdgeIdx(rafcstab, p_IsuperdiagEdgesIdx)
    call afcstab_getbase_IsubdiagEdge(rafcstab, p_IsubdiagEdges)
    call afcstab_getbase_IsubdiagEdgeIdx(rafcstab, p_IsubdiagEdgesIdx)
    call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
    call lsyssc_getbase_double(rjacobian, p_Jac)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! Assembled extended Jacobian matrix?
    if (present(bextendedSparsity)) then
      bisExtended = bextendedSparsity
    else
      bisExtended = .true.
    end if


    ! What kind of matrix format are we?
    select case(rjacobian%cmatrixFormat)
    case(LSYSSC_MATRIX7)
      !-------------------------------------------------------------------------
      ! Matrix format 7
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_Symm(&
          p_IsuperdiagEdgesIdx, p_IedgeList,&
          p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
          p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kld,&
          p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm,&
          dscale, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
          rafcstab%NNVEDGE, bisExtended, .true., p_Ksep, p_Jac)

      ! Free storage
      call storage_free(h_Ksep)

    case(LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Matrix format 9
      !-------------------------------------------------------------------------
      
      ! Set pointers
      call lsyssc_getbase_Kld(rjacobian, p_Kld)
      call lsyssc_getbase_Kcol(rjacobian, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rjacobian, p_Kdiagonal)
      
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rjacobian%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rjacobian%NEQ+1)
      call lalg_vectorAddScalarInt(p_Ksep, 1)
      
      call doJacobianMat79_Symm(&
          p_IsuperdiagEdgesIdx, p_IedgeList,&
          p_IsubdiagEdgesIdx, p_IsubdiagEdges,&
          p_DcoefficientsAtEdge, p_Kld, p_Kcol, p_Kdiagonal,&
          p_Dx, p_Dflux, p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm,&
          dscale, hstep, rafcstab%NEQ, rafcstab%NEDGE,&
          rafcstab%NNVEDGE, bisExtended, .false., p_Ksep, p_Jac)

      ! Free storage
      call storage_free(h_Ksep)
      
    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsc_buildJacobianSymmScalar')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routine follow
    
    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat7(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,l
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k)+1, Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat7


    !**************************************************************    
    ! Adjust the diagonal separator.
    ! The separator is initialied by the column separator (increased
    ! by one if this is necessary for matrix format 7).
    ! Based on the matrix structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: k

      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      integer :: ild,l
      
      
      ! Loop over all entries of the k-th row
      do ild = Kld(k), Kld(k+1)-1
        
        ! Get the column number
        l = Kcol(ild)
        
        ! Move separator to next position
        Ksep(l) = Ksep(l)+1
      end do
    end subroutine adjustKsepMat9


    !**************************************************************
    ! Assemble the Jacobian matrix for symmetric flux limiting
    subroutine doJacobianMat79_Symm(IsuperdiagEdgesIdx,&
        IedgeList, IsubdiagEdgesIdx, IsubdiagEdges,&
        DcoefficientsAtEdge, Kld, Kcol, Kdiagonal, Dx, Dflux, Dpp, Dpm,&
        Dqp, Dqm, Drp, Drm, dscale, hstep, NEQ, NEDGE, NNVEDGE,&
        bisExtended, bisMat7, Ksep, Jac)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dflux,Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), intent(in) :: dscale,hstep
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IsuperdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdgesIdx
      integer, dimension(:), intent(in) :: IsubdiagEdges
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE,NNVEDGE
      logical, intent(in) :: bisExtended,bisMat7

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(2,0:NNVEDGE) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Drploc,Drmloc,Dfluxloc
      integer, dimension(5,NNVEDGE) :: Kloc
      integer :: ild,iedge,i,j,k,l,iloc,nloc
      

      ! Loop over all columns of the Jacobian matrix
      do k = 1, NEQ
        
        ! Assemble nodal coefficients P and Q for node k and all vertices 
        ! surrounding node k. Note that it suffices to initialize only
        ! those quantities which belong to node k. All other quantities
        ! will be overwritten in the update procedure below
        Dpploc(:,0) = 0; Dpmloc(:,0) = 0
        Dqploc(:,0) = 0; Dqmloc(:,0) = 0
        
        ! Initialize local counter
        iloc = 0

        ! Loop over all subdiagonal edges
        do ild = IsubdiagEdgesIdx(k), IsubdiagEdgesIdx(k+1)-1
          
          ! Get edge number
          iedge = IsubdiagEdges(ild)
          
          ! Increase local counter
          iloc = iloc+1

          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
          ! Update local coefficients
          call updateJacobianMat79_Symm(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm,&
              hstep, iedge, i, j, iloc, k,&
              Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Loop over all superdiagonal edges
        do iedge = IsuperdiagEdgesIdx(k), IsuperdiagEdgesIdx(k+1)-1

          ! Increase local counter
          iloc = iloc+1
             
          ! Determine indices
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Update local coefficients
          call updateJacobianMat79_Symm(&
              DcoefficientsAtEdge, Dx, Dpp, Dpm, Dqp, Dqm,&
              hstep, iedge, i, j, iloc, k,&
              Dpploc, Dpmloc, Dqploc, Dqmloc, Dfluxloc, Kloc)
        end do

        ! Save total number of local neighbors
        nloc = iloc
        

        ! Compute nodal correction factors for node k and all other
        ! nodes l_1,l_2,...,l_|k| which are direct neighbors to k
        Drploc(:,0:nloc) = gfsc_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = gfsc_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


        ! Now we have all required information, the local fluxes, the
        ! nodal correction factors, etc. for assembling the k-th
        ! column of the Jacobian matrix. Hence, loop over all direct
        ! neighbors of node k (stored during coefficient assembly)
        do iloc = 1, nloc
          
          ! Get the global node number of the node l opposite to k
          l = Kloc(1,iloc)
          
          ! Loop over all subdiagonal edges
          do ild = IsubdiagEdgesIdx(l), IsubdiagEdgesIdx(l+1)-1
            
            ! Get edge number
            iedge = IsubdiagEdges(ild)
            
            call assembleJacobianMat79_Symm(&
                IedgeList, Kld, Kcol, Dflux, Drp, Drm,&
                Kloc, Drploc, Drmloc, Dfluxloc, dscale, hstep,&
                iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
          
          ! Loop over all superdiagonal edges
          do iedge = IsuperdiagEdgesIdx(l), IsuperdiagEdgesIdx(l+1)-1
            
            call assembleJacobianMat79_Symm(&
                IedgeList, Kld, Kcol, Dflux, Drp, Drm,&
                Kloc, Drploc, Drmloc, Dfluxloc, dscale, hstep,&
                iedge, iloc, k, l, bisExtended, Ksep, Jac)
          end do
        end do
        
        ! Adjust the diagonal separator
        if (bisMat7) then
          call adjustKsepMat7(Kld, Kcol, k, Ksep)
        else
          call adjustKsepMat9(Kld, Kcol, k, Ksep)
        end if
      end do   ! end-of k-loop
    end subroutine doJacobianMat79_Symm

    
    !**************************************************************
    ! Update the local coefficients for symmetric flux limiting
    subroutine updateJacobianMat79_Symm(DcoefficientsAtEdge, Dx, Dpp,&
        Dpm, Dqp, Dqm, hstep, iedge, i, j, iloc, k, Dpploc, Dpmloc, Dqploc,&
        Dqmloc, Dfluxloc, Kloc)

      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx,Dpp,Dpm,Dqp,Dqm
      real(DP), intent(in) :: hstep
      integer, intent(in) :: iedge,i,j,k,iloc
      
      ! We actually know, that all local quantities start at index zero
      real(DP), dimension(:,0:), intent(inout) :: Dpploc,Dpmloc,Dqploc,Dqmloc,Dfluxloc
      integer, dimension(:,:), intent(inout) :: Kloc

      ! local variables
      real(DP) :: d_ij,f_ij,s_ij,diff,hstep_ik,hstep_jk,dsign
      integer :: iperturb

      
      !------------------------------------------------------------
      ! (1) unperturbed values: Retrieve the global Ps and Qs and
      !     copy their content to the local ones. Moreover,
      !     eliminate the contribution of the edge IJ for the
      !     unperturbed solution values Dx_i and Dx_j.
      !------------------------------------------------------------
      ! Determine coefficients
      d_ij = DcoefficientsAtEdge(1,iedge)
      s_ij = DcoefficientsAtEdge(2,iedge)

      ! Determine solution difference
      diff = Dx(i)-Dx(j)

      if (i .eq. k) then
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = j

        ! Compute signed perturbation parameters
        hstep_ik = hstep; hstep_jk = 0.0_DP
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        Dpploc(:,iloc) = Dpp(j)-max(0.0_DP, -f_ij)
        Dpmloc(:,iloc) = Dpm(j)-max(0.0_DP, -f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        Dqploc(:,iloc) = Dqp(j)-max(0.0_DP, -f_ij)
        Dqmloc(:,iloc) = Dqm(j)-min(0.0_DP, -f_ij)
        
      else
        
        ! Store global node number of the opposite node
        Kloc(1,iloc) = i
        
        ! Compute signed perturbation parameters
        hstep_ik = 0.0_DP; hstep_jk = hstep
        
        ! Compute raw antidiffusve flux
        f_ij = d_ij*diff
        
        ! Update sums of raw antidiffusive fluxes
        Dpploc(:,iloc) = Dpp(i)-max(0.0_DP, f_ij)
        Dpmloc(:,iloc) = Dpm(i)-min(0.0_DP, f_ij)
        
        ! Compute admissible edge contribution
        f_ij = -s_ij*diff
        
        ! Update upper/lower bounds
        Dqploc(:,iloc) = Dqp(i)-max(0.0_DP, f_ij)
        Dqmloc(:,iloc) = Dqm(i)-min(0.0_DP, f_ij)
      end if
      
      !------------------------------------------------------------
      ! (2) perturbed values: Now, the local Ps and Qs still
      !     require the contribution of the perturbed solution
      !     values u +/- h*e_k, whereby e_k denotes the k-th unit
      !     vector and h stands for the perturbation step length
      !------------------------------------------------------------
      
      !------------------------------------------------------------
      ! (3) perform the perturbation for "+/-h*e_k"
      !------------------------------------------------------------
        
      do iperturb = 1, 2
        
        ! Compute correct sign of perturbation
        dsign = -2*iperturb+3  
        
        ! Save local node numbers
        Kloc(2*iperturb:2*iperturb+1,iloc) = (/i,j/)
        
        if (i .eq. k) then
          
          ! Compute raw antidiffusve flux
          f_ij = d_ij*(diff+dsign*(hstep_ik-hstep_jk))
          Dfluxloc(iperturb,iloc) = f_ij

          ! Update sums of raw antidiffusive fluxes
          Dpploc(iperturb,0)    = Dpploc(iperturb,0)+max(0.0_DP, f_ij)
          Dpmloc(iperturb,0)    = Dpmloc(iperturb,0)+min(0.0_DP, f_ij)
          Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, -f_ij)
          Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))

          ! Update upper/lower bounds
          Dqploc(iperturb,0)    = Dqploc(iperturb,0)+max(0.0_DP, f_ij)
          Dqmloc(iperturb,0)    = Dqmloc(iperturb,0)+min(0.0_DP, f_ij)
          Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, -f_ij)
          Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, -f_ij)
          
        else
          
          ! Compute raw antidiffusve flux
          f_ij = d_ij*(diff+dsign*(hstep_ik-hstep_jk))
          Dfluxloc(iperturb,iloc) = f_ij

          ! Update sums of raw antidiffusive fluxes
          Dpploc(iperturb,iloc) = Dpploc(iperturb,iloc)+max(0.0_DP, f_ij)
          Dpmloc(iperturb,iloc) = Dpmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          Dpploc(iperturb,0)    = Dpploc(iperturb,0)+max(0.0_DP, -f_ij)
          Dpmloc(iperturb,0)    = Dpmloc(iperturb,0)+min(0.0_DP, -f_ij)

          ! Compute admissible edge contribution
          f_ij = -s_ij*(diff+dsign*(hstep_ik-hstep_jk))
          
          ! Update upper/lower bounds
          Dqploc(iperturb,iloc) = Dqploc(iperturb,iloc)+max(0.0_DP, f_ij)
          Dqmloc(iperturb,iloc) = Dqmloc(iperturb,iloc)+min(0.0_DP, f_ij)
          Dqploc(iperturb,0)    = Dqploc(iperturb,0)+max(0.0_DP, -f_ij)
          Dqmloc(iperturb,0)    = Dqmloc(iperturb,0)+min(0.0_DP, -f_ij)
        end if
      end do
    end subroutine updateJacobianMat79_Symm

    
    !**************************************************************
    ! Assemble the given column of the Jacobian for symmetric flux limiting
    subroutine assembleJacobianMat79_Symm(IedgeList, Kdiagonal,&
        Kcol, Dflux, Drp, Drm, Kloc, Drploc, Drmloc, Dfluxloc, dscale,&
        hstep, iedge, iloc, k, l, bisExtended, Ksep, Jac)

      real(DP), dimension(:,0:), intent(in) :: Drploc,Drmloc,Dfluxloc
      real(DP), dimension(:), intent(in) :: Drp,Drm,Dflux
      real(DP), intent(in) :: dscale,hstep
      integer, dimension(:,:), intent(in)  :: IedgeList,Kloc
      integer, dimension(:), intent(in) :: Kdiagonal,Kcol
      integer, intent(in) :: iedge,iloc,k,l
      logical, intent(in) :: bisExtended

      real(DP), dimension(:), intent(inout) :: Jac
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP) :: f_ij
      integer :: ik,jk,i,j,m,iperturb
      
      ! Get global node number for edge IJ and the 
      ! number of the node m which is not l
      i = IedgeList(1,iedge)
      j = IedgeList(2,iedge)
      m = (i+j)-l
      
      ! We need to find out, which kind of edge is processed
      if (m .eq. k) then

        !------------------------------------------------------------
        ! 1. Case: primary edge
        !------------------------------------------------------------
        ! The current edge connects the perturbed node k with its
        ! direct neighbor l. Hence, all required information can be
        ! extracted from the local arrays and no global data
        ! retrieval has to be performed.
        !
        ! (a) The edge orientation needs to be adjusted for each
        !     perturbation direction
        
        do iperturb = 1, 2
          
          ! Retrieve precomputed flux
          f_ij = Dfluxloc(iperturb,iloc)

          ! Adjust edge orientation
          i = Kloc(2*iperturb,iloc)
          j = Kloc(2*iperturb+1,iloc)
          
          ! Which node is located upwind?
          if (i .eq. k) then
            
            ! Get corresponding matrix indices
            ik = Kdiagonal(i); jk = Ksep(j)
            
            ! Limit Dflux 
            if (f_ij > 0.0_DP) then
              f_ij = dscale*min(Drploc(iperturb,0), Drmloc(iperturb,iloc))*f_ij
            else
              f_ij = dscale*min(Drmloc(iperturb,0), Drploc(iperturb,iloc))*f_ij
            end if
            
          else
            
            ! Get corresponding matrix indices
            jk = Kdiagonal(j); ik = Ksep(i)
            
            ! Limit flux
            if (f_ij > 0.0_DP) then
              f_ij = dscale*min(Drploc(iperturb,iloc), Drmloc(iperturb,0))*f_ij
            else
              f_ij = dscale*min(Drmloc(iperturb,iloc), Drploc(iperturb,0))*f_ij
            end if
            
          end if
          
          ! Adopt sign for perturbation direction
          f_ij = -(iperturb-1.5_DP)*f_ij/hstep

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij
        end do
        
      elseif (bisExtended) then
        
        !------------------------------------------------------------
        ! 2. Case: secondary edge
        !------------------------------------------------------------
        ! The current edge connects two nodes l and m which both are
        ! not equal to the perturbed vertex k. Thus, the influence of
        ! the solution perturbation can only be due to a change in
        ! the correction factors alpha_ij. Moreover, for upwind
        ! -biased flux limiting techniques only the nodal correction
        ! factors for the upwind node i is used. Hence, it suffices
        ! to check if node i corresponds to the modified vertex l.
        ! Interestingly enough, some edge LM which connects two
        ! direct neighbors of the perturbed vertex k is only
        ! processed once due to the fact that either l or (!) m
        ! corresponds to the upwind node.

        if (i .eq. l) then

          if (Dflux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*dscale*(min(Drploc(1,iloc), Drm(j))-&
                                  min(Drploc(2,iloc), Drm(j)))*Dflux(iedge)/hstep
          else
            f_ij = 0.5_DP*dscale*(min(Drmloc(1,iloc), Drp(j))-&
                                  min(Drmloc(2,iloc), Drp(j)))*Dflux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        else

          if (Dflux(iedge) > 0.0_DP) then
            f_ij = 0.5_DP*dscale*(min(Drp(i), Drmloc(1,iloc))-&
                                  min(Drp(i), Drmloc(2,iloc)))*Dflux(iedge)/hstep
          else
            f_ij = 0.5_DP*dscale*(min(Drm(i), Drploc(1,iloc))-&
                                  min(Drm(i), Drploc(2,iloc)))*Dflux(iedge)/hstep
          end if
          
          ! Get corresponding matrix indices
          ik = Ksep(i); jk = Ksep(j)

          ! Apply perturbed antidiffusive contribution
          Jac(ik) = Jac(ik)-f_ij
          Jac(jk) = Jac(jk)+f_ij

        end if

      end if
    end subroutine assembleJacobianMat79_Symm

  end subroutine gfsc_buildJacobianSymmScalar

  ! ****************************************************************************
  ! Here, some private auxiliary subroutine follow
  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine gfsc_combineFluxesDble(NEDGE, dscale, Dflux1, Dflux2, Dalpha)

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
!</input>

!<inputoutput>
    ! Second flux
    real(DP), dimension(:), intent(inout) :: Dflux2
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: iedge
    
    if (present(Dalpha)) then
      
      if (dscale .eq. 1.0_DP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge)&
                        + Dalpha(iedge) * Dflux1(iedge)
        end do
        !$omp end parallel do
        
      elseif (dscale .eq. -1.0_DP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge)&
                        - Dalpha(iedge) * Dflux1(iedge)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
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
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge) + Dflux1(iedge)
        end do
        !$omp end parallel do
        
      elseif (dscale .eq. -1.0_DP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge) - Dflux1(iedge)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(iedge) = Dflux2(iedge) + dscale * Dflux1(iedge)
        end do
        !$omp end parallel do
        
      end if
      
    end if

  end subroutine gfsc_combineFluxesDble

  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine gfsc_combineFluxesSngl(NEDGE, fscale, Fflux1, Fflux2, Falpha)

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
!</input>

!<inputoutput>
    ! Second flux
    real(SP), dimension(:), intent(inout) :: Fflux2
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: iedge
    
    if (present(Falpha)) then
      
      if (fscale .eq. 1.0_SP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge)&
                        + Falpha(iedge) * Fflux1(iedge)
        end do
        !$omp end parallel do
        
      elseif (fscale .eq. -1.0_SP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge)&
                        - Falpha(iedge) * Fflux1(iedge)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
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
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge) + Fflux1(iedge)
        end do
        !$omp end parallel do
        
      elseif (fscale .eq. -1.0_SP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge) - Fflux1(iedge)
        end do
        !$omp end parallel do
        
      else
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(iedge) = Fflux2(iedge) + fscale * Fflux1(iedge)
        end do
        !$omp end parallel do
        
      end if
      
    end if

  end subroutine gfsc_combineFluxesSngl

  !*****************************************************************************

!<function>
  
  elemental function gfsc_limitUnboundedDble(p, q, dval) result(r)

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
  end function gfsc_limitUnboundedDble
  
  !*****************************************************************************

!<function>
  
  elemental function gfsc_limitUnboundedSngl(p, q, fval) result(r)

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
  end function gfsc_limitUnboundedSngl

  !*****************************************************************************

!<function>
  
  elemental function gfsc_limitBoundedDble(p, q, dval, dbound) result(r)

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
  end function gfsc_limitBoundedDble

  !*****************************************************************************

!<function>
  
  elemental function gfsc_limitBoundedSngl(p, q, fval, fbound) result(r)

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
  end function gfsc_limitBoundedSngl

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_upwindOrientationDble(Dcoefficients, IedgeList,&
      DcoeffsAtEdge, ipos, jpos)

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
!</input>

!<inputoutput>
    ! Coefficients used to determine the upwind direction
    ! DIMENSION(1:IPOS:JPOS:NPOS,1:NEDGE),
    ! where IPOS and JPOS are the positions where the data
    ! for node I and J is stored, respectively
    real(DP), dimension(:,:), intent(inout) :: Dcoefficients
    
    ! Additional coefficients at the edges
    ! DIMENSION(1:N2,1:2,1:NEDGE)
    real(DP), dimension(:,:,:), intent(inout) :: DcoeffsAtEdge

    ! List of edges
    ! DIMENSION(1:N3,1:NEDGE)
    integer, dimension(:,:), intent(inout) :: IedgeList   
!</inputoutput>
!</subroutine

    ! local variables
    real(DP) :: daux
    integer :: iedge,iaux,naux,nedge
    
    nedge = size(IedgeList,2)
    naux  = size(DcoeffsAtEdge,1)
    
    !$omp parallel do default(shared) private(iaux,daux)&
    !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
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
        
        ! Swap edge data data_ij <-> data_ji
        daux = Dcoefficients(ipos,iedge)
        Dcoefficients(ipos,iedge) = Dcoefficients(jpos,iedge)
        Dcoefficients(jpos,iedge) = daux
        
        ! Swap edgewise coefficients Data_ij <-> Data_ji
        do iaux = 1, naux
          daux = DcoeffsAtEdge(iaux,1,iedge)
          DcoeffsAtEdge(iaux,1,iedge) = DcoeffsAtEdge(iaux,2,iedge)
          DcoeffsAtEdge(iaux,2,iedge) = daux
        end do
      end if
    end do
    !$omp end parallel do
    
  end subroutine gfsc_upwindOrientationDble

  !*****************************************************************************

!<subroutine>

  subroutine gfsc_upwindOrientationSngl(Fcoefficients, IedgeList,&
      FcoeffsAtEdge, ipos, jpos)

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
!</input>

!<inputoutput>
    ! Coefficients used to determine the upwind direction
    ! DIMENSION(1:IPOS:JPOS:NPOS,1:NEDGE),
    ! where IPOS and JPOS are the positions where the data
    ! for node I and J is stored, respectively
    real(SP), dimension(:,:), intent(inout) :: Fcoefficients
    
    ! Additional coefficients at the edges
    ! DIMENSION(1:N2,1:2,1:NEDGE)
    real(SP), dimension(:,:,:), intent(inout) :: FcoeffsAtEdge

    ! List of edges
    ! DIMENSION(1:N3,1:NEDGE)
    integer, dimension(:,:), intent(inout) :: IedgeList   
!</inputoutput>
!</subroutine

    ! local variables
    real(SP) :: faux
    integer :: iedge,iaux,naux,nedge
    
    nedge = size(IedgeList,2)
    naux  = size(FcoeffsAtEdge,1)
    
    !$omp parallel do default(shared) private(iaux,faux)&
    !$omp if (NEDGE > GFSC_NEDGEMIN_OMP)
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
        
        ! Swap edge data data_ij <-> data_ji
        faux = Fcoefficients(ipos,iedge)
        Fcoefficients(ipos,iedge) = Fcoefficients(jpos,iedge)
        Fcoefficients(jpos,iedge) = faux
        
        ! Swap edgewise coefficients Data_ij <-> Data_ji
        do iaux = 1, naux
          faux = FcoeffsAtEdge(iaux,1,iedge)
          FcoeffsAtEdge(iaux,1,iedge) = FcoeffsAtEdge(iaux,2,iedge)
          FcoeffsAtEdge(iaux,2,iedge) = faux
        end do
      end if
    end do
    !$omp end parallel do
    
  end subroutine gfsc_upwindOrientationSngl

end module groupfemscalar

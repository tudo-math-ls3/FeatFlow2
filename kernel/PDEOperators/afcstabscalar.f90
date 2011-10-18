!##############################################################################
!# ****************************************************************************
!# <name> afcstabscalar </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the algebraic
!# flux correction methodology proposed by Kuzmin, Moeller and Turek
!# in a series of publications. As a starting point for scalar
!# conservation laws, the reader is referred to the book chapter
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
!# afcstabbase. The initialisation as a scalar stabilisation
!# structure is done by the routine afcsc_initStabilisation.
!#
!# The following routines are available:
!#
!# 1.) afcsc_initStabilisation = afcsc_initStabByMatrix
!#                               afcsc_initStabByGroupFEMSet
!#     -> Initialises the stabilisation structure
!#
!# 2.) afcsc_buildVectorFCT = afcsc_buildVectorFCTScalar /
!#                            afcsc_buildVectorFCTBlock /
!#     -> Assembles the vector for AFC stabilisation of FCT type
!#
!# 3.) afcsc_buildVectorTVD = afcsc_buildVectorTVDScalar /
!#                            afcsc_buildVectorTVDBlock
!#     -> Assembles the vector for AFC stabilisation of TVD type
!#
!# 4.) afcsc_buildVectorGP = afcsc_buildVectorGPScalar /
!#                           afcsc_buildVectorGPBlock
!#     -> Assembles the vector for AFC stabilisation of general-purpose type
!#
!# 5.) afcsc_buildVectorSymm = afcsc_buildVectorSymmScalar /
!#                             afcsc_buildVectorSymmBlock
!#     -> Assembles the vector for AFC stabilisation of symmetric type
!#
!# 6.) afcsc_buildVectorLPT = afcsc_buildVecLPTScalar /
!#                            afcsc_buildVecLPTBlock
!#      -> Assembles the vector for stabilisation by means of
!#         linearity-preserving flux correction
!#
!# 7.) afcsc_buildFluxFCT = afcsc_buildFluxFCTScalar /
! #                         afcsc_buildFluxFCTBlock
!#     -> Assembles the raw antidiffusive flux for AFC stabilisation of FCT type
!#
!# 8.) afcsc_buildFluxLPT = afcsc_buildFluxLPTScalar /
!#                          afcsc_buildFluxLPTBlock
!#     -> Assembles the raw antidiffusive flux for the
!#        linearity-preserving stabilisation
!#
!# 9.) afcsc_buildJacobianFCT = afcsc_buildJacLinearFCTScalar /
!#                              afcsc_buildJacLinearFCTBlock /
!#                              afcsc_buildJacobianFCTScalar /
!#                              afcsc_buildJacobianFCTBlock
!#     -> Assembles the Jacobian matrix for the stabilisation part of FCT
!#        type; For the first two routines, the velocity is assumed
!#        to be linear which simplifies the evaluation of the
!#        Jacobian matrix significantly.  For the second two
!#        routines, the velocity can be arbitrary.
!#
!# 10.) afcsc_buildJacobianTVD = afcsc_buildJacLinearTVDScalar /
!#                               afcsc_buildJacLinearTVDBlock /
!#                               afcsc_buildJacobianTVDScalar /
!#                               afcsc_buildJacobianTVDBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of TVD
!#         type; For the first two routines, the velocity is assumed
!#         to be linear which simplifies the evaluation of the
!#         Jacobian matrix significantly.  For the second two
!#         routines, the velocity can be arbitrary.
!#
!# 11.) afcsc_buildJacobianGP = afcsc_buildJacLinearGPScalar /
!#                              afcsc_buildJacLinearGPBlock /
!#                              afcsc_buildJacobianGPScalar /
!#                              afcsc_buildJacobianGPBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of
!#         general purpose limiter. For the first two routines, the
!#         velocity is assumed to be linear which simplifies the
!#         evaluation of the Jacobian matrix significantly. For the
!#         second two routines, the velocity can be arbitrary.
!#
!# 12.) afcsc_buildJacobianSymm = afcsc_buildJacobianSymmScalar /
!#                                afcsc_buildJacobianSymmBlock
!#      -> Assembles the Jacobian matrix for the stabilisation part of
!#         symmetric flux limiting.
!#
!# 13.) afcsc_initPerfConfig
!#       -> Initialises the global performance configuration
!#
!# </purpose>
!##############################################################################

module afcstabscalar

  use afcstabbase
  use basicgeometry
  use collection
  use fsystem
  use genoutput
  use groupfembase
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use perfconfig
  use spatialdiscretisation
  use storage

  implicit none

  private
  
  public :: afcsc_initStabilisation
  public :: afcsc_initPerfConfig

  public :: afcsc_buildVectorFCT
  public :: afcsc_buildVectorTVD
  public :: afcsc_buildVectorGP
  public :: afcsc_buildVectorSymm
  public :: afcsc_buildVectorLPT

  public :: afcsc_buildFluxFCT
  public :: afcsc_buildFluxLPT

  public :: afcsc_buildJacobianFCT
  public :: afcsc_buildJacobianTVD
  public :: afcsc_buildJacobianGP
  public :: afcsc_buildJacobianSymm

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of nodes to handle simultaneously when building matrices
#ifndef AFCSC_NEQSIM
  integer, parameter, public :: AFCSC_NEQSIM = 128
#endif

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of edges to handle simultaneously when building matrices
#ifndef AFCSC_NEDGESIM
  integer, parameter, public :: AFCSC_NEDGESIM = 64
#endif

!</constantblock>

!</constants>

  !*****************************************************************************
  
  ! global performance configuration
  type(t_perfconfig), target, save :: afcsc_perfconfig

  !*****************************************************************************
  
  interface afcsc_initStabilisation
    module procedure afcsc_initStabByMatrix
    module procedure afcsc_initStabByGroupFEMSet
  end interface

  interface afcsc_buildVectorFCT
    module procedure afcsc_buildVectorFCTScalar
    module procedure afcsc_buildVectorFCTBlock
  end interface

  interface afcsc_buildVectorTVD
    module procedure afcsc_buildVectorTVDScalar
    module procedure afcsc_buildVectorTVDBlock
  end interface

  interface afcsc_buildVectorGP
    module procedure afcsc_buildVectorGPScalar
    module procedure afcsc_buildVectorGPBlock
  end interface

  interface afcsc_buildVectorSymm
    module procedure afcsc_buildVectorSymmScalar
    module procedure afcsc_buildVectorSymmBlock
  end interface

  interface afcsc_buildVectorLPT
    module procedure afcsc_buildVecLPTScalar
    module procedure afcsc_buildVecLPTBlock
  end interface

  interface afcsc_buildFluxFCT
    module procedure afcsc_buildFluxFCTScalar
    module procedure afcsc_buildFluxFCTBlock
  end interface

  interface afcsc_buildFluxLPT
    module procedure afcsc_buildFluxLPTScalar
    module procedure afcsc_buildFluxLPTBlock
  end interface

  interface afcsc_buildJacobianFCT
    module procedure afcsc_buildJacLinearFCTScalar
    module procedure afcsc_buildJacLinearFCTBlock
    module procedure afcsc_buildJacobianFCTScalar
    module procedure afcsc_buildJacobianFCTBlock
  end interface

  interface afcsc_buildJacobianTVD
    module procedure afcsc_buildJacLinearTVDScalar
    module procedure afcsc_buildJacLinearTVDBlock
    module procedure afcsc_buildJacobianTVDScalar
    module procedure afcsc_buildJacobianTVDBlock
  end interface

  interface afcsc_buildJacobianGP
    module procedure afcsc_buildJacLinearGPScalar
    module procedure afcsc_buildJacLinearGPBlock
    module procedure afcsc_buildJacobianGPScalar
    module procedure afcsc_buildJacobianGPBlock
  end interface
  
  interface afcsc_buildJacobianSymm
    module procedure afcsc_buildJacobianSymmScalar
    module procedure afcsc_buildJacobianSymmBlock
  end interface

contains

  !****************************************************************************

!<subroutine>

  subroutine afcsc_initPerfConfig(rperfconfig)

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
      afcsc_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(afcsc_perfconfig)
      afcsc_perfconfig%NEQSIM   = AFCSC_NEQSIM
      afcsc_perfconfig%NEDGESIM = AFCSC_NEDGESIM
    end if
  
  end subroutine afcsc_initPerfConfig

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_initStabByMatrix(rmatrix, rafcstab, rblockDiscretisation)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
!</description>

!<input>
    ! The template matrix
    type(t_matrixScalar), intent(in) :: rmatrix

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
    rafcstab%NVARtransformed = rmatrix%NVAR
    rafcstab%NVAR            = rmatrix%NVAR
    rafcstab%NEQ             = rmatrix%NEQ
    rafcstab%NEDGE           = (rmatrix%NA-rmatrix%NEQ)/2
    rafcstab%NNVEDGE         = 0
    
    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED
    
    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)
      
    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,6)
      call afcstab_genEdgeList(rmatrix, rafcstab)

      !-------------------------------------------------------------------------
      
    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE)
      
      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,6)
      call afcstab_genEdgeList(rmatrix, rafcstab)
      
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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,6)
      call afcstab_genEdgeList(rmatrix, rafcstab)

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
      
      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,6)
      call afcstab_genEdgeList(rmatrix, rafcstab)

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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,4)
      call afcstab_genEdgeList(rmatrix, rafcstab)

      ! Handle for DcoefficientsAtEdge: (/d_ij,s_ij/)
      call afcstab_allocCoeffsAtEdge(rafcstab,2)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the edgewise vector for the fluxes
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_MASS,&
          AFCSTAB_LINLPT_MASS)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,6)
      call afcstab_genEdgeList(rmatrix, rafcstab)

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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,6)
      call afcstab_genEdgeList(rmatrix, rafcstab)

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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      call afcstab_allocEdgeStructure(rafcstab,6)
      call afcstab_genEdgeList(rmatrix, rafcstab)

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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByMatrix')
      call sys_halt()
    end select

  end subroutine afcsc_initStabByMatrix
  
  !*****************************************************************************

!<subroutine>

  subroutine afcsc_initStabByGroupFEMSet(rgroupFEMSet, rafcstab, rblockDiscretisation)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The group finite element set
    ! is used to determine the number of equations and the number of
    ! edges. Common data structures are shared with the group FE set.
!</description>

!<input>
    ! The group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

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
    rafcstab%NVARtransformed = rgroupFEMSet%NVAR
    rafcstab%NVAR            = rgroupFEMSet%NVAR
    rafcstab%NEQ             = rgroupFEMSet%NEQ
    rafcstab%NEDGE           = rgroupFEMSet%NEDGE
    rafcstab%NNVEDGE         = 0
    
    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED
    
    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)
      
    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND,&
          AFCSTAB_DMP)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if

      !-------------------------------------------------------------------------
      
    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE)
      
      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if
      
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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if

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
      
      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if

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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if

      ! Handle for DcoefficientsAtEdge: (/d_ij,s_ij/)
      call afcstab_allocCoeffsAtEdge(rafcstab,2)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the edgewise vector for the fluxes
      call afcstab_allocFlux(rafcstab)

      !-------------------------------------------------------------------------

    case (AFCSTAB_NLINLPT_MASS,&
          AFCSTAB_LINLPT_MASS)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if

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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if

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

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
      if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
        rafcstab%h_IedgeListIdx   = rgroupFEMSet%h_IedgeListIdx
        rafcstab%h_IedgeList      = rgroupFEMSet%h_IedgeList
        rafcstab%iduplicationFlag = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
        rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
            iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
      else
        call output_line('Group finite element set does not provide edge structure',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
        call sys_halt()
      end if

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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_initStabByGroupFEMSet')
      call sys_halt()
    end select

  end subroutine afcsc_initStabByGroupFEMSet

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorFCTBlock(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-FCT type. Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
    
    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry
    !</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks .ne. 1 .or. ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTBlock')
      call sys_halt()

    else
      
      call afcsc_buildVectorFCTScalar(rafcstab, rmatrix, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1), rperfconfig)
      
    end if
    
  end subroutine afcsc_buildVectorFCTBlock

  ! *****************************************************************************
  
!<subroutine>
  
  subroutine afcsc_buildVectorFCTScalar(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-FCT type. The idea of flux corrected transport can be
    ! traced back to the early SHASTA algorithm by Boris and Bock in
    ! the early 1970s. Zalesak suggested a fully multi-dimensional
    ! generalisation of this approach and paved the way for a large
    ! family of FCT algorithms.
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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
    
    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_DfluxPrel
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
      call sys_halt()
    end if

    ! Clear vector?
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
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
          call sys_halt()
        end if
        
        ! Check if stabilisation provides edge-based structure
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide edge structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix positions
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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

  end subroutine afcsc_buildVectorFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorTVDBlock(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-TVD type.  Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the corresponding
    ! scalar routine is called. Otherwise, an error is thrown.
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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry    
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks   .ne. 1 .or.&
        ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDBlock')
      call sys_halt()

    else

      call afcsc_buildVectorTVDScalar(rafcstab, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1), rperfconfig)
      
    end if

  end subroutine afcsc_buildVectorTVDBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorTVDScalar(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-TVD type.
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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
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
    
    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq.0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq.0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq.0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          
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
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
            call sys_halt()
          end if
          
          ! Assemble antidiffusive fluxes, sums of antidiffusive
          ! increments and the upper and lower bounds simultaneously
          
          !$omp parallel default(shared)&
          !$omp private(i,j,f_ij,fm_ij,fp_ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          
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
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        
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
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        ! Apply the nodal limiter
        
        !$omp parallel sections
        !$omp section
        Drp = afcstab_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)

        !$omp section
        Drm = afcstab_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)
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
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorTVDScalar')
          call sys_halt()
        end if

        !$omp parallel default(shared)&
        !$omp private(i,j,d_ij,l_ji,diff,f_ij)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
        
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
    
  end subroutine afcsc_buildVectorTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorGPBlock(rafcstab, rmatrix,&
      rx, rx0, theta, dscale, bclear, ioperationSpec, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! of FEM-GP type.  Note that this routine serves as a wrapper for
    ! block vectors. If there is only one block, then the corresponding
    ! scalar routine is called.  Otherwise, an error is thrown.
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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry   
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if ((rx%nblocks   .ne. 1) .or.&
        (rx0%nblocks  .ne. 1) .or.&
        (ry%nblocks .ne. 1)) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorGPBlock')
      call sys_halt()

    else

      call afcsc_buildVectorGPScalar(rafcstab, rmatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), theta, dscale, bclear, ioperationSpec,&
          ry%RvectorBlock(1), rperfconfig)
      
    end if
  end subroutine afcsc_buildVectorGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorGPScalar(rafcstab, rmatrix,&
      rx, rx0, theta, dscale, bclear, ioperationSpec, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! using the general purpose limiter.
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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_MC,p_Dx,p_Dx0,p_Dy,p_Dflux,p_Dflux0
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    
    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if
    
    ! Check if stabilisation is prepared
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorGPScalar')
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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
      Drp = afcstab_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)
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
    
  end subroutine afcsc_buildVectorGPScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorSymmBlock(rafcstab, rx, dscale, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! by means of symmetric flux limiting. Note that this routine
    ! serves as a wrapper for block vectors. If there is only one
    ! block, then the corresponding scalar routine is
    ! called. Otherwise, an error is thrown.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks .ne. 1 .or. ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorSymmBlock')
      call sys_halt()

    else

      call afcsc_buildVectorSymmScalar(rafcstab,&
          rx%RvectorBlock(1), dscale, ry%RvectorBlock(1), rperfconfig)

    end if
  end subroutine afcsc_buildVectorSymmBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVectorSymmScalar(rafcstab, rx, dscale, ry, rperfconfig)

!<description>
    ! This subroutine assembles the vector and applies stabilisation
    ! by means of symmetric flux limiting.
    !
    ! Yet, there is no publication available. This routine is based on
    ! private communication with D. Kuzmin.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dx,p_Dy,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    
    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if
    
    ! Check if stabilisation is prepared
    if ((rafcstab%cafcstabType .ne. AFCSTAB_SYMMETRIC) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)      .eq. 0)    .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)    .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorSymmScalar')
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      
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
      Drp = afcstab_limit(Dpp, Dqp, 0.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 0.0_DP, 1.0_DP)
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

  end subroutine afcsc_buildVectorSymmScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVecLPTBlock(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rmatrix, rperfconfig)

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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry    
!</inputoutput>
!</subroutine>

    ! Check if block vectors contain exactly one block
    if (rx%nblocks .ne. 1 .or.&
        ry%nblocks .ne. 1) then

      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTBlock')
      call sys_halt()

    else

      call afcsc_buildVecLPTScalar(rafcstab, rx%RvectorBlock(1),&
          dscale, bclear, ioperationSpec, ry%RvectorBlock(1),&
          rmatrix, rperfconfig)
      
    end if

  end subroutine afcsc_buildVecLPTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildVecLPTScalar(rafcstab, rx, dscale,&
      bclear, ioperationSpec, ry, rmatrix, rperfconfig)

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
    
    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry    
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy,p_Dq
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if

    ! Check if stabilisation is prepeared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
      call sys_halt()
    end if

    ! Clear vector?
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
    ! 6) Apply the limited antidifusive fluxes
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVectorFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildVecLPTScalar')
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)
          
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

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

  end subroutine afcsc_buildVecLPTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildFluxFCTBlock(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, fcb_calcFluxFCTSc_sim,&
      rgroupFEMSet, rmatrix, rxTimeDeriv, rxPredictor, rcollection, rperfconfig)

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

    ! OPTIONAL: callback functions to compute antidiffusive fluxes
    include 'intf_calcFluxFCTSc_sim.inc'
    optional :: fcb_calcFluxFCTSc_sim

    ! OPTIONAL: group finite element set
    type(t_groupFEMSet), intent(in), optional :: rgroupFEMSet

    ! OPTIONAL: Consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorBlock), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorBlock), intent(in), optional :: rxPredictor

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    integer :: nblocks

    ! Check if block vector(s) contains exactly one block
    nblocks = rx%nblocks
    if (present(rxTimeDeriv)) nblocks = max(nblocks, rxTimeDeriv%nblocks)
    if (present(rxPredictor)) nblocks = max(nblocks, rxPredictor%nblocks)

    if (nblocks .ne. 1) then
      call output_line('Vector(s) must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Call subroutine for scalar vectors
    if (present(rxTimeDeriv)) then
      if (present(rxPredictor)) then
        ! ... both approximate time derivative and predictor are present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rxTimeDeriv%RvectorBlock(1), rxPredictor%RvectorBlock(1), &
            rcollection, rperfconfig)
      else
        ! ... only the approximate time derivative is present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rxTimeDeriv=rxTimeDeriv%RvectorBlock(1),&
            rcollection=rcollection, rperfconfig=rperfconfig)
      end if
    else
      if (present(rxPredictor)) then
        ! ... only the predictor is present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rxPredictor=rxPredictor%RvectorBlock(1),&
            rcollection=rcollection, rperfconfig=rperfconfig)
      else
        ! ... neither the approximate time derivative nor the predictor is present
        call afcsc_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
            theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
            fcb_calcFluxFCTSc_sim, rgroupFEMSet, rmatrix,&
            rcollection=rcollection, rperfconfig=rperfconfig)
      end if
    end if
    
  end subroutine afcsc_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>


  subroutine afcsc_buildFluxFCTScalar(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, fcb_calcFluxFCTSc_sim,&
      rgroupFEMSet, rmatrix, rxTimeDeriv, rxPredictor, rcollection, rperfconfig)
    
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

    ! OPTIONAL: callback functions to compute antidiffusive fluxes
    include 'intf_calcFluxFCTSc_sim.inc'
    optional :: fcb_calcFluxFCTSc_sim

    ! OPTIONAL: group finite element set
    type(t_groupFEMSet), intent(in), optional :: rgroupFEMSet

    ! OPTIONAL: consistent mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorScalar), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorScalar), intent(in), optional :: rxPredictor

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dmatrix,p_Dx
    real(DP), dimension(:), pointer :: p_DxTimeDeriv, p_DxPredictor
    real(DP), dimension(:), pointer :: p_Dflux0,p_Dflux,p_DfluxPrel,p_Dalpha
    real(DP), dimension(:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    logical :: buseCallback

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
      call sys_halt()
    end if
    
    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge data structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
        call sys_halt()
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call lsyssc_getbase_double(rx, p_Dx)

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)
      
    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_IMPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE)

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)
    
      ! Use callback routine?
      if (present(fcb_calcFluxFCTSc_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        buseCallback = .true.
      else
        
        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if
        
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        buseCallback = .false.       
      end if

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
      
      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_EXPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble explicit part of raw-antidiffive fluxes
        !-----------------------------------------------------------------------
                
        if (theta .ne. 1.0_DP) then
          ! Assemble the explicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij}^n = (1-\theta)\Delta t d_{ij}^n(u_i^n-u_j^n) $$
          if (buseCallback) then
            call doFluxesByCallbackDble(p_IedgeList, rafcstab%NEDGE,&
                p_DcoeffsAtEdge, p_Dx, dscale*(1.0_DP-theta),&
                bclear, p_Dflux0)
          else
            call doFluxesByCoeffsDble(p_IedgeList, rafcstab%NEDGE,&
                p_Dcoefficients, p_Dx, dscale*(1.0_DP-theta),&
                bclear, p_Dflux0)
          end if
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
          if (buseCallback) then
            call doFluxesByCallbackDble(p_IedgeList, rafcstab%NEDGE,&
                p_DcoeffsAtEdge, p_Dx, dscale, .true., p_DfluxPrel)
          else
            call doFluxesByCoeffsDble(p_IedgeList, rafcstab%NEDGE,&
                p_Dcoefficients, p_Dx, dscale, .true., p_DfluxPrel)
          end if

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
              if (buseCallback) then
                call doFluxesByCallbackDble(p_IedgeList, rafcstab%NEDGE,&
                    p_DcoeffsAtEdge, p_DxPredictor, dscale, .true., p_DfluxPrel)
              else
                call doFluxesByCoeffsDble(p_IedgeList, rafcstab%NEDGE,&
                    p_Dcoefficients, p_DxPredictor, dscale, .true., p_DfluxPrel)
              end if
            else
              call output_line('Invalid type of prelimiting!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
              call sys_halt()
            end if
          else
            call output_line('Fluxes for prelimiting cannot be assembled without predictor!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
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
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
            call sys_halt()
          end if
          
          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
            call sys_halt()
          end if
          
          ! Set pointer
          call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
          
          ! Subtract amount of rejected antidiffusion
          call afcstab_combineFluxes(rafcstab%NEDGE, -1.0_DP, p_Dflux, p_Dflux0, p_Dalpha)
        end if
        
        !-----------------------------------------------------------------------
        
        if (theta .ne. 0.0_DP) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij} = \theta\Delta t d_{ij}(u_i-u_j) $$
          if (buseCallback) then
            call doFluxesByCallbackDble(p_IedgeList, rafcstab%NEDGE,&
                p_DcoeffsAtEdge, p_Dx, dscale*theta, bclear, p_Dflux)
          else
            call doFluxesByCoeffsDble(p_IedgeList, rafcstab%NEDGE,&
                p_Dcoefficients, p_Dx, dscale*theta, bclear, p_Dflux)
          end if
        end if
        
        if (bquickAssembly) then
          ! We may check of either the implicit or explicit part are
          ! missing so that some redundant computations may be skipped
          if (theta .ne. 1.0_DP) then
            ! The explicit part of the raw-antidiffusive fluxes exists
            if (theta .ne. 0.0_DP) then
              ! The implicit part of the raw-antidiffusive fluxes
              ! exists; so combine them both into common fluxes
              call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_DP, p_Dflux0, p_Dflux)
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
          call afcstab_combineFluxes(rafcstab%NEDGE, 1.0_DP, p_Dflux0, p_Dflux)
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
      
      ! Set pointer
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      
      ! Use callback routine?
      if (present(fcb_calcFluxFCTSc_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        buseCallback = .true.
      else
        
        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
          call sys_halt()
        end if
        
        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        buseCallback = .false.
      end if

      ! Assemble spatial part of raw-antidiffusive fluxes
      if (buseCallback) then
        call doFluxesByCallbackDble(p_IedgeList, rafcstab%NEDGE,&
            p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Dflux)
      else
        call doFluxesByCoeffsDble(p_IedgeList, rafcstab%NEDGE,&
            p_Dcoefficients, p_Dx, dscale, bclear, p_Dflux)
      end if
      
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
        call sys_halt()
      end if
      
      
    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxFCTScalar')
      call sys_halt()
    end select

  contains
    
    ! Here, the working routines follow

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the edge-by-edge array DcoefficientsAtEdge
    
    subroutine doFluxesByCoeffsDble(IedgeList, NEDGE,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)
      
      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear
      
      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dflux
      
      ! local variables
      integer :: iedge,i,j
      
      if (dscale .eq. 0.0_DP) then
        
        if (bclear) call lalg_clearVector(Dflux, NEDGE)
        
      elseif (dscale .eq. 1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if
        
      elseif (dscale .eq. -1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        end if
        
      else
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Compute raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if
        
      end if
      
    end subroutine doFluxesByCoeffsDble

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix stored in Dmatrix
    
    subroutine doFluxesByMatrixDble(IedgeList, NEDGE,&
        Dmatrix, Dx, dscale, bclear, Dflux)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dmatrix, Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear
      
      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dflux
    
      ! local variables
      integer :: iedge,i,j,ij
      
      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
    ! Assemble raw antidiffusive fluxes with aid of callback function
    
    subroutine doFluxesByCallbackDble(IedgeList, NEDGE,&
        DcoeffsAtEdge, Dx, dscale, bclear, Dflux)
      
      ! input parameters
      real(DP), dimension(:), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dflux

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtEdge
      real(DP), dimension(:), pointer :: DfluxAtEdge
      
      ! local variables
      integer :: idx,iedge,IEDGEset,IEDGEmax


      if (dscale .eq. 0.0_DP) then
        
        if (bclear) call lalg_clearVector(Dflux, NEDGE)

      elseif (bclear) then

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,idx,iedge,IEDGEmax)
        
        ! Allocate temporal memory
        allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))
        
        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, p_rperfconfig%NEDGESIM
          
          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most  edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+p_rperfconfig%NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dflux(IEDGEset:IEDGEmax), rcollection)
        end do
        !$omp end do
        
        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        !$omp end parallel

      else   ! bclear = .false.

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,DfluxAtEdge,idx,iedge,IEDGEmax)
        
        ! Allocate temporal memory
        allocate(DdataAtEdge(2,p_rperfconfig%NEDGESIM))
        allocate(DfluxAtEdge(p_rperfconfig%NEDGESIM))
        
        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, p_rperfconfig%NEDGESIM
          
          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most  edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+p_rperfconfig%NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(1,idx) = Dx(IedgeList(1,iedge))
            DdataAtEdge(2,idx) = Dx(IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSc_sim(&
              DdataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxAtEdge(1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Add antidiffusive fluxes            
            Dflux(iedge) = Dflux(iedge) + DfluxAtEdge(idx)
          end do
        end do
        !$omp end do
        
        ! Deallocate temporal memory
        deallocate(DdataAtEdge, DfluxAtEdge)
        !$omp end parallel

      end if

    end subroutine doFluxesByCallbackDble
    
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
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE
        
        ! Determine indices
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        
        ! Compute solution difference; in contrast to the literature,
        ! we compute the solution difference $u_i-u_j$ and check if
        ! $f_{ij}(u_i-u_j)<0$ in the prelimiting step.
        Dflux(iedge) = Dx(i)-Dx(j)
      end do
      !$omp end parallel do

    end subroutine doDifferencesDble
       
  end subroutine afcsc_buildFluxFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildFluxLPTBlock(rafcstab, rx, dscale, bclear,&
      ioperationSpec, rmatrix, rperfconfig)

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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    if (rx%nblocks .ne. 1) then
      call output_line('Vector must not contain more than one block!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTBlock')
      call sys_halt()
    end if

    ! Call subroutine for scalar vectors
    call afcsc_buildFluxLPTScalar(rafcstab, rx%RvectorBlock(1),&
        dscale, bclear, ioperationSpec, rmatrix, rperfconfig)
    
  end subroutine afcsc_buildFluxLPTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildFluxLPTScalar(rafcstab, rx, dscale, bclear,&
      ioperationSpec, rmatrix, rperfconfig)

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

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
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
    
    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig
    
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsc_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge data structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
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
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if

      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
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
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if
        
      else
        call output_line('Stabilisation does not provide (oriented) edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
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
                OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
            call sys_halt()
          end if
        end if

      else
        call output_line('Stabilisation does not provide edge data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
        call sys_halt()
      end if
      
    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildFluxLPTScalar')
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

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
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        end if
        
      elseif (dscale .eq. -1.0_DP) then
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Update raw antidiffusive flux
            Dflux(iedge) = Dflux(iedge)&
                + DcoefficientsAtEdge(1,iedge) * (Dx(j)-Dx(i))
          end do
          !$omp end parallel do
        end if
        
      else
        
        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
            ! Store raw antidiffusive flux
            Dflux(iedge) = dscale * DcoefficientsAtEdge(1,iedge) * (Dx(i)-Dx(j))
          end do
          !$omp end parallel do
        else
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE
            
            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)
            
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
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Update bounds
          Dq(i) = Dq(i) + DboundsAtEdge(1,iedge)*DcoefficientsAtEdge(1,iedge)
          Dq(j) = Dq(j) + DboundsAtEdge(2,iedge)*DcoefficientsAtEdge(1,iedge)
        end do
        !$omp end do
      end do
      !$omp end parallel
      
    end subroutine doBoundsByCoeffDble

  end subroutine afcsc_buildFluxLPTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearFCTBlock(rx, theta, tstep, hstep,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTBlock')
      call sys_halt()

    else

      call afcsc_buildJacLinearFCTScalar(&
          rx%RvectorBlock(1), theta, tstep, hstep, bclear,&
          rafcstab, rjacobian, rmatrix)

    end if
  end subroutine afcsc_buildJacLinearFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearFCTScalar(rx, theta, tstep, hstep,&
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
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)   .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTScalar')
        call sys_halt()
      end select
      
    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearFCTScalar')
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

  end subroutine afcsc_buildJacLinearFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearTVDBlock(rx, tstep, hstep,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearTVDBlock')
      call sys_halt()

    else

      call afcsc_buildJacLinearTVDScalar(rx%RvectorBlock(1), tstep,&
          hstep, bclear, rafcstab, rjacobian, bextendedSparsity)

    end if
  end subroutine afcsc_buildJacLinearTVDBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearTVDScalar(rx, tstep, hstep,&
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
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearTVDScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearTVDScalar')
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
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    pure subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

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

  end subroutine afcsc_buildJacLinearTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearGPBlock(rmatrix, rx,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearGPBlock')
      call sys_halt()

    else

      call afcsc_buildJacLinearGPScalar(&
          rmatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), theta, tstep, hstep,&
          bclear, rafcstab, rjacobian, bextendedSparsity)
      
    end if
  end subroutine afcsc_buildJacLinearGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacLinearGPScalar(rmatrix, rx,&
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
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearGPScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacLinearGPScalar')
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
    ! Based on the matric structure given by Kld/Kcol, the separator
    ! is moved to the given column k. For efficiency reasons, only
    ! those entries are considered which are present in column k.
    pure subroutine adjustKsepMat9(Kld, Kcol, k, Ksep)
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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


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
  end subroutine afcsc_buildJacLinearGPScalar

  !*****************************************************************************
  
!<subroutine>

  subroutine afcsc_buildJacobianFCTBlock(RcoeffMatrices, rx,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTBlock')
      call sys_halt()

    else

      call afcsc_buildJacobianFCTScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrixSc_sim,&
          theta, tstep, hstep, bclear, rafcstab, rjacobian,&
          rmatrix, rcollection)

    end if
  end subroutine afcsc_buildJacobianFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianFCTScalar(RcoeffMatrices, rx,&
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
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)   .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid type of AFC stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianFCTScalar')
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

  end subroutine afcsc_buildJacobianFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianTVDBlock(RcoeffMatrices, rx,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDBlock')
      call sys_halt()

    else

      call afcsc_buildJacobianTVDScalar(&
          RcoeffMatrices, rx%RvectorBlock(1), fcb_calcMatrixSc_sim,&
          tstep, hstep, bclear, rafcstab, rjacobian,&
          bextendedSparsity, rcollection)

    end if
  end subroutine afcsc_buildJacobianTVDBlock
  
  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianTVDScalar(RcoeffMatrices, rx,&
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
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianTVDScalar')
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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)

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
    
  end subroutine afcsc_buildJacobianTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianGPBlock(RcoeffMatrices, rmatrix,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianGPBlock')
      call sys_halt()
      
    else
      
      call afcsc_buildJacobianGPScalar(&
          RcoeffMatrices, rmatrix, rx%RvectorBlock(1),&
          rx0%RvectorBlock(1), fcb_calcMatrixSc_sim, theta, tstep, hstep,&
          bclear, rafcstab, rjacobian,bextendedSparsity, rcollection)

    end if
  end subroutine afcsc_buildJacobianGPBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianGPScalar(RcoeffMatrices, rmatrix,&
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
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS)    .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)      .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianGPScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianGPScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianGPScalar')
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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


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
    
  end subroutine afcsc_buildJacobianGPScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianSymmBlock(rx, dscale, hstep, bclear,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianSymmBlock')
      call sys_halt()

    else

      call afcsc_buildJacobianSymmScalar(rx%RvectorBlock(1), dscale,&
          hstep, bclear, rafcstab, rjacobian, bextendedSparsity)

    end if
  end subroutine afcsc_buildJacobianSymmBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsc_buildJacobianSymmScalar(rx, dscale, hstep, bclear,&
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
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)     .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)   .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)     .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianSymmScalar')
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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsc_buildJacobianSymmScalar')
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
        Drploc(:,0:nloc) = afcstab_limit(Dpploc(:,0:nloc), Dqploc(:,0:nloc), 0.0_DP, 1.0_DP)
        Drmloc(:,0:nloc) = afcstab_limit(Dpmloc(:,0:nloc), Dqmloc(:,0:nloc), 0.0_DP, 1.0_DP)


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

  end subroutine afcsc_buildJacobianSymmScalar

end module afcstabscalar

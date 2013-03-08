!##############################################################################
!# ****************************************************************************
!# <name> afcstabsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides the basic routines for applying the algebraic
!# flux correction methodology proposed by Kuzmin, Moeller and Turek
!# in a series of publications. As a starting point for systems of
!# conservation laws, the reader is referred to the book chapter
!#
!#     D. Kuzmin and M. Moeller, Algebraic flux correction
!#     II. Compressible Euler Equations, In: D. Kuzmin et al. (eds),
!#     Flux-Corrected Transport: Principles, Algorithms, and
!#     Applications, Springer, 2005, 207-250.
!#
!# A more detailed description of the algorithms is given in the
!# comments of the subroutine implementing the corresponding
!# discretisation schemes. All methods are based on the stabilisation
!# structure t_afcstab which is defined in the underlying module
!# afcstabbase. The initialisation as a scalar stabilisation
!# structure is done by the routine afcsys_initStabilisation.
!#
!# The following routines are available:
!#
!# 1.) afcsys_initStabilisation = afcsys_initStabByMatrixSc /
!#                                afcsys_initStabByMatrixBl /
!#                                afcsys_initStabByGroupFEMSetSc /
!#                                afcsys_initStabByGroupFEMSetBl
!#     -> Initialises the stabilisation structure
!#
!# 2.) afcsys_buildVectorFCT = afcsys_buildVectorFCTScalar /
!#                             afcsys_buildVectorFCTBlock
!#     -> Assembles the vector for AFC stabilisation of FCT type
!#
!# 3.) afcsys_buildVectorTVD = afcsys_buildVectorTVDScalar /
!#                             afcsys_buildVectorTVDBlock
!#     -> Assembles the vector for AFC stabilisation of TVD type
!#
!# 4.) afcsys_buildFluxFCT = afcsys_buildFluxFCTScalar /
!#                           afcsys_buildFluxFCTBlock
!#     -> Assembles the raw antidiffusive flux for FEM-FCT stabilisation
!#
!# 5.) afcsys_failsafeFCT = afcsys_failsafeFCTScalar /
!#                          afcsys_failsafeFCTBlock
!#     -> perform failsafe limiting of FCT type
!#
!# 6.) afcsys_initPerfConfig
!#       -> Initialises the global performance configuration
!#
!#
!# The following auxiliary routines are available:
!#
!# 1.) minmodDP / minmodSP
!#     -> Computes minmod(a,b)
!#
!# </purpose>
!##############################################################################

module afcstabsystem

!$use omp_lib
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

  public :: afcsys_initStabilisation
  public :: afcsys_initPerfConfig

  public :: afcsys_buildVectorFCT
  public :: afcsys_buildVectorTVD

  public :: afcsys_buildFluxFCT
  public :: afcsys_failsafeFCT

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of nodes to handle simultaneously when building matrices
#ifndef AFCSYS_NEQSIM
  integer, parameter, public :: AFCSYS_NEQSIM = 128
#endif

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of edges to handle simultaneously when building matrices
#ifndef AFCSYS_NEDGESIM
  integer, parameter, public :: AFCSYS_NEDGESIM = 64
#endif
!</constantblock>

!<constantblock description="Global constants for directional splitting">

  ! unit vector in X-direction in 1D
  real(DP), dimension(NDIM1D) :: XDir1D = (/ 1.0_DP /)

  ! unit vector in X-direction in 2D
  real(DP), dimension(NDIM2D) :: XDir2D = (/ 1.0_DP, 0.0_DP /)

  ! unit vector in Y-direction in 2D
  real(DP), dimension(NDIM2D) :: YDir2D = (/ 0.0_DP, 1.0_DP /)

  ! unit vector in X-direction in 3D
  real(DP), dimension(NDIM3D) :: XDir3D = (/ 1.0_DP, 0.0_DP, 0.0_DP /)

  ! unit vector in Y-direction in 3D
  real(DP), dimension(NDIM3D) :: YDir3D = (/ 0.0_DP, 1.0_DP, 0.0_DP /)

  ! unit vector in Z-direction in 3D
  real(DP), dimension(NDIM3D) :: ZDir3D = (/ 0.0_DP, 0.0_DP, 1.0_DP /)

!</constantblock>

!</constants>

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: afcsys_perfconfig

  ! ****************************************************************************

  interface afcsys_initStabilisation
    module procedure afcsys_initStabByMatrixSc
    module procedure afcsys_initStabByMatrixBl
    module procedure afcsys_initStabByGroupFEMSetSc
    module procedure afcsys_initStabByGroupFEMSetBl
  end interface

  interface afcsys_buildVectorFCT
    module procedure afcsys_buildVectorFCTScalar
    module procedure afcsys_buildVectorFCTBlock
  end interface

  interface afcsys_buildVectorTVD
    module procedure afcsys_buildVectorTVDScalar
    module procedure afcsys_buildVectorTVDBlock
  end interface

  interface afcsys_buildFluxFCT
    module procedure afcsys_buildFluxFCTScalar
    module procedure afcsys_buildFluxFCTBlock
  end interface

  interface afcsys_failsafeFCT
    module procedure afcsys_failsafeFCTScalar
    module procedure afcsys_failsafeFCTBlock
  end interface

contains

  !****************************************************************************

!<subroutine>

  subroutine afcsys_initPerfConfig(rperfconfig)

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
      afcsys_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(afcsys_perfconfig)
      afcsys_perfconfig%NEQSIM   = AFCSYS_NEQSIM
      afcsys_perfconfig%NEDGESIM = AFCSYS_NEDGESIM
    end if

  end subroutine afcsys_initPerfConfig

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_initStabByMatrixBl(rmatrix, rafcstab,&
      rblockDiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as block matrix. If this matrix
    ! contains only one block, then the scalar counterpart of this
    ! subroutine is called with the corresponding scalar submatrix.
!</description>

!<input>
    ! template block matrix
    type(t_matrixBlock), intent(in) :: rmatrix

    ! OPTIONAL: block discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_blockDiscretisation), intent(in), optional :: rblockDiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Check if block matrix has only one block
    if ((rmatrix%nblocksPerCol .eq. 1) .and.&
        (rmatrix%nblocksPerRow .eq. 1)) then
      if (present(rblockDiscretisation)) then
        call afcsys_initStabByMatrixSc(rmatrix%RmatrixBlock(1,1), rafcstab,&
            rblockDiscretisation%RspatialDiscr(1), NVARtransformed)
      else
        call afcsys_initStabByMatrixSc(rmatrix%RmatrixBlock(1,1), rafcstab,&
            NVARtransformed=NVARtransformed)
      end if
      return
    end if

    ! Check that number of columns equations number of rows
    if (rmatrix%nblocksPerCol .ne.&
        rmatrix%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByMatrixBl')
      call sys_halt()
    end if

    ! Check if matrix exhibits group structure
    if (rmatrix%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByMatrixBl')
      call sys_halt()
    end if

    ! Set atomic data from first block
    rafcstab%NVAR    = rmatrix%nblocksPerCol
    rafcstab%NEQ     = rmatrix%RmatrixBlock(1,1)%NEQ
    rafcstab%NEDGE   = (rmatrix%RmatrixBlock(1,1)%NA-&
                        rmatrix%RmatrixBlock(1,1)%NEQ)/2
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rafcstab%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    call afcstab_allocEdgeStructure(rafcstab,6)
    call afcstab_genEdgeList(rmatrix%RmatrixBlock(1,1), rafcstab)

    ! Allocate internal data structures and
    call afcstab_allocInternalData(rafcstab, .false.,&
        rblockDiscretisation=rblockDiscretisation)

  end subroutine afcsys_initStabByMatrixBl

  ! ****************************************************************************

!<subroutine>
  subroutine afcsys_initStabByMatrixSc(rmatrix, rafcstab,&
      rdiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as scalar matrix. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Set atomic data
    rafcstab%NVAR    = rmatrix%NVAR
    rafcstab%NEQ     = rmatrix%NEQ
    rafcstab%NEDGE   = (rmatrix%NA-rmatrix%NEQ)/2
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rmatrix%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    call afcstab_allocEdgeStructure(rafcstab,6)
    call afcstab_genEdgeList(rmatrix, rafcstab)

    ! Allocate internal data structures
    call afcstab_allocInternalData(rafcstab, .false.,&
        rdiscretisation=rdiscretisation)

  end subroutine afcsys_initStabByMatrixSc

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_initStabByGroupFEMSetBl(rgroupFEMSet, rafcstab,&
      rblockDiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The group finite element set
    ! is used to determine the number of equations and the number of
    ! edges. Common data structures are shared with the group FE set.
!</description>

!<input>
    ! The group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! Block discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_blockDiscretisation), intent(in) :: rblockDiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Check if block discretisation has only one block
    if (rblockDiscretisation%ncomponents .eq. 1) then
      call afcsys_initStabByGroupFEMSetSc(rgroupFEMSet, rafcstab,&
          rblockDiscretisation%RspatialDiscr(1), NVARtransformed)
      return
    end if

    ! Set atomic data
    rafcstab%NVAR    = rgroupFEMSet%NVAR
    rafcstab%NEQ     = rgroupFEMSet%NEQ
    rafcstab%NEDGE   = rgroupFEMSet%NEDGE
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rafcstab%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

    ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji,ii,jj/)
    if (iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST) .ne. 0) then
      rafcstab%h_IedgeListIdx     = rgroupFEMSet%h_IedgeListIdx
      rafcstab%h_IedgeList        = rgroupFEMSet%h_IedgeList
      rafcstab%iduplicationFlag   = ior(rafcstab%iduplicationFlag,&
                                        AFCSTAB_SHARE_EDGELIST)
      rafcstab%istabilisationSpec = ior(rafcstab%istabilisationSpec,&
          iand(rgroupFEMSet%isetSpec, GFEM_HAS_EDGELIST))
    else
      call output_line('Group finite element set does not provide edge structure',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByGroupFEMSetBl')
      call sys_halt()
    end if

    ! Allocate internal data structures
    call afcstab_allocInternalData(rafcstab, .false.,&
        rblockDiscretisation=rblockDiscretisation)

  end subroutine afcsys_initStabByGroupFEMSetBl

!*****************************************************************************

!<subroutine>

  subroutine afcsys_initStabByGroupFEMSetSc(rgroupFEMSet, rafcstab,&
      rdiscretisation, NVARtransformed)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a system stabilisation. The group finite element set
    ! is used to determine the number of equations and the number of
    ! edges. Common data structures are shared with the group FE set.
!</description>

!<input>
    ! The group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Set atomic data
    rafcstab%NVAR    = rgroupFEMSet%NVAR
    rafcstab%NEQ     = rgroupFEMSet%NEQ
    rafcstab%NEDGE   = rgroupFEMSet%NEDGE
    rafcstab%NNVEDGE = 0

    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rafcstab%NVAR
    end if

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

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
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_initStabByGroupFEMSetSc')
      call sys_halt()
    end if

    ! Allocate internal data structures
    call afcstab_allocInternalData(rafcstab, .false.,&
        rdiscretisation=rdiscretisation)

  end subroutine afcsys_initStabByGroupFEMSetSc

  ! ****************************************************************************

!<subroutine>

  subroutine afcsys_buildVectorFCTBlock(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry, NVARtransformed,&
      fcb_calcFluxTransformation_sim, fcb_calcDiffTransformation_sim,&
      fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
      fcb_limitEdgewise, fcb_calcCorrection, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector for nonlinear FEM-FCT
    ! schemes. If the vectors contain only one block, then the scalar
    ! counterpart of this routine is called with the scalar
    ! subvectors. Consider the documentation of subroutine
    ! 'afcsys_buildVectorFCTScalar' for further details.
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
    ! combination of different AFCSTAB_FCT_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVARtransformed is taken from the stabilisation structure
    integer, intent(in), optional :: NVARtransformed

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_calcFluxTransformation_sim.inc'
    optional :: fcb_calcFluxTransformation_sim

    include 'intf_calcDiffTransformation_sim.inc'
    optional :: fcb_calcDiffTransformation_sim

    ! OPTIONAL: callback functions to overwrite the standard operations
    include 'intf_calcADIncrements.inc'
    optional :: fcb_calcADIncrements

    include 'intf_calcBounds.inc'
    optional :: fcb_calcBounds

    include 'intf_limitNodal.inc'
    optional :: fcb_limitNodal

    include 'intf_limitEdgewise.inc'
    optional :: fcb_limitEdgewise

    include 'intf_calcCorrection.inc'
    optional :: fcb_calcCorrection

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_DfluxPrel
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: nvariable

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1)) then
      call afcsys_buildVectorFCTScalar(&
          rafcstab, rmatrix, rx%RvectorBlock(1), dscale, bclear,&
          ioperationSpec, ry%RvectorBlock(1), NVARtransformed,&
          fcb_calcFluxTransformation_sim, fcb_calcDiffTransformation_sim,&
          fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
          fcb_limitEdgewise, fcb_calcCorrection, rcollection, rperfconfig)
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBLock')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)

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
    call lsysbl_getbase_double(rx, p_Dx)
    call lsysbl_getbase_double(ry, p_Dy)

    !---------------------------------------------------------------------------
    ! The nonlinear FEM-FCT algorithm is split into the following
    ! steps which can be skipped and performed externally by the user:
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
    ! 6) Apply the limited antidifusive fluxes to the vector
    !
    !    Step 6) may be split into the following substeps
    !
    !    6.1) Compute the edgewise correction factors based on the pre-
    !         computed raw-antidiffusive fluxes.
    !
    !    6.2) Compute the raw antidiffusive fluxes for a different set of
    !         variables and limit them by the precomputed correction factors.
    !-------------------------------------------------------------------------

    ! Determine number of transformed variables (if any)
    if (present(NVARtransformed)) then
      nvariable = NVARtransformed
    else
      nvariable = rafcstab%NVARtransformed
    end if

    if (iand(ioperationSpec, AFCSTAB_FCTALGO_INITALPHA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Initialise the edgewise correction factors by unity
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
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides edge-based structure
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide edge structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
          call sys_halt()
        end if

        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

        if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
          ! Perform standard prelimiting
          call doStdPrelimitDP(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
          ! Perform minmod prelimiting
          call doMinModPrelimitDP(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        else
          call output_line('Invalid type of prelimiting!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

        ! Compute sums of antidiffusive increments
        ! based on the prelimiting fluxes
        if (present(fcb_calcADIncrements)) then
          ! User-defined callback routine
          call fcb_calcADIncrements(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_DfluxPrel, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim, rcollection=rcollection)
        else
          if (present(fcb_calcFluxTransformation_sim)) then
            ! Standard routine with flux transformation
            call doADIncrementsTransformedDP(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Standard routine without flux transformation
            call doADIncrementsDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
          end if
        end if

      else

        ! Compute sums of antidiffusive increments
        ! based on the raw-antidiffusive fluxes
        if (present(fcb_calcADIncrements)) then
          ! User-defined callback routine
          call fcb_calcADIncrements(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim, rcollection=rcollection)
        else
          if (present(fcb_calcFluxTransformation_sim)) then
            ! Compute antidiffusive incrementswith flux transformation
            call doADIncrementsTransformedDP(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Compute antidiffusive increments without flux transformation
            call doADIncrementsDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
          end if
        end if

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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Compute bounds
      if (present(fcb_calcBounds)) then
        ! User-supplied callback routine
        call fcb_calcBounds(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
            rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dqp, p_Dqm,&
            fcb_calcDiffTransformation_sim, rcollection)
      elseif (present(fcb_calcDiffTransformation_sim)) then
        ! Standard routine with difference transformation
        call doBoundsTransformedDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            nvariable, p_Dx, p_Dqp, p_Dqm)
      else
        ! Standard routine without difference transformation
        call doBoundsDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dqp, p_Dqm)
      end if

      ! Set specifiers
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITNODAL) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute nodal correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides antidiffusive increments and local bounds
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADINCREMENTS) .eq. 0) .or.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODEBOUNDS)   .eq. 0)) then
        call output_line('Stabilisation does not provide increments and/or bounds!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Set additional pointers
      call lsyssc_getbase_double(rmatrix, p_ML)

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, rcollection)
      elseif (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        ! Standard routine without constraints
        call doLimitNodalDP(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        ! Standard routine with constraints
        call doLimitNodalConstrainedDP(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      end if

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides nodal correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provides nodal correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Compute edgewise correction factors
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Drp, p_Drm, fcb_calcFluxTransformation_sim, p_DfluxPrel, rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseConstrTransfDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_DfluxPrel, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseConstrainedDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_DfluxPrel, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      else

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Drp, p_Drm, fcb_calcFluxTransformation_sim, rcollection=rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseTransformedDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if
      end if

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_CORRECT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! Correct antidiffusive fluxes and apply them
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides edgewise correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provides edgewise correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTBlock')
        call sys_halt()
      end if

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then

        ! Set pointer
        call lsyssc_getbase_double(rmatrix, p_ML)

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NEQ,&
              rafcstab%NVAR, dscale, p_Dx, p_Dalpha, p_Dflux, p_Dy, p_ML, rcollection)
        else
          ! Standard routine
          call doCorrectScaleByMassDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
              p_ML, p_Dalpha, p_Dflux, p_Dy)
        end if

      else

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NEQ,&
              rafcstab%NVAR, dscale, p_Dx, p_Dalpha, p_Dflux, p_Dy,&
              rcollection=rcollection)
        else
          ! Standard routine
          call doCorrectDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
              p_Dalpha, p_Dflux, p_Dy)
        end if
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes the standard way, as
    ! suggested by Boris and Book in their first FCT algorithm

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doStdPrelimitDP(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,DfluxPrel
      integer, intent(in) :: NEDGE,NVAR

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      integer :: iedge,ivar

      ! Loop over all edges
      !$omp parallel do default(shared) private(ivar)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      edgeloop: do iedge = 1, NEDGE

        ! Check if the antidiffusive flux is directed down the gradient
        !   $F_{ij}*(u_i-u_j) < 0$
        ! and if its magnitude is larger than an absolute tolerance
        !  $ |F_{ij}| > tol$
        ! In this case, cancel the flux completely.
        do ivar = 1, NVAR
          if ((Dflux(ivar,iedge)*DfluxPrel(ivar,iedge) .lt. 0.0_DP) .and.&
              abs(Dflux(ivar,iedge)) .gt. AFCSTAB_PRELIMABS) then
            Dalpha(iedge) = 0.0_DP
            cycle edgeloop
          end if
        end do
      end do edgeloop
      !$omp end parallel do

    end subroutine doStdPrelimitDP

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes using minmod limiter

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doMinModPrelimitDP(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,DfluxPrel
      integer, intent(in) :: NEDGE,NVAR

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      integer :: iedge,ivar

      ! Loop over all edges
      !$omp parallel do default(shared) private(ivar)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      edgeloop: do iedge = 1, NEDGE

        do ivar = 1,NVAR

          ! Check if the magnitude of the antidiffusive flux is larger
          ! than an absolute tolerance; otherwise no prelimiting is done
          if (abs(Dflux(ivar,iedge)) .gt. AFCSTAB_PRELIMABS) then
            ! Check if the antidiffusive flux is directed down the gradient
            !   $F_{ij}*fp_ij < 0$
            if (Dflux(ivar,iedge)*DfluxPrel(ivar,iedge) .lt. 0.0_DP) then
              ! Then, cancel the antidiffusive flux completely
              Dalpha(iedge) = 0.0_DP
              cycle edgeloop
            elseif (abs(Dflux(ivar,iedge)) .gt. abs(DfluxPrel(ivar,iedge))) then
              ! Check if the magnitude of the raw antidiffusive flux
              ! exceeds the magnitude of the prelimiting flux
              !   $|F_{ij}| > |fp_ij|$
              ! then set the correction factor as follows
              Dalpha(iedge) = min(Dalpha(iedge),&
                                  DfluxPrel(ivar,iedge)/Dflux(ivar,iedge))
            end if
          end if
        end do
      end do edgeloop
      !$omp end parallel do

    end subroutine doMinModPrelimitDP

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and prelimiting

    subroutine doADIncrementsDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! output parameters
      real(DP), dimension(NVAR,NEQ), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j


      !$omp parallel default(shared) private(i,j,F_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
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
          F_ij = Dalpha(iedge) * Dflux(:,iedge)

          ! Compute the sums of antidiffusive increments
          Dpp(:,i) = Dpp(:,i) + max(0.0_DP, F_ij)
          Dpp(:,j) = Dpp(:,j) + max(0.0_DP,-F_ij)
          Dpm(:,i) = Dpm(:,i) + min(0.0_DP, F_ij)
          Dpm(:,j) = Dpm(:,j) + min(0.0_DP,-F_ij)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsDP

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-defined
    ! set of variables prior to computing the sums

    subroutine doADIncrementsTransformedDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dpp,Dpm

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j


      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,DtransformedFluxesAtEdge,&
      !$omp         IEDGEmax,i,idx,iedge,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
            DfluxesAtEdge(:,idx) = Dalpha(iedge)*Dflux(:,iedge)
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vectors
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute the sums of positive/negative antidiffusive increments
            Dpp(:,i) = Dpp(:,i) + max(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
            Dpp(:,j) = Dpp(:,j) + max(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
            Dpm(:,i) = Dpm(:,i) + min(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
            Dpm(:,j) = Dpm(:,j) + min(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      deallocate(DtransformedFluxesAtEdge)
      !$omp end parallel

    end subroutine doADIncrementsTransformedDP

    !**************************************************************
    ! Assemble the local bounds from the predicted solution without
    ! transformation

    subroutine doBoundsDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! output parameters
      real(DP), dimension(NVAR,NEQ), intent(out) :: Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Diff
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,Diff)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)

          ! Compute solution difference
          Diff = Dx(j,:)-Dx(i,:)

          ! Compute the distance to a local extremum
          ! of the predicted solution
          Dqp(:,i) = max(Dqp(:,i), Diff)
          Dqp(:,j) = max(Dqp(:,j),-Diff)
          Dqm(:,i) = min(Dqm(:,i), Diff)
          Dqm(:,j) = min(Dqm(:,j),-Diff)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doBoundsDP

    !**************************************************************
    ! Assemble the local bounds from the predicted solution which is
    ! transformed to a user-defined set of variables

    subroutine doBoundsTransformedDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dqp,Dqm

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: DtransformedDataAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DtransformedDataAtEdge,idx,IEDGEmax,i,j,iedge)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedDataAtEdge(NVARtransformed,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute transformed differences
          call fcb_calcDiffTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DtransformedDataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute the distance to a local extremum of the predicted solution
            Dqp(:,i) = max(Dqp(:,i), DtransformedDataAtEdge(:,idx))
            Dqp(:,j) = max(Dqp(:,j),-DtransformedDataAtEdge(:,idx))
            Dqm(:,i) = min(Dqm(:,i), DtransformedDataAtEdge(:,idx))
            Dqm(:,j) = min(Dqm(:,j),-DtransformedDataAtEdge(:,idx))
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedDataAtEdge)
      !$omp end parallel

    end subroutine doBoundsTransformedDP

    !**************************************************************
    ! Compute the nodal correction factors without constraints

    subroutine doLimitNodalDP(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      real(DP) :: daux
      integer :: ieq,ivar

      if (dscale .eq. 0.0_DP) then

        ! Clear R`s
        !$omp parallel sections
        !$omp section
        call lalg_clearVector(Drp)
        !$omp section
        call lalg_clearVector(Drm)
        !$omp end parallel sections

      else

        !$omp parallel sections default(shared) private(daux,ieq,ivar)

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpp(ivar,ieq) .gt. AFCSTAB_EPSABS/dscale) then
              Drp(ivar,ieq) = daux*Dqp(ivar,ieq)/Dpp(ivar,ieq)
            else
              Drp(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpm(ivar,ieq) .lt. -AFCSTAB_EPSABS/dscale) then
              Drm(ivar,ieq) = daux*Dqm(ivar,ieq)/Dpm(ivar,ieq)
            else
              Drm(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp end parallel sections

      end if

    end subroutine doLimitNodalDP

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrainedDP(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      real(DP) :: daux
      integer :: ieq,ivar

      if (dscale .eq. 0.0_DP) then

        ! Clear R`s
        !$omp parallel sections
        !$omp section
        call lalg_clearVector(Drp)
        !$omp section
        call lalg_clearVector(Drm)
        !$omp end parallel sections

      else

        !$omp parallel sections default(shared) private(daux,ieq,ivar)

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpp(ivar,ieq) .gt. AFCSTAB_EPSABS/dscale) then
              Drp(ivar,ieq) = min(1.0_DP, daux*Dqp(ivar,ieq)/Dpp(ivar,ieq))
            else
              Drp(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpm(ivar,ieq) .lt. -AFCSTAB_EPSABS/dscale) then
              Drm(ivar,ieq) = min(1.0_DP, daux*Dqm(ivar,ieq)/Dpm(ivar,ieq))
            else
              Drm(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp end parallel sections

      end if

    end subroutine doLimitNodalConstrainedDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewiseDP(IedgeList,&
        NEDGE, NEQ, NVAR, Dflux, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F_ij,R_ij
      integer :: iedge,i,j

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,F_ij,R_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)

        ! Get precomputed raw antidiffusive fluxes
        F_ij = Dflux(:,iedge)

        ! Compute nodal correction factors
        where (F_ij .gt. AFCSTAB_EPSABS)
          R_ij = min(Drp(:,i),Drm(:,j))
        elsewhere (F_ij .lt. -AFCSTAB_EPSABS)
          R_ij = min(Drp(:,j),Drm(:,i))
        elsewhere
          R_ij = 1.0_DP
        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(R_ij)
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    ! which are transformed to a user-defined set of variables
    ! priori to computing the correction factors

    subroutine doLimitEdgewiseTransformedDP(IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Drp, Drm, Dalpha)

      ! input  parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DtransformedFluxesAtEdge,&
      !$omp         IEDGEmax,R_ij,R_ji,i,idx,iedge,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))

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
          DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Compute nodal correction factors for fluxes into node i
          where (DtransformedFluxesAtEdge(:,1,idx) .gt. AFCSTAB_EPSABS)
            R_ij = Drp(:,i)
          elsewhere (DtransformedFluxesAtEdge(:,1,idx) .lt. -AFCSTAB_EPSABS)
            R_ij = Drm(:,i)
          elsewhere
            R_ij = 1.0_DP
          end where

          ! Compute nodal correction factors for fluxes into node j
          where (DtransformedFluxesAtEdge(:,2,idx) .gt. AFCSTAB_EPSABS)
            R_ji = Drp(:,j)
          elsewhere (DtransformedFluxesAtEdge(:,2,idx) .lt. -AFCSTAB_EPSABS)
            R_ji = Drm(:,j)
          elsewhere
            R_ji = 1.0_DP
          end where

          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
        end do
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxesAtEdge)
      !$omp end parallel

    end subroutine doLimitEdgewiseTransformedDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrainedDP(IedgeList,&
        NEDGE, NEQ, NVAR, Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F1_ij,F2_ij,R_ij
      integer :: iedge,i,j

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,F1_ij,F2_ij,R_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)

        ! Get precomputed raw antidiffusive fluxes
        F1_ij = Dflux1(:,iedge)
        F2_ij = Dflux2(:,iedge)

        ! Compute nodal correction factors
        where (F1_ij*F2_ij .le. 0.0_DP)
          R_ij = 0.0_DP
        elsewhere
          where (F1_ij .ge. 0.0_DP)
            R_ij = min(1.0_DP, F1_ij/F2_ij*min(Drp(:,i),Drm(:,j)))
          elsewhere
            R_ij = min(1.0_DP, F1_ij/F2_ij*min(Drp(:,j),Drm(:,i)))
          end where
        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(R_ij)
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseConstrainedDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes which are transformed
    ! to a user-defined set of variables priori to computing the
    ! correction factors

    subroutine doLimitEdgewiseConstrTransfDP(IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes1AtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes2AtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DtransformedFluxes1AtEdge,&
      !$omp         DtransformedFluxes2AtEdge,IEDGEmax,R_ij,R_ji,i,idx,iedge,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxes1AtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxes2AtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))

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
          DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
          DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            Dflux1(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
            DtransformedFluxes1AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            Dflux2(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
            DtransformedFluxes2AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Compute nodal correction factors
          where (DtransformedFluxes1AtEdge(:,1,idx)*&
                 DtransformedFluxes2AtEdge(:,1,idx) .le. 0.0_DP)
            R_ij = 0.0_DP
          elsewhere
            R_ij = min(1.0_DP, DtransformedFluxes1AtEdge(:,1,idx)/&
                               DtransformedFluxes2AtEdge(:,1,idx)*&
                         merge(Drp(:,i), Drm(:,i),&
                               DtransformedFluxes1AtEdge(:,1,idx) .ge. 0.0_DP))
          end where

          where (DtransformedFluxes1AtEdge(:,2,idx)*&
                 DtransformedFluxes2AtEdge(:,2,idx) .le. 0.0_DP)
            R_ji = 0.0_DP
          elsewhere
            R_ji = min(1.0_DP, DtransformedFluxes1AtEdge(:,2,idx)/&
                               DtransformedFluxes2AtEdge(:,2,idx)*&
                         merge(Drp(:,j), Drm(:,j),&
                               DtransformedFluxes1AtEdge(:,2,idx) .ge. 0.0_DP))
          end where

          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
        end do
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxes1AtEdge)
      deallocate(DtransformedFluxes2AtEdge)
      !$omp end parallel

    end subroutine doLimitEdgewiseConstrTransfDP

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrectDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, dscale, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,F_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

          ! Apply limited antidiffusive fluxes
          Dy(i,:) = Dy(i,:) + F_ij
          Dy(j,:) = Dy(j,:) - F_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectDP

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix

    subroutine doCorrectScaleByMassDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, dscale, ML, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha,ML
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,F_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

          ! Apply limited antidiffusive fluxes
          Dy(i,:) = Dy(i,:) + F_ij/ML(i)
          Dy(j,:) = Dy(j,:) - F_ij/ML(j)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectScaleByMassDP

  end subroutine afcsys_buildVectorFCTBlock

  ! ****************************************************************************

!<subroutine>

  subroutine afcsys_buildVectorFCTScalar(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry, NVARtransformed,&
      fcb_calcFluxTransformation_sim, fcb_calcDiffTransformation_sim,&
      fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
      fcb_limitEdgewise, fcb_calcCorrection, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector for nonlinear FEM-FCT
    ! schemes. Note that the vectors are required as scalar vectors
    ! which are stored in the interleave format. The idea of flux
    ! corrected transport can be traced back to the early SHASTA
    ! algorithm by Boris and Bock in the early 1970s. Zalesak
    ! suggested a fully multi-dimensional generalisation of this
    ! approach and paved the way for a large family of FCT algorithms.
    !
    ! This subroutine provides different nonlinear FEM-FCT algorithms:
    !
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
    ! combination of different AFCSTAB_FCT_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVARtransformed is taken from the stabilisation structure
    integer, intent(in), optional :: NVARtransformed

    ! OPTIONAL: callback function to compute variable transformation
    include 'intf_calcFluxTransformation_sim.inc'
    optional :: fcb_calcFluxTransformation_sim

    include 'intf_calcDiffTransformation_sim.inc'
    optional :: fcb_calcDiffTransformation_sim

    ! OPTIONAL: callback functions to overwrite the standard operations
    include 'intf_calcADIncrements.inc'
    optional :: fcb_calcADIncrements

    include 'intf_calcBounds.inc'
    optional :: fcb_calcBounds

    include 'intf_limitNodal.inc'
    optional :: fcb_limitNodal

    include 'intf_limitEdgewise.inc'
    optional :: fcb_limitEdgewise

    include 'intf_calcCorrection.inc'
    optional :: fcb_calcCorrection

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_DfluxPrel
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: nvariable

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
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
    ! The nonlinear FEM-FCT algorithm is split into the following
    ! steps which can be skipped and performed externally by the user:
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
    ! 6) Apply the limited antidifusive fluxes to the vector
    !
    !    Step 6) may be split into the following substeps
    !
    !    6.1) Compute the edgewise correction factors based on the pre-
    !         computed raw-antidiffusive fluxes.
    !
    !    6.2) Compute the raw antidiffusive fluxes for a different set of
    !         variables and limit them by the precomputed correction factors.
    !-------------------------------------------------------------------------

    ! Determine number of transformed variables (if any)
    if (present(NVARtransformed)) then
      nvariable = NVARtransformed
    else
      nvariable = rafcstab%NVARtransformed
    end if

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
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
          call sys_halt()
        end if

        ! Check if stabilisation provides edge-based structure
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide edge structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
          call sys_halt()
        end if

        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

        if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
          ! Perform standard prelimiting
          call doStdPrelimitDP(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
          ! Perform minmod prelimiting
          call doMinModPrelimitDP(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        else
          call output_line('Invalid type of prelimiting!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

        ! Compute sums of antidiffusive increments
        ! based on the prelimiting fluxes
        if (present(fcb_calcADIncrements)) then
          ! User-defined callback routine
          call fcb_calcADIncrements(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_DfluxPrel, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim, rcollection=rcollection)
        else
          if (present(fcb_calcFluxTransformation_sim)) then
            ! Standard routine with flux transformation
            call doADIncrementsTransformedDP(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Standard routine without flux transformation
            call doADIncrementsDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
          end if
        end if

      else

        ! Compute sums of antidiffusive increments
        ! based on the raw-antidiffusive fluxes
        if (present(fcb_calcADIncrements)) then
          ! User-defined callback routine
          call fcb_calcADIncrements(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Dpp, p_Dpm, fcb_calcFluxTransformation_sim, rcollection=rcollection)
        else
          if (present(fcb_calcFluxTransformation_sim)) then
            ! Compute antidiffusive incrementswith flux transformation
            call doADIncrementsTransformedDP(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Compute antidiffusive increments without flux transformation
            call doADIncrementsDP(p_IedgeListIdx, p_IedgeList,&
                rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
          end if
        end if

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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Compute bounds
      if (present(fcb_calcBounds)) then
        ! User-supplied callback routine
        call fcb_calcBounds(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
            rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dqp, p_Dqm,&
            fcb_calcDiffTransformation_sim, rcollection)
      elseif (present(fcb_calcDiffTransformation_sim)) then
        ! Standard routine with difference transformation
        call doBoundsTransformedDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEQ, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
            p_Dx, p_Dqp, p_Dqm)
      else
        ! Standard routine without difference transformation
        call doBoundsDP(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_Dqp, p_Dqm)
      end if

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
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Set additional pointers
      call lsyssc_getbase_double(rmatrix, p_ML)

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, rcollection)
      elseif (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        ! Standard routine without constraints
        call doLimitNodalDP(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        ! Standard routine with constraints
        call doLimitNodalConstrainedDP(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      end if

      ! Set specifier
      rafcstab%istabilisationSpec =&
          ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)
    end if


    if (iand(ioperationSpec, AFCSTAB_FCTALGO_LIMITEDGE) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 7) Compute edgewise correction factors
      !-------------------------------------------------------------------------

      ! Check if stabilisation provides nodal correction factors
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER) .eq. 0) then
        call output_line('Stabilisation does not provides nodal correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Compute edgewise correction factors
      if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Drp, p_Drm, fcb_calcFluxTransformation_sim, p_DfluxPrel, rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseConstrTransfDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_DfluxPrel, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseConstrainedDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_DfluxPrel, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if

      else

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Drp, p_Drm, fcb_calcFluxTransformation_sim, rcollection=rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseTransformedDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseDP(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux, p_Drp, p_Drm, p_Dalpha)
        end if
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
        call output_line('Stabilisation does not provides edgewise correction factors!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorFCTScalar')
        call sys_halt()
      end if

      ! Apply antidiffusive fluxes
      if (iand(ioperationSpec, AFCSTAB_FCTALGO_SCALEBYMASS) .ne. 0) then

        ! Set pointer
        call lsyssc_getbase_double(rmatrix, p_ML)

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NVAR,&
              rafcstab%NEQ, dscale, p_Dx, p_Dalpha, p_Dflux, p_Dy, p_ML, rcollection)
        else
          ! Standard routine
          call doCorrectScaleByMassDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
              p_ML, p_Dalpha, p_Dflux, p_Dy)
        end if

      else

        if (present(fcb_calcCorrection)) then
          ! User-supplied callback routine
          call fcb_calcCorrection(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NVAR,&
              rafcstab%NEQ, dscale, p_Dx, p_Dalpha, p_Dflux, p_Dy,&
              rcollection=rcollection)
        else
          ! Standard routine
          call doCorrectDP(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
              p_Dalpha, p_Dflux, p_Dy)
        end if
      end if
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes the standard way, as
    ! suggested by Boris and Book in their first FCT algorithm

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doStdPrelimitDP(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,DfluxPrel
      integer, intent(in) :: NEDGE,NVAR

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      integer :: iedge,ivar

      ! Loop over all edges
      !$omp parallel do default(shared) private(ivar)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      edgeloop: do iedge = 1, NEDGE

        ! Check if the antidiffusive flux is directed down the gradient
        !   $f_ij*(u_i-u_j) < 0$
        ! and if its magnitude is larger than an absolute tolerance
        !  $ |f_ij| > tol$
        ! In this case, cancel the flux completely.
        do ivar = 1, NVAR
          if ((Dflux(ivar,iedge)*DfluxPrel(ivar,iedge) .lt. 0.0_DP) .and.&
              abs(Dflux(ivar,iedge)) .gt. AFCSTAB_PRELIMABS) then
            Dalpha(iedge) = 0.0_DP
            cycle edgeloop
          end if
        end do
      end do edgeloop
      !$omp end parallel do

    end subroutine doStdPrelimitDP

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes using minmod limiter

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doMinModPrelimitDP(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux,DfluxPrel
      integer, intent(in) :: NEDGE,NVAR

      ! On input: the edge-wise correction factor from previous
      !           multiplicative correction steps
      ! On exit:  the edge-wise correction factor with prelimiting
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      integer :: iedge,ivar

      ! Loop over all edges
      !$omp parallel do default(shared) private(ivar)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      edgeloop: do iedge = 1, NEDGE

        do ivar = 1,NVAR

          ! Check if the magnitude of the antidiffusive flux is larger
          ! than an absolute tolerance; otherwise no prelimiting is done
          if (abs(Dflux(ivar,iedge)) .gt. AFCSTAB_PRELIMABS) then
            ! Check if the antidiffusive flux is directed down the gradient
            !   $f_ij*fp_ij < 0$
            if (Dflux(ivar,iedge)*DfluxPrel(ivar,iedge) .lt. 0.0_DP) then
              ! Then, cancel the antidiffusive flux completely
              Dalpha(iedge) = 0.0_DP
              cycle edgeloop
            elseif (abs(Dflux(ivar,iedge)) .gt. abs(DfluxPrel(ivar,iedge))) then
              ! Check if the magnitude of the raw antidiffusive flux
              ! exceeds the magnitude of the prelimiting flux
              !   $|f_ij| > |fp_ij|$
              ! then set the correction factor as follows
              Dalpha(iedge) = min(Dalpha(iedge),&
                                  DfluxPrel(ivar,iedge)/Dflux(ivar,iedge))
            end if
          end if
        end do
      end do edgeloop
      !$omp end parallel do

    end subroutine doMinModPrelimitDP

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and prelimiting

    subroutine doADIncrementsDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! output parameters
      real(DP), dimension(NVAR,NEQ), intent(out) :: Dpp,Dpm

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j


      !$omp parallel default(shared) private(i,j,F_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
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
          F_ij = Dalpha(iedge) * Dflux(:,iedge)

          ! Compute the sums of antidiffusive increments
          Dpp(:,i) = Dpp(:,i) + max(0.0_DP, F_ij)
          Dpp(:,j) = Dpp(:,j) + max(0.0_DP,-F_ij)
          Dpm(:,i) = Dpm(:,i) + min(0.0_DP, F_ij)
          Dpm(:,j) = Dpm(:,j) + min(0.0_DP,-F_ij)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doADIncrementsDP

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-defined
    ! set of variables prior to computing the sums

    subroutine doADIncrementsTransformedDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Dalpha, Dpp, Dpm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(:), intent(in) :: Dalpha
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dpp,Dpm

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DfluxesAtEdge
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,DtransformedFluxesAtEdge,&
      !$omp         IEDGEmax,i,idx,iedge,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
            DfluxesAtEdge(:,idx) = Dalpha(iedge)*Dflux(:,iedge)
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DfluxesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vectors
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute the sums of positive/negative antidiffusive increments
            Dpp(:,i) = Dpp(:,i) + max(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
            Dpp(:,j) = Dpp(:,j) + max(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
            Dpm(:,i) = Dpm(:,i) + min(0.0_DP, DtransformedFluxesAtEdge(:,1,idx))
            Dpm(:,j) = Dpm(:,j) + min(0.0_DP, DtransformedFluxesAtEdge(:,2,idx))
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      deallocate(DtransformedFluxesAtEdge)
      !$omp end parallel

    end subroutine doADIncrementsTransformedDP

    !**************************************************************
    ! Assemble the local bounds from the predicted solution
    ! without transformation

    subroutine doBoundsDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE, NEQ,NVAR

      ! output parameters
      real(DP), dimension(NVAR,NEQ), intent(out) :: Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Diff
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,Diff)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)

          ! Compute solution difference
          Diff = Dx(:,j)-Dx(:,i)

          ! Compute the distance to a local extremum
          ! of the predicted solution
          Dqp(:,i) = max(Dqp(:,i), Diff)
          Dqp(:,j) = max(Dqp(:,j),-Diff)
          Dqm(:,i) = min(Dqm(:,i), Diff)
          Dqm(:,j) = min(Dqm(:,j),-Diff)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doBoundsDP

    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! which is transformed to a user-defined set of variables

    subroutine doBoundsTransformedDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! output parameters
      real(DP), dimension(NVARtransformed,NEQ), intent(out) :: Dqp,Dqm

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: DtransformedDataAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DtransformedDataAtEdge,idx,IEDGEmax,i,j,iedge)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedDataAtEdge(NVARtransformed,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute transformed differences
          call fcb_calcDiffTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DtransformedDataAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute the distance to a local extremum of the predicted solution
            Dqp(:,i) = max(Dqp(:,i), DtransformedDataAtEdge(:,idx))
            Dqp(:,j) = max(Dqp(:,j),-DtransformedDataAtEdge(:,idx))
            Dqm(:,i) = min(Dqm(:,i), DtransformedDataAtEdge(:,idx))
            Dqm(:,j) = min(Dqm(:,j),-DtransformedDataAtEdge(:,idx))
          end do
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedDataAtEdge)
      !$omp end parallel

    end subroutine doBoundsTransformedDP

    !**************************************************************
    ! Compute nodal correction factors without constraints

    subroutine doLimitNodalDP(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      real(DP) :: daux
      integer :: ieq,ivar

      if (dscale .eq. 0.0_DP) then

        ! Clear R`s
        !$omp parallel sections
        !$omp section
        call lalg_clearVector(Drp)
        !$omp section
        call lalg_clearVector(Drm)
        !$omp end parallel sections

      else

        !$omp parallel sections default(shared) private(daux,ieq,ivar)

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpp(ivar,ieq) .gt. AFCSTAB_EPSABS/dscale) then
              Drp(ivar,ieq) = daux*Dqp(ivar,ieq)/Dpp(ivar,ieq)
            else
              Drp(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpm(ivar,ieq) .lt. -AFCSTAB_EPSABS/dscale) then
              Drm(ivar,ieq) = daux*Dqm(ivar,ieq)/Dpm(ivar,ieq)
            else
              Drm(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp end parallel sections

      end if

    end subroutine doLimitNodalDP

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrainedDP(NEQ, NVAR, dscale,&
        ML, Dpp, Dpm, Dqp, Dqm, Drp, Drm)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dpp,Dpm,Dqp,Dqm
      real(DP), dimension(:), intent(in) :: ML
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Drp,Drm

      ! local variables
      real(DP) :: daux
      integer :: ieq,ivar

      if (dscale .eq. 0.0_DP) then

        ! Clear R`s
        !$omp parallel sections
        !$omp section
        call lalg_clearVector(Drp)
        !$omp section
        call lalg_clearVector(Drm)
        !$omp end parallel sections

      else

        !$omp parallel sections default(shared) private(daux,ieq,ivar)

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpp(ivar,ieq) .gt. AFCSTAB_EPSABS/dscale) then
              Drp(ivar,ieq) = min(1.0_DP, daux*Dqp(ivar,ieq)/Dpp(ivar,ieq))
            else
              Drp(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp section

        !$omp parallel do default(shared) private(daux,ieq,ivar)
        do ieq = 1, NEQ
          daux = ML(ieq)/dscale
          do ivar = 1, NVAR
            if (Dpm(ivar,ieq) .lt. -AFCSTAB_EPSABS/dscale) then
              Drm(ivar,ieq) = min(1.0_DP, daux*Dqm(ivar,ieq)/Dpm(ivar,ieq))
            else
              Drm(ivar,ieq) = 1.0_DP
            end if
          end do
        end do
        !$omp end parallel do

        !$omp end parallel sections

      end if

    end subroutine doLimitNodalConstrainedDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewiseDP(IedgeList,&
        NEDGE, NEQ, NVAR, Dflux, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F_ij,R_ij
      integer :: iedge,i,j

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,F_ij,R_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)

        ! Get precomputed raw antidiffusive fluxes
        F_ij = Dflux(:,iedge)

        ! Compute nodal correction factors
        where (F_ij .gt. AFCSTAB_EPSABS)
          R_ij = min(Drp(:,i),Drm(:,j))
        elsewhere (F_ij .lt. -AFCSTAB_EPSABS)
          R_ij = min(Drp(:,j),Drm(:,i))
        elsewhere
          R_ij = 1.0_DP
        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(R_ij)
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    ! which are transformed to a user-defined set of variables
    ! priori to computing the correction factors

    subroutine doLimitEdgewiseTransformedDP(IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge


      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DtransformedFluxesAtEdge,&
      !$omp         IEDGEmax,R_ij,R_ji,i,idx,iedge,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))

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
          DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
            DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Compute nodal correction factors for fluxes into node i
          where (DtransformedFluxesAtEdge(:,1,idx) .gt. AFCSTAB_EPSABS)
            R_ij = Drp(:,i)
          elsewhere (DtransformedFluxesAtEdge(:,1,idx) .lt. -AFCSTAB_EPSABS)
            R_ij = Drm(:,i)
          elsewhere
            R_ij = 1.0_DP
          end where

          ! Compute nodal correction factors for fluxes into node j
          where (DtransformedFluxesAtEdge(:,2,idx) .gt. AFCSTAB_EPSABS)
            R_ji = Drp(:,j)
          elsewhere (DtransformedFluxesAtEdge(:,2,idx) .lt. -AFCSTAB_EPSABS)
            R_ji = Drm(:,j)
          elsewhere
            R_ji = 1.0_DP
          end where

          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij,R_ji))
        end do
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxesAtEdge)
      !$omp end parallel

    end subroutine doLimitEdgewiseTransformedDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrainedDP(IedgeList,&
        NEDGE, NEQ, NVAR, Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVAR,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! local variables
      real(DP), dimension(NVAR) :: F1_ij,F2_ij,R_ij
      integer :: iedge,i,j

      ! Loop over all edges
      !$omp parallel do default(shared) private(i,j,F1_ij,F2_ij,R_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Get node numbers
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)

        ! Get precomputed raw antidiffusive fluxes
        F1_ij = Dflux1(:,iedge)
        F2_ij = Dflux2(:,iedge)

        ! Compute nodal correction factors
        where (F1_ij*F2_ij .le. 0.0_DP)
          R_ij = 0.0_DP
        elsewhere
          where (F1_ij .ge. 0.0_DP)
            R_ij = min(1.0_DP, F1_ij/F2_ij*min(Drp(:,i),Drm(:,j)))
          elsewhere
            R_ij = min(1.0_DP, F1_ij/F2_ij*min(Drp(:,j),Drm(:,i)))
          end where
        end where

        ! Compute multiplicative correction factor
        Dalpha(iedge) = Dalpha(iedge) * minval(R_ij)
      end do
      !$omp end parallel do

    end subroutine doLimitEdgewiseConstrainedDP

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes which are transformed
    ! to a user-defined set of variables priori to computing the
    ! correction factors

    subroutine doLimitEdgewiseConstrTransfDP(IedgeList,&
        NEDGE, NEQ, NVAR, NVARtransformed, Dx, Dflux1, Dflux2, Drp, Drm, Dalpha)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1,Dflux2
      real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR,NVARtransformed

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dalpha

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes1AtEdge
      real(DP), dimension(:,:,:), pointer :: DtransformedFluxes2AtEdge

      ! local variables
      real(DP), dimension(NVARtransformed) :: R_ij,R_ji
      integer :: idx,IEDGEset,IEDGEmax,i,j,iedge

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DtransformedFluxes1AtEdge,&
      !$omp         DtransformedFluxes2AtEdge,IEDGEmax,R_ij,R_ji,i,idx,iedge,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxes1AtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))
      allocate(DtransformedFluxes2AtEdge(NVARtransformed,2,p_rperfconfig%NEDGESIM))

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
          DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
          DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
        end do

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            Dflux1(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
            DtransformedFluxes1AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Use callback function to compute transformed fluxes
        call fcb_calcFluxTransformation_sim(&
            DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            Dflux2(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
            DtransformedFluxes2AtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
            rcollection)

        ! Loop through all edges in the current set
        ! and scatter the entries to the global vector
        do idx = 1, IEDGEmax-IEDGEset+1

          ! Get actual edge number
          iedge = idx+IEDGEset-1

          ! Get position of nodes
          i = IedgeList(1,iedge)
          j = IedgeList(2,iedge)

          ! Compute nodal correction factors
          where (DtransformedFluxes1AtEdge(:,1,idx)*&
                 DtransformedFluxes2AtEdge(:,1,idx) .le. 0.0_DP)
            R_ij = 0.0_DP
          elsewhere
            R_ij = min(1.0_DP, DtransformedFluxes1AtEdge(:,1,idx)/&
                               DtransformedFluxes2AtEdge(:,1,idx)*&
                         merge(Drp(:,i), Drm(:,i),&
                               DtransformedFluxes1AtEdge(:,1,idx) .ge. 0.0_DP))
          end where

          where (DtransformedFluxes1AtEdge(:,2,idx)*&
                 DtransformedFluxes2AtEdge(:,2,idx) .le. 0.0_DP)
            R_ji = 0.0_DP
          elsewhere
            R_ji = min(1.0_DP, DtransformedFluxes1AtEdge(:,2,idx)/&
                               DtransformedFluxes2AtEdge(:,2,idx)*&
                         merge(Drp(:,j), Drm(:,j),&
                               DtransformedFluxes1AtEdge(:,2,idx) .ge. 0.0_DP))
          end where

          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) * minval(min(R_ij, R_ji))
        end do
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DtransformedFluxes1AtEdge)
      deallocate(DtransformedFluxes2AtEdge)
      !$omp end parallel

    end subroutine doLimitEdgewiseConstrTransfDP

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrectDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, dscale, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,F_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

          ! Apply limited antidiffusive fluxes
          Dy(:,i) = Dy(:,i) + F_ij
          Dy(:,j) = Dy(:,j) - F_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectDP

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix

    subroutine doCorrectScaleByMassDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, dscale, ML, Dalpha, Dflux, Dy)

      ! input parameters
      real(DP), dimension(:), intent(in) :: Dalpha,ML
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j

      !$omp parallel default(shared) private(i,j,F_ij)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)

          ! Apply limited antidiffusive fluxes
          Dy(:,i) = Dy(:,i) + F_ij/ML(i)
          Dy(:,j) = Dy(:,j) - F_ij/ML(j)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel
    end subroutine doCorrectScaleByMassDP

  end subroutine afcsys_buildVectorFCTScalar

  ! ****************************************************************************

!<subroutine>

  subroutine afcsys_buildVectorTVDBlock(rafcstab, rgroupFEMSet, rx, ndim,&
      fcb_calcFluxSys_sim, fcb_calcCharacteristics_sim, dscale, bclear, ry,&
      rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector for FEM-TVD schemes. If
    ! the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
    ! group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! number of spatial dimensions
    integer, intent(in) :: ndim

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback function to compute local fluxes
    include 'intf_calcFluxSys_sim.inc'

    ! callback function to compute local characteristics
    include 'intf_calcCharacteristics_sim.inc'

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
    !</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call afcsys_buildVectorTVDScalar(rafcstab, rgroupFEMSet, rx%RvectorBlock(1),&
          ndim, fcb_calcFluxSys_sim, fcb_calcCharacteristics_sim, dscale, bclear,&
          ry%RvectorBlock(1), rcollection, rperfconfig)
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST).eq. 0) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDBlock')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)

    ! Set pointers
    call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsysbl_getbase_double(rx, p_Dx)
    call lsysbl_getbase_double(ry, p_Dy)

    ! How many dimensions do we have?
    select case(ndim)
    case (NDIM1D)
      call doLimitTVD_1D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM2D)
      call doLimitTVD_2D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM3D)
      call doLimitTVD_3D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D

    subroutine doLimitTVD_1D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
            Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

        ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_1D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D

    subroutine doLimitTVD_2D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
            Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_2D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D

    subroutine doLimitTVD_3D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(i,:) = Dy(i,:)+DfluxesAtEdge(:,1,idx)
            Dy(j,:) = Dy(j,:)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(3, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Z-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(3, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_3D


    !**************************************************************
    ! Assemble the upper and lower bounds Q and the sums of
    ! antidiffusive contributions P for a given set of edges

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doBoundsAndIncrements_sim(idirection, NVAR,&
        DcoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IedgeList, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dpp,Dpm,Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      integer :: idx,ivar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,i,ivar,j)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)

        ! Loop over all characteristic variables
        do ivar = 1, NVAR
          ! Set node orientation
          if (Daux1(ivar) .gt. 0) then
            i = IedgeList(2,idx)
            j = IedgeList(1,idx)
            Dflux(ivar) = -Dflux(ivar)
          else
            i = IedgeList(1,idx)
            j = IedgeList(2,idx)
          end if

          ! Assemble P's and Q's
          if (Dflux(ivar) .gt. 0) then
            Dpp(ivar,i) = Dpp(ivar,i)+Dflux(ivar)
            Dqm(ivar,i) = Dqm(ivar,i)-Dflux(ivar)
            Dqp(ivar,j) = Dqp(ivar,j)+Dflux(ivar)
          else
            Dpm(ivar,i) = Dpm(ivar,i)+Dflux(ivar)
            Dqp(ivar,i) = Dqp(ivar,i)-Dflux(ivar)
            Dqm(ivar,j) = Dqm(ivar,j)+Dflux(ivar)
          end if
        end do

      end do
      !$omp end parallel do

    end subroutine doBoundsAndIncrements_sim


    !**************************************************************
    ! Limit the antidiffusive fluxes and apply them to the vector

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doLimitADFluxes_sim(idirection, NVAR, dscale, Drp, Drm,&
        DcoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IedgeList, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      real(DP), dimension(:,:), intent(in) :: DrighteigenvectorsAtEdge,Drp,Drm
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      real(DP) :: daux
      integer :: idx,ivar,jvar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,daux,i,ivar,j,jvar)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)
        Daux2 = abs(Daux1)*DcharVariablesAtEdge(:,idx)

        ! Get position of nodes
        i = IedgeList(1,idx)
        j = IedgeList(2,idx)

        ! Loop over all characteristic variables
        ! and limit characteristic fluxes
        do ivar = 1, NVAR

          if (Daux1(ivar) .lt. 0) then
            if (Dflux(ivar) .gt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,i)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,i)*Dflux(ivar)
            end if
          else
            if (Dflux(ivar) .lt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,j)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,j)*Dflux(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        do ivar = 1, NVAR
          daux = 0.0_DP
          do jvar = 1, NVAR
            daux = daux+DrighteigenvectorsAtEdge(NVAR*(jvar-1)+ivar,idx)*Daux2(jvar)
          end do
          Dflux(ivar) = dscale*daux
        end do

        ! Apply limited fluxes to global vector
        Dy(i,:) = Dy(i,:)+Dflux
        Dy(j,:) = Dy(j,:)-Dflux

      end do

      !$omp end parallel do

    end subroutine doLimitADFluxes_sim

  end subroutine afcsys_buildVectorTVDBlock

  ! ****************************************************************************

!<subroutine>

  subroutine afcsys_buildVectorTVDScalar(rafcstab, rgroupFEMSet, rx, ndim,&
      fcb_calcFluxSys_sim, fcb_calcCharacteristics_sim, dscale, bclear, ry,&
      rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the vector for FEM-TVD schemes.
!</description>

!<input>
    ! group finite element set
    type(t_groupFEMSet), intent(in) :: rgroupFEMSet

    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! number of spatial dimensions
    integer, intent(in) :: ndim

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! callback function to compute local fluxes
    include 'intf_calcFluxSys_sim.inc'

    ! callback function to compute local characteristics
    include 'intf_calcCharacteristics_sim.inc'

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! destination vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:,:,:), pointer :: p_DcoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST) .eq. 0) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildVectorTVDScalar')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call gfem_getbase_DcoeffsAtEdge(rgroupFEMset, p_DcoeffsAtEdge)
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call lsyssc_getbase_double(rafcstab%p_rvectorPp, p_Dpp)
    call lsyssc_getbase_double(rafcstab%p_rvectorPm, p_Dpm)
    call lsyssc_getbase_double(rafcstab%p_rvectorQp, p_Dqp)
    call lsyssc_getbase_double(rafcstab%p_rvectorQm, p_Dqm)
    call lsyssc_getbase_double(rafcstab%p_rvectorRp, p_Drp)
    call lsyssc_getbase_double(rafcstab%p_rvectorRm, p_Drm)
    call lsyssc_getbase_double(rx, p_Dx)
    call lsyssc_getbase_double(ry, p_Dy)

    ! How many dimensions do we have?
    select case(ndim)
    case (NDIM1D)
      call doLimitTVD_1D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM2D)
      call doLimitTVD_2D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM3D)
      call doLimitTVD_3D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DcoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D

    subroutine doLimitTVD_1D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameter
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
            Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir1D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_1D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D

    subroutine doLimitTVD_2D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
            Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir2D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_2D


    !**************************************************************
    ! Assemble vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D

    subroutine doLimitTVD_3D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DcoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dpp,Dpm,Dqp,Dqm,Drp,Drm
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge,DfluxesAtEdge
      real(DP), dimension(:,:), pointer :: DcharVariablesAtEdge
      real(DP), dimension(:,:), pointer :: DeigenvaluesAtEdge
      real(DP), dimension(:,:), pointer :: DrighteigenvectorsAtEdge

      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DcharVariablesAtEdge,DdataAtEdge,DeigenvaluesAtEdge,&
      !$omp         DfluxesAtEdge,DrighteigenvectorsAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,p_rperfconfig%NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,p_rperfconfig%NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Assemble high-order Galerkin fluxes
          !-----------------------------------------------------------------------

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Update the global vector
            Dy(:,i) = Dy(:,i)+DfluxesAtEdge(:,1,idx)
            Dy(:,j) = Dy(:,j)+DfluxesAtEdge(:,2,idx)
          end do

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(1, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate some temporal memory
      deallocate(DfluxesAtEdge)

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (X-direction)
      !-------------------------------------------------------------------------

      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,p_rperfconfig%NEDGESIM))

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (X-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the X-direction
          call fcb_calcCharacteristics_sim(XDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(1, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(2, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Y-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Y-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Y-direction
          call fcb_calcCharacteristics_sim(YDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(2, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)

          !-----------------------------------------------------------------------
          ! Assemble artificial viscosities and antidiffusive fluxes (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Assemble the upper and lower bounds Q and the sums of
          ! antidiffusive contributions P for the set of edges
          call doBoundsAndIncrements_sim(3, NVAR,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              Dpp, Dpm, Dqp, Dqm)
        end do
        !$omp end do

      end do ! igroup

      !-------------------------------------------------------------------------
      ! Compute nodal correction factors (Z-direction)
      !-------------------------------------------------------------------------

      !$omp single
      Drp = afcstab_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = afcstab_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, p_rperfconfig%NEDGESIM

          ! We always handle  edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most  edges simultaneously.

          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+p_rperfconfig%NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          !-----------------------------------------------------------------------
          ! Apply artificial viscosities and limited antidiffusion (Z-direction)
          !-----------------------------------------------------------------------

          ! Use callback function to compute the characteristic variables
          ! and corresponding eigenvalues along the Z-direction
          call fcb_calcCharacteristics_sim(ZDir3D,&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              IEDGEmax-IEDGEset+1,&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              rcollection=rcollection)

          ! Apply limited characteristic fluxes to global vector
          call doLimitADFluxes_sim(3, NVAR, dscale, Drp, Drm,&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              DcharVariablesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DeigenvaluesAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              DrightEigenvectorsAtEdge(:,1:IEDGEmax-IEDGEset+1),&
              IedgeList(:,IEDGEset:IEDGEmax), Dy)
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcharVariablesAtEdge)
      deallocate(DeigenvaluesAtEdge)
      deallocate(DrighteigenvectorsAtEdge)
      !$omp end parallel

    end subroutine doLimitTVD_3D

    !**************************************************************
    ! Assemble the upper and lower bounds Q and the sums of
    ! antidiffusive contributions P for a given set of edges

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doBoundsAndIncrements_sim(idirection, NVAR,&
        DcoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IedgeList, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dpp,Dpm,Dqp,Dqm

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      integer :: idx,ivar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,i,ivar,j)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)

        ! Loop over all characteristic variables
        do ivar = 1, NVAR
          ! Set node orientation
          if (Daux1(ivar) .gt. 0) then
            i = IedgeList(2,idx)
            j = IedgeList(1,idx)
            Dflux(ivar) = -Dflux(ivar)
          else
            i = IedgeList(1,idx)
            j = IedgeList(2,idx)
          end if

          ! Assemble P's and Q's
          if (Dflux(ivar) .gt. 0) then
            Dpp(ivar,i) = Dpp(ivar,i)+Dflux(ivar)
            Dqm(ivar,i) = Dqm(ivar,i)-Dflux(ivar)
            Dqp(ivar,j) = Dqp(ivar,j)+Dflux(ivar)
          else
            Dpm(ivar,i) = Dpm(ivar,i)+Dflux(ivar)
            Dqp(ivar,i) = Dqp(ivar,i)-Dflux(ivar)
            Dqm(ivar,j) = Dqm(ivar,j)+Dflux(ivar)
          end if
        end do

      end do
      !$omp end parallel do

    end subroutine doBoundsAndIncrements_sim

    !**************************************************************
    ! Limit the antidiffusive fluxes and apply them to the vector

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doLimitADFluxes_sim(idirection, NVAR, dscale, Drp, Drm,&
        DcoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IedgeList, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), dimension(:,:), intent(in) :: DcharVariablesAtEdge,DeigenvaluesAtEdge
      real(DP), dimension(:,:), intent(in) :: DrighteigenvectorsAtEdge,Drp,Drm
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: idirection, NVAR

      ! input/output parameters
      real(DP), dimension(:,:), intent(inout) :: Dy

      ! local variables
      real(DP), dimension(NVAR) :: Daux1,Daux2,Dflux
      real(DP) :: daux
      integer :: idx,ivar,jvar,i,j

      ! Loop over all edges in the set
      !$omp parallel do default(shared)&
      !$omp private(Daux1,Daux2,Dflux,daux,i,ivar,j,jvar)
      do idx = 1, size(DcharVariablesAtEdge,2)

        ! Compute unidirectional antidiffusive fluxes
        Daux1 = (DcoeffsAtEdge(idirection,2,idx)-&
                 DcoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DcoeffsAtEdge(idirection,1,idx)+&
                 DcoeffsAtEdge(idirection,2,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Dflux = -max(0.0_DP, min(abs(Daux1)-Daux2,&
                                 2.0_DP*abs(Daux1)))*DcharVariablesAtEdge(:,idx)
        Daux2 = abs(Daux1)*DcharVariablesAtEdge(:,idx)

        ! Get position of nodes
        i = IedgeList(1,idx)
        j = IedgeList(2,idx)

        ! Loop over all characteristic variables
        ! and limit characteristic fluxes
        do ivar = 1, NVAR

          if (Daux1(ivar) .lt. 0) then
            if (Dflux(ivar) .gt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,i)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,i)*Dflux(ivar)
            end if
          else
            if (Dflux(ivar) .lt. 0) then
              Daux2(ivar) = Daux2(ivar)+Drp(ivar,j)*Dflux(ivar)
            else
              Daux2(ivar) = Daux2(ivar)+Drm(ivar,j)*Dflux(ivar)
            end if
          end if
        end do

        ! Transform back into conservative variables
        do ivar = 1, NVAR
          daux = 0.0_DP
          do jvar = 1, NVAR
            daux = daux+DrighteigenvectorsAtEdge(NVAR*(jvar-1)+ivar,idx)*Daux2(jvar)
          end do
          Dflux(ivar) = dscale*daux
        end do

        ! Apply limited fluxes to global vector
        Dy(:,i) = Dy(:,i)+Dflux
        Dy(:,j) = Dy(:,j)-Dflux

      end do
      !$omp end parallel do

    end subroutine doLimitADFluxes_sim

  end subroutine afcsys_buildVectorTVDScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_buildFluxFCTBlock(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, fcb_calcFluxFCTSys_sim,&
      rgroupFEMSet, rmatrix, rxTimeDeriv, rxPredictor, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! algebraic flux correction of FCT-type with or without the
    ! contribution of the consistent mass matrix. If the vectors
    ! contain only one block, then the scalar counterpart of this
    ! routine is called with the scalar subvectors.
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling parameter
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
    include 'intf_calcFluxFCTSys_sim.inc'
    optional :: fcb_calcFluxFCTSys_sim

    ! OPTIONAL: group finite element set
    type(t_groupFEMSet), intent(in), optional :: rgroupFEMSet

    ! OPTIONAL: mass matrix
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
    ! Stabilisation structure
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
    integer :: nblocks
    logical :: buseCallback

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if block vector(s) contains exactly one block
    nblocks = rx%nblocks
    if (present(rxTimeDeriv)) nblocks = max(nblocks, rxTimeDeriv%nblocks)
    if (present(rxPredictor)) nblocks = max(nblocks, rxPredictor%nblocks)

    if (nblocks .eq. 1) then
      ! Call subroutine for scalar vectors
      if (present(rxTimeDeriv)) then
        if (present(rxPredictor)) then
          ! ... both approximate time derivative and predictor are present
          call afcsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              fcb_calcFluxFCTSys_sim, rgroupFEMSet, rmatrix,&
              rxTimeDeriv%RvectorBlock(1), rxPredictor%RvectorBlock(1),&
              rcollection, rperfconfig)
        else
          ! ... only the approximate time derivative is present
          call afcsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              fcb_calcFluxFCTSys_sim, rgroupFEMSet, rmatrix,&
              rxTimeDeriv%RvectorBlock(1),&
              rcollection=rcollection, rperfconfig=rperfconfig)
        end if
      else
        if (present(rxPredictor)) then
          ! ... only the predictor is present
          call afcsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              theta, tstep, dscale, bclear,  bquickAssembly, ioperationSpec,&
              fcb_calcFluxFCTSys_sim, rgroupFEMSet, rmatrix,&
              rxPredictor=rxPredictor%RvectorBlock(1),&
              rcollection=rcollection, rperfconfig=rperfconfig)
        else
          ! ... neither the approximate time derivative nor the predictor is present
          call afcsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              fcb_calcFluxFCTSys_sim, rgroupFEMSet, rmatrix,&
              rcollection=rcollection, rperfconfig=rperfconfig)
        end if
      end if

      ! That`s it
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    !---------------------------------------------------------------------------

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix (if present)
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
        call sys_halt()
      end if
    end if

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call lsysbl_getbase_double(rx, p_Dx)

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_NLINFCT_IMPLICIT)

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)

      ! Use callback routine?
      if (present(fcb_calcFluxFCTSys_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
          call sys_halt()
        end if
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        buseCallback = .true.
      else

        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
          call sys_halt()
        end if

        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        buseCallback = .false.
      end if

      !-------------------------------------------------------------------------
      ! Classical, iterative and semi-implicit nonlinear FEM-FCT algorithm
      ! The raw antidiffusive flux for all algorithms can be assembled
      ! in the same way. The only difference is that the amount of rejected
      ! antidiffusion is subtracted from the initial fluxes in subsequent
      ! iterations if the iterative FEM-FCT algorithm is applied.
      ! Moreover the initial flux is without mass contribution is stored
      ! separately for the semi-implicit FEM-FCT algorithm since it is
      ! used to constrain the raw antidiffusive fluxes in each iteration.
      !-------------------------------------------------------------------------

      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_EXPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble explicit part of raw-antidiffive fluxes
        !-----------------------------------------------------------------------

        if (theta .ne. 1.0_DP) then
          ! Assemble the explicit part of the raw-antidiffusive fluxes
          ! $$ F_{ij}^n = (1-\theta)\Delta t D_{ij}^n(U_i^n-U_j^n) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale*(1.0_DP-theta),&
                bclear, p_Dflux0)
          else
            call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale*(1.0_DP-theta),&
                bclear, p_Dflux0)
          end if
        elseif (.not.bquickAssembly .and. bclear) then
          ! Clear the explicit part of the raw-antidiffusive fluxes
          ! $$ F_{ij}^n = 0 $$
          call lalg_clearVector(p_Dflux0, rafcstab%NEDGE)
          ! if bquickAssembly = TRUE then this step can be skipped
        end if

        !-----------------------------------------------------------------------

        ! Check for special treatment
        if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

          ! Set pointers
          call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

          ! We have to store the raw-antidiffusive fluxes based on the
          ! initial solution without contribution of the consistent
          ! mass matrix and without scaling by the implicitness parameter
          ! $$ F_{ij} = \Delta t D_{ij}^n(U_i^n-U_j^n) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale, .true., p_DfluxPrel)
          else
            call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale, .true., p_DfluxPrel)
          end if

        elseif (rafcstab%cprelimitingType .ne. AFCSTAB_PRELIMITING_NONE) then

          ! We have to assemble the raw-antidiffusive fluxes for
          ! prelimiting separately based on the low-order predictor
          if (present(rxPredictor)) then

            ! Set pointers
            call lsysbl_getbase_double(rxPredictor, p_DxPredictor)
            call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

            if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
              ! Compute solution difference for standard prelimiting
              call doDifferencesDP(p_IedgeList, rafcstab%NEDGE,&
                  rafcstab%NEQ, rafcstab%NVAR, p_DxPredictor, p_DfluxPrel)
            elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
              ! Compute fluxes for minmod prelimiting
              if (buseCallback) then
                call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                    rafcstab%NVAR, p_DcoeffsAtEdge, p_DxPredictor, dscale,&
                    .true., p_DfluxPrel)
              else
                call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                    rafcstab%NVAR, p_Dcoefficients, p_DxPredictor, dscale,&
                    .true., p_DfluxPrel)
              end if
            else
              call output_line('Invalid type of prelimiting!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
              call sys_halt()
            end if
          else
            call output_line('Fluxes for prelimiting cannot be assembled without predictor!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
            call sys_halt()
          end if
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then

          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the explicit part of the mass-antidiffusive fluxes
          ! $$ F_{ij}^n := F_{ij}^n - M_{ij}(U_i^n-U_j^n) $$
          call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, -dscale/tstep,&
              .false., p_Dflux0)
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
                OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
            call sys_halt()
          end if

          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
            call sys_halt()
          end if

          ! Set pointer
          call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)

          ! Subtract amount of rejected antidiffusion
          call afcstab_combineFluxesDP(rafcstab%NVAR, rafcstab%NEDGE,&
              -1.0_DP, p_Dflux, p_Dflux0, p_Dalpha)
        end if

        !-----------------------------------------------------------------------

        if (theta .ne. 0.0_DP) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! $$ F_{ij} = \theta\Delta t D_{ij}(U_i-U_j) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale*theta,&
                bclear, p_Dflux)
          else
            call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale*theta,&
                bclear, p_Dflux)
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
              call afcstab_combineFluxesDP(rafcstab%NVAR, rafcstab%NEDGE,&
                  1.0_DP, p_Dflux0, p_Dflux)
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
          call afcstab_combineFluxesDP(rafcstab%NVAR, rafcstab%NEDGE,&
              1.0_DP, p_Dflux0, p_Dflux)
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then

          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the implicit part of the mass-antidiffusive fluxes
          ! $$ F_{ij}^m := F_{ij}^m + M_{ij}(U_i^m-U_j^m) $$
          call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, dscale/tstep,&
              .false., p_Dflux)

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
      if (present(fcb_calcFluxFCTSys_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
          call sys_halt()
        end if
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        buseCallback = .true.
      else

        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
          call sys_halt()
        end if

        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        buseCallback = .false.
      end if

      ! Assemble spatial part of raw-antidiffusive fluxes
      if (buseCallback) then
        call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Dflux)
      else
        call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale, bclear, p_Dflux)
      end if

      !-------------------------------------------------------------------------

      if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
        ! Compute fluxes for standard prelimiting based on the
        ! low-order solution which serves as predictor
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        call doDifferencesDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_DfluxPrel)

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
        call lsysbl_getbase_double(rxTimeDeriv, p_DxTimeDeriv)

        ! Apply mass antidiffusion to antidiffusive fluxes based on
        ! the approximation to the time derivative
        call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_DxTimeDeriv,&
            dscale, .false., p_Dflux)
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
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)

        ! Assemble mass-antidiffusive fluxes based on the solution
        call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, dscale, .true., p_Dflux)

        ! Set specifiers for raw antidiffusive fluxes
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)

      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
        call sys_halt()
      end if


    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTBlock')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the edge-by-edge array DcoefficientsAtEdge

    subroutine doFluxesByCoeffsDP(IedgeList, NEDGE, NEQ, NVAR,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux

      ! local variables
      integer :: iedge,i,j

      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute raw antidiffusive flux
            Dflux(:,iedge) = DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(i,:)-Dx(j,:))
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
            Dflux(:,iedge) = Dflux(:,iedge)&
                + DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(i,:)-Dx(j,:))
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
            Dflux(:,iedge) = DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(j,:)-Dx(i,:))
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
            Dflux(:,iedge) = Dflux(:,iedge)&
                + DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(j,:)-Dx(i,:))
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
            Dflux(:,iedge) = dscale * DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(i,:)-Dx(j,:))
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
            Dflux(:,iedge) = Dflux(:,iedge)&
                + dscale * DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(i,:)-Dx(j,:))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesByCoeffsDP

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without
    ! contribution of the consistent mass matrix.

    subroutine doFluxesByCallbackDP(IedgeList, NEDGE, NEQ, NVAR,&
        DcoeffsAtEdge, Dx, dscale, bclear, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      logical, intent(in) :: bclear

      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: DfluxAtEdge

      ! local variables
      integer :: idx,iedge,IEDGEset,IEDGEmax


      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux)

      elseif (bclear) then

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,idx,iedge,IEDGEmax)

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))

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
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dflux(:,IEDGEset:IEDGEmax), rcollection)
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        !$omp end parallel

      else   ! bclear = .false.

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,DfluxAtEdge,idx,iedge,IEDGEmax)

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
        allocate(DfluxAtEdge(NVAR,p_rperfconfig%NEDGESIM))

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
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxAtEdge(:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Add antidiffusive fluxes
            Dflux(:,iedge) = Dflux(:,iedge) + DfluxAtEdge(:,idx)
          end do
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge, DfluxAtEdge)
        !$omp end parallel

      end if

    end subroutine doFluxesByCallbackDP

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix Dmatrix

    subroutine doFluxesByMatrixDP(IedgeList, NEDGE, NEQ, NVAR,&
        Dmatrix, Dx, dscale, bclear, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: Dmatrix
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      logical, intent(in) :: bclear

      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux

      ! local variables
      integer :: iedge,ij,i,j

      if (dscale .eq. 0.0_DP) then

        ! Do we have to clear the vector?
        if (bclear) call lalg_clearVector(Dflux)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dmatrix(ij) * (Dx(i,:)-Dx(j,:))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dflux(:,iedge)&
                           + Dmatrix(ij) * (Dx(i,:)-Dx(j,:))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dmatrix(ij) * (Dx(j,:)-Dx(i,:))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dflux(:,iedge)&
                           + Dmatrix(ij) * (Dx(j,:)-Dx(i,:))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = dscale * Dmatrix(ij) * (Dx(i,:)-Dx(j,:))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dflux(:,iedge)&
                           + dscale * Dmatrix(ij) * (Dx(i,:)-Dx(j,:))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesByMatrixDP

    !**************************************************************
    ! Assemble fluxes for classical prelimiting.

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doDifferencesDP(IedgeList, NEDGE, NEQ, NVAR, Dx, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      integer :: iedge,i,j

      !$omp parallel do default(shared) private(i,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Determine indices
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)

        ! Compute solution difference; in contrast to the literature,
        ! we compute the solution difference $u_i-u_j$ and check if
        ! $F_{ij}(U_i-U_j)<0$ in the prelimiting step.
        Dflux(:,iedge) = Dx(i,:)-Dx(j,:)
      end do
      !$omp end parallel do

    end subroutine doDifferencesDP

  end subroutine afcsys_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_buildFluxFCTScalar(rafcstab, rx, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, fcb_calcFluxFCTSys_sim,&
      rgroupFEMSet, rmatrix, rxTimeDeriv, rxPredictor, rcollection, rperfconfig)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! FEM-FCT schemes. Note that the vectors are required as scalar
    ! vectors which are stored in the interleave format.
!</description>

!<input>
    ! solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling parameter
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
    include 'intf_calcFluxFCTSys_sim.inc'
    optional :: fcb_calcFluxFCTSys_sim

    ! OPTIONAL: group finite element set
    type(t_groupFEMSet), intent(in), optional :: rgroupFEMSet

    ! OPTIONAL: mass matrix
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
    ! Stabilisation structure
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
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix (if present)
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
        call sys_halt()
      end if
    end if

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call lsyssc_getbase_double(rx, p_Dx)

    ! What kind of stabilisation are we?
    select case(rafcstab%cafcstabType)

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_NLINFCT_IMPLICIT)

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)

      ! Use callback routine?
      if (present(fcb_calcFluxFCTSys_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
          call sys_halt()
        end if
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        buseCallback = .true.
      else

        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
          call sys_halt()
        end if

        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        buseCallback = .false.
      end if

      !-------------------------------------------------------------------------
      ! Classical, iterative and semi-implicit nonlinear FEM-FCT algorithm
      ! The raw antidiffusive flux for all algorithms can be assembled
      ! in the same way. The only difference is that the amount of rejected
      ! antidiffusion is subtracted from the initial fluxes in subsequent
      ! iterations if the iterative FEM-FCT algorithm is applied.
      ! Moreover the initial flux is without mass contribution is stored
      ! separately for the semi-implicit FEM-FCT algorithm since it is
      ! used to constrain the raw antidiffusive fluxes in each iteration.
      !-------------------------------------------------------------------------

      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_EXPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble explicit part of raw-antidiffive fluxes
        !-----------------------------------------------------------------------

        if (theta .ne. 1.0_DP) then
          ! Assemble the explicit part of the raw-antidiffusive fluxes
          ! $$ F_{ij}^n = (1-\theta)\Delta t D_{ij}^n(U_i^n-U_j^n) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale*(1.0_DP-theta),&
                bclear, p_Dflux0)
          else
            call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale*(1.0_DP-theta),&
                bclear, p_Dflux0)
          end if
        elseif (.not.bquickAssembly .and. bclear) then
          ! Clear the explicit part of the raw-antidiffusive fluxes
          ! $$ F_{ij}^n = 0 $$
          call lalg_clearVector(p_Dflux0, rafcstab%NEDGE)
          ! if bquickAssembly = TRUE then this step can be skipped
        end if

        !-----------------------------------------------------------------------

        ! Check for special treatment
        if (rafcstab%cafcstabType .eq. AFCSTAB_NLINFCT_IMPLICIT) then

          ! Set pointers
          call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

          ! We have to store the raw-antidiffusive fluxes based on the
          ! initial solution without contribution of the consistent
          ! mass matrix and without scaling by the implicitness parameter
          ! $$ F_{ij} = \Delta t D_{ij}^n(U_i^n-U_j^n) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale, .true., p_DfluxPrel)
          else
            call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale, .true., p_DfluxPrel)
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
              call doDifferencesDP(p_IedgeList, rafcstab%NEDGE,&
                  rafcstab%NEQ, rafcstab%NVAR, p_DxPredictor, p_DfluxPrel)
            elseif (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_MINMOD) then
              ! Compute fluxes for minmod prelimiting
              if (buseCallback) then
                call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                    rafcstab%NVAR, p_DcoeffsAtEdge, p_DxPredictor, dscale,&
                    .true., p_DfluxPrel)
              else
                call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                    rafcstab%NVAR, p_Dcoefficients, p_DxPredictor, dscale,&
                    .true., p_DfluxPrel)
              end if
            else
              call output_line('Invalid type of prelimiting!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
              call sys_halt()
            end if
          else
            call output_line('Fluxes for prelimiting cannot be assembled without predictor!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
            call sys_halt()
          end if
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then

          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the explicit part of the mass-antidiffusive fluxes
          ! $$ F_{ij}^n := F_{ij}^n - M_{ij}(U_i^n-U_j^n) $$
          call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, -dscale/tstep,&
              .false., p_Dflux0)
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
                OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
            call sys_halt()
          end if

          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
            call sys_halt()
          end if

          ! Set pointer
          call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)

          ! Subtract amount of rejected antidiffusion
          call afcstab_combineFluxesDP(rafcstab%NVAR, rafcstab%NEDGE,&
              -1.0_DP, p_Dflux, p_Dflux0, p_Dalpha)
        end if

        !-----------------------------------------------------------------------

        if (theta .ne. 0.0_DP) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! $$ F_{ij} = \theta\Delta t D_{ij}(U_i-U_j) $$
          if (buseCallback) then
            call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale*theta,&
                bclear, p_Dflux)
          else
            call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale*theta,&
                bclear, p_Dflux)
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
              call afcstab_combineFluxesDP(rafcstab%NVAR, rafcstab%NEDGE,&
                  1.0_DP, p_Dflux0, p_Dflux)
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
          call afcstab_combineFluxesDP(rafcstab%NVAR, rafcstab%NEDGE,&
              1.0_DP, p_Dflux0, p_Dflux)
        end if

        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then

          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the implicit part of the mass-antidiffusive fluxes
          ! $$ F_{ij}^m := F_{ij}^m + M_{ij}(U_i^m-U_j^m) $$
          call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
              rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, dscale/tstep,&
              .false., p_Dflux)

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
      if (present(fcb_calcFluxFCTSys_sim) .and. present(rgroupFEMSet)) then

        ! Check if group finite element set and stabilisation structure are compatible
        if ((rafcstab%NEQ   .ne. rgroupFEMSet%NEQ) .or.&
            (rafcstab%NEDGE .ne. rgroupFEMSet%NEDGE)) then
          call output_line('Stabilisation and group finite element set are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
          call sys_halt()
        end if
        call gfem_getbase_DcoeffsAtEdge(rgroupFEMSet, p_DcoeffsAtEdge)
        buseCallback = .true.
      else

        ! Check if stabilisation provides edge data
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEVALUES) .eq. 0) then
          call output_line('Stabilisation does not provide edge data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
          call sys_halt()
        end if

        call afcstab_getbase_DcoeffsAtEdge(rafcstab, p_Dcoefficients)
        buseCallback = .false.
      end if

      ! Assemble spatial part of raw-antidiffusive fluxes
      if (buseCallback) then
        call doFluxesByCallbackDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DcoeffsAtEdge, p_Dx, dscale, bclear, p_Dflux)
      else
        call doFluxesByCoeffsDP(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_Dcoefficients, p_Dx, dscale, bclear, p_Dflux)
      end if

      !-------------------------------------------------------------------------

      if (rafcstab%cprelimitingType .eq. AFCSTAB_PRELIMITING_STD) then
        ! Compute fluxes for standard prelimiting based on the
        ! low-order solution which serves as predictor
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        call doDifferencesDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_DfluxPrel)

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
        call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_DxTimeDeriv,&
            dscale, .false., p_Dflux)
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
        call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
        call lsyssc_getbase_double(rmatrix, p_Dmatrix)

        ! Clear vector and assemble antidiffusive fluxes
        call lalg_clearVector(p_Dflux)
        call doFluxesByMatrixDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, dscale, .true., p_Dflux)

        ! Set specifiers for raw antidiffusive fluxes
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)

      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
        call sys_halt()
      end if


    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_buildFluxFCTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the edge-by-edge array DcoefficientsAtEdge

    subroutine doFluxesByCoeffsDP(IedgeList, NEDGE, NEQ, NVAR,&
        DcoefficientsAtEdge, Dx, dscale, bclear, Dflux)

      ! input parameters
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux

      ! local variables
      integer :: iedge,i,j

      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux)

      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute raw antidiffusive flux
            Dflux(:,iedge) = DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(:,i)-Dx(:,j))
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
            Dflux(:,iedge) = Dflux(:,iedge)&
                + DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(:,i)-Dx(:,j))
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
            Dflux(:,iedge) = DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(:,j)-Dx(:,i))
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
            Dflux(:,iedge) = Dflux(:,iedge)&
                + DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(:,j)-Dx(:,i))
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
            Dflux(:,iedge) = dscale * DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(:,i)-Dx(:,j))
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
            Dflux(:,iedge) = Dflux(:,iedge)&
                + dscale * DcoefficientsAtEdge(1:NVAR,iedge) * (Dx(:,i)-Dx(:,j))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesByCoeffsDP

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix stored in Dmatrix

    subroutine doFluxesByMatrixDP(IedgeList, NEDGE, NEQ, NVAR,&
        Dmatrix, Dx, dscale, bclear, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:), intent(in) :: Dmatrix
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      logical, intent(in) :: bclear

      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux

      ! local variables
      integer :: iedge,ij,i,j


      if (dscale .eq. 0.0_DP) then

        ! Do we have to clear the vector?
        if (bclear) call lalg_clearVector(Dflux)


      elseif (dscale .eq. 1.0_DP) then

        if (bclear) then
          !$omp parallel do default(shared) private(i,j,ij)&
          !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
          do iedge = 1, NEDGE

            ! Determine indices
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)
            ij = IedgeList(3,iedge)

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dmatrix(ij) * (Dx(:,i)-Dx(:,j))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dflux(:,iedge)&
                           + Dmatrix(ij) * (Dx(:,i)-Dx(:,j))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dmatrix(ij) * (Dx(:,j)-Dx(:,i))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dflux(:,iedge)&
                           + Dmatrix(ij) * (Dx(:,j)-Dx(:,i))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = dscale * Dmatrix(ij) * (Dx(:,j)-Dx(:,i))
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

            ! Compute the raw antidiffusives fluxes
            Dflux(:,iedge) = Dflux(:,iedge)&
                           + dscale * Dmatrix(ij) * (Dx(:,j)-Dx(:,i))
          end do
          !$omp end parallel do
        end if

      end if

    end subroutine doFluxesByMatrixDP

    !**************************************************************
    ! Assemble raw antidiffusive fluxes with aid of callback function

    subroutine doFluxesByCallbackDP(IedgeList, NEDGE, NEQ, NVAR,&
        DcoeffsAtEdge, Dx, dscale, bclear, Dflux)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DcoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      logical, intent(in) :: bclear

      ! input/output parameters
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:), pointer :: DfluxAtEdge

      ! local variables
      integer :: idx,iedge,IEDGEset,IEDGEmax


      if (dscale .eq. 0.0_DP) then

        if (bclear) call lalg_clearVector(Dflux)

      elseif (bclear) then

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,idx,iedge,IEDGEmax)

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))

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
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              Dflux(:,IEDGEset:IEDGEmax), rcollection)
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        !$omp end parallel

      else   ! bclear = .false.

        !$omp parallel default(shared)&
        !$omp private(DdataAtEdge,DfluxAtEdge,idx,iedge,IEDGEmax)

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,p_rperfconfig%NEDGESIM))
        allocate(DfluxAtEdge(NVAR,p_rperfconfig%NEDGESIM))

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
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute internodal fluxes
          call fcb_calcFluxFCTSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DcoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DfluxAtEdge(:,1:IEDGEmax-IEDGEset+1), rcollection)

          ! Loop through all edges in the current set
          do idx = 1, IEDGEmax-IEDGEset+1

            ! Get actual edge number
            iedge = idx+IEDGEset-1

            ! Add antidiffusive fluxes
            Dflux(:,iedge) = Dflux(:,iedge) + DfluxAtEdge(:,idx)
          end do
        end do
        !$omp end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge, DfluxAtEdge)
        !$omp end parallel

      end if

    end subroutine doFluxesByCallbackDP

    !**************************************************************
    ! Assemble fluxes for classical prelimiting.

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doDifferencesDP(IedgeList, NEDGE, NEQ, NVAR, Dx, Dflux)

      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR

      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      integer :: iedge,i,j

      !$omp parallel do default(shared) private(i,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)
      do iedge = 1, NEDGE

        ! Determine indices
        i  = IedgeList(1,iedge)
        j  = IedgeList(2,iedge)

        ! Compute solution difference; in contrast to the literature,
        ! we compute the solution difference $u_i-u_j$ and check if
        ! $F_{ij}(U_i-U_j)<0$ in the prelimiting step.
        Dflux(:,iedge) = Dx(:,i)-Dx(:,j)
      end do
      !$omp end parallel do

    end subroutine doDifferencesDP

  end subroutine afcsys_buildFluxFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_failsafeFCTBlock(rafcstab, rmatrix, rx, dscale, dtol,&
      ioperationSpec, bisAccepted, dfactor, nsteps, CvariableNames,&
      fcb_extractVariable, fcb_calcFailsafe, rxBackup, rvectorCorr,&
      rvectorTmp, rcollection, rperfconfig)

!<description>
    ! This subroutine performs failsafe flux limiting as described in
    ! the paper by Kuzmin, Moeller, Shadid, and Shashkov: "Failsafe
    ! flux limiting and constrained data projection for equations of
    ! gas dynamics" Journal of Computational Physics, vol. 229,
    ! Nov. 2010, p. 8766-8779.
    !
    ! Note that vector rx must provide the low-order solution and the
    ! stabilisation structure rafcstab must provide the edgewise
    ! limiting factors and the raw-antidiffusive fluxes. That is, the
    ! standard FCT algorithm must have been evoked before except for
    ! the last correction step which applied the limited antidiffusive
    ! fluxes to the solution vector.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! Lumped mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Tolerance parameter
    real(DP), intent(in) :: dtol

    ! Operation specification tag. This is a bitfield coming from an
    ! OR combination of different AFCSTAB_FAILSAFE_xxxx constants and
    ! specifies which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: Failsafe correction factor
    ! If not present then the correction factor is computed as follows
    ! dfactor = istep/nsteps
    ! If nsteps is also not present, then dfactor=0 is used.
    real(DP), intent(in), optional :: dfactor

    ! OPTIONAL: Number of failsafe correction steps
    ! If not present then a single step is performed using the value
    ! dfactor as failsafe correction factor.
    integer, intent(in), optional :: nsteps

    ! OPTIONAL: Array containing the names of control variables
    character(len=*), dimension(:), intent(in), optional :: CvariableNames

    ! OPTIONAL: Callback function to extract variables
    include 'intf_extractVariable.inc'
    optional :: fcb_extractVariable

    ! OPTIONAL: Callback function to calculate the failsafe
    ! correction. If present, then most tasks of this subroutine are
    ! skipped and deligated to the user-defined callback function
    include 'intf_calcFailsafe.inc'
    optional :: fcb_calcFailsafe

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Solution vector
    ! On input:  This vector must provide the low-order solution to
    !            which failsafe correction resorts in the worst case
    ! On Output: This vector is the corrected solution vector after
    !            failsafe correction has been applied. By explicitly
    !            unspecifying AFCSTAB_FAILSAFEALGO_CORRECT this vector
    !            will be the original low-order solution without any
    !            failsafe correction being applied to it.
    type(t_vectorBlock), intent(inout) :: rx

    ! OPTIONAL: Collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: Auxiliary vector storing a backup of rx
    type(t_vectorBlock), intent(inout), target, optional :: rxBackup

    ! OPTIONAL: Scalar vector of length equal to the number of edges
    ! which contains the correction factors resulting from the
    ! failsafe limiter. If not present, then the failsafe correction
    ! is directly applied to the correction factors provided by the
    ! stabilisation structure rafcstab.
    type(t_vectorScalar), intent(inout), target, optional :: rvectorCorr

    ! OPTIONAL: Auxiliary block vector storing internal data. If not
    ! present, then internal memory is allocated.
    type(t_vectorBlock), intent(inout), target, optional :: rvectorTmp
!</inputoutput>

!<output>
    ! Flag which is TRUE if the computed solution vector
    ! is accepted by the failsafe correction procedure
    logical, intent(out) :: bisAccepted
!</output>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rxBackup,p_rvectorTmp
    type(t_vectorScalar), pointer :: p_rvectorCorr
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dbeta,p_Dx,p_DxBackup
    real(DP), dimension(:), pointer :: p_Ddata,p_DdataLbound,p_DdataUbound
    real(DP), dimension(:), pointer :: p_DdataMatrix,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    real(DP) :: dcorr
    integer :: ivariable,nvariables,istep
    logical :: bextractVariables

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    ! Check if block vector contains only one block.
    if (rx%nblocks .eq. 1) then
      if (present(rxBackup)) then
        call afcsys_failsafeFCTScalar(&
            rafcstab, rmatrix, rx%RvectorBlock(1), dscale, dtol,&
            ioperationSpec, bisAccepted, dfactor, nsteps,&
            CvariableNames, fcb_extractVariable, fcb_calcFailsafe,&
            rxBackup%RvectorBlock(1),&
            rvectorCorr, rvectorTmp, rcollection, rperfconfig)
      else
        call afcsys_failsafeFCTScalar(&
            rafcstab, rmatrix, rx%RvectorBlock(1), dscale, dtol,&
            ioperationSpec, bisAccepted, dfactor, nsteps,&
            CvariableNames, fcb_extractVariable, fcb_calcFailsafe,&
            rvectorCorr=rvectorCorr, rvectorTmp=rvectorTmp,&
            rcollection=rcollection, rperfconfig=rperfconfig)
      end if
      return
    end if

    ! Set pointer to performance configuration
    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)     .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_failsafeFCTBlock')
      call sys_halt()
    else
      ! Set pointers
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rmatrix, p_DdataMatrix)
      call lsysbl_getbase_double(rx, p_Dx)
    end if

    ! Determine strategy for failsafe correction
    if (present(fcb_calcFailsafe)) then

      ! Call user-defined callback function
      if (present(rvectorCorr)) then
        call lsyssc_getbase_double(rvectorCorr, p_Dbeta)
        call fcb_calcFailsafe(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NEQ,&
            rafcstab%NVAR, ioperationSpec, dscale, dtol, p_DdataMatrix,&
            p_Dx, p_Dalpha, p_Dflux, bisAccepted, p_Dbeta, rcollection)
      else
        call fcb_calcFailsafe(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NEQ,&
            rafcstab%NVAR, ioperationSpec, dscale, dtol, p_DdataMatrix,&
            p_Dx, p_Dalpha, p_Dflux, bisAccepted, rcollection=rcollection)
      end if

      ! That's it return
      return
    end if


    ! Failsafe correction is performed internally.
    if (present(CvariableNames) .and. present(fcb_extractVariable)) then
      ! Failsafe correction is performed in terms of user-defined
      ! variables which are extracted from the given solution vector
      nvariables = size(CvariableNames)
      bextractVariables = .true.
    else
      ! Failsafe correction is performed in terms of the solution
      ! variables so that no variable extraction is required
      nvariables = rx%nblocks
      bextractVariables = .false.
    end if

    ! Set pointer to given vector rvectorCorr or create new one
    if (present(rvectorCorr)) then
      p_rvectorCorr => rvectorCorr
    else
      allocate(p_rvectorCorr)
      call lsyssc_createVector(p_rvectorCorr, rafcstab%NEDGE, .false.)
    end if
    call lsyssc_getbase_double(p_rvectorCorr, p_Dbeta)

    ! Set pointer to given vector rxBackup or create new one
    if (present(rxBackup)) then
      p_rxBackup => rxBackup
    else
      allocate(p_rxBackup)
      call lsysbl_createVectorBlock(rx, p_rxBackup, .false.)
    end if
    call lsysbl_getbase_double(p_rxBackup, p_DxBackup)
    call lalg_copyVector(p_Dx, p_DxBackup)

    ! Set pointer to given vector rvectorTmp or create new one (if
    ! required). The dimension of the temporal vector depends on the
    ! failsafe correction strategy. If the vector is provided then
    ! it is tacidly assumed that it has the correct dimension
    if (present(rvectorTmp)) then
      p_rvectorTmp => rvectorTmp
    else
      allocate(p_rvectorTmp)
      if (bextractVariables) then
        call lsysbl_createVectorBlock(p_rvectorTmp, rafcstab%NEQ,&
            3*nvariables, .false.)
      else
        call lsysbl_createVectorBlock(p_rvectorTmp, rafcstab%NEQ,&
            2*nvariables, .false.)
      end if
    end if

    ! Set pointers to subvectors
    if (bextractVariables) then
      call lsysbl_getbase_double(p_rvectorTmp, p_Ddata, 1, nvariables)
      call lsysbl_getbase_double(&
          p_rvectorTmp, p_DdataLbound, nvariables+1, 2*nvariables)
      call lsysbl_getbase_double(&
          p_rvectorTmp, p_DdataUbound, 2*nvariables+1, 3*nvariables)
    else
      call lsysbl_getbase_double(rx, p_Ddata)
      call lsysbl_getbase_double(p_rvectorTmp, p_DdataLbound, 1, nvariables)
      call lsysbl_getbase_double(&
          p_rvectorTmp, p_DdataUbound, nvariables+1, 2*nvariables)
    end if


    !---------------------------------------------------------------------------
    ! The failsafe FCT algorithm is split into the following steps
    ! which can be skipped and performed externally by the user.
    !
    ! 1) Initialise the edgewise correction factors by unity
    !
    ! 2) Initialise the upper and lower bounds (if required)
    !
    ! 3) Perform failsafe correction
    !
    ! 4) Reject solution of required
    !---------------------------------------------------------------------------

    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_INITBETA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 1) Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------
      call lalg_setVector(p_Dbeta, 1.0_DP)
    end if


    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 2) Initialise the upper and lower bounds
      !-------------------------------------------------------------------------
      if (bextractVariables) then
        ! Extract variables by user-defined callback function
        do ivariable = 1, nvariables
          call fcb_extractVariable(rx, trim(CvariableNames(ivariable)),&
              p_rvectorTmp%RvectorBlock(ivariable))
        end do
      end if

      ! Compute the upper and lower solution bounds
      call doBoundsDP(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEDGE, rafcstab%NEQ, nvariables,&
          p_Ddata, p_DdataLbound, p_DdataUbound)
    end if


    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_LIMIT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 3) Perform failsafe limiting
      !-------------------------------------------------------------------------

      if (present(nsteps)) then

        ! Perform prescribed number of failsafe steps
        failsafe: do istep = 1, nsteps

          ! Determine correction factor for this step
          dcorr = 1.0_DP-real(istep,DP)/real(nsteps,DP)

          ! Apply the corrected fluxes to the low-order solution
          call doCorrectScaleByMass(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
              p_DdataMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Dx)

          ! Recompute the control variables (if required)
          if (bextractVariables) then
            ! Extract variables by user-defined callback function
            do ivariable = 1, nvariables
              call fcb_extractVariable(rx, trim(CvariableNames(ivariable)),&
                  p_rvectorTmp%RvectorBlock(ivariable))
            end do
          end if

          ! Compute failsafe correction factors
          call doFailsafeLimitDP(p_IedgeList, rafcstab%NEDGE,&
              rafcstab%NEQ, nvariables, p_Ddata, p_DdataLbound,&
              p_DdataUbound, dcorr, dtol, p_Dbeta, bisAccepted)

          if (bisAccepted) exit failsafe

          ! Solution is not acceptable. Another failsafe correction
          ! step starting from the low-order solution is performed
          call lalg_copyVector(p_DxBackup, p_Dx)
        end do failsafe

        ! If failsafe correction did not lead to an acceptable
        ! solution then we have to recompute it using zero as failsafe
        ! correction factor which was set in the very last step
        if (.not.bisAccepted)&
            call doCorrectScaleByMass(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
            p_DdataMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Dx)

      else ! nsteps not present

        ! Perform exactly one failsafe correction step with
        if (present(dfactor)) then
          dcorr = dfactor
        else
          dcorr = 0.0_DP
        end if

        ! Apply the corrected fluxes to the low-order solution
        call doCorrectScaleByMass(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
            p_DdataMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Dx)

        ! Recompute the control variables (if required)
        if (bextractVariables) then
          ! Extract variables by user-defined callback function
          do ivariable = 1, nvariables
            call fcb_extractVariable(rx, trim(CvariableNames(ivariable)),&
                p_rvectorTmp%RvectorBlock(ivariable))
          end do
        end if

        ! Compute failsafe correction factor
        call doFailsafeLimitDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Ddata, p_DdataLbound,&
            p_DdataUbound, dcorr, dtol, p_Dbeta, bisAccepted)
      end if
    end if

    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_CORRECT) .eq. 0) then
      !-------------------------------------------------------------------------
      ! 4) Apply failsafe correction. In the current implementation
      !    the correciton term has already been applied to the
      !    low-order solution. Therefe, we have to overwrite the
      !    solution vector rx by the low-order backup rxBackup if
      !    this routine is enforced NOT to apply the correction.
      !-------------------------------------------------------------------------
      call lalg_copyVector(p_DxBackup, p_Dx)
    end if


    ! Release temporal vectors
    if (.not.present(rvectorCorr)) then
      call lsyssc_releaseVector(p_rvectorCorr); deallocate(p_rvectorCorr)
    end if

    if (.not.present(rxBackup)) then
      call lsysbl_releaseVector(p_rxBackup); deallocate(p_rxBackup)
    end if

    if (.not.present(rvectorTmp)) then
      call lsysbl_releaseVector(p_rvectorTmp); deallocate(p_rvectorTmp)
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Compute the local upper and lower bounds based on the double
    ! values array Dx evaluated at the neighbouring nodes

    subroutine doBoundsDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVARfailsafe, Dx, Dlbound, Dubound)

      ! input parameters
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), dimension(NEQ,NVARfailsafe), intent(in) :: Dx
      integer, intent(in) :: NEDGE,NEQ,NVARfailsafe

      ! output parameters
      real(DP), dimension(NEQ,NVARfailsafe), intent(out) :: Dlbound, Dubound

      ! local variables
      integer :: i,iedge,igroup,j

      !$omp parallel sections
      !$omp section
      call lalg_copyVector(Dx, Dlbound)
      !$omp section
      call lalg_copyVector(Dx, Dubound)
      !$omp end parallel sections

      !$omp parallel default(shared) private(i,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)

          ! Compute minimum/maximum value of neighboring nodes
          Dlbound(i,:) = min(Dlbound(i,:), Dx(j,:))
          Dlbound(j,:) = min(Dlbound(j,:), Dx(i,:))
          Dubound(i,:) = max(Dubound(i,:), Dx(j,:))
          Dubound(j,:) = max(Dubound(j,:), Dx(i,:))
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doBoundsDP

    !**************************************************************
    ! Premultiply the correction factor Dalpha by the failsafe factor
    ! Dbeta and limit the raw antidiffusive fluxes Dflux by the
    ! resulting net correction factor Dcorr = Dalpha*Dbeta. Apply the
    ! corrected antidiffusive fluxes to the low-order solution Dx and
    ! scale each entry by the entry of the lumped mass matrix.

    subroutine doCorrectScaleByMass(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, dscale, ML, Dalpha, Dbeta, Dflux, Dx)

      ! input parameters
      real(DP), dimension(:), intent(in) :: ML,Dalpha,Dbeta
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dx

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j

      if (dscale .eq. 0.0_DP) then
        ! Do nothing
        return

      elseif (dscale .eq. 1.0_DP) then

        !$omp parallel default(shared) private(i,j,F_ij)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1

          ! Do nothing for empty groups
          if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

          ! Loop over all edges
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)

            ! Compute portion of corrected antidiffusive flux
            F_ij = Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)

            ! Remove flux from solution
            Dx(i,:) = Dx(i,:) + F_ij/ML(i)
            Dx(j,:) = Dx(j,:) - F_ij/ML(j)
          end do
          !$omp end do

        end do ! igroup
        !$omp end parallel

      else ! dscale /= 1.0

        !$omp parallel default(shared) private(i,j,F_ij)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1

          ! Do nothing for empty groups
          if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

          ! Loop over all edges
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)

            ! Compute portion of corrected antidiffusive flux
            F_ij = dscale * Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)

            ! Remove flux from solution
            Dx(i,:) = Dx(i,:) + F_ij/ML(i)
            Dx(j,:) = Dx(j,:) - F_ij/ML(j)
          end do
          !$omp end do

        end do ! igroup
        !$omp end parallel

      end if

    end subroutine doCorrectScaleByMass

    !**************************************************************
    ! Compute the edgewise failsafe correction factors

    subroutine doFailsafeLimitDP(IedgeList, NEDGE, NEQ, NVARfailsafe,&
        Dx, Dlbound, Dubound, dcorr, dtol, Dbeta, baccept)

      ! input parameters
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), dimension(NEQ,NVARfailsafe), intent(in) :: Dx,Dlbound,Dubound
      real(DP), intent(in) :: dcorr,dtol
      integer, intent(in) :: NEDGE,NEQ,NVARfailsafe

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dbeta

      ! output parameters
      logical, intent(out) :: baccept

      ! local variables
      integer :: iedge,i,j,ivar

      ! Initialisation
      baccept = .true.

      ! Loop over all variables
      !$omp parallel do default(shared) private(i,j,iedge)&
      !$omp reduction(.and.:baccept) schedule(static,1)
      do ivar = 1, NVARfailsafe

        do iedge = 1, NEDGE

          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)

          ! Check if solution exceeds
          if ((Dx(i,ivar) .lt. Dlbound(i,ivar)-dtol) .or.&
              (Dx(j,ivar) .lt. Dlbound(j,ivar)-dtol) .or.&
              (Dx(i,ivar) .gt. Dubound(i,ivar)+dtol) .or.&
              (Dx(j,ivar) .gt. Dubound(j,ivar)+dtol)) then
            Dbeta(iedge) = dcorr
            baccept = .false.
          end if
        end do
      end do
      !$omp end parallel do

    end subroutine doFailsafeLimitDP

  end subroutine afcsys_failsafeFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine afcsys_failsafeFCTScalar(rafcstab, rmatrix, rx, dscale, dtol,&
      ioperationSpec, bisAccepted, dfactor, nsteps, CvariableNames,&
      fcb_extractVariable, fcb_calcFailsafe, rxBackup, rvectorCorr,&
      rvectorTmp, rcollection, rperfconfig)

!<description>
    ! This subroutine performs failsafe flux limiting as described in
    ! the paper by Kuzmin, Moeller, Shadid, and Shashkov: "Failsafe
    ! flux limiting and constrained data projection for equations of
    ! gas dynamics" Journal of Computational Physics, vol. 229,
    ! Nov. 2010, p. 8766-8779.
!</description>

!<input>
    ! Stabilisation structure
    type(t_afcstab), intent(in) :: rafcstab

    ! Lumped mass matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Tolerance parameter
    real(DP), intent(in) :: dtol

    ! Operation specification tag. This is a bitfield coming from an
    ! OR combination of different AFCSTAB_FAILSAFE_xxxx constants and
    ! specifies which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: Failsafe correction factor
    ! If not present then the correction factor is computed as follows
    ! dfactor = istep/nsteps
    ! If nsteps is also not present, then dfactor=0 is used.
    real(DP), intent(in), optional :: dfactor

    ! OPTIONAL: Number of failsafe correction steps
    ! If not present then a single step is performed using the value
    ! dfactor as failsafe correction factor.
    integer, intent(in), optional :: nsteps

    ! OPTIONAL: Array containing the names of control variables
    character(len=*), dimension(:), intent(in), optional :: CvariableNames

    ! OPTIONAL: Callback function to extract variables
    include 'intf_extractVariable.inc'
    optional :: fcb_extractVariable

    ! OPTIONAL: Callback function to calculate the failsafe correction
    include 'intf_calcFailsafe.inc'
    optional :: fcb_calcFailsafe

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Solution vector
    ! On input:  This vector must provide the low-order solution to
    !            which failsafe correction resorts in the worst case
    ! On Output: This vector is the corrected solution vector after
    !            failsafe correction has been applied. By explicitly
    !            unspecifying AFCSTAB_FAILSAFEALGO_CORRECT this vector
    !            will be the original low-order solution without any
    !            failsafe correction being applied to it.
    type(t_vectorScalar), intent(inout) :: rx

    ! OPTIONAL: Collection structure
    type(t_collection), intent(inout), optional :: rcollection

    ! OPTIONAL: Auxiliary vector storing a backup of rx
    type(t_vectorScalar), intent(inout), target, optional :: rxBackup

    ! OPTIONAL: Scalar vector of length equal to the number of edges
    ! which contains the correction factors resulting from the
    ! failsafe limiter. If not present, then the failsafe correction
    ! is directly applied to the correction factors provided by the
    ! stabilisation structure rafcstab.
    type(t_vectorScalar), intent(inout), target, optional :: rvectorCorr

    ! OPTIONAL: Auxiliary block vector storing internal data. If not
    ! present, then internal memory is allocated.
    type(t_vectorBlock), intent(inout), target, optional :: rvectorTmp
!</inputoutput>

!<output>
    ! Flag which is TRUE if the computed solution vector
    ! is accepted by the failsafe correction procedure
    logical, intent(out) :: bisAccepted
!</output>
!</subroutine>

    ! local variables
    type(t_blockDiscretisation) :: rblockDiscr
    type(t_vectorBlock) :: rxBlock
    type(t_vectorBlock), pointer :: p_rvectorTmp
    type(t_vectorScalar), pointer :: p_rvectorCorr,p_rxBackup
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dbeta,p_Dx,p_DxBackup
    real(DP), dimension(:), pointer :: p_Ddata,p_DdataLbound,p_DdataUbound
    real(DP), dimension(:), pointer :: p_DdataMatrix,p_Dflux
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    real(DP) :: dcorr
    integer :: ivariable,nvariables,istep,ivar,ieq
    logical :: bextractVariables

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => afcsys_perfconfig
    end if

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIST)        .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)     .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'afcsys_failsafeFCTScalar')
      call sys_halt()
    else
      ! Set pointers
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rmatrix, p_DdataMatrix)
      call lsyssc_getbase_double(rx, p_Dx)
    end if

    ! Determine strategy for failsafe correction
    if (present(fcb_calcFailsafe)) then

      ! Call user-defined callback function
      if (present(rvectorCorr)) then
        call lsyssc_getbase_double(rvectorCorr, p_Dbeta)
        call fcb_calcFailsafe(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NVAR,&
            rafcstab%NEQ, ioperationSpec, dscale, dtol, p_DdataMatrix,&
            p_Dx, p_Dalpha, p_Dflux, bisAccepted, p_Dbeta, rcollection)
      else
        call fcb_calcFailsafe(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, rafcstab%NVAR,&
            rafcstab%NEQ, ioperationSpec, dscale, dtol, p_DdataMatrix,&
            p_Dx, p_Dalpha, p_Dflux, bisAccepted, rcollection=rcollection)
      end if

      ! That's it return
      return
    end if


    ! Failsafe correction is performed internally.
    if (present(CvariableNames) .and. present(fcb_extractVariable)) then
      ! Failsafe correction is performed in terms of user-defined
      ! variables which are extracted from the given solution vector
      nvariables = size(CvariableNames)
      bextractVariables = .true.

      ! Create temporal 1-block vector required as argument for the
      ! callback function which extracts variables from the solution
      if (associated(rx%p_rspatialDiscr)) then
        call spdiscr_createBlockDiscrInd(rx%p_rspatialDiscr, rblockDiscr)
        call lsysbl_createVecFromScalar(rx, rxBlock, rblockDiscr)
      else
        call lsysbl_createVecFromScalar(rx, rxBlock)
      end if
    else
      ! Failsafe correction is performed in terms of the solution
      ! variables so that no variable extraction is required
      nvariables = rx%NVAR
      bextractVariables = .false.
    end if

    ! Set pointer to given vector rvectorCorr or create new one
    if (present(rvectorCorr)) then
      p_rvectorCorr => rvectorCorr
    else
      allocate(p_rvectorCorr)
      call lsyssc_createVector(p_rvectorCorr, rafcstab%NEDGE, .false.)
    end if
    call lsyssc_getbase_double(p_rvectorCorr, p_Dbeta)

    ! Set pointer to given vector rxBackup or create new one
    if (present(rxBackup)) then
      p_rxBackup => rxBackup
      call lsyssc_getbase_double(p_rxBackup, p_DxBackup)
      call lalg_copyVector(p_Dx, p_DxBackup)
    else
      allocate(p_rxBackup)
      call lsyssc_copyVector(rx, p_rxBackup)
      call lsyssc_getbase_double(p_rxBackup, p_DxBackup)
    end if


    ! Set pointer to given vector rvectorTmp or create new one (if
    ! required). The dimension of the temporal vector depends on the
    ! failsafe correction strategy. If the vector is provided then
    ! it is tacidly assumed that it has the correct dimension
    if (present(rvectorTmp)) then
      p_rvectorTmp => rvectorTmp
    else
      allocate(p_rvectorTmp)
      call lsysbl_createVectorBlock(p_rvectorTmp, rafcstab%NEQ,&
          3*nvariables, .false.)
    end if

    ! Set pointers to subvectors
    call lsysbl_getbase_double(p_rvectorTmp, p_Ddata, 1, nvariables)
    call lsysbl_getbase_double(&
        p_rvectorTmp, p_DdataLbound, nvariables+1, 2*nvariables)
    call lsysbl_getbase_double(&
        p_rvectorTmp, p_DdataUbound, 2*nvariables+1, 3*nvariables)


    !---------------------------------------------------------------------------
    ! The failsafe FCT algorithm is split into the following steps
    ! which can be skipped and performed externally by the user.
    !
    ! 1) Initialise the edgewise correction factors by unity
    !
    ! 2) Initialise the upper and lower bounds (if required)
    !
    ! 3) Perform failsafe correction
    !
    ! 4) Reject solution of required
    !---------------------------------------------------------------------------

    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_INITBETA) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 1) Initialise the edgewise correction factors by unity
      !-------------------------------------------------------------------------
      call lalg_setVector(p_Dbeta, 1.0_DP)
    end if


    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_BOUNDS) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 2) Initialise the upper and lower bounds
      !-------------------------------------------------------------------------
      if (bextractVariables) then
        ! Extract variables by user-defined callback function
        do ivariable = 1, nvariables
          call fcb_extractVariable(rxBlock, trim(CvariableNames(ivariable)),&
              p_rvectorTmp%RvectorBlock(ivariable))
        end do
      else
        ! The solution vector is store in interleaved format whereas
        ! the upper and lower bounds are stored in block format.
        ! Therefore, we have to convert the solution vector.
        do ivar = 1,rafcstab%NVAR
          do ieq = 1,rafcstab%NEQ
            p_Ddata((ivar-1)*rafcstab%NVAR+ieq) = p_Dx((ieq-1)*rafcstab%NEQ+ivar)
          end do
        end do
      end if

      ! Compute the upper and lower solution bounds
      call doBoundsDP(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEDGE, rafcstab%NEQ, nvariables,&
          p_Ddata, p_DdataLbound, p_DdataUbound)
    end if


    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_LIMIT) .ne. 0) then
      !-------------------------------------------------------------------------
      ! 3) Perform failsafe limiting
      !-------------------------------------------------------------------------

      if (present(nsteps)) then

        ! Perform prescribed number of failsafe steps
        failsafe: do istep = 1, nsteps

          ! Determine correction factor for this step
          dcorr = 1.0_DP-real(istep,DP)/real(nsteps,DP)

          ! Apply the corrected fluxes to the low-order solution
          call doCorrectScaleByMass(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
              p_DdataMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Dx)

          ! Recompute the control variables (if required)
          if (bextractVariables) then
            ! Extract variables by user-defined callback function
            do ivariable = 1, nvariables
              call fcb_extractVariable(rxBlock, trim(CvariableNames(ivariable)),&
                  p_rvectorTmp%RvectorBlock(ivariable))
            end do
          end if

          ! Compute failsafe correction factors
          call doFailsafeLimitDP(p_IedgeList, rafcstab%NEDGE,&
              rafcstab%NEQ, nvariables, p_Ddata, p_DdataLbound,&
              p_DdataUbound, dcorr, dtol, p_Dbeta, bisAccepted)

          if (bisAccepted) exit failsafe

          ! Solution is not acceptable. Another failsafe correction
          ! step starting from the low-order solution is performed
          call lalg_copyVector(p_DxBackup, p_Dx)
        end do failsafe

        ! If failsafe correction did not lead to an acceptable
        ! solution then we have to recompute it using zero as failsafe
        ! correction factor which was set in the very last step
        if (.not.bisAccepted)&
            call doCorrectScaleByMass(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
            p_DdataMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Dx)

      else ! nsteps not present

        ! Perform exactly one failsafe correction step with
        if (present(dfactor)) then
          dcorr = dfactor
        else
          dcorr = 0.0_DP
        end if

        ! Apply the corrected fluxes to the low-order solution
        call doCorrectScaleByMass(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, dscale,&
            p_DdataMatrix, p_Dalpha, p_Dbeta, p_Dflux, p_Dx)

        ! Recompute the control variables (if required)
        if (bextractVariables) then
          ! Extract variables by user-defined callback function
          do ivariable = 1, nvariables
            call fcb_extractVariable(rxBlock, trim(CvariableNames(ivariable)),&
                p_rvectorTmp%RvectorBlock(ivariable))
          end do
        end if

        ! Compute failsafe correction factor
        call doFailsafeLimitDP(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Ddata, p_DdataLbound,&
            p_DdataUbound, dcorr, dtol, p_Dbeta, bisAccepted)
      end if
    end if

    if (iand(ioperationSpec, AFCSTAB_FAILSAFEALGO_CORRECT) .eq. 0) then
      !-------------------------------------------------------------------------
      ! 4) Apply failsafe correction. In the current implementation
      !    the correciton term has already been applied to the
      !    low-order solution. Therefe, we have to overwrite the
      !    solution vector rx by the low-order backup rxBackup if
      !    this routine is enforced NOT to apply the correction.
      !-------------------------------------------------------------------------
      call lalg_copyVector(p_DxBackup, p_Dx)
    end if


    ! Release temporal vectors
    if (.not.present(rvectorCorr)) then
      call lsyssc_releaseVector(p_rvectorCorr); deallocate(p_rvectorCorr)
    end if

    if (.not.present(rxBackup)) then
      call lsyssc_releaseVector(p_rxBackup); deallocate(p_rxBackup)
    end if

    if (.not.present(rvectorTmp)) then
      call lsysbl_releaseVector(p_rvectorTmp); deallocate(p_rvectorTmp)
    end if

    ! Release temporal 1-block vector
    call lsysbl_releaseVector(rxBlock)
    if (associated(rx%p_rspatialDiscr))&
        call spdiscr_releaseBlockDiscr(rblockDiscr)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Compute the local upper and lower bounds based on the double
    ! values array Dx evaluated at the neighbouring nodes

    subroutine doBoundsDP(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVARfailsafe, Dx, Dlbound, Dubound)

      ! input parameters
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), dimension(NEQ,NVARfailsafe), intent(in) :: Dx
      integer, intent(in) :: NEDGE,NEQ,NVARfailsafe

      ! output parameters
      real(DP), dimension(NEQ,NVARfailsafe), intent(out) :: Dlbound, Dubound

      ! local variables
      integer :: i,iedge,igroup,j

      !$omp parallel sections
      !$omp section
      call lalg_copyVector(Dx, Dlbound)
      !$omp section
      call lalg_copyVector(Dx, Dubound)
      !$omp end parallel sections

      !$omp parallel default(shared) private(i,j)&
      !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronise memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over all edges
        !$omp do
        do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

          ! Get node numbers
          i  = IedgeList(1,iedge)
          j  = IedgeList(2,iedge)

          ! Compute minimum/maximum value of neighboring nodes
          Dlbound(i,:) = min(Dlbound(i,:), Dx(j,:))
          Dlbound(j,:) = min(Dlbound(j,:), Dx(i,:))
          Dubound(i,:) = max(Dubound(i,:), Dx(j,:))
          Dubound(j,:) = max(Dubound(j,:), Dx(i,:))
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doBoundsDP

    !**************************************************************
    ! Premultiply the correction factor Dalpha by the failsafe factor
    ! Dbeta and limit the raw antidiffusive fluxes Dflux by the
    ! resulting net correction factor Dcorr = Dalpha*Dbeta. Apply the
    ! corrected antidiffusive fluxes to the low-order solution Dx and
    ! scale each entry by the entry of the lumped mass matrix.

    subroutine doCorrectScaleByMass(IedgeListIdx, IedgeList,&
        NEDGE, NEQ, NVAR, dscale, ML, Dalpha, Dbeta, Dflux, Dx)

      ! input parameters
      real(DP), dimension(:), intent(in) :: ML,Dalpha,Dbeta
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEDGE,NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dx

      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: i,iedge,igroup,j

      if (dscale .eq. 0.0_DP) then
        ! Do nothing
        return

      elseif (dscale .eq. 1.0_DP) then

        !$omp parallel default(shared) private(i,j,F_ij)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1

          ! Do nothing for empty groups
          if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

          ! Loop over all edges
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)

            ! Compute portion of corrected antidiffusive flux
            F_ij = Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)

            ! Remove flux from solution
            Dx(:,i) = Dx(:,i) + F_ij/ML(i)
            Dx(:,j) = Dx(:,j) - F_ij/ML(j)
          end do
          !$omp end do

        end do ! igroup
        !$omp end parallel

      else ! dscale /= 1.0

        !$omp parallel default(shared) private(i,j,F_ij)&
        !$omp if (NEDGE > p_rperfconfig%NEDGEMIN_OMP)

        ! Loop over the edge groups and process all edges of one group
        ! in parallel without the need to synchronise memory access
        do igroup = 1, size(IedgeListIdx)-1

          ! Do nothing for empty groups
          if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

          ! Loop over all edges
          !$omp do
          do iedge = IedgeListIdx(igroup), IedgeListIdx(igroup+1)-1

            ! Get node numbers
            i  = IedgeList(1,iedge)
            j  = IedgeList(2,iedge)

            ! Compute portion of corrected antidiffusive flux
            F_ij = dscale * Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)

            ! Remove flux from solution
            Dx(:,i) = Dx(:,i) + F_ij/ML(i)
            Dx(:,j) = Dx(:,j) - F_ij/ML(j)
          end do
          !$omp end do

        end do ! igroup
        !$omp end parallel

      end if

    end subroutine doCorrectScaleByMass

    !**************************************************************
    ! Compute the edgewise failsafe correction factors

    subroutine doFailsafeLimitDP(IedgeList, NEDGE, NEQ, NVARfailsafe,&
        Dx, Dlbound, Dubound, dcorr, dtol, Dbeta, baccept)

      ! input parameters
      integer, dimension(:,:), intent(in) :: IedgeList
      real(DP), dimension(NEQ,NVARfailsafe), intent(in) :: Dx,Dlbound,Dubound
      real(DP), intent(in) :: dcorr,dtol
      integer, intent(in) :: NEDGE,NEQ,NVARfailsafe

      ! input/output parameters
      real(DP), dimension(:), intent(inout) :: Dbeta

      ! output parameters
      logical, intent(out) :: baccept

      ! local variables
      integer :: iedge,i,j,ivar

      ! Initialisation
      baccept = .true.

      ! Loop over all variables
      !$omp parallel do default(shared) private(i,j,iedge)&
      !$omp reduction(.and.:baccept) schedule(static,1)
      do ivar = 1, NVARfailsafe

        do iedge = 1, NEDGE

          ! Get node numbers
          i  = IedgeList(1, iedge)
          j  = IedgeList(2, iedge)

          ! Check if solution exceeds
          if ((Dx(i,ivar) .lt. Dlbound(i,ivar)-dtol) .or.&
              (Dx(j,ivar) .lt. Dlbound(j,ivar)-dtol) .or.&
              (Dx(i,ivar) .gt. Dubound(i,ivar)+dtol) .or.&
              (Dx(j,ivar) .gt. Dubound(j,ivar)+dtol)) then
            Dbeta(iedge) = dcorr
            baccept = .false.
          end if
        end do
      end do
      !$omp end parallel do

    end subroutine doFailsafeLimitDP

  end subroutine afcsys_failsafeFCTScalar

  !*****************************************************************************

!<function>

  elemental function minmodDP(da,db)

!<description>
    ! This function computes minmod(a,b) following the implementation
    ! suggested by N. Robidoux, M. Gong, J. Cupitt, A. Turcott, and
    ! K. Martinez in "CPU, SMP and GPU Implementations of Nohalo Level
    ! 1, a Fast Co-Convex Antialiasing Image Resampler", ACM 2009
!</description>

!<input>
    real(DP), intent(in) :: da,db
!</input>

!<result>
    real(DP) :: minmodDP
!</result>

    ! local variables
    real(DP) :: dsa,dsb

    dsa = sign(1._DP,da)
    dsb = sign(1._DP,db)

    minmodDP = 0.5_DP * (dsa+dsb) * min(da*dsa, db*dsb)

  end function minmodDP

  !*****************************************************************************

!<function>

  elemental function minmodSP(fa,fb)

!<description>
    ! This function computes minmod(a,b) following the implementation
    ! suggested by N. Robidoux, M. Gong, J. Cupitt, A. Turcott, and
    ! K. Martinez in "CPU, SMP and GPU Implementations of Nohalo Level
    ! 1, a Fast Co-Convex Antialiasing Image Resampler", ACM 2009
!</description>

!<input>
    real(SP), intent(in) :: fa,fb
!</input>

!<result>
    real(SP) :: minmodSP
!</result>

    ! local variables
    real(SP) :: fsa,fsb

    fsa = sign(1._SP,fa)
    fsb = sign(1._SP,fb)

    minmodSP = 0.5_SP * (fsa+fsb) * min(fa*fsa, fb*fsb)

  end function minmodSP

end module afcstabsystem

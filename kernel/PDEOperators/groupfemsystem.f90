!##############################################################################
!# ****************************************************************************
!# <name> groupfemsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines for applying the
!# group-finite element formulation to systems of conservation laws.
!# The technique was proposed by C.A.J. Fletcher in:
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
!# starting point for systems of conservation laws, the reader is
!# referred to the book chapter
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
!# afcstabilisation. The initialisation as a system stabilisation
!# structure is done by the routine gfsys_initStabilisation.
!#
!# There are three types of routines. The gfsys_buildDivOperator
!# routines can be used to assemble the discrete divergence operators
!# resulting from the standard Galerkin finite element discretisation
!# plus some discretely defined artificial viscosities. This technique
!# represents a generalisation of the discrete upwinding approach
!# which has been used to construct upwind finite element scheme for
!# scalar conservation laws (see module groupfemscalar for details).
!#
!# The second type of routines is given by
!# gfsys_buildDivVectorXXX. They can be used to update/initialise the
!# divergence term applying some sort of algebraic flux
!# correction. Importantly, the family of AFC schemes gives rise to
!# nonlinear algebraic equations that need to be solved
!# iteratively. Thus, it is usefull to build the compensating
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
!# 1.) gfsys_initStabilisation = gfsys_initStabilisationScalar /
!#                               gfsys_initStabilisationBlock
!#     -> initialize the stabilisation structure
!#
!# 2.) gfsys_buildDivOperator = gfsys_buildDivOperatorScalar /
!#                              gfsys_buildDivOperatorBlock
!#     -> assembles the global operator that results from the discretisation
!#        of the divergence term $div(F)$ by means of the Galerkin method
!#        and some sort of artificial dissipation (if required)
!#
!# 3.) gfsys_buildDivVector = gfsys_buildDivVectorScalar /
!#                            gfsys_buildDivVectorBlock
!#     -> assembles the divergence vector
!#
!# 4.) gfsys_buildDivVectorTVD = gfsys_buildDivVecTVDScalar /
!#                               gfsys_buildDivVecTVDBlock
!#     -> assembles the divergence term for FEM-TVD stabilisation
!#
!# 5.) gfsys_buildDivVectorFCT = gfsys_buildDivVecFCTScalar /
!#                               gfsys_buildDivVecFCTBlock
!#     -> assembles the divergence term for FEM-FCT stabilisation
!#
!# 6.) gfsys_buildFluxFCT = gfsys_buildFluxFCTScalar /
!#                          gfsys_buildFluxFCTBlock
!#     -> assembles the raw antidiffusive flux for FEM-FCT stabilisation
!#
!# 7.) gfsys_failsafeFCT = gfsys_failsafeFCTScalar /
!#                         gfsys_failsafeFCTBlock
!#     -> perform failsafe limiting of FCT type
!#
!# </purpose>
!##############################################################################

module groupfemsystem

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

  public :: gfsys_initStabilisation
  public :: gfsys_buildDivOperator
  public :: gfsys_buildDivVector
  public :: gfsys_buildDivVectorTVD
  public :: gfsys_buildDivVectorFCT
  public :: gfsys_buildFluxFCT
  public :: gfsys_failsafeFCT

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of nodes to handle simultaneously when building matrices
#ifndef GFSYS_NEQSIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEQSIM = 128
#else
  integer, public            :: GFSYS_NEQSIM = 128
#endif
#endif

  ! Number of edges to handle simultaneously when building matrices
#ifndef GFSYS_NEDGESIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEDGESIM = 64
#else
  integer, public            :: GFSYS_NEDGESIM = 64
#endif
#endif
  
  ! Minimum number of nodes for OpenMP parallelisation: If the number of
  ! nodes is below this value, then no parallelisation is performed.
#ifndef GFSYS_NEQMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEQMIN_OMP = 1000
#else
  integer, public            :: GFSYS_NEQMIN_OMP = 1000
#endif
#endif

  ! Minimum number of edges for OpenMP parallelisation: If the number of
  ! edges is below this value, then no parallelisation is performed.
#ifndef GFSYS_NEDGEMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: GFSYS_NEDGEMIN_OMP = 1000
#else
  integer, public            :: GFSYS_NEDGEMIN_OMP = 1000
#endif
#endif
!</constantblock>

!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
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

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface gfsys_initStabilisation
    module procedure gfsys_initStabilisationScalar
    module procedure gfsys_initStabilisationBlock
  end interface

  interface gfsys_buildDivOperator
     module procedure gfsys_buildDivOperatorScalar
     module procedure gfsys_buildDivOperatorBlock
  end interface

  interface gfsys_buildDivVector
    module procedure gfsys_buildDivVectorScalar
    module procedure gfsys_buildDivVectorBlock
  end interface

  interface gfsys_buildDivVectorTVD
    module procedure gfsys_buildDivVecTVDScalar
    module procedure gfsys_buildDivVecTVDBlock
  end interface

  interface gfsys_buildDivVectorFCT
    module procedure gfsys_buildDivVecFCTScalar
    module procedure gfsys_buildDivVecFCTBlock
  end interface

  interface gfsys_buildFluxFCT
    module procedure gfsys_buildFluxFCTScalar
    module procedure gfsys_buildFluxFCTBlock
  end interface

  interface gfsys_failsafeFCT
    module procedure gfsys_failsafeFCTScalar
    module procedure gfsys_failsafeFCTBlock
  end interface

  interface gfsys_combineFluxes
    module procedure gfsys_combineFluxesDble
    module procedure gfsys_combineFluxesSngl
  end interface

  interface gfsys_limit
    module procedure gfsys_limitUnboundedDble
    module procedure gfsys_limitUnboundedSngl
    module procedure gfsys_limitBoundedDble
    module procedure gfsys_limitBoundedSngl
  end interface

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

contains

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_initStabilisationBlock(rmatrixBlockTemplate,&
      rafcstab, NVARtransformed, rblockDiscretisation)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as block matrix. If this matrix
    ! contains only one block, then the scalar counterpart of this
    ! subroutine is called with the corresponding scalar submatrix.
!</description>

!<input>
    ! template block matrix
    type(t_matrixBlock), intent(in) :: rmatrixBlockTemplate

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed

    ! OPTIONAL: block discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_blockDiscretisation), intent(in), optional :: rblockDiscretisation
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>


    ! Check if block matrix has only one block
    if ((rmatrixBlockTemplate%nblocksPerCol .eq. 1) .and.&
        (rmatrixBlockTemplate%nblocksPerRow .eq. 1)) then
      if (present(rblockDiscretisation)) then
        call gfsys_initStabilisationScalar(&
            rmatrixBlockTemplate%RmatrixBlock(1,1), rafcstab, NVARtransformed,&
            rblockDiscretisation%RspatialDiscr(1))
      else
        call gfsys_initStabilisationScalar(&
            rmatrixBlockTemplate%RmatrixBlock(1,1), rafcstab, NVARtransformed)
      end if
      return
    end if

    ! Check that number of columns equans number of rows
    if (rmatrixBlockTemplate%nblocksPerCol .ne.&
        rmatrixBlockTemplate%nblocksPerRow) then
      call output_line('Block matrix must have equal number of columns and rows!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      call sys_halt()
    end if

    ! Check if matrix exhibits group structure
    if (rmatrixBlockTemplate%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      call sys_halt()
    end if

    ! Set atomic data from first block
    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rmatrixBlockTemplate%nblocksPerCol
    end if
    rafcstab%NVAR  = rmatrixBlockTemplate%nblocksPerCol
    rafcstab%NEQ   = rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ
    rafcstab%NEDGE = (rmatrixBlockTemplate%RmatrixBlock(1,1)%NA-&
                      rmatrixBlockTemplate%RmatrixBlock(1,1)%NEQ)/2
    rafcstab%NNVEDGE = 0

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)


    case (AFCSTAB_TVD)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)
      

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_NLINFCT_IMPLICIT)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the 3 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux0(rafcstab)
      call afcstab_allocFlux(rafcstab)

      ! We need the edgewise vector for the prelimited fluxes
      if ((rafcstab%ctypePrelimiting      .ne. AFCSTAB_PRELIMITING_NONE).or.&
          (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT)) then
        call afcstab_allocFluxPrel(rafcstab)
      end if

      ! We need the nodal block vector for the low-order predictor
      allocate(rafcstab%p_rvectorPredictor)
      if (present(rblockDiscretisation)) then
        call lsysbl_createVectorBlock(rblockDiscretisation,&
            rafcstab%p_rvectorPredictor, .false., rafcstab%cdataType)
      else
        call lsysbl_createVectorBlock(rafcstab%p_rvectorPredictor,&
            rafcstab%NEQ, rafcstab%NVAR, .false., rafcstab%cdataType)
      end if


    case (AFCSTAB_LINFCT,&
          AFCSTAB_LINFCT_MASS)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)

      ! We need the edgewise vector if raw antidiffusive fluxes should
      ! be prelimited
      if (rafcstab%ctypePrelimiting .ne. AFCSTAB_PRELIMITING_NONE) then
        call afcstab_allocFluxPrel(rafcstab)
      end if

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationBlock')
      call sys_halt()
    end select

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

  end subroutine gfsys_initStabilisationBlock

  ! *****************************************************************************

!<subroutine>
  subroutine gfsys_initStabilisationScalar(rmatrixTemplate, rafcstab,&
      NVARtransformed, rspatialDiscretisation)

!<description>
    ! This subroutine initialises the discrete stabilisation structure
    ! for use as a scalar stabilisation. The template matrix is used
    ! to determine the number of equations and the number of edges.
    !
    ! Note that the matrix is required as scalar matrix. It can be
    ! stored in interleave format.
!</description>

!<input>
    ! template matrix
    type(t_matrixScalar), intent(in) :: rmatrixTemplate

    ! OPTIONAL: number of transformed variables
    ! If not present, then the number of variables
    ! NVAR is taken from the template matrix
    integer, intent(in), optional :: NVARtransformed

    ! OPTIONAL: spatial discretisation structure which is used to
    ! create auxiliary vectors, e.g., for the predictor
    type(t_spatialDiscretisation), intent(in), optional :: rspatialDiscretisation
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rvectorTmp


    ! Set atomic data
    if (present(NVARtransformed)) then
      rafcstab%NVARtransformed = NVARtransformed
    else
      rafcstab%NVARtransformed = rmatrixTemplate%NVAR
    end if
    rafcstab%NVAR  = rmatrixTemplate%NVAR
    rafcstab%NEQ   = rmatrixTemplate%NEQ
    rafcstab%NEDGE = (rmatrixTemplate%NA-rmatrixTemplate%NEQ)/2
    rafcstab%NNVEDGE = 0

    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_GALERKIN,&
          AFCSTAB_UPWIND)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)


    case (AFCSTAB_TVD)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      
    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_NLINFCT_IMPLICIT)

      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)

      ! We need the 3 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux0(rafcstab)
      call afcstab_allocFlux(rafcstab)

      ! We need the edgewise vector for the prelimited fluxes
      if ((rafcstab%ctypePrelimiting      .ne. AFCSTAB_PRELIMITING_NONE) .or.&
          (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT)) then
        call afcstab_allocFluxPrel(rafcstab)
      end if

      ! We need the nodal block vector for the low-order predictor
      allocate(rafcstab%p_rvectorPredictor)
      if (present(rspatialDiscretisation)) then
        call lsyssc_createVecByDiscr(rspatialDiscretisation,&
            rvectorTmp, rafcstab%NVAR, .false., rafcstab%cdataType)
      else
        call lsyssc_createVector(rvectorTmp,&
            rafcstab%NEQ, rafcstab%NVAR, .false., rafcstab%cdataType)
      end if
      
      ! Convert into 1-block vector
      call lsysbl_convertVecFromScalar(rvectorTmp, rafcstab%p_rvectorPredictor)
      call lsyssc_releaseVector(rvectorTmp)
      

    case (AFCSTAB_LINFCT,&
          AFCSTAB_LINFCT_MASS)
      
      ! Handle for IedgeListIdx and IedgeList: (/i,j,ij,ji/)
      call afcstab_allocEdgeStructure(rafcstab,4)

      ! We need the 6 nodal vectors P, Q and R each for '+' and '-'
      call afcstab_allocVectorsPQR(rafcstab)
      
      ! We need the 2 edgewise vectors for the correction factors and the fluxes
      call afcstab_allocAlpha(rafcstab)
      call afcstab_allocFlux(rafcstab)
      
      ! We need the edgewise vector if raw antidiffusive fluxes should
      ! be prelimited
      if (rafcstab%ctypePrelimiting .ne. AFCSTAB_PRELIMITING_NONE) then
        call afcstab_allocFluxPrel(rafcstab)
      end if
      
      ! We need the nodal block vector for the low-order predictor
      allocate(rafcstab%p_rvectorPredictor)
      if (present(rspatialDiscretisation)) then
        call lsyssc_createVecByDiscr(rspatialDiscretisation,&
            rvectorTmp, rafcstab%NVAR, .false., ST_DOUBLE)
      else
        call lsyssc_createVector(rvectorTmp,&
            rafcstab%NEQ, rafcstab%NVAR, .false., ST_DOUBLE)
      end if
      
      ! Convert into 1-block vector
      call lsysbl_convertVecFromScalar(rvectorTmp, rafcstab%p_rvectorPredictor)
      call lsyssc_releaseVector(rvectorTmp)
      
      
    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_initStabilisationScalar')
      call sys_halt()
    end select

    ! Set specifier
    rafcstab%istabilisationSpec = AFCSTAB_INITIALISED

  end subroutine gfsys_initStabilisationScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivOperatorBlock(rafcstab, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, dscale,&
      bclear, rdivOp, rcollection, fcb_calcDivOperator)

!<description>
    ! This subroutine assembles the discrete divergence operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !   <tex> $$ \nabla\cdot{\bf F}(u) $$ </tex>
    !
    ! where ${\bf f}(U)$ is a user-defined flux function for the
    ! multi-component field $U$.
    !
    ! This routine can be used to apply the following discretisations:
    !
    ! (1) the standard Galerkin finite element method
    !     which will be referred to as high-order approximation
    !
    ! (2) discrete upwinding for hyperbolic systems which results from
    !     the conservative elimination of negative eigenvalues from
    !     the Galerkin operator. This technique is for instance
    !     described in the reference:
    !
    !     D. Kuzmin and M. Moeller, Algebraic flux correction II. Compressible
    !     Euler equations, In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 207-250.
    !
    ! Note that this routine is designed for block matrices/vectors.
    ! If there is only one block, then the corresponding scalar routine
    ! is called. Otherwise, the global operator is treated as block matrix.
    ! This block matrix has to be in group structure, that is, the structure
    ! of subblock(1,1) will serve as template for all other submatrices.
!</description>

!<input>
    ! The solution vector
    type(t_vectorBlock), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute local matrices
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'

    ! OPTIONAL: callback  function to overwrite the standard operation
    include 'intf_calcDivOperator.inc'
    optional :: fcb_calcDivOperator
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The divergence operator
    type(t_matrixBlock), intent(inout) :: rdivOp

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_array), dimension(:,:), allocatable  :: rarray
    real(DP), dimension(:), pointer :: p_Dx
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx,p_Kdiagonal
    logical :: bisFullMatrix


    ! Check if block vector contains only one block, if global
    ! operator is stored in interleave format and no user-defined
    ! callback function is provided.
    if (.not.present(fcb_calcDivOperator) .and.&
        (rx%nblocks .eq. 1)               .and.&
        (rdivOp%nblocksPerCol .eq. 1) .and.&
        (rdivOp%nblocksPerRow .eq. 1)) then
      call gfsys_buildDivOperatorScalar(rafcstab, rx%RvectorBlock(1),&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim,&
          dscale, bclear, rdivOp%RmatrixBlock(1,1), rcollection,&
          fcb_calcDivOperator)
      return
    end if

    ! Check if block matrix exhibits group structure
    if (rdivOp%imatrixSpec .ne. LSYSBS_MSPEC_GROUPMATRIX) then
      call output_line('Block matrix must have group structure!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivOperator)) then
      ! Call used-defined assembly
      call fcb_calcDivOperator(rafcstab, rx, rdivOp, dscale, bclear,&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, rcollection)
    else
      ! Allocate temporal memory
      allocate(rarray(rx%nblocks,rx%nblocks))
      
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtNode(rafcstab, p_DmatrixCoeffsAtNode)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call afcstab_getbase_array(rdivOp, rarray, bisFullMatrix)
      call lsysbl_getbase_double(rx, p_Dx)
      
      ! What kind of matrix are we?
      select case(rdivOp%RmatrixBlock(1,1)%cmatrixFormat)
      case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)
        !-----------------------------------------------------------------------
        ! Matrix format 7 and 9
        !-----------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rdivOp%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX7) then
          call lsyssc_getbase_Kld(rdivOp%RmatrixBlock(1,1), p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rdivOp%RmatrixBlock(1,1), p_Kdiagonal)
        end if
        
        ! What type of matrix are we?
        if (bisFullMatrix) then
          
          call doOperatorMat79(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rx%nblocks, p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge,&
              p_Dx, dscale, bclear, rarray)
          
        else   ! bisFullMatrix == no
          
          call doOperatorMat79Diag(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rx%nblocks, p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge,&
              p_Dx, dscale, bclear, rarray)
          
        end if   ! bisFullMatrix
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorBlock')
        call sys_halt()
      end select
      
      ! Deallocate temporal memory
      deallocate(rarray)
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble block-diagonal divergence operator K
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79Diag(Kdiagonal, IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx, dscale, bclear, K)
      
      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: K

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IverticesAtNode
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,IEQmax,IEQset,idx,igroup
      integer :: i,iedge,ii,ij,ivar,ji,jj

      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtEdge,DcoefficientsAtNode,DdataAtEdge,&
      !$omp         DdataAtNode,IEDGEmax,IEQmax,IverticesAtNode,i,idx,&
      !$omp         iedge,ii,ij,ivar,ji,jj)

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IverticesAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(NVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSYS_NEQSIM

        ! We always handle GFSYS_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSYS_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSYS_NEQSIM)
        
        ! Loop through all equations in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IEQmax-IEQset+1

          ! Get actual equation number
          i = idx+IEQset-1
          
          ! Get position of diagonal entry
          ii = Kdiagonal(i)
          
          ! Fill auxiliary arrays
          IverticesAtNode(1,idx) = i
          IverticesAtNode(2,idx) = ii
          DdataAtNode(:,idx)     = Dx(i,:)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IverticesAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)
        
        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IverticesAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              K(ivar,ivar)%p_Ddata(ii) = DcoefficientsAtNode(ivar,1,idx)
            end do
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IverticesAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              K(ivar,ivar)%p_Ddata(ii) = K(ivar,ivar)%p_Ddata(ii)+&
                  DcoefficientsAtNode(ivar,1,idx)
            end do
          end do

        end if
      end do
      !$omp end do
      
      ! Deallocate temporal memory
      deallocate(IverticesAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcoefficientsAtEdge(NVAR,3,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
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
              do ivar = 1, NVAR
                K(ivar,ivar)%p_Ddata(ii) = K(ivar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(jj) = K(ivar,ivar)%p_Ddata(jj)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(ij) = &
                    DcoefficientsAtEdge(ivar,2,idx) + DcoefficientsAtEdge(ivar,1,idx) 
                K(ivar,ivar)%p_Ddata(ji) = &
                    DcoefficientsAtEdge(ivar,3,idx) + DcoefficientsAtEdge(ivar,1,idx) 
              end do
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
              do ivar = 1, NVAR
                K(ivar,ivar)%p_Ddata(ii) = K(ivar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(jj) = K(ivar,ivar)%p_Ddata(jj)-&
                    DcoefficientsAtEdge(ivar,1,idx)
                K(ivar,ivar)%p_Ddata(ij) = K(ivar,ivar)%p_Ddata(ij)+&
                    DcoefficientsAtEdge(ivar,2,idx) + DcoefficientsAtEdge(ivar,1,idx) 
                K(ivar,ivar)%p_Ddata(ji) = K(ivar,ivar)%p_Ddata(ji)+&
                    DcoefficientsAtEdge(ivar,3,idx) + DcoefficientsAtEdge(ivar,1,idx) 
              end do
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)
      !$omp end parallel

    end subroutine doOperatorMat79Diag

    
    !**************************************************************
    ! Assemble divergence operator K in 1D
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79(Kdiagonal, IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx, dscale, bclear, K)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx,Kdiagonal
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      type(t_array), dimension(:,:), intent(inout) :: K

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IverticesAtNode
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,IEQmax,IEQset,idx,igroup
      integer :: i,iedge,ii,ij,ijpos,ivar,ji,jj,jvar

      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtEdge,DcoefficientsAtNode,DdataAtEdge,&
      !$omp         DdataAtNode,IEDGEmax,IEQmax,IverticesAtNode,i,idx,&
      !$omp         iedge,ii,ij,ijpos,ivar,jvar,ji,jj)

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IverticesAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(NVAR*NVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSYS_NEQSIM

        ! We always handle GFSYS_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSYS_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSYS_NEQSIM)
        
        ! Loop through all equations in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IEQmax-IEQset+1

          ! Get actual equation number
          i = idx+IEQset-1

          ! Get position of diagonal entry
          ii = Kdiagonal(i)
          
          ! Fill auxiliary arrays
          IverticesAtNode(1,idx) = i
          IverticesAtNode(2,idx) = ii
          DdataAtNode(:,idx)     = Dx(i,:)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IverticesAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IverticesAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              do jvar = 1, NVAR
                ijpos = NVAR*(ivar-1)+jvar
                K(jvar,ivar)%p_Ddata(ii) = DcoefficientsAtNode(ijpos,1,idx)
              end do
            end do
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1
            
            ! Get position of diagonal entry
            ii = IverticesAtNode(2,idx)
            
            ! Update the diagonal coefficient
            do ivar = 1, NVAR
              do jvar = 1, NVAR
                ijpos = NVAR*(ivar-1)+jvar
                K(jvar,ivar)%p_Ddata(ii) = K(jvar,ivar)%p_Ddata(ii)+&
                    DcoefficientsAtNode(ijpos,1,idx)
              end do
            end do
          end do

        end if
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(IverticesAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcoefficientsAtEdge(NVAR*NVAR,3,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
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
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  K(jvar,ivar)%p_Ddata(ii) = K(jvar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(jj) = K(jvar,ivar)%p_Ddata(jj)-&
                      DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(ij) = &
                      DcoefficientsAtEdge(ijpos,2,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                  K(jvar,ivar)%p_Ddata(ji) = &
                      DcoefficientsAtEdge(ijpos,3,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                end do
              end do
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
              do ivar = 1, NVAR
                do jvar = 1, NVAR
                  ijpos = NVAR*(ivar-1)+jvar
                  K(jvar,ivar)%p_Ddata(ii) = K(jvar,ivar)%p_Ddata(ii)-&
                    DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(jj) = K(jvar,ivar)%p_Ddata(jj)-&
                      DcoefficientsAtEdge(ijpos,1,idx)
                  K(jvar,ivar)%p_Ddata(ij) = K(jvar,ivar)%p_Ddata(ij)+&
                      DcoefficientsAtEdge(ijpos,2,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                  K(jvar,ivar)%p_Ddata(ji) = K(jvar,ivar)%p_Ddata(ji)+&
                      DcoefficientsAtEdge(ijpos,3,idx) + DcoefficientsAtEdge(ijpos,1,idx) 
                end do
              end do
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)
      !$omp end parallel

    end subroutine doOperatorMat79

  end subroutine gfsys_buildDivOperatorBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivOperatorScalar(rafcstab, rx,&
      fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, dscale,&
      bclear, rdivOp, rcollection, fcb_calcDivOperator)

!<description>
    ! This subroutine assembles the discrete divergence operator which results
    ! from the group finite element formulation of the continuous problem
    !
    !   <tex> $$ \nabla\cdot{\bf F}(u) $$ </tex>
    !
    ! where ${\bf f}(U)$ is a user-defined flux function for the
    ! multi-component field $U$.
    !
    ! This routine can be used to apply the following discretisations:
    !
    ! (1) the standard Galerkin finite element method
    !     which will be referred to as high-order approximation
    !
    ! (2) discrete upwinding for hyperbolic systems which results from
    !     the conservative elimination of negative eigenvalues from
    !     the Galerkin operator. This technique is for instance
    !     described in the reference:
    !
    !     D. Kuzmin and M. Moeller, Algebraic flux correction II. Compressible
    !     Euler equations, In: D. Kuzmin et al. (eds), Flux-Corrected
    !     Transport: Principles,  Algorithms, and Applications,
    !     Springer, 2005, 207-250.
    !
    ! Note that this routine requires the scalar matrices/vectors
    ! are stored in the interleave format.
!</description>

!<input>
    ! The solution vector
    type(t_vectorScalar), intent(in) :: rx

    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Switch for matrix assembly
    ! TRUE  : clear matrix before assembly
    ! FLASE : assemble matrix in an additive way
    logical, intent(in) :: bclear

    ! Callback functions to compute local matrices
    include 'intf_calcMatrixDiagSys_sim.inc'
    include 'intf_calcMatrixSys_sim.inc'

    ! OPTIONAL: callback  function to overwrite the standard operation
    include 'intf_calcDivOperator.inc'
    optional :: fcb_calcDivOperator
!</input>

!<inputoutput>
    ! The stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! The divergence operator
    type(t_matrixScalar), intent(inout) :: rdivOp

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock
    type(t_matrixBlock) :: rdivOpBlock
    real(DP), dimension(:), pointer :: p_DivOp,p_Dx
    real(DP), dimension(:,:), pointer :: p_DmatrixCoeffsAtNode
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kdiagonal,p_IedgeListIdx


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivOperator)) then
      ! Create auxiliary 1-block vector and matrix
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createMatFromScalar(rdivOp, rdivOpBlock)
      
      ! Call user-defined assembly
      call fcb_calcDivOperator(rafcstab, rxBlock, rdivOpBlock, dscale, bclear,&
          fcb_calcMatrixDiagSys_sim, fcb_calcMatrixSys_sim, rcollection)
      
      ! Release auxiliary 1-block vector and matrix
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseMatrix(rdivOpBlock)
    else
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtNode(rafcstab, p_DmatrixCoeffsAtNode)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call lsyssc_getbase_double(rdivOp, p_DivOp)
      call lsyssc_getbase_double(rx, p_Dx)
      
      ! What kind of matrix are we?
      select case(rdivOp%cmatrixFormat)
      case(LSYSSC_MATRIX7INTL, LSYSSC_MATRIX9INTL)
        !-------------------------------------------------------------------------
        ! Matrix format 7 and 9 interleaved
        !-------------------------------------------------------------------------
        
        ! Set diagonal pointer
        if (rdivOp%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
          call lsyssc_getbase_Kld(rdivOp, p_Kdiagonal)
        else
          call lsyssc_getbase_Kdiagonal(rdivOp, p_Kdiagonal)
        end if
        
        ! What type of matrix are we?
        select case(rdivOp%cinterleavematrixFormat)
          
        case (LSYSSC_MATRIX1)
          call doOperatorMat79(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rdivOp%NA, rx%NVAR, rx%NVAR*rx%NVAR,&
              p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, bclear, p_DivOp)
          
        case (LSYSSC_MATRIXD)
          call doOperatorMat79(p_Kdiagonal, p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEQ, rdivOp%NA, rx%NVAR, rx%NVAR,&
              p_DmatrixCoeffsAtNode, p_DmatrixCoeffsAtEdge, p_Dx, dscale, bclear, p_DivOp)
          
        case DEFAULT
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
          call sys_halt()
        end select
        
      case DEFAULT
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivOperatorScalar')
        call sys_halt()
      end select
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence operator K
    ! All matrices are stored in matrix format 7 and 9

    subroutine doOperatorMat79(Kdiagonal, IedgeListIdx, IedgeList,&
        NEQ, NA, NVAR, MVAR, DmatrixCoeffsAtNode, DmatrixCoeffsAtEdge, Dx,&
        dscale, bclear, K)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:), intent(in) :: DmatrixCoeffsAtNode
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      logical, intent(in) :: bclear
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: Kdiagonal,IedgeListIdx
      integer, intent(in) :: NEQ,NA,NVAR,MVAR

      ! input/output parameters
      real(DP), dimension(MVAR,NA), intent(inout) :: K

      ! auxiliary arrays
      real(DP), dimension(:,:), pointer :: DdataAtNode
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtNode
      real(DP), dimension(:,:,:), pointer :: DcoefficientsAtEdge
      integer, dimension(:,:), pointer  :: IverticesAtNode
      
      ! local variables
      integer :: igroup,idx,IEQset,IEQmax,IEDGEset,IEDGEmax
      integer :: i,ii,jj,ij,ji,iedge

      
      !$omp parallel default(shared)&
      !$omp private(DcoefficientsAtEdge,DcoefficientsAtNode,DdataAtEdge,&
      !$omp         DdataAtNode,IEDGEmax,IEQmax,IverticesAtNode,i,idx,&
      !$omp         iedge,ii,ij,ji,jj)

      !-------------------------------------------------------------------------
      ! Assemble diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(IverticesAtNode(2,GFSYS_NEQSIM))
      allocate(DdataAtNode(NVAR,GFSYS_NEQSIM))
      allocate(DcoefficientsAtNode(MVAR,1,GFSYS_NEQSIM))

      ! Loop over the equations
      !$omp do schedule(static,1)
      do IEQset = 1, NEQ, GFSYS_NEQSIM

        ! We always handle GFSYS_NEQSIM equations simultaneously.
        ! How many equations have we actually here?
        ! Get the maximum equation number, such that we handle 
        ! at most GFSYS_NEQSIM equations simultaneously.
        
        IEQmax = min(NEQ, IEQset-1+GFSYS_NEQSIM)
        
        ! Loop through all equations in the current set
        ! and prepare the auxiliary arrays
        do idx = 1, IEQmax-IEQset+1

          ! Get actual equation number
          i = idx+IEQset-1

          ! Get position of diagonal entry
          ii = Kdiagonal(i)
          
          ! Fill auxiliary arrays
          IverticesAtNode(1,idx) = i
          IverticesAtNode(2,idx) = ii
          DdataAtNode(:,idx)     = Dx(:,i)
        end do

        ! Use callback function to compute diagonal entries
        call fcb_calcMatrixDiagSys_sim(&
            DdataAtNode(:,1:IEQmax-IEQset+1),&
            DmatrixCoeffsAtNode(:,IEQset:IEQmax),&
            IverticesAtNode(:,1:IEQmax-IEQset+1),&
            dscale, IEQmax-IEQset+1,&
            DcoefficientsAtNode(:,:,1:IEQmax-IEQset+1), rcollection)

        ! Loop through all equations in the current set
        ! and scatter the entries to the global matrix
        if (bclear) then
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IverticesAtNode(2,idx)
            
            ! Update the diagonal coefficient
            K(:,ii) = DcoefficientsAtNode(:,1,idx)
          end do

        else   ! do not clear matrix
          do idx = 1, IEQmax-IEQset+1

            ! Get position of diagonal entry
            ii = IverticesAtNode(2,idx)
            
            ! Update the diagonal coefficient
            K(:,ii) = K(:,ii) + DcoefficientsAtNode(:,1,idx)
          end do

        end if
      end do
      !$omp end do

      ! Deallocate temporal memory
      deallocate(IverticesAtNode)
      deallocate(DdataAtNode)
      deallocate(DcoefficientsAtNode)

      !-------------------------------------------------------------------------
      ! Assemble off-diagonal entries
      !-------------------------------------------------------------------------

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcoefficientsAtEdge(MVAR,3,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1,IEDGEset-1+GFSYS_NEDGESIM)

          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do
          
          ! Use callback function to compute off-diagonal entries
          call fcb_calcMatrixSys_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
              IedgeList(:,IEDGEset:IEDGEmax),&
              dscale, IEDGEmax-IEDGEset+1,&
              DcoefficientsAtEdge(:,:,1:IEDGEmax-IEDGEset+1), rcollection)
          
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
              K(:,ii) = K(:,ii) - DcoefficientsAtEdge(:,1,idx)
              K(:,jj) = K(:,jj) - DcoefficientsAtEdge(:,1,idx)
              K(:,ij) = DcoefficientsAtEdge(:,2,idx) + DcoefficientsAtEdge(:,1,idx) 
              K(:,ji) = DcoefficientsAtEdge(:,3,idx) + DcoefficientsAtEdge(:,1,idx) 
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
              K(:,ii) = K(:,ii) - DcoefficientsAtEdge(:,1,idx)
              K(:,jj) = K(:,jj) - DcoefficientsAtEdge(:,1,idx)
              K(:,ij) = K(:,ij) + DcoefficientsAtEdge(:,2,idx) + DcoefficientsAtEdge(:,1,idx) 
              K(:,ji) = K(:,ji) + DcoefficientsAtEdge(:,3,idx) + DcoefficientsAtEdge(:,1,idx) 
            end do

          end if
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DcoefficientsAtEdge)
      !$omp end parallel

    end subroutine doOperatorMat79

  end subroutine gfsys_buildDivOperatorScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVectorBlock(rafcstab, rx, fcb_calcFlux_sim,&
      dscale, bclear, ry, rcollection, fcb_calcDivVector)

!<description>
    ! This subroutine assembles the divergence vector for block vectors.
    ! If the vector contains only one block, then the scalar
    ! counterpart of this routine is called with the scalar subvector.
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

    ! callback function to compute local fluxes
    include 'intf_calcFlux_sim.inc'

    ! OPTIONAL: callback  function to overwrite the standard operation
    include 'intf_calcDivVector.inc'
    optional :: fcb_calcDivVector
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx


    ! Check if block vectors contain only one block and no
    ! user-defined callback function is provided.
    if (.not.present(fcb_calcDivVector) .and.&
        (rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call gfsys_buildDivVectorScalar(rafcstab, rx%RvectorBlock(1),&
          fcb_calcFlux_sim, dscale, bclear, ry%RvectorBlock(1), rcollection)
      return
    end if

    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorBlock')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivVector)) then
      ! Call used-defined assembly
      call fcb_calcDivVector(rafcstab, rx, ry, dscale, bclear,&
          fcb_calcFlux_sim, rcollection)
    else
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call lsysbl_getbase_double(rx, p_Dx)
      call lsysbl_getbase_double(ry, p_Dy)
      
      ! Clear vector?
      if (bclear) call lsysbl_clearVector(ry)

      ! Assemble the divergence vector
      call doDivVector(p_IedgeListIdx, p_IedgeList, rafcstab%NEQ,&
          rx%nblocks, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dy)
    end if

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector
    
    subroutine doDivVector(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j


      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      !$omp end parallel
      
    end subroutine doDivVector
    
  end subroutine gfsys_buildDivVectorBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVectorScalar(rafcstab, rx, fcb_calcFlux_sim,&
      dscale, bclear, ry, rcollection, fcb_calcDivVector)

!<description>
    ! This subroutine assembles the divergence vector. Note that the
    ! vectors are required as scalar vectors which are stored in the
    ! interleave format.
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

    ! callback functions to compute local fluxes
    include 'intf_calcFlux_sim.inc'

    ! OPTIONAL: callback function to overwrite the standard operation
    include 'intf_calcDivVector.inc'
    optional :: fcb_calcDivVector
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rxBlock,ryBlock
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx


    ! Check if stabilisation has been initialised
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVectorScalar')
      call sys_halt()
    end if

    ! Check if user-defined assembly is provided
    if (present(fcb_calcDivVector)) then
      ! Create auxiliary 1-block vectors
      call lsysbl_createVecFromScalar(rx, rxBlock)
      call lsysbl_createVecFromScalar(ry, ryBlock)

      ! Call user-defined assembly
      call fcb_calcDivVector(rafcstab, rxBlock, ryBlock, dscale, bclear,&
          fcb_calcFlux_sim, rcollection)

      ! Release auxiliary 1-block vectors
      call lsysbl_releaseVector(rxBlock)
      call lsysbl_releaseVector(ryBlock)
    else
      ! Set pointers
      call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
      call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
      call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
      call lsyssc_getbase_double(rx, p_Dx)
      call lsyssc_getbase_double(ry, p_Dy)

      ! Clear vector?
      if (bclear) call lsyssc_clearVector(ry)
      
      ! Assemble the divergence vector
      call doDivVector(p_IedgeListIdx, p_IedgeList, rafcstab%NEQ,&
          rx%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, p_Dy)
    end if
    
  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector

    subroutine doDivVector(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
      real(DP), intent(in) :: dscale
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, dimension(:), intent(in) :: IedgeListIdx
      integer, intent(in) :: NEQ,NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dy

      ! auxiliary arrays
      real(DP), dimension(:,:,:), pointer :: DdataAtEdge
      real(DP), dimension(:,:,:), pointer :: DfluxesAtEdge
      
      ! local variables
      integer :: IEDGEmax,IEDGEset,i,idx,iedge,igroup,j

      !$omp parallel default(shared)&
      !$omp private(DdataAtEdge,DfluxesAtEdge,IEDGEmax,i,idx,iedge,j)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
        end do
        !$omp end do

      end do ! igroup

      ! Deallocate temporal memory
      deallocate(DdataAtEdge)
      deallocate(DfluxesAtEdge)
      !$omp end parallel
      
    end subroutine doDivVector

  end subroutine gfsys_buildDivVectorScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecTVDBlock(rafcstab, rx, ndim, fcb_calcFlux_sim,&
      fcb_calcCharacteristics_sim, dscale, bclear, ry, rcollection)

!<description>
    ! This subroutine assembles the divergence vector for FEM-TVD schemes.
    ! If the vectors contain only one block, then the scalar counterpart
    ! of this routine is called with the scalar subvectors.
!</description>

!<input>
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
    include 'intf_calcFlux_sim.inc'

    ! callback function to compute local characteristics
    include 'intf_calcCharacteristics_sim.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
    !</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx


    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1) ) then
      call gfsys_buildDivVecTVDScalar(rafcstab, rx%RvectorBlock(1), ndim,&
          fcb_calcFlux_sim, fcb_calcCharacteristics_sim, dscale, bclear,&
          ry%RvectorBlock(1), rcollection)
      return
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDBlock')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(ry)


    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
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
          rafcstab%NEQ, rx%nblocks, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM2D)
      call doLimitTVD_2D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM3D)
      call doLimitTVD_3D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%nblocks, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D

    subroutine doLimitTVD_1D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)

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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single
      
      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
        
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D

    subroutine doLimitTVD_2D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)

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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))
      
      !$omp single
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
        
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
        
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D

    subroutine doLimitTVD_3D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)

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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))

      !$omp single
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
        
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)

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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IedgeList, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
        Daux1 = (DmatrixCoeffsAtEdge(idirection,2,idx)-&
                 DmatrixCoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DmatrixCoeffsAtEdge(idirection,1,idx)+&
                 DmatrixCoeffsAtEdge(idirection,2,idx))*&
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
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IedgeList, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
        Daux1 = (DmatrixCoeffsAtEdge(idirection,2,idx)-&
                 DmatrixCoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DmatrixCoeffsAtEdge(idirection,1,idx)+&
                 DmatrixCoeffsAtEdge(idirection,2,idx))*&
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

  end subroutine gfsys_buildDivVecTVDBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecTVDScalar(rafcstab, rx, ndim, fcb_calcFlux_sim,&
      fcb_calcCharacteristics_sim, dscale, bclear, ry, rcollection)

!<description>
    ! This subroutine assembles the divergence vector for FEM-TVD schemes
!</description>

!<input>
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
    include 'intf_calcFlux_sim.inc'

    ! callback function to compute local characteristics
    include 'intf_calcCharacteristics_sim.inc'
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx


    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecTVDScalar')
      call sys_halt()
    end if

    ! Clear vector?
    if (bclear) call lsyssc_clearVector(ry)

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_IedgeListIdx(rafcstab, p_IedgeListIdx)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
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
          rafcstab%NEQ, rx%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM2D)
      call doLimitTVD_2D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    case (NDIM3D)
      call doLimitTVD_3D(p_IedgeListIdx, p_IedgeList,&
          rafcstab%NEQ, rx%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
          p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, p_Dy)
    end select

    ! Set specifiers for Ps, Qs and Rs
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_NODELIMITER)

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 1D

    subroutine doLimitTVD_1D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single
      
      ! Allocate some temporal memory
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 2D

    subroutine doLimitTVD_2D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))

      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))

      !$omp single
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)

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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
    ! Assemble divergence vector for low-order operator plus
    ! algebraic flux correction of TVD-type in 3D

    subroutine doLimitTVD_3D(IedgeListIdx, IedgeList,&
        NEQ, NVAR, DmatrixCoeffsAtEdge, Dx, dscale,&
        Dpp, Dpm, Dqp, Dqm, Drp, Drm, Dy)

      ! input parameters
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DcharVariablesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DeigenvaluesAtEdge(NVAR,GFSYS_NEDGESIM))
      
      ! Clear P's and Q's (X-direction)
      !$omp single
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)

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
          call fcb_calcFlux_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      allocate(DrightEigenvectorsAtEdge(NVAR*NVAR,GFSYS_NEDGESIM))

      !$omp single
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      
      ! Clear P's and Q's (Y-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)

      ! Clear P's and Q's (Z-direction)
      call lalg_clearVector(Dpp)
      call lalg_clearVector(Dpm)
      call lalg_clearVector(Dqp)
      call lalg_clearVector(Dqm)
      !$omp end single

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
      Drp = gfsys_limit(Dpp, Dqp, 1.0_DP, 1.0_DP)
      Drm = gfsys_limit(Dpm, Dqm, 1.0_DP, 1.0_DP)
      !$omp end single
      
      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
        
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge,&
        DeigenvaluesAtEdge, IedgeList, Dpp, Dpm, Dqp, Dqm)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
        Daux1 = (DmatrixCoeffsAtEdge(idirection,2,idx)-&
                 DmatrixCoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DmatrixCoeffsAtEdge(idirection,1,idx)+&
                 DmatrixCoeffsAtEdge(idirection,2,idx))*&
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
        DmatrixCoeffsAtEdge, DcharVariablesAtEdge, DeigenvaluesAtEdge,&
        DrightEigenvectorsAtEdge, IedgeList, Dy)

      ! input parameters
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
        Daux1 = (DmatrixCoeffsAtEdge(idirection,2,idx)-&
                 DmatrixCoeffsAtEdge(idirection,1,idx))*&
                 DeigenvaluesAtEdge(:,idx)/2.0_DP
        Daux2 = (DmatrixCoeffsAtEdge(idirection,1,idx)+&
                 DmatrixCoeffsAtEdge(idirection,2,idx))*&
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

  end subroutine gfsys_buildDivVecTVDScalar

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecFCTBlock(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry, NVARtransformed,&
      fcb_calcFluxTransformation_sim, fcb_calcDiffTransformation_sim,&
      fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
      fcb_limitEdgewise, fcb_calcCorrection, rcollection)

!<description>
    ! This subroutine assembles the divergence vector for nonlinear
    ! FEM-FCT schemes.  If the vectors contain only one block, then
    ! the scalar counterpart of this routine is called with the scalar
    ! subvectors. Consider the documentation of subroutine
    ! 'gfsys_buildDivVecFCTScalar' for further details.
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
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorBlock), intent(inout) :: ry

    ! OPTIONAL collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_Dflux0,p_DfluxPrel
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: nvariable

    ! Check if block vectors contain only one block.
    if ((rx%nblocks .eq. 1) .and. (ry%nblocks .eq. 1)) then
      call gfsys_buildDivVecFCTScalar(&
          rafcstab, rmatrix, rx%RvectorBlock(1), dscale, bclear,&
          ioperationSpec, ry%RvectorBlock(1), NVARtransformed,&
          fcb_calcFluxTransformation_sim, fcb_calcDiffTransformation_sim,&
          fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
          fcb_limitEdgewise, fcb_calcCorrection, rcollection)
      return
    end if

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBLock')
      call sys_halt()
    end if

    ! Clear divergence vector?
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
    ! 6) Apply the limited antidifusive fluxes to the divergence
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
      if (rafcstab%ctypePrelimiting .ne. AFCSTAB_PRELIMITING_NONE) then
        
        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
          call sys_halt()
        end if
        
        ! Check if stabilisation provides edge-based structure
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide edge structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
          call sys_halt()
        end if
        
        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        
        if (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_STD) then
          ! Perform standard prelimiting
          call doStdPrelimitDble(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        elseif (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_MINMOD) then
          ! Perform minmod prelimiting
          call doMinModPrelimitDble(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        else
          call output_line('Invalid type of prelimiting!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        
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
            call doADIncrementsTransformedDble(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Standard routine without flux transformation
            call doADIncrementsDble(p_IedgeListIdx, p_IedgeList,&
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
            call doADIncrementsTransformedDble(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Compute antidiffusive increments without flux transformation
            call doADIncrementsDble(p_IedgeListIdx, p_IedgeList,&
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
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
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
        call doBoundsTransformedDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
            nvariable, p_Dx, p_Dqp, p_Dqm)
      else
        ! Standard routine without difference transformation
        call doBoundsDble(p_IedgeListIdx, p_IedgeList,&
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Set additional pointers
      call lsyssc_getbase_double(rmatrix, p_ML)

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, rcollection)
      elseif (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        ! Standard routine without constraints
        call doLimitNodalDble(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        ! Standard routine with constraints
        call doLimitNodalConstrainedDble(rafcstab%NEQ, nvariable, dscale, p_ML,&
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Compute edgewise correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux0)

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Drp, p_Drm, fcb_calcFluxTransformation_sim, p_Dflux0, rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseConstrTransfDble(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseConstrainedDble(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
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
          call doLimitEdgewiseTransformedDble(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseDble(p_IedgeList,&
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTBlock')
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
          call doCorrectScaleByMassDble(p_IedgeListIdx, p_IedgeList,&
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
          call doCorrectDble(p_IedgeListIdx, p_IedgeList,&
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
    subroutine doStdPrelimitDble(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)
      
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doStdPrelimitDble

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes using minmod limiter

#ifndef USE_OPENMP
    pure&
#endif
    subroutine doMinModPrelimitDble(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)
      
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
      
    end subroutine doMinModPrelimitDble

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and prelimiting

    subroutine doADIncrementsDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

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

    end subroutine doADIncrementsDble

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-defined
    ! set of variables prior to computing the sums

    subroutine doADIncrementsTransformedDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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

    end subroutine doADIncrementsTransformedDble

    !**************************************************************
    ! Assemble the local bounds from the predicted solution without
    ! transformation

    subroutine doBoundsDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
      
      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
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

    end subroutine doBoundsDble

    !**************************************************************
    ! Assemble the local bounds from the predicted solution which is
    ! transformed to a user-defined set of variables

    subroutine doBoundsTransformedDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
      
      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections
      
      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedDataAtEdge(NVARtransformed,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1
        
        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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

    end subroutine doBoundsTransformedDble

    !**************************************************************
    ! Compute the nodal correction factors without constraints

    subroutine doLimitNodalDble(NEQ, NVAR, dscale,&
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

    end subroutine doLimitNodalDble

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrainedDble(NEQ, NVAR, dscale,&
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
      
    end subroutine doLimitNodalConstrainedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewiseDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doLimitEdgewiseDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    ! which are transformed to a user-defined set of variables
    ! priori to computing the correction factors

    subroutine doLimitEdgewiseTransformedDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Loop over the edges
      !$omp do schedule(static,1)
      do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM

        ! We always handle GFSYS_NEDGESIM edges simultaneously.
        ! How many edges have we actually here?
        ! Get the maximum edge number, such that we handle 
        ! at most GFSYS_NEDGESIM edges simultaneously.
        
        IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)

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
     
    end subroutine doLimitEdgewiseTransformedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrainedDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doLimitEdgewiseConstrainedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes which are transformed
    ! to a user-defined set of variables priori to computing the
    ! correction factors

    subroutine doLimitEdgewiseConstrTransfDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes1AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes2AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      
      ! Loop over the edges
      !$omp do schedule(static,1)
      do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM

        ! We always handle GFSYS_NEDGESIM edges simultaneously.
        ! How many edges have we actually here?
        ! Get the maximum edge number, such that we handle 
        ! at most GFSYS_NEDGESIM edges simultaneously.
        
        IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)

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
      
    end subroutine doLimitEdgewiseConstrTransfDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrectDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge) 
 
          ! Apply limited antidiffusive fluxes
          Dy(i,:) = Dy(i,:) + F_ij
          Dy(j,:) = Dy(j,:) - F_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix

    subroutine doCorrectScaleByMassDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge) 

          ! Apply limited antidiffusive fluxes
          Dy(i,:) = Dy(i,:) + F_ij/ML(i)
          Dy(j,:) = Dy(j,:) - F_ij/ML(j)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel

    end subroutine doCorrectScaleByMassDble

  end subroutine gfsys_buildDivVecFCTBlock

  ! *****************************************************************************

!<subroutine>

  subroutine gfsys_buildDivVecFCTScalar(rafcstab, rmatrix, rx,&
      dscale, bclear, ioperationSpec, ry, NVARtransformed,&
      fcb_calcFluxTransformation_sim, fcb_calcDiffTransformation_sim,&
      fcb_calcADIncrements, fcb_calcBounds, fcb_limitNodal,&
      fcb_limitEdgewise, fcb_calcCorrection, rcollection)

!<description>
    ! This subroutine assembles the divergence vector for nonlinear
    ! FEM-FCT schemes.  Note that the vectors are required as scalar
    ! vectors which are stored in the interleave format. The idea of
    ! flux corrected transport can be traced back to the early SHASTA
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
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab

    ! divergence vector
    type(t_vectorScalar), intent(inout) :: ry

    ! OPTIONAL collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_ML,p_Dx,p_Dy
    real(DP), dimension(:), pointer :: p_Dpp,p_Dpm,p_Dqp,p_Dqm,p_Drp,p_Drm
    real(DP), dimension(:), pointer :: p_Dalpha,p_Dflux,p_Dflux0,p_DfluxPrel
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: nvariable


    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
      call sys_halt()
    end if

    ! Clear divergence vector?
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
    ! 6) Apply the limited antidifusive fluxes to the divergence
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
      if (rafcstab%ctypePrelimiting .ne. AFCSTAB_PRELIMITING_NONE) then
        
        ! Check if stabilisation provides raw antidiffusive fluxes
        if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
          call output_line('Stabilisation does not provide antidiffusive fluxes!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
          call sys_halt()
        end if
        
        ! Check if stabilisation provides edge-based structure
        if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
            (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
          call output_line('Stabilisation does not provide edge structure!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
          call sys_halt()
        end if
        
        ! Set additional pointer
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        
        if (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_STD) then
          ! Perform standard prelimiting
          call doStdPrelimitDble(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        elseif (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_MINMOD) then
          ! Perform minmod prelimiting
          call doMinModPrelimitDble(rafcstab%NEDGE, rafcstab%NVAR,&
              p_Dflux, p_DfluxPrel, p_Dalpha)
        else
          call output_line('Invalid type of prelimiting!',&
              OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Special treatment for semi-implicit FEM-FCT algorithm
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        
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
            call doADIncrementsTransformedDble(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_DfluxPrel, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Standard routine without flux transformation
            call doADIncrementsDble(p_IedgeListIdx, p_IedgeList,&
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
            call doADIncrementsTransformedDble(p_IedgeListIdx,&
                p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
                nvariable, p_Dx, p_Dflux, p_Dalpha, p_Dpp, p_Dpm)
          else
            ! Compute antidiffusive increments without flux transformation
            call doADIncrementsDble(p_IedgeListIdx, p_IedgeList,&
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
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
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
        call doBoundsTransformedDble(p_IedgeListIdx, p_IedgeList,&
            rafcstab%NEQ, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
            p_Dx, p_Dqp, p_Dqm)
      else
        ! Standard routine without difference transformation
        call doBoundsDble(p_IedgeListIdx, p_IedgeList,&
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if
      
      ! Set additional pointers
      call lsyssc_getbase_double(rmatrix, p_ML)

      ! Compute nodal correction factors
      if (present(fcb_limitNodal)) then
        ! User-supplied callback routine
        call fcb_limitNodal(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm, rcollection)
      elseif (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then
        ! Standard routine without constraints
        call doLimitNodalDble(rafcstab%NEQ, nvariable, dscale, p_ML,&
            p_Dpp, p_Dpm, p_Dqp, p_Dqm, p_Drp, p_Drm)
      else
        ! Standard routine with constraints
        call doLimitNodalConstrainedDble(rafcstab%NEQ, nvariable, dscale, p_ML,&
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if
      
      ! Compute edgewise correction factors
      if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then

        ! Special treatment for semi-implicit FEM-FCT algorithm
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_Dflux0)

        if (present(fcb_limitEdgewise)) then
          ! User-supplied callback routine
          call fcb_limitEdgewise(p_IedgeListIdx, p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              rafcstab%NVAR, rafcstab%NEQ, p_Dx, p_Dflux, p_Dalpha,&
              p_Drp, p_Drm, fcb_calcFluxTransformation_sim, p_Dflux0, rcollection)
        elseif (present(fcb_calcFluxTransformation_sim)) then
          ! Standard routine with flux transformation
          call doLimitEdgewiseConstrTransfDble(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseConstrainedDble(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR,&
              p_Dflux0, p_Dflux, p_Drp, p_Drm, p_Dalpha)
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
          call doLimitEdgewiseTransformedDble(p_IedgeList,&
              rafcstab%NEDGE, rafcstab%NEQ, rafcstab%NVAR, nvariable,&
              p_Dx, p_Dflux, p_Drp, p_Drm, p_Dalpha)
        else
          ! Standard routine without flux transformation
          call doLimitEdgewiseDble(p_IedgeList,&
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
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides raw antidiffusive fluxes
      if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
        call output_line('Stabilisation does not provide antidiffusive fluxes!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
        call sys_halt()
      end if

      ! Check if stabilisation provides edge-based structure
      if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
          (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
        call output_line('Stabilisation does not provide edge structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildDivVecFCTScalar')
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
          call doCorrectScaleByMassDble(p_IedgeListIdx, p_IedgeList,&
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
          call doCorrectDble(p_IedgeListIdx, p_IedgeList,&
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
    subroutine doStdPrelimitDble(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)
      
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doStdPrelimitDble

    !**************************************************************
    ! Prelimit the raw antidiffusive fluxes using minmod limiter
    
#ifndef USE_OPENMP
    pure&
#endif
    subroutine doMinModPrelimitDble(NEDGE, NVAR, Dflux, DfluxPrel, Dalpha)
      
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
      
    end subroutine doMinModPrelimitDble

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes without transformation and prelimiting

    subroutine doADIncrementsDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

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

    end subroutine doADIncrementsDble

    !**************************************************************
    ! Assemble the sums of antidiffusive increments for the given
    ! antidiffusive fluxes which are transformed to a user-defined
    ! set of variables prior to computing the sums

    subroutine doADIncrementsTransformedDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

      ! Clear P`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dpp)
      !$omp section
      call lalg_clearVector(Dpm)
      !$omp end sections

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DfluxesAtEdge(NVAR,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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

    end subroutine doADIncrementsTransformedDble
    
    !**************************************************************
    ! Assemble the local bounds from the predicted solution
    ! without transformation

    subroutine doBoundsDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
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

    end subroutine doBoundsDble

    !**************************************************************
    ! Assemble local bounds from the predicted solution
    ! which is transformed to a user-defined set of variables

    subroutine doBoundsTransformedDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
      
      ! Clear Q`s
      !$omp sections
      !$omp section
      call lalg_clearVector(Dqp)
      !$omp section
      call lalg_clearVector(Dqm)
      !$omp end sections

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedDataAtEdge(NVARtransformed,GFSYS_NEDGESIM))

      ! Loop over the edge groups and process all edges of one group
      ! in parallel without the need to synchronize memory access
      do igroup = 1, size(IedgeListIdx)-1

        ! Do nothing for empty groups
        if (IedgeListIdx(igroup+1)-IedgeListIdx(igroup) .le. 0) cycle

        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = IedgeListIdx(igroup),&
                      IedgeListIdx(igroup+1)-1, GFSYS_NEDGESIM

          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(IedgeListIdx(igroup+1)-1, IEDGEset-1+GFSYS_NEDGESIM)
          
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

    end subroutine doBoundsTransformedDble

    !**************************************************************
    ! Compute nodal correction factors without constraints

    subroutine doLimitNodalDble(NEQ, NVAR, dscale,&
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
      
    end subroutine doLimitNodalDble

    !**************************************************************
    ! Compute nodal correction factors with constraints

    subroutine doLimitNodalConstrainedDble(NEQ, NVAR, dscale,&
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

    end subroutine doLimitNodalConstrainedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes

    subroutine doLimitEdgewiseDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doLimitEdgewiseDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of antidiffusive fluxes
    ! which are transformed to a user-defined set of variables
    ! priori to computing the correction factors

    subroutine doLimitEdgewiseTransformedDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

      ! Loop over the edges
      !$omp do schedule(static,1)
      do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM

        ! We always handle GFSYS_NEDGESIM edges simultaneously.
        ! How many edges have we actually here?
        ! Get the maximum edge number, such that we handle 
        ! at most GFSYS_NEDGESIM edges simultaneously.
        
        IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)

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
      
    end subroutine doLimitEdgewiseTransformedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes

    subroutine doLimitEdgewiseConstrainedDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doLimitEdgewiseConstrainedDble

    !**************************************************************
    ! Compute edgewise correction factors based on the precomputed
    ! nodal correction factors and the sign of a pair of explicit
    ! and implicit raw antidiffusive fluxes which are transformed
    ! to a user-defined set of variables priori to computing the
    ! correction factors

    subroutine doLimitEdgewiseConstrTransfDble(IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

      ! Allocate temporal memory
      allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes1AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      allocate(DtransformedFluxes2AtEdge(NVARtransformed,2,GFSYS_NEDGESIM))
      
      ! Loop over the edges
      !$omp do schedule(static,1)
      do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM

        ! We always handle GFSYS_NEDGESIM edges simultaneously.
        ! How many edges have we actually here?
        ! Get the maximum edge number, such that we handle 
        ! at most GFSYS_NEDGESIM edges simultaneously.
        
        IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)

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
      
    end subroutine doLimitEdgewiseConstrTransfDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them

    subroutine doCorrectDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)
          
          ! Apply limited antidiffusive fluxes
          Dy(:,i) = Dy(:,i) + F_ij
          Dy(:,j) = Dy(:,j) - F_ij
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel
      
    end subroutine doCorrectDble

    !**************************************************************
    ! Correct the antidiffusive fluxes and apply them
    ! scaled by the inverse of the lumped mass matrix

    subroutine doCorrectScaleByMassDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)

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
          F_ij = dscale * Dalpha(iedge) * Dflux(:,iedge)
          
          ! Apply limited antidiffusive fluxes
          Dy(:,i) = Dy(:,i) + F_ij/ML(i)
          Dy(:,j) = Dy(:,j) - F_ij/ML(j)
        end do
        !$omp end do

      end do ! igroup
      !$omp end parallel
    end subroutine doCorrectScaleByMassDble

  end subroutine gfsys_buildDivVecFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildFluxFCTBlock(rafcstab, rx, fcb_calcFluxFCT_sim,&
      theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
      rmatrix, rxTimeDeriv, rxPredictor, rcollection)

!<description>
    ! This subroutine assembles the raw antidiffusive fluxes for
    ! FEM-FCT schemes.  If the vectors contain only one block, then
    ! the scalar counterpart of this routine is called with the scalar
    ! subvectors.
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

    ! Callback functions to compute antidiffusive fluxes
    include 'intf_calcFluxFCT_sim.inc'

    ! OPTIONAL: mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorBlock), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorBlock), intent(in), optional :: rxPredictor
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
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    integer :: nblocks
    

    ! Check if block vector(s) contains exactly one block
    nblocks = rx%nblocks
    if (present(rxTimeDeriv)) nblocks = max(nblocks, rxTimeDeriv%nblocks)
    if (present(rxPredictor)) nblocks = max(nblocks, rxPredictor%nblocks)
    
    if (nblocks .eq. 1) then
      ! Call subroutine for scalar vectors
      if (present(rxTimeDeriv)) then
        if (present(rxPredictor)) then
          ! ... both approximate time derivative and predictor are present
          call gfsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              fcb_calcFluxFCT_sim, theta, tstep, dscale, bclear,&
              bquickAssembly, ioperationSpec, rmatrix,&
              rxTimeDeriv%RvectorBlock(1), rxPredictor%RvectorBlock(1),&
              rcollection=rcollection)
        else
          ! ... only the approximate time derivative is present
          call gfsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              fcb_calcFluxFCT_sim, theta, tstep, dscale, bclear,&
              bquickAssembly, ioperationSpec, rmatrix,&
              rxTimeDeriv%RvectorBlock(1), rcollection=rcollection)
        end if
      else
        if (present(rxPredictor)) then
          ! ... only the predictor is present
          call gfsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              fcb_calcFluxFCT_sim, theta, tstep, dscale, bclear,&
              bquickAssembly, ioperationSpec, rmatrix,&
              rxPredictor=rxPredictor%RvectorBlock(1), rcollection=rcollection)
        else
          ! ... neither the approximate time derivative nor the predictor is present
          call gfsys_buildFluxFCTScalar(rafcstab, rx%RvectorBlock(1),&
              fcb_calcFluxFCT_sim, theta, tstep, dscale, bclear,&
              bquickAssembly, ioperationSpec, rmatrix, rcollection=rcollection)
        end if
      end if

      ! That`s it
      return
    end if

    !---------------------------------------------------------------------------
    
    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0) .and.&
        (rafcstab%ctypeAFCstabilisation .ne. AFCSTAB_LINFCT_MASS)) then
      call output_line('Stabilisation does not provide edge data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix (if present)
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ) .or.&
          (rafcstab%cmatrixFormat .ne. rmatrix%cmatrixFormat)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
        call sys_halt()
      end if
    end if

    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call lsysbl_getbase_double(rx, p_Dx)
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_NLINFCT_IMPLICIT)

      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)

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
          ! $$ f_{ij}^n = (1-\theta)\Delta t d_{ij}^n(u_i^n-u_j^n) $$
          call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale*(1.0_DP-theta),&
            bclear, p_Dflux0)
        elseif (.not.bquickAssembly .and. bclear) then
          ! Clear the explicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij}^n = 0 $$
          call lalg_clearVector(p_Dflux0, rafcstab%NEDGE)
          ! if bquickAssembly = TRUE then this step can be skipped
        end if

        !-----------------------------------------------------------------------

        ! Check for special treatment
        if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)
          call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

          ! We have to store the raw-antidiffusive fluxes based on the
          ! initial solution without contribution of the consistent
          ! mass matrix and without scaling by the implicitness parameter
          ! $$ f_{ij} = \Delta t d_{ij}^n(u_i^n-u_j^n) $$
          call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
              .true., p_DfluxPrel)

        elseif (rafcstab%ctypePrelimiting .ne. AFCSTAB_PRELIMITING_NONE) then
          
          ! We have to assemble the raw-antidiffusive fluxes for
          ! prelimiting separately based on the low-order predictor
          if (present(rxPredictor)) then
            
            ! Set pointers
            call lsysbl_getbase_double(rxPredictor, p_DxPredictor)
            call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
            
            if (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_STD) then
              ! Compute solution difference for standard prelimiting
              call doDifferencesDble(p_IedgeList, rafcstab%NEDGE,&
                  rafcstab%NEQ, rafcstab%NVAR, p_DxPredictor, p_DfluxPrel)
            elseif (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_MINMOD) then
              ! Compute fluxes for minmod prelimiting
              call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                  rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_DxPredictor, dscale,&
                  .true., p_DfluxPrel)
            else
              call output_line('Invalid type of prelimiting!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
              call sys_halt()
            end if
          else
            call output_line('Fluxes for prelimiting cannot be assembled without predictor!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
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
              rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, -dscale/tstep,&
              .false., p_Dflux0)
        end if

      end if
      

      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_IMPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble implicit part of raw-antidiffusive fluxes
        !-----------------------------------------------------------------------
        
        if ((rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_ITERATIVE) .and.&
            iand(ioperationSpec, AFCSTAB_FCTFLUX_REJECTED) .ne. 0) then
          !---------------------------------------------------------------------
          ! Apply the rejected antidiffusive fluxes from the previous limiting
          ! step to the implicit part of the raw-antidiffusive fluxes
          ! --------------------------------------------------------------------
          
          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
            call output_line('Stabilisation does not provide correction factors!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
            call sys_halt()
          end if
          
          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
            call sys_halt()
          end if
          
          ! Set pointer
          call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
          
          ! Subtract amount of rejected antidiffusion
          call gfsys_combineFluxesDble(rafcstab%NVAR, rafcstab%NEDGE, -1.0_DP,&
              p_Dflux, p_Dflux0, p_Dalpha)
        end if

        !-----------------------------------------------------------------------
        
        if (theta .ne. 0.0_DP) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij} = \theta\Delta t d_{ij}(u_i-u_j) $$
          call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale*theta,&
              bclear, p_Dflux)
        end if
        
        if (bquickAssembly) then
          ! We may check of either the implicit or explicit part are
          ! missing so that some redundant computations may be skipped
          if (theta .ne. 1.0_DP) then
            ! The explicit part of the raw-antidiffusive fluxes exists
            if (theta .ne. 0.0_DP) then
              ! The implicit part of the raw-antidiffusive fluxes
              ! exists; so combine them both into common fluxes
              call gfsys_combineFluxesDble(rafcstab%NVAR, rafcstab%NEDGE,&
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
          call gfsys_combineFluxesDble(rafcstab%NVAR, rafcstab%NEDGE,&
              1.0_DP, p_Dflux0, p_Dflux)
        end if
        
        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the implicit part of the mass-antidiffusive fluxes
          ! $$ f_{ij}^m := f_{ij}^m + m_{ij}(u_i^m-u_j^m) $$
          call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
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

      ! Assemble spatial part of raw-antidiffusive fluxes
      call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
          rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, bclear, p_Dflux)

      !-------------------------------------------------------------------------

      if (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_STD) then
        ! Compute fluxes for standard prelimiting based on the
        ! low-order solution which serves as predictor
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        call doDifferencesDble(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_DfluxPrel)
        
      elseif (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_MINMOD) then
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
        call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
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
        call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, dscale, .true., p_Dflux)

        ! Set specifiers for raw antidiffusive fluxes
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)
        
      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
        call sys_halt()
      end if
      

    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTBlock')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without
    ! contribution of the consistent mass matrix.

    subroutine doFluxesDble(IedgeList, NEDGE, NEQ, NVAR,&
        DmatrixCoeffsAtEdge, Dx, dscale, bclear, Dflux)

      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        
        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFluxFCT_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DfluxAtEdge(NVAR,GFSYS_NEDGESIM))
        
        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFluxFCT_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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

    end subroutine doFluxesDble
    
    !**************************************************************
    ! Assemble fluxes for classical prelimiting.

      subroutine doDifferencesDble(IedgeList, NEDGE, NEQ, NVAR, Dx, Dflux)
      
      real(DP), dimension(NEQ,NVAR), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      
      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      integer :: iedge,i,j

      !$omp parallel do default(shared) private(i,j)&
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doDifferencesDble

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix Dmatrix

    subroutine doFluxesByMatrixDble(IedgeList, NEDGE, NEQ, NVAR,&
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doFluxesByMatrixDble

  end subroutine gfsys_buildFluxFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_buildFluxFCTScalar(rafcstab, rx, fcb_calcFluxFCT_sim,&
      theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
      rmatrix, rxTimeDeriv, rxPredictor, rcollection)

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

    ! Callback functions to compute antidiffusive fluxes
    include 'intf_calcFluxFCT_sim.inc'

    ! OPTIONAL: mass matrix
    type(t_matrixScalar), intent(in), optional :: rmatrix

    ! OPTIONAL: approximate time derivative of vector rx
    type(t_vectorScalar), intent(in), optional :: rxTimeDeriv

    ! OPTIONAL: low-order predictor of vector rx
    ! This vector is required to assemble the fluxes for prelimiting
    ! in some variants of the FCT algorithm.
    type(t_vectorScalar), intent(in), optional :: rxPredictor
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
    real(DP), dimension(:,:,:), pointer :: p_DmatrixCoeffsAtEdge
    integer, dimension(:,:), pointer :: p_IedgeList
    

    ! Check if stabilisation is prepared
    if (iand(rafcstab%istabilisationSpec, AFCSTAB_INITIALISED) .eq. 0) then
      call output_line('Stabilisation has not been initialised!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation provides edge-based data structures
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_MATRIXCOEFFS)  .eq. 0) .and.&
        (rafcstab%ctypeAFCstabilisation .ne. AFCSTAB_LINFCT_MASS)) then
      call output_line('Stabilisation does not provide edge-based data structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end if

    ! Check if stabilisation is compatible with matrix (if present)
    if (present(rmatrix)) then
      if ((rafcstab%NEQ       .ne. rmatrix%NEQ) .or.&
          (rafcstab%NEDGE * 2 .ne. rmatrix%NA-rmatrix%NEQ) .or.&
          (rafcstab%cmatrixFormat .ne. rmatrix%cmatrixFormat)) then
        call output_line('Matrix is not compatible with stabilisation structure!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
        call sys_halt()
      end if
    end if
    
    ! Set pointers
    call afcstab_getbase_IedgeList(rafcstab, p_IedgeList)
    call afcstab_getbase_DmatCoeffAtEdge(rafcstab, p_DmatrixCoeffsAtEdge)
    call lsyssc_getbase_double(rx, p_Dx)
    
    ! What kind of stabilisation are we?
    select case(rafcstab%ctypeAFCstabilisation)

    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_NLINFCT_IMPLICIT)
      
      ! Set pointers
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(rafcstab%p_rvectorFlux0, p_Dflux0)

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
          ! $$ f_{ij}^n = (1-\theta)\Delta t d_{ij}^n(u_i^n-u_j^n) $$
          call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
            rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale*(1.0_DP-theta),&
            bclear, p_Dflux0)
        elseif (.not.bquickAssembly .and. bclear) then
          ! Clear the explicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij}^n = 0 $$
          call lalg_clearVector(p_Dflux0, rafcstab%NEDGE)
          ! if bquickAssembly = TRUE then this step can be skipped
        end if

        !-----------------------------------------------------------------------

        ! Check for special treatment
        if (rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_IMPLICIT) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)
          call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)

          ! We have to store the raw-antidiffusive fluxes based on the
          ! initial solution without contribution of the consistent
          ! mass matrix and without scaling by the implicitness parameter
          ! $$ f_{ij} = \Delta t d_{ij}^n(u_i^n-u_j^n) $$
          call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale,&
              .true., p_DfluxPrel)

        elseif (rafcstab%ctypePrelimiting .ne. AFCSTAB_PRELIMITING_NONE) then
          
          ! We have to assemble the raw-antidiffusive fluxes for
          ! prelimiting separately based on the low-order predictor
          if (present(rxPredictor)) then
            
            ! Set pointers
            call lsyssc_getbase_double(rxPredictor, p_DxPredictor)
            call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
            
            if (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_STD) then
              ! Compute solution difference for standard prelimiting
              call doDifferencesDble(p_IedgeList, rafcstab%NEDGE,&
                  rafcstab%NEQ, rafcstab%NVAR, p_DxPredictor, p_DfluxPrel)
            elseif (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_MINMOD) then
              ! Compute fluxes for minmod prelimiting
              call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
                  rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_DxPredictor, dscale,&
                  .true., p_DfluxPrel)
            else
              call output_line('Invalid type of prelimiting!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
              call sys_halt()
            end if
          else
            call output_line('Fluxes for prelimiting cannot be assembled without predictor!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
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
              rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, -dscale/tstep,&
              .false., p_Dflux0)
        end if

      end if
      

      if (iand(ioperationSpec, AFCSTAB_FCTFLUX_IMPLICIT) .ne. 0) then
        !-----------------------------------------------------------------------
        ! Assemble implicit part of raw-antidiffusive fluxes
        !-----------------------------------------------------------------------
        
        if ((rafcstab%ctypeAFCstabilisation .eq. AFCSTAB_NLINFCT_ITERATIVE) .and.&
            iand(ioperationSpec, AFCSTAB_FCTFLUX_REJECTED) .ne. 0) then
          !---------------------------------------------------------------------
          ! Apply the rejected antidiffusive fluxes from the previous limiting
          ! step to the implicit part of the raw-antidiffusive fluxes
          ! --------------------------------------------------------------------
          
          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER) .eq. 0) then
            call output_line('Stabilisation does not provide correction factors!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
            call sys_halt()
          end if
          
          ! Check if stabilisation provides raw antidiffusive fluxes
          if (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES) .eq. 0) then
            call output_line('Stabilisation does not provide antidiffusive fluxes!',&
                OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
            call sys_halt()
          end if
          
          ! Set pointer
          call lsyssc_getbase_double(rafcstab%p_rvectorAlpha, p_Dalpha)
          
          ! Subtract amount of rejected antidiffusion
          call gfsys_combineFluxesDble(rafcstab%NVAR, rafcstab%NEDGE, -1.0_DP,&
              p_Dflux, p_Dflux0, p_Dalpha)
        end if

        !-----------------------------------------------------------------------
        
        if (theta .ne. 0.0_DP) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! $$ f_{ij} = \theta\Delta t d_{ij}(u_i-u_j) $$
          call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
              rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale*theta,&
              bclear, p_Dflux)
        end if
        
        if (bquickAssembly) then
          ! We may check of either the implicit or explicit part are
          ! missing so that some redundant computations may be skipped
          if (theta .ne. 1.0_DP) then
            ! The explicit part of the raw-antidiffusive fluxes exists
            if (theta .ne. 0.0_DP) then
              ! The implicit part of the raw-antidiffusive fluxes
              ! exists; so combine them both into common fluxes
              call gfsys_combineFluxesDble(rafcstab%NVAR, rafcstab%NEDGE,&
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
          call gfsys_combineFluxesDble(rafcstab%NVAR, rafcstab%NEDGE,&
              1.0_DP, p_Dflux0, p_Dflux)
        end if
        
        !-----------------------------------------------------------------------

        ! Do we have to include mass-antidiffuion?
        if (present(rmatrix)) then
          
          ! Set pointers
          call lsyssc_getbase_double(rmatrix, p_Dmatrix)

          ! Assemble the implicit part of the mass-antidiffusive fluxes
          ! $$ f_{ij}^m := f_{ij}^m + m_{ij}(u_i^m-u_j^m) $$
          call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
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

      ! Assemble spatial part of raw-antidiffusive fluxes
      call doFluxesDble(p_IedgeList, rafcstab%NEDGE, rafcstab%NEQ,&
          rafcstab%NVAR, p_DmatrixCoeffsAtEdge, p_Dx, dscale, bclear, p_Dflux)

      !-------------------------------------------------------------------------

      if (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_STD) then
        ! Compute fluxes for standard prelimiting based on the
        ! low-order solution which serves as predictor
        call lsyssc_getbase_double(rafcstab%p_rvectorFluxPrel, p_DfluxPrel)
        call doDifferencesDble(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dx, p_DfluxPrel)
        
      elseif (rafcstab%ctypePrelimiting .eq. AFCSTAB_PRELIMITING_MINMOD) then
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
        call doFluxesByMatrixDble(p_IedgeList, rafcstab%NEDGE,&
            rafcstab%NEQ, rafcstab%NVAR, p_Dmatrix, p_Dx, dscale, .true., p_Dflux)
        
        ! Set specifiers for raw antidiffusive fluxes
        rafcstab%istabilisationSpec =&
            ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)

      else
        call output_line('Unable to compute mass antidiffusion without mass matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
        call sys_halt()
      end if

      
    case default
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_buildFluxFCTScalar')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Assemble raw antidiffusive fluxes without
    ! contribution of the consistent mass matrix

    subroutine doFluxesDble(IedgeList, NEDGE, NEQ, NVAR,&
        DmatrixCoeffsAtEdge, Dx, dscale, bclear, Dflux)
      
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      real(DP), dimension(:,:,:), intent(in) :: DmatrixCoeffsAtEdge
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
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        
        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFluxFCT_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DfluxAtEdge(NVAR,GFSYS_NEDGESIM))
        
        ! Loop over the edges
        !$omp do schedule(static,1)
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle 
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
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
          call fcb_calcFluxFCT_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              DmatrixCoeffsAtEdge(:,:,IEDGEset:IEDGEmax),&
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

    end subroutine doFluxesDble

    !**************************************************************
    ! Assemble fluxes for classical prelimiting.

#ifndef USE_OPENMP
    pure&
#endif
        subroutine doDifferencesDble(IedgeList, NEDGE, NEQ, NVAR, Dx, Dflux)
      
      real(DP), dimension(NVAR,NEQ), intent(in) :: Dx
      integer, dimension(:,:), intent(in) :: IedgeList
      integer, intent(in) :: NEDGE,NEQ,NVAR
      
      real(DP), dimension(NVAR,NEDGE), intent(out) :: Dflux

      ! local variables
      integer :: iedge,i,j

      !$omp parallel do default(shared) private(i,j)&
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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

    end subroutine doDifferencesDble

    !**************************************************************
    ! Assemble raw antidiffusive fluxes using the coefficients
    ! supplied by the CSR-matrix Dmatrix

    subroutine doFluxesByMatrixDble(IedgeList, NEDGE, NEQ, NVAR,&
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
          !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
      
    end subroutine doFluxesByMatrixDble

  end subroutine gfsys_buildFluxFCTScalar

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_failsafeFCTBlock(rafcstab, rmatrix, rx, dscale, dtol,&
      ioperationSpec, bisAccepted, dfactor, nsteps, CvariableNames,&
      fcb_extractVariable, fcb_calcFailsafe, rxBackup, rvectorCorr,&
      rvectorTmp, rcollection)

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

    ! Check if block vector contains only one block.
    if (rx%nblocks .eq. 1) then
      if (present(rxBackup)) then
        call gfsys_failsafeFCTScalar(&
            rafcstab, rmatrix, rx%RvectorBlock(1), dscale, dtol,&
            ioperationSpec, bisAccepted, dfactor, nsteps,&
            CvariableNames, fcb_extractVariable, fcb_calcFailsafe,&
            rxBackup%RvectorBlock(1),&
            rvectorCorr, rvectorTmp, rcollection)
      else
        call gfsys_failsafeFCTScalar(&
            rafcstab, rmatrix, rx%RvectorBlock(1), dscale, dtol,&
            ioperationSpec, bisAccepted, dfactor, nsteps,&
            CvariableNames, fcb_extractVariable, fcb_calcFailsafe,&
            rvectorCorr=rvectorCorr, rvectorTmp=rvectorTmp,&
            rcollection=rcollection)
      end if
      return
    end if
    
    
    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)     .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_failsafeFCTBlock')
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
      call doBoundsDble(p_IedgeListIdx, p_IedgeList,&
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
          call doFailsafeLimitDble(p_IedgeList, rafcstab%NEDGE,&
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
        call doFailsafeLimitDble(p_IedgeList, rafcstab%NEDGE,&
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
    
    subroutine doBoundsDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
   
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
          
          ! Compute minimum/maximum value of neighboring nodes
          Dlbound(i,:) = min(Dlbound(i,:), Dx(j,:))
          Dlbound(j,:) = min(Dlbound(j,:), Dx(i,:))
          Dubound(i,:) = max(Dubound(i,:), Dx(j,:))
          Dubound(j,:) = max(Dubound(j,:), Dx(i,:))
        end do
        !$omp end do
        
      end do ! igroup
      !$omp end parallel
      
    end subroutine doBoundsDble

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
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        
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
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        
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
    
    subroutine doFailsafeLimitDble(IedgeList, NEDGE, NEQ, NVARfailsafe,&
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

    end subroutine doFailsafeLimitDble
      
  end subroutine gfsys_failsafeFCTBlock

  !*****************************************************************************

!<subroutine>

  subroutine gfsys_failsafeFCTScalar(rafcstab, rmatrix, rx, dscale, dtol,&
      ioperationSpec, bisAccepted, dfactor, nsteps, CvariableNames,&
      fcb_extractVariable, fcb_calcFailsafe, rxBackup, rvectorCorr,&
      rvectorTmp, rcollection)

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

    ! Check if stabilisation provides edge-based structure
    if ((iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGESTRUCTURE)   .eq. 0) .and.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_ADFLUXES)        .eq. 0) .or.&
        (iand(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)     .eq. 0)) then
      call output_line('Stabilisation does not provide required structures!',&
          OU_CLASS_ERROR,OU_MODE_STD,'gfsys_failsafeFCTScalar')
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
      call doBoundsDble(p_IedgeListIdx, p_IedgeList,&
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
          call doFailsafeLimitDble(p_IedgeList, rafcstab%NEDGE,&
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
        call doFailsafeLimitDble(p_IedgeList, rafcstab%NEDGE,&
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
    
    subroutine doBoundsDble(IedgeListIdx, IedgeList,&
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
      !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
   
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
          
          ! Compute minimum/maximum value of neighboring nodes
          Dlbound(i,:) = min(Dlbound(i,:), Dx(j,:))
          Dlbound(j,:) = min(Dlbound(j,:), Dx(i,:))
          Dubound(i,:) = max(Dubound(i,:), Dx(j,:))
          Dubound(j,:) = max(Dubound(j,:), Dx(i,:))
        end do
        !$omp end do
        
      end do ! igroup
      !$omp end parallel
      
    end subroutine doBoundsDble

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
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        
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
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        
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
    
    subroutine doFailsafeLimitDble(IedgeList, NEDGE, NEQ, NVARfailsafe,&
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

    end subroutine doFailsafeLimitDble

  end subroutine gfsys_failsafeFCTScalar

  ! ****************************************************************************
  ! Here, some private auxiliary subroutine follow
  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
    pure&
#endif
    subroutine gfsys_combineFluxesDble(NVAR, NEDGE, dscale, Dflux1, Dflux2, Dalpha)

!<description>
    ! This subroutine combines the two fluxes:
    ! Dflux2 := Dflux2 + dscale * Dalpha * Dflux1
!</description>

!<input>
    ! First flux
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux1

    ! Individual scaling factor for each entry of flux1
    real(DP), dimension(:), intent(in), optional :: Dalpha

    ! Global scaling factor for all entries of flux1
    real(DP), intent(in) :: dscale

    ! Number of entries/edges
    integer, intent(in) :: NEDGE

    ! Number of variables
    integer, intent(in) :: NVAR
!</input>

!<inputoutput>
    ! Second flux
    real(DP), dimension(NVAR,NEDGE), intent(inout) :: Dflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge
    
    if (present(Dalpha)) then

      if (dscale .eq. 1.0_DP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          + Dalpha(iedge) * Dflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (dscale .eq. -1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          - Dalpha(iedge) * Dflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          + Dflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (dscale .eq. -1.0_DP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          - Dflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Dflux2(:,iedge) = Dflux2(:,iedge)&
                          + dscale * Dflux1(:,iedge)
        end do
        !$omp end parallel do
        
      end if

    end if

  end subroutine gfsys_combineFluxesDble

  ! ****************************************************************************

!<subroutine>

#ifndef USE_OPENMP
    pure&
#endif
    subroutine gfsys_combineFluxesSngl(NVAR, NEDGE, fscale, Fflux1, Fflux2, Falpha)

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
!</input>

!<inputoutput>
    ! Second flux
    real(SP), dimension(NVAR,NEDGE), intent(inout) :: Fflux2
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: iedge
    
    if (present(Falpha)) then

      if (fscale .eq. 1.0_SP) then
        
        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          + Falpha(iedge) * Fflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (fscale .eq. -1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          - Falpha(iedge) * Fflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
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
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          + Fflux1(:,iedge)
        end do
        !$omp end parallel do

      elseif (fscale .eq. -1.0_SP) then

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          - Fflux1(:,iedge)
        end do
        !$omp end parallel do

      else

        ! Loop over all edges
        !$omp parallel do default(shared)&
        !$omp if (NEDGE > GFSYS_NEDGEMIN_OMP)
        do iedge = 1, NEDGE
          Fflux2(:,iedge) = Fflux2(:,iedge)&
                          + fscale * Fflux1(:,iedge)
        end do
        !$omp end parallel do
        
      end if

    end if

  end subroutine gfsys_combineFluxesSngl

  !*****************************************************************************

!<function>
  
  elemental function gfsys_limitUnboundedDble(p, q, dval) result(r)

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
  end function gfsys_limitUnboundedDble
  
  !*****************************************************************************

!<function>
  
  elemental function gfsys_limitUnboundedSngl(p, q, fval) result(r)

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
  end function gfsys_limitUnboundedSngl

  !*****************************************************************************

!<function>
  
  elemental function gfsys_limitBoundedDble(p, q, dval, dbound) result(r)

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
  end function gfsys_limitBoundedDble

  !*****************************************************************************

!<function>
  
  elemental function gfsys_limitBoundedSngl(p, q, fval, fbound) result(r)

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
  end function gfsys_limitBoundedSngl

end module groupfemsystem
